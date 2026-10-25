"""The circular (spherical) latent factor model and its MCMC sampler.

Model (Yu & Rodriguez 2021, eqs. 4-9), for legislators i = 1..I and roll
calls j = 1..J with angles on (-pi, pi]:

    y_ij | .     ~ Bernoulli(G_kappa_j(e_ij)),
    e_ij         = rho(zeta_j, beta_i)^2 - rho(psi_j, beta_i)^2,
    rho(a, b)    = arccos(cos(a - b))      (geodesic distance on the circle)
    beta_i       ~ vonMises(mu, omega),     omega ~ Gamma(omega_shape, omega_rate)
    psi_j, zeta_j ~ Uniform(circle)          ("Yea" / "Nay" positions)
    kappa_j      ~ Gamma(kappa_shape, rate=xi),  xi ~ Gamma(xi_shape, xi_rate)

Missing votes are treated as missing completely at random and simply drop out
of the likelihood (the R code imputed them every 50 iterations, which targets
the same posterior for the parameters).

Each sweep reproduces the original Gibbs scheme:

    kappa_j  random-walk MH on log scale, proposal sd tuned during burn-in
    xi       conjugate Gibbs draw
    beta_i   geodesic HMC, all legislators in parallel
    psi_j    geodesic HMC, all roll calls in parallel
    zeta_j   geodesic HMC, all roll calls in parallel
    omega    random-walk MH on log scale
    delta    exact draw of a common rotation of all angles (new; see below)

The conditionals within each GHMC block are independent across coordinates,
so the block is one vectorised GHMC call (see ``ghmc.py``).  The leapfrog
gradients are computed in ``grad_dtype`` (float32 by default) while the
accept/reject step uses float64, which keeps the chain exact (see ``ghmc.py``).
"""

from __future__ import annotations

import time
from dataclasses import asdict, dataclass, field
from functools import partial
from typing import NamedTuple

import jax
import jax.numpy as jnp
import numpy as np
from jax.scipy.special import i0e

from .ghmc import angle_to_circle, circle_to_angle, ghmc_step, sample_von_mises
from .special import KAPPA_MAX, log_link_cdf

jax.config.update("jax_enable_x64", True)


@dataclass(frozen=True)
class Priors:
    """Hyperparameters; defaults are the paper's (Section 4.3)."""

    mu: float = 0.0  # prior mean direction of the ideal points (tau_beta)
    omega_shape: float = 1.0  # omega ~ Gamma(shape, rate): exponential, mean 10
    omega_rate: float = 0.1
    kappa_shape: float = 1.0  # kappa_j ~ Gamma(shape, rate=xi): exponential, mean 1/xi
    xi_shape: float = 1.0  # xi = 1/lambda ~ Gamma(shape, rate): exponential, mean 1/25
    xi_rate: float = 25.0


@dataclass(frozen=True)
class SamplerConfig:
    n_samples: int = 20_000
    burnin: int = 30_000
    thin: int = 1
    jitter_every: int = 50  # redraw step sizes / leapfrog counts this often
    beta_step: tuple[float, float] = (0.01, 0.04)  # uniform range for beta step sizes
    tau_step: tuple[float, float] = (0.01, 0.105)  # ... for psi and zeta
    beta_max_leap: int = 10  # leapfrog steps ~ U{1..max_leap}
    tau_max_leap: int = 10
    kappa_prop_sd: float = 0.5  # initial log-scale RW sd for kappa_j
    omega_prop_sd: float = 0.1
    # During burn-in, rescale each coordinate's step-size range so its GHMC
    # acceptance lands in [0.6, 0.9] (the paper tuned ranges by hand to hit
    # this window).  Frozen after burn-in.
    adapt_steps: bool = True
    grad_dtype: str = "float32"  # "float64" for exact-precision gradients
    # Exact Gibbs draw of a common rotation of (beta, psi, zeta), to which the
    # likelihood is invariant; without it the whole configuration drifts very
    # slowly around the circle (only the von Mises prior on beta pins it).
    # Off by default: it is exact, but on the 116th House every chain with it
    # converged more slowly in log-likelihood than every chain without it.
    rotation_move: bool = False
    seed: int = 0
    print_every: int = 1000  # iterations; 0 to silence


class Data(NamedTuple):
    sign: jax.Array  # (I, J) +1 for Yea, -1 for Nay, 0 if missing
    mask: jax.Array  # (I, J) bool, True where observed


class State(NamedTuple):
    beta: jax.Array  # (I,)
    yes: jax.Array  # (J,) psi
    no: jax.Array  # (J,) zeta
    kappa: jax.Array  # (J,)
    xi: jax.Array  # ()   1 / lambda
    omega: jax.Array  # ()
    ll: jax.Array  # (I, J) log-likelihood matrix at the current state (0 if missing)
    kappa_sd: jax.Array  # (J,) tuning
    beta_scale: jax.Array  # (I,) step-size multipliers (adapted in burn-in)
    yes_scale: jax.Array  # (J,)
    no_scale: jax.Array  # (J,)


class Jitter(NamedTuple):
    beta_eps: jax.Array
    beta_leap: jax.Array
    yes_eps: jax.Array
    yes_leap: jax.Array
    no_eps: jax.Array
    no_leap: jax.Array


class SweepStats(NamedTuple):
    acc_kappa: jax.Array
    acc_beta: jax.Array
    acc_yes: jax.Array
    acc_no: jax.Array
    acc_omega: jax.Array
    loglik: jax.Array


def wrap(a):
    """Map angles to [-pi, pi)."""
    return jnp.mod(a + jnp.pi, 2 * jnp.pi) - jnp.pi


def loglik_matrix(beta, yes, no, kappa, data: Data):
    """(I, J) matrix of log P(y_ij | params); zero where y_ij is missing."""
    dt = beta.dtype
    d_yes = wrap(yes[None, :] - beta[:, None])
    d_no = wrap(no[None, :] - beta[:, None])
    u = data.sign.astype(dt) * (d_no * d_no - d_yes * d_yes)
    ll = log_link_cdf(u, kappa[None, :].astype(dt))
    return jnp.where(data.mask, ll, 0.0)


def log_prior_omega(omega, beta, priors: Priors):
    """log p(omega | beta) up to a constant (von Mises normaliser included)."""
    n = beta.shape[0]
    log_i0 = jnp.log(i0e(omega)) + omega
    return (
        -n * log_i0
        + omega * jnp.sum(jnp.cos(beta - priors.mu))
        + (priors.omega_shape - 1.0) * jnp.log(omega)
        - priors.omega_rate * omega
    )


def _sweep(key, state: State, jitter: Jitter, data: Data, priors: Priors, cfg: SamplerConfig) -> tuple[State, SweepStats]:
    gdt = jnp.dtype(cfg.grad_dtype)
    data_g = Data(data.sign.astype(gdt), data.mask)
    k_kappa, k_kappa_u, k_xi, k_beta, k_yes, k_no, k_om, k_om_u, k_rot = jax.random.split(key, 9)
    beta, yes, no, kappa, xi, omega, ll = state.beta, state.yes, state.no, state.kappa, state.xi, state.omega, state.ll

    # --- kappa_j: log-normal random walk MH --------------------------------
    log_k = jnp.log(kappa)
    log_kp = log_k + state.kappa_sd * jax.random.normal(k_kappa, kappa.shape)
    kappa_p = jnp.exp(log_kp)
    ll_p = loglik_matrix(beta, yes, no, kappa_p, data)
    log_r = (
        ll_p.sum(0)
        - ll.sum(0)
        + priors.kappa_shape * (log_kp - log_k)  # (shape-1) log-prior term + Jacobian
        - xi * (kappa_p - kappa)
    )
    # The likelihood code is validated for kappa <= KAPPA_MAX; proposals beyond
    # it are rejected: a truncation of the prior at a point where it has mass
    # ~exp(-KAPPA_MAX * xi), i.e. nothing for any plausible xi.
    acc_kappa = (jnp.log(jax.random.uniform(k_kappa_u, kappa.shape)) < log_r) & (kappa_p <= KAPPA_MAX)
    kappa = jnp.where(acc_kappa, kappa_p, kappa)
    ll = jnp.where(acc_kappa[None, :], ll_p, ll)

    # --- xi = 1/lambda: conjugate Gamma ----------------------------------------
    J = kappa.shape[0]
    xi = jax.random.gamma(k_xi, priors.kappa_shape * J + priors.xi_shape) / (priors.xi_rate + kappa.sum())

    # --- beta_i: GHMC, all legislators at once ----------------------------------
    # Each block's log_prob also returns its slice of the log-likelihood matrix
    # (rows for beta, columns for psi/zeta) so the matrix is carried through the
    # sweep; only the proposals need fresh float64 evaluations.
    def lp_beta(x, yes, no, kappa, omega, data):
        th = circle_to_angle(x)
        return loglik_matrix(th, yes, no, kappa, data).sum(1) + omega * jnp.cos(th - priors.mu)

    def lp_beta_aux(x):
        th = circle_to_angle(x)
        m = loglik_matrix(th, yes, no, kappa, data)
        return m.sum(1) + omega * jnp.cos(th - priors.mu), m

    def grad_of(lp, *args):
        """Leapfrog gradient of lp(x, *args, data) evaluated in grad_dtype."""
        cast = [a.astype(gdt) for a in args]
        return lambda x: jax.grad(lambda x: lp(x, *cast, data_g).sum())(x.astype(gdt))

    x, info_b = ghmc_step(
        k_beta,
        angle_to_circle(beta),
        lp_beta_aux,
        jitter.beta_eps * state.beta_scale,
        jitter.beta_leap,
        cfg.beta_max_leap,
        grad_fn=grad_of(lp_beta, yes, no, kappa, omega),
        log_prob_x=(ll.sum(1) + omega * jnp.cos(beta - priors.mu), ll),
        has_aux=True,
    )
    beta = circle_to_angle(x)
    ll = info_b.aux

    # --- psi_j (Yea positions): GHMC ----------------------------------------------
    def lp_yes(x, beta, no, kappa, data):
        return loglik_matrix(beta, circle_to_angle(x), no, kappa, data).sum(0)

    def lp_yes_aux(x):
        m = loglik_matrix(beta, circle_to_angle(x), no, kappa, data)
        return m.sum(0), m.T

    x, info_y = ghmc_step(
        k_yes,
        angle_to_circle(yes),
        lp_yes_aux,
        jitter.yes_eps * state.yes_scale,
        jitter.yes_leap,
        cfg.tau_max_leap,
        grad_fn=grad_of(lp_yes, beta, no, kappa),
        log_prob_x=(ll.sum(0), ll.T),
        has_aux=True,
    )
    yes = circle_to_angle(x)
    ll = info_y.aux.T

    # --- zeta_j (Nay positions): GHMC ---------------------------------------------
    def lp_no(x, beta, yes, kappa, data):
        return loglik_matrix(beta, yes, circle_to_angle(x), kappa, data).sum(0)

    def lp_no_aux(x):
        m = loglik_matrix(beta, yes, circle_to_angle(x), kappa, data)
        return m.sum(0), m.T

    x, info_n = ghmc_step(
        k_no,
        angle_to_circle(no),
        lp_no_aux,
        jitter.no_eps * state.no_scale,
        jitter.no_leap,
        cfg.tau_max_leap,
        grad_fn=grad_of(lp_no, beta, yes, kappa),
        log_prob_x=(info_y.log_prob, ll.T),
        has_aux=True,
    )
    no = circle_to_angle(x)
    ll = info_n.aux.T

    # --- omega: log-normal random walk MH ----------------------------------------
    log_om = jnp.log(omega)
    log_om_p = log_om + cfg.omega_prop_sd * jax.random.normal(k_om)
    om_p = jnp.exp(log_om_p)
    log_r = log_prior_omega(om_p, beta, priors) - log_prior_omega(omega, beta, priors) + (log_om_p - log_om)
    acc_omega = jnp.log(jax.random.uniform(k_om_u)) < log_r
    omega = jnp.where(acc_omega, om_p, omega)

    # --- common rotation delta of all angles (group move, Liu & Wu 1999) ------
    # The likelihood depends on angle differences only, so p(delta | rest) is
    # proportional to prod_i exp(omega cos(beta_i + delta - mu)) =
    # vonMises(-phi, omega * R), with (R, phi) the resultant of beta - mu.
    if cfg.rotation_move:
        C, S = jnp.sum(jnp.cos(beta - priors.mu)), jnp.sum(jnp.sin(beta - priors.mu))
        delta = sample_von_mises(k_rot, -jnp.arctan2(S, C), omega * jnp.hypot(C, S))
        beta, yes, no = wrap(beta + delta), wrap(yes + delta), wrap(no + delta)
    new_state = state._replace(beta=beta, yes=yes, no=no, kappa=kappa, xi=xi, omega=omega, ll=ll)
    stats = SweepStats(acc_kappa, info_b.accepted, info_y.accepted, info_n.accepted, acc_omega, ll.sum())
    return new_state, stats


class WaicAccum(NamedTuple):
    """Running sums for pointwise WAIC over the observed entries."""

    lse: jax.Array  # log sum_s exp(ll_ij^s)
    s1: jax.Array  # sum_s ll_ij^s
    s2: jax.Array  # sum_s (ll_ij^s)^2
    n: jax.Array


def _accumulate(acc: WaicAccum, ll):
    return WaicAccum(jnp.logaddexp(acc.lse, ll), acc.s1 + ll, acc.s2 + ll * ll, acc.n + 1)


@dataclass
class Posterior:
    """Posterior draws (angles in (-pi, pi]) plus diagnostics."""

    beta: np.ndarray  # (S, I)
    yes: np.ndarray  # (S, J)  psi
    no: np.ndarray  # (S, J)   zeta
    kappa: np.ndarray  # (S, J)
    omega: np.ndarray  # (S,)
    xi: np.ndarray  # (S,)
    ll_legislator: np.ndarray  # (S, I) per-legislator log-likelihood (grouped WAIC)
    loglik_trace: np.ndarray  # (burnin + S*thin,) joint log-likelihood, every iteration
    accept: dict  # mean acceptance rates after burn-in
    waic_pointwise: dict
    priors: dict
    config: dict
    seconds: float
    tuning: dict = field(default_factory=dict)

    def waic_grouped(self) -> dict:
        """Grouped (by legislator) WAIC as in the paper's Table 4."""
        ll = self.ll_legislator.astype(np.float64)
        S = ll.shape[0]
        m = ll.max(0)
        lppd = np.sum(m + np.log(np.exp(ll - m).sum(0)) - np.log(S))
        p_waic = np.sum(ll.var(0, ddof=1))
        return {"lppd": lppd, "p_waic": p_waic, "elpd_waic": lppd - p_waic, "waic_deviance": -2 * (lppd - p_waic)}

    def save(self, path):
        arrays = {k: v for k, v in asdict(self).items() if isinstance(v, np.ndarray)}
        meta = {k: v for k, v in asdict(self).items() if not isinstance(v, np.ndarray)}
        np.savez_compressed(path, **arrays, meta=np.array(repr(meta)))

    @classmethod
    def load(cls, path):
        import ast

        z = np.load(path, allow_pickle=False)
        meta = ast.literal_eval(str(z["meta"]))
        return cls(**{k: z[k] for k in z.files if k != "meta"}, **meta)


def initial_state(I, J, key, priors: Priors, cfg: SamplerConfig, party=None) -> dict:
    """Initial values following the original code.

    With party labels, Democrats start in (-pi/2, 0) and Republicans in
    (0, pi/2) (others at 0), which also fixes the reflection of the circle so
    that liberals sit at negative angles.
    """
    k1, k2, k3 = jax.random.split(key, 3)
    if party is not None:
        party = np.asarray(party)
        u = np.asarray(jax.random.uniform(k1, (I,), minval=0, maxval=np.pi / 2))
        beta = np.where(party == "D", -u, np.where(party == "R", u, 0.0))
    else:
        beta = np.asarray(jax.random.uniform(k1, (I,), minval=-np.pi / 2, maxval=np.pi / 2))
    return dict(
        beta=beta,
        yes=np.asarray(jax.random.uniform(k2, (J,), minval=-np.pi, maxval=np.pi)),
        no=np.asarray(jax.random.uniform(k3, (J,), minval=-np.pi, maxval=np.pi)),
        kappa=np.full(J, 0.1),
        xi=priors.xi_shape / priors.xi_rate,
        omega=priors.omega_shape / priors.omega_rate,
    )


def _adapt_scale(scale, acc_rate, low=0.6, high=0.9):
    scale = jnp.where(acc_rate < low, scale * 0.9, scale)
    scale = jnp.where(acc_rate > high, scale * 1.1, scale)
    return jnp.clip(scale, 0.05, 20.0)


def sample(
    votes,
    priors: Priors = Priors(),
    config: SamplerConfig = SamplerConfig(),
    party=None,
    init: dict | None = None,
) -> Posterior:
    """Run the MCMC sampler on an (I, J) vote matrix (1 Yea / 0 Nay / NaN missing).

    Args:
      party: optional length-I array of "D"/"R"/... labels for initialisation.
      init: optional dict of initial values (beta, yes, no, kappa, xi, omega)
        overriding the defaults, e.g. to continue a previous run.
    """
    votes = np.asarray(votes, dtype=float)
    I, J = votes.shape
    mask = ~np.isnan(votes)
    data = Data(jnp.asarray(np.where(mask, 2 * np.nan_to_num(votes) - 1, 0.0)), jnp.asarray(mask))
    cfg = config
    key = jax.random.key(cfg.seed)
    key, k_init = jax.random.split(key)
    values = initial_state(I, J, k_init, priors, cfg, party)
    if init:
        values.update({k: v for k, v in init.items() if k in values})
    values = {k: jnp.asarray(v, jnp.float64) for k, v in values.items()}
    state = State(
        **values,
        ll=loglik_matrix(values["beta"], values["yes"], values["no"], values["kappa"], data),
        kappa_sd=jnp.full(J, cfg.kappa_prop_sd),
        beta_scale=jnp.ones(I),
        yes_scale=jnp.ones(J),
        no_scale=jnp.ones(J),
    )

    block = cfg.jitter_every
    n_iter = cfg.burnin + cfg.n_samples * cfg.thin
    n_blocks = -(-n_iter // block)

    @jax.jit
    def draw_jitter(key):
        ks = jax.random.split(key, 6)
        u = lambda k, n, r: jax.random.uniform(k, (n,), minval=r[0], maxval=r[1])
        return Jitter(
            u(ks[0], I, cfg.beta_step),
            jax.random.randint(ks[1], (I,), 1, cfg.beta_max_leap + 1),
            u(ks[2], J, cfg.tau_step),
            jax.random.randint(ks[3], (J,), 1, cfg.tau_max_leap + 1),
            u(ks[4], J, cfg.tau_step),
            jax.random.randint(ks[5], (J,), 1, cfg.tau_max_leap + 1),
        )

    @partial(jax.jit, static_argnames=("n", "record"), donate_argnums=(1, 4))
    def run_block(key, state, jitter, data, waic, n, record):
        def body(carry, key):
            state, waic = carry
            state, stats = _sweep(key, state, jitter, data, priors, cfg)
            if record:
                waic = _accumulate(waic, state.ll)
                draws = (state.beta, state.yes, state.no, state.kappa, state.omega, state.xi, state.ll.sum(1))
            else:
                draws = ()
            return (state, waic), (stats, draws)

        (state, waic), (stats, draws) = jax.lax.scan(body, (state, waic), jax.random.split(key, n))
        mean = lambda a: a.astype(jnp.float64).mean(0)
        summary = SweepStats(*(mean(s) for s in stats[:-1]), stats.loglik)
        return state, waic, summary, draws

    waic = WaicAccum(jnp.full((I, J), -jnp.inf), jnp.zeros((I, J)), jnp.zeros((I, J)), jnp.zeros((), jnp.int32))
    S = cfg.n_samples
    out = {
        "beta": np.empty((S, I)),
        "yes": np.empty((S, J), np.float32),
        "no": np.empty((S, J), np.float32),
        "kappa": np.empty((S, J), np.float32),
        "omega": np.empty(S),
        "xi": np.empty(S),
        "ll_legislator": np.empty((S, I), np.float32),
    }
    trace = np.empty(n_iter)
    acc_sum = {k: 0.0 for k in ("kappa", "beta", "yes", "no", "omega")}
    n_post_iter = 0
    filled = 0
    it = 0
    t0 = time.perf_counter()
    last_print, last_time = 0, t0
    for b in range(n_blocks):
        n = min(block, n_iter - it)
        key, k_jit, k_run = jax.random.split(key, 3)
        jitter = draw_jitter(k_jit)
        record = it >= cfg.burnin
        if not record and it + n > cfg.burnin:  # keep blocks aligned with the burn-in boundary
            n = cfg.burnin - it
        state, waic, summary, draws = run_block(k_run, state, jitter, data, waic, n, record)
        trace[it : it + n] = np.asarray(summary.loglik)
        if record:
            idx = np.arange(it - cfg.burnin, it - cfg.burnin + n)
            keep = (idx % cfg.thin) == cfg.thin - 1
            m = int(keep.sum())
            for name, arr in zip(out, draws):
                out[name][filled : filled + m] = np.asarray(arr)[keep]
            filled += m
            n_post_iter += n
            for name, v in zip(acc_sum, summary[:5]):
                acc_sum[name] += float(np.mean(np.asarray(v))) * n
        elif it > 0 or b > 0:
            # Burn-in tuning (as in the R code: kappa proposal sd towards 30-60%).
            ks = summary.acc_kappa
            sd = jnp.where(ks <= 0.3, state.kappa_sd * 0.95, state.kappa_sd)
            sd = jnp.where(ks >= 0.6, sd * 1.05, sd)
            state = state._replace(kappa_sd=sd)
            if cfg.adapt_steps:
                state = state._replace(
                    beta_scale=_adapt_scale(state.beta_scale, summary.acc_beta),
                    yes_scale=_adapt_scale(state.yes_scale, summary.acc_yes),
                    no_scale=_adapt_scale(state.no_scale, summary.acc_no),
                )
        it += n
        if cfg.print_every and (it - last_print >= cfg.print_every or it == n_iter):
            now = time.perf_counter()
            rate = (it - last_print) / (now - last_time)
            last_print, last_time = it, now
            print(
                f"iter {it:>7}/{n_iter}  loglik {trace[it - 1]:>12.1f}  "
                f"acc beta {float(summary.acc_beta.mean()):.2f} psi {float(summary.acc_yes.mean()):.2f} "
                f"zeta {float(summary.acc_no.mean()):.2f} kappa {float(summary.acc_kappa.mean()):.2f} "
                f"omega {float(summary.acc_omega):.2f}  omega={float(state.omega):.2f}  "
                f"{rate:.1f} it/s",
                flush=True,
            )

    seconds = time.perf_counter() - t0
    lse, s1, s2, n_rec = (np.asarray(a) for a in waic)
    n_rec = int(n_rec)
    obs = mask
    lppd = np.sum(lse[obs] - np.log(n_rec))
    p_waic = np.sum((s2[obs] - s1[obs] ** 2 / n_rec) / (n_rec - 1)) if n_rec > 1 else np.nan
    return Posterior(
        **out,
        loglik_trace=trace,
        accept={k: v / max(n_post_iter, 1) for k, v in acc_sum.items()},
        waic_pointwise={"lppd": float(lppd), "p_waic": float(p_waic), "elpd_waic": float(lppd - p_waic)},
        priors=asdict(priors),
        config=asdict(cfg),
        seconds=seconds,
        tuning={
            "kappa_sd": np.asarray(state.kappa_sd).tolist(),
            "beta_scale_mean": float(state.beta_scale.mean()),
            "yes_scale_mean": float(state.yes_scale.mean()),
            "no_scale_mean": float(state.no_scale.mean()),
        },
    )
