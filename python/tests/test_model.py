"""Model-level correctness tests: likelihood and the full Gibbs sweep."""

import jax
import jax.numpy as jnp
import numpy as np
import pytest
from scipy import stats

from circfactor.model import Data, Jitter, Priors, SamplerConfig, State, _sweep, loglik_matrix
from circfactor.special import log_link_cdf


def test_loglik_matrix_matches_direct_formula():
    """Paper eq. (4)-(5): P(yea) = BetaCDF((e + pi^2) / (2 pi^2); kappa, kappa)."""
    rng = np.random.default_rng(0)
    I, J = 7, 9
    beta, yes, no = (rng.uniform(-np.pi, np.pi, n) for n in (I, J, J))
    kappa = rng.exponential(5.0, J) + 0.05
    y = rng.integers(0, 2, (I, J)).astype(float)
    y[rng.random((I, J)) < 0.2] = np.nan
    mask = ~np.isnan(y)
    data = Data(jnp.asarray(np.where(mask, 2 * np.nan_to_num(y) - 1, 0.0)), jnp.asarray(mask))

    rho = lambda a, b: np.arccos(np.cos(a - b))  # geodesic distance, as in the paper
    e = rho(no[None, :], beta[:, None]) ** 2 - rho(yes[None, :], beta[:, None]) ** 2
    p_yea = stats.beta.cdf((e + np.pi**2) / (2 * np.pi**2), kappa, kappa)
    expected = np.where(mask, np.where(y == 1, np.log(p_yea), np.log1p(-p_yea)), 0.0)

    got = np.asarray(loglik_matrix(*(jnp.asarray(a) for a in (beta, yes, no, kappa)), data))
    np.testing.assert_allclose(got, expected, rtol=1e-10, atol=1e-12)


def _draw_prior(key, I, J, priors: Priors):
    """Exact draws of all parameters from the prior."""
    k = jax.random.split(key, 6)
    omega = jax.random.gamma(k[0], priors.omega_shape) / priors.omega_rate
    xi = jax.random.gamma(k[1], priors.xi_shape) / priors.xi_rate
    kappa = jax.random.gamma(k[2], priors.kappa_shape, (J,)) / xi
    # von Mises via the circular normal draws in jax.random? Not available:
    # use rejection-free inverse CDF on a fine grid instead.
    grid = jnp.linspace(-jnp.pi, jnp.pi, 4097)
    logp = omega * jnp.cos(grid - priors.mu)
    cdf = jnp.cumsum(jnp.exp(logp - logp.max()))
    cdf = cdf / cdf[-1]
    beta = jnp.interp(jax.random.uniform(k[3], (I,)), cdf, grid)
    yes = jax.random.uniform(k[4], (J,), minval=-jnp.pi, maxval=jnp.pi)
    no = jax.random.uniform(k[5], (J,), minval=-jnp.pi, maxval=jnp.pi)
    return beta, yes, no, kappa, xi, omega


def _simulate_votes(key, beta, yes, no, kappa, mask):
    """y_ij ~ Bernoulli(G_kappa_j(e_ij)), returned in the sampler's Data format."""
    d_yes = jnp.mod(yes[None, :] - beta[:, None] + jnp.pi, 2 * jnp.pi) - jnp.pi
    d_no = jnp.mod(no[None, :] - beta[:, None] + jnp.pi, 2 * jnp.pi) - jnp.pi
    p = jnp.exp(log_link_cdf(d_no**2 - d_yes**2, kappa[None, :]))
    y = jax.random.uniform(key, p.shape) < p
    return Data(jnp.where(mask, jnp.where(y, 1.0, -1.0), 0.0), mask)


@pytest.mark.parametrize("grad_dtype,rotation", [("float32", False), ("float64", False), ("float32", True)])
def test_geweke_joint_distribution(grad_dtype, rotation):
    """Geweke (2004) successive-conditional test of the whole Gibbs sweep.

    Each of N independent chains starts at an exact prior draw with data
    simulated from it, then alternates (i) one MCMC sweep given y and (ii) a
    fresh y | parameters.  The joint law of (parameters, y) is invariant under
    both steps iff the sweep leaves the posterior invariant, so after T
    rounds the parameters must still be exact prior draws.  Any error in an
    acceptance ratio, Jacobian, conjugate update or bookkeeping of the
    carried log-likelihood matrix shows up as drift away from the prior.
    """
    I, J, N, T = 6, 8, 2000, 40
    # Priors chosen so kappa stays moderate (the paper's xi prior makes kappa
    # extremely heavy tailed a priori, which no small data set can tame).
    priors = Priors(mu=0.4, omega_shape=2.0, omega_rate=0.5, kappa_shape=2.0, xi_shape=6.0, xi_rate=20.0)
    cfg = SamplerConfig(omega_prop_sd=0.5, grad_dtype=grad_dtype, rotation_move=rotation)
    rng = np.random.default_rng(1)
    mask = jnp.asarray(rng.random((I, J)) > 0.1)

    def init(key):
        k1, k2 = jax.random.split(key)
        beta, yes, no, kappa, xi, omega = _draw_prior(k1, I, J, priors)
        data = _simulate_votes(k2, beta, yes, no, kappa, mask)
        state = State(
            beta, yes, no, kappa, xi, omega,
            ll=loglik_matrix(beta, yes, no, kappa, data),
            kappa_sd=jnp.full(J, 0.7), beta_scale=jnp.ones(I), yes_scale=jnp.ones(J), no_scale=jnp.ones(J),
        )
        return state, data

    def round_(carry, key):
        state, data = carry
        k_jit, k_sweep, k_data = jax.random.split(key, 3)
        kj = jax.random.split(k_jit, 6)
        # Step sizes large enough that the MH corrections matter.
        jitter = Jitter(
            jax.random.uniform(kj[0], (I,), minval=0.05, maxval=0.4), jax.random.randint(kj[1], (I,), 1, 11),
            jax.random.uniform(kj[2], (J,), minval=0.05, maxval=0.4), jax.random.randint(kj[3], (J,), 1, 11),
            jax.random.uniform(kj[4], (J,), minval=0.05, maxval=0.4), jax.random.randint(kj[5], (J,), 1, 11),
        )
        # Several sweeps per data refresh, so later sweeps rely on the
        # log-likelihood matrix carried in the state.
        for k in jax.random.split(k_sweep, 3):
            state, stats = _sweep(k, state, jitter, data, priors, cfg)
        data = _simulate_votes(k_data, state.beta, state.yes, state.no, state.kappa, mask)
        state = state._replace(ll=loglik_matrix(state.beta, state.yes, state.no, state.kappa, data))
        return (state, data), (stats.acc_beta.mean(), stats.acc_yes.mean(), stats.acc_kappa.mean())

    @jax.jit
    def run(key):
        k0, k1 = jax.random.split(key)
        carry = init(k0)
        (state, _), acc = jax.lax.scan(round_, carry, jax.random.split(k1, T))
        return state, acc

    states, acc = jax.vmap(run)(jax.random.split(jax.random.key(42), N))
    acc = [float(a.mean()) for a in acc]
    assert all(0.1 < a < 0.99 for a in acc), acc  # the MH steps really are rejecting

    # Fresh prior draws for comparison (two-sample tests).
    ref = jax.vmap(lambda k: _draw_prior(k, I, J, priors))(jax.random.split(jax.random.key(7), 20000))
    names = ["beta", "yes", "no", "kappa", "xi", "omega"]
    got = [states.beta, states.yes, states.no, states.kappa, states.xi, states.omega]
    pvals = {}
    for name, g, r in zip(names, got, ref):
        g, r = np.asarray(g), np.asarray(r)
        if g.ndim == 2:  # pool coordinates of the first chain row-wise: use coordinate 0 and the mean
            for c in (0, g.shape[1] - 1):
                pvals[f"{name}[{c}]"] = stats.ks_2samp(g[:, c], r[:, c]).pvalue
        else:
            pvals[name] = stats.ks_2samp(g, r).pvalue
    # Exact one-sample checks where the marginal is known in closed form.
    pvals["xi_exact"] = stats.kstest(np.asarray(states.xi), stats.gamma(priors.xi_shape, scale=1 / priors.xi_rate).cdf).pvalue
    pvals["omega_exact"] = stats.kstest(np.asarray(states.omega), stats.gamma(priors.omega_shape, scale=1 / priors.omega_rate).cdf).pvalue
    pvals["yes_uniform"] = stats.kstest(np.asarray(states.yes).ravel()[::3], stats.uniform(-np.pi, 2 * np.pi).cdf).pvalue
    pvals["no_uniform"] = stats.kstest(np.asarray(states.no).ravel()[::3], stats.uniform(-np.pi, 2 * np.pi).cdf).pvalue
    bad = {k: v for k, v in pvals.items() if v < 1e-3}
    assert not bad, (bad, pvals)


def test_carried_loglik_matrix_is_consistent():
    """After every sweep the carried matrix must equal a fresh evaluation."""
    I, J = 9, 11
    priors = Priors()
    cfg = SamplerConfig()
    rng = np.random.default_rng(3)
    mask = jnp.asarray(rng.random((I, J)) > 0.15)
    beta, yes, no, kappa, xi, omega = _draw_prior(jax.random.key(0), I, J, Priors(xi_shape=6.0, xi_rate=20.0))
    data = _simulate_votes(jax.random.key(1), beta, yes, no, kappa, mask)
    state = State(beta, yes, no, kappa, xi, omega, loglik_matrix(beta, yes, no, kappa, data),
                  jnp.full(J, 0.7), jnp.ones(I), jnp.ones(J), jnp.ones(J))
    jitter = Jitter(jnp.full(I, 0.3), jnp.full(I, 5), jnp.full(J, 0.3), jnp.full(J, 5), jnp.full(J, 0.3), jnp.full(J, 5))
    sweep = jax.jit(lambda k, s: _sweep(k, s, jitter, data, priors, cfg))
    n_rej = 0
    for k in jax.random.split(jax.random.key(2), 30):
        state, stats = sweep(k, state)
        n_rej += int((~stats.acc_yes).sum() + (~stats.acc_no).sum() + (~stats.acc_beta).sum())
        fresh = loglik_matrix(state.beta, state.yes, state.no, state.kappa, data)
        np.testing.assert_allclose(np.asarray(state.ll), np.asarray(fresh), rtol=1e-12, atol=1e-12)
    assert n_rej > 0  # rejections happened, so the accept masks were exercised
