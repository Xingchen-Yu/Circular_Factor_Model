"""Validate the geodesic HMC kernel against von Mises targets.

Each test runs thousands of independent GHMC chains in one vectorised call and
compares the chains' final states (i.i.d. given enough iterations) with the
exact distribution, using Kolmogorov-Smirnov tests and moment checks.
"""

import jax
import jax.numpy as jnp
import numpy as np
import pytest
from scipy import stats
from scipy.special import i0e, i1e

from circfactor.ghmc import angle_to_circle, circle_to_angle, ghmc_step, sample_von_mises

N_CHAINS = 5000
P_MIN = 1e-3  # per-test KS threshold; ~50 KS tests in this file


def wrap(a):
    return (np.asarray(a) + np.pi) % (2 * np.pi) - np.pi


def run_chains(key, x0, log_prob, n_iter, eps_range, max_leap, grad_fn=None):
    """Run n_iter GHMC steps with per-chain jittered (eps, L) every iteration."""
    B = x0.shape[0]

    @jax.jit
    def body(carry, key):
        x, n_acc = carry
        k_eps, k_leap, k_step = jax.random.split(key, 3)
        eps = jax.random.uniform(k_eps, (B,), minval=eps_range[0], maxval=eps_range[1])
        leap = jax.random.randint(k_leap, (B,), 1, max_leap + 1)
        x, info = ghmc_step(k_step, x, log_prob, eps, leap, max_leap, grad_fn=grad_fn)
        return (x, n_acc + info.accepted), info.accept_prob.mean()

    (x, n_acc), acc = jax.lax.scan(body, (x0, jnp.zeros(B)), jax.random.split(key, n_iter))
    return x, np.asarray(acc)


def vm_log_prob(mu, omega):
    direction = jnp.array([np.cos(mu), np.sin(mu)])
    return lambda x: omega * (x @ direction)


def uniform_start(key, B):
    return angle_to_circle(jax.random.uniform(key, (B,), minval=-np.pi, maxval=np.pi))


def check_von_mises(theta, mu, omega):
    """KS test of angles against vonMises(mu, omega) plus a resultant-length check."""
    delta = wrap(theta - mu)
    if omega == 0:
        cdf = lambda d: (d + np.pi) / (2 * np.pi)
    else:
        cdf = lambda d: stats.vonmises.cdf(d, omega)
    p = stats.kstest(delta, cdf).pvalue
    assert p > P_MIN, f"KS p={p:.2e} for mu={mu}, omega={omega}"
    # E[cos(theta - mu)] = I1(omega)/I0(omega); compare within 5 standard errors.
    c = np.cos(delta)
    expected = 0.0 if omega == 0 else i1e(omega) / i0e(omega)
    se = c.std() / np.sqrt(len(c))
    assert abs(c.mean() - expected) < 5 * se + 1e-12
    # Circular mean direction should match mu (skip when nearly uniform).
    if omega >= 1:
        s = np.sin(delta)
        assert abs(s.mean()) < 5 * s.std() / np.sqrt(len(s))
    return p


MUS = [0.0, 2.0, np.pi - 0.05, -np.pi / 2]
OMEGAS = [0.0, 0.1, 1.0, 4.0, 16.0, 64.0, 256.0]


@pytest.mark.parametrize("omega", OMEGAS)
@pytest.mark.parametrize("mu", MUS)
def test_ghmc_recovers_von_mises(mu, omega):
    """Tuned step sizes: eps scaled with the target's width 1/sqrt(omega)."""
    key = jax.random.key(hash((round(mu, 3), omega)) % 2**31)
    k0, k1 = jax.random.split(key)
    scale = 1.0 / np.sqrt(1.0 + omega)
    x, acc = run_chains(k1, uniform_start(k0, N_CHAINS), vm_log_prob(mu, omega), 300, (0.2 * scale, 0.6 * scale), 10)
    check_von_mises(np.asarray(circle_to_angle(x)), mu, omega)
    assert acc[-50:].mean() > 0.6


@pytest.mark.parametrize("omega", [0.0, 1.0, 10.0, 50.0])
def test_ghmc_paper_jitter_ranges(omega):
    """The paper's tuning for tau: eps ~ U(0.01, 0.105), L ~ U{1..10}."""
    key = jax.random.key(int(omega) + 7)
    k0, k1 = jax.random.split(key)
    x, acc = run_chains(k1, uniform_start(k0, N_CHAINS), vm_log_prob(1.0, omega), 1500, (0.01, 0.105), 10)
    check_von_mises(np.asarray(circle_to_angle(x)), 1.0, omega)
    assert acc[-100:].mean() > 0.9


@pytest.mark.parametrize("omega", [0.5, 5.0, 40.0])
def test_ghmc_exact_with_wrong_gradient(omega):
    """A wrong (but deterministic) gradient must not bias the chain.

    This is what licenses computing leapfrog gradients in float32 while the
    accept/reject step uses the float64 density.
    """
    mu = -2.5
    log_prob = vm_log_prob(mu, omega)
    direction = jnp.array([np.cos(mu), np.sin(mu)])
    # 0.6 x true gradient plus a smooth, non-gradient-consistent distortion.
    bad_grad = lambda x: 0.6 * omega * direction + 0.8 * jnp.sin(3 * x[:, ::-1])

    key = jax.random.key(int(omega * 10))
    k0, k1, k2 = jax.random.split(key, 3)
    scale = 1.0 / np.sqrt(1.0 + omega)
    # From exact draws: stationarity.
    exact = mu + stats.vonmises.rvs(omega, size=N_CHAINS, random_state=int(omega * 10))
    x, _ = run_chains(k1, angle_to_circle(jnp.asarray(exact)), log_prob, 30, (0.2 * scale, 0.6 * scale), 10, grad_fn=bad_grad)
    check_von_mises(np.asarray(circle_to_angle(x)), mu, omega)
    # From uniform: convergence.
    x, acc = run_chains(k2, uniform_start(k0, N_CHAINS), log_prob, 400, (0.2 * scale, 0.6 * scale), 10, grad_fn=bad_grad)
    check_von_mises(np.asarray(circle_to_angle(x)), mu, omega)


@pytest.mark.parametrize("kappa,n", [(0.5, 20), (2.0, 50), (8.0, 10), (30.0, 200)])
def test_posterior_of_mean_direction(kappa, n):
    """Recover the location of simulated von Mises data.

    With x_1..x_n ~ vonMises(mu, kappa) (kappa known) and a uniform prior on
    mu, the posterior is vonMises(atan2(S, C), kappa * R), where (C, S) =
    sum(cos x, sin x) and R = |(C, S)|.  For each of 40 simulated datasets we
    sample the posterior with GHMC, test it against the closed form, and check
    the coverage of 90% credible intervals for the true mu.
    """
    rng = np.random.default_rng(int(kappa * 100) + n)
    n_rep, n_chain = 40, 2000
    mu_true = rng.uniform(-np.pi, np.pi, n_rep)
    data = mu_true[:, None] + stats.vonmises.rvs(kappa, size=(n_rep, n), random_state=rng)
    resultant = np.stack([np.cos(data).sum(1), np.sin(data).sum(1)], axis=1)  # (rep, 2)

    # One batch: n_rep * n_chain chains, each targeting its dataset's posterior.
    eta = jnp.asarray(np.repeat(kappa * resultant, n_chain, axis=0))  # (B, 2)
    log_prob = lambda x: jnp.sum(eta * x, axis=-1)
    post_conc = kappa * np.linalg.norm(resultant, axis=1)
    step = np.repeat(0.5 / np.sqrt(1 + post_conc), n_chain)
    key = jax.random.key(n)
    k0, k1 = jax.random.split(key)
    B = n_rep * n_chain

    @jax.jit
    def body(x, key):
        k_eps, k_leap, k_step = jax.random.split(key, 3)
        eps = step * jax.random.uniform(k_eps, (B,), minval=0.4, maxval=1.2)
        leap = jax.random.randint(k_leap, (B,), 1, 11)
        x, _ = ghmc_step(k_step, x, log_prob, eps, leap, 10)
        return x, None

    x, _ = jax.lax.scan(body, uniform_start(k0, B), jax.random.split(k1, 300))
    theta = np.asarray(circle_to_angle(x)).reshape(n_rep, n_chain)

    post_mean = np.arctan2(resultant[:, 1], resultant[:, 0])
    pvals = [stats.kstest(wrap(theta[r] - post_mean[r]), lambda d, c=post_conc[r]: stats.vonmises.cdf(d, c)).pvalue for r in range(n_rep)]
    # Combine the 40 tests: p-values must look uniform.
    assert stats.kstest(pvals, "uniform").pvalue > P_MIN
    assert min(pvals) > P_MIN / n_rep

    # Coverage of central 90% intervals around the posterior mean direction.
    dev = np.abs(wrap(theta - post_mean[:, None]))
    half_width = np.quantile(dev, 0.9, axis=1)
    covered = np.abs(wrap(mu_true - post_mean)) <= half_width
    assert covered.mean() >= 0.75  # binomial(40, 0.9): P(<30) ~ 0.002


@pytest.mark.parametrize("kappa", [0.0, 2.0, 25.0])
def test_ghmc_von_mises_fisher_s2(kappa):
    """The kernel is dimension-generic: von Mises-Fisher on S^2.

    For vMF(mu, kappa) on S^2, t = mu.x has density proportional to
    exp(kappa t) on [-1, 1].
    """
    mu = jnp.array([0.0, 0.6, 0.8])
    log_prob = lambda x: kappa * (x @ mu)
    key = jax.random.key(int(kappa) + 3)
    k0, k1 = jax.random.split(key)
    x0 = jax.random.normal(k0, (N_CHAINS, 3))
    x0 = x0 / jnp.linalg.norm(x0, axis=1, keepdims=True)
    scale = 1.0 / np.sqrt(1.0 + kappa)
    x, _ = run_chains(k1, x0, log_prob, 300, (0.2 * scale, 0.6 * scale), 10)
    t = np.asarray(x @ mu)
    if kappa == 0:
        cdf = lambda t: (t + 1) / 2
    else:
        cdf = lambda t: np.expm1(kappa * (t + 1)) / np.expm1(2 * kappa)
    assert stats.kstest(t, cdf).pvalue > P_MIN
    assert np.allclose(np.linalg.norm(np.asarray(x), axis=1), 1.0, atol=1e-10)


@pytest.mark.parametrize("omega", [2.0, 20.0, 200.0])
def test_ghmc_large_steps(omega):
    """Coarse steps make the integrator error large (acceptance well below 1),
    so this checks that the Metropolis correction itself is right."""
    mu = 0.7
    key = jax.random.key(int(omega) + 11)
    k0, k1 = jax.random.split(key)
    scale = 1.0 / np.sqrt(omega)
    x, acc = run_chains(k1, uniform_start(k0, N_CHAINS), vm_log_prob(mu, omega), 400, (1.2 * scale, 2.2 * scale), 3)
    assert 0.2 < acc[-100:].mean() < 0.85
    check_von_mises(np.asarray(circle_to_angle(x)), mu, omega)


@pytest.mark.parametrize("kappa", [0.0, 0.3, 2.0, 50.0, 5000.0])
def test_exact_von_mises_sampler(kappa):
    mu = 2.9
    draws = jax.jit(jax.vmap(lambda k: sample_von_mises(k, mu, kappa)))(jax.random.split(jax.random.key(int(kappa) + 5), 20000))
    check_von_mises(np.asarray(draws), mu, kappa)
