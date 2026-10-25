import jax
import jax.numpy as jnp
import mpmath as mp
import numpy as np
import pytest
from scipy import special as sp

from circfactor.special import PI2, log_link_cdf

KAPPAS = [0.001, 0.01, 0.1, 0.5, 1.0, 3.0, 10.0, 50.0, 200.0, 1000.0, 5000.0, 10000.0]


def _reference(kappa, u):
    """log G_kappa(u) and its derivative at the exact float64 values of u.

    s is formed with the same float64 pi^2 the model uses.  scipy is used
    wherever it does not underflow (mpmath is slow near s = 1/2 for large
    kappa); mpmath covers the tails.
    """
    mp.mp.dps = 40
    pi2 = mp.mpf(PI2)
    out, grad = np.empty_like(u), np.empty_like(u)
    for i, ui in enumerate(u):
        s_lo = (PI2 + ui) / (2 * PI2)
        s_hi = (PI2 - ui) / (2 * PI2)  # 1 - s without cancellation
        small = min(s_lo, s_hi)
        p_small = sp.betainc(kappa, kappa, small)
        s = (pi2 + mp.mpf(float(ui))) / (2 * pi2)
        if kappa > 5000 and p_small < 1e-250:
            out[i] = grad[i] = np.nan  # mpmath far too slow here; skipped below
            continue
        if kappa <= 50 or p_small < 1e-250:
            lower = mp.betainc(kappa, kappa, 0, s, regularized=True)
            logcdf = mp.log(lower) if lower < 0.5 else mp.log1p(-mp.betainc(kappa, kappa, 0, 1 - s, regularized=True))
        else:
            logcdf = np.log(p_small) if s_lo <= 0.5 else np.log1p(-p_small)
        out[i] = float(logcdf)
        grad[i] = float(mp.exp((kappa - 1) * (mp.log(s) + mp.log(1 - s)) - mp.log(mp.beta(kappa, kappa)) - logcdf) / (2 * pi2))
    return out, grad


def _inputs(kappa, rng, edge):
    # Uniform over the range, plus points within a few link sds of the centre
    # (for large kappa almost everything else is in the far tail).
    sd = PI2 / np.sqrt(2 * kappa + 1)
    parts = [rng.uniform(-PI2, PI2, 300), np.clip(rng.normal(0, 3 * sd, 200), -PI2 * (1 - 1e-9), PI2 * (1 - 1e-9)), [0.0, 1e-4, -1e-4]]
    # Deep tails on both sides.
    parts += [-PI2 * (1 - 10 ** rng.uniform(np.log10(edge), -1, 60)), PI2 * (1 - 10 ** rng.uniform(np.log10(edge), -1, 60))]
    return np.concatenate(parts)


def _rel(a, b):
    return np.abs(a - b) / np.maximum(1.0, np.abs(b))


@pytest.mark.parametrize("kappa", KAPPAS)
def test_log_link_cdf_float64(kappa):
    u = _inputs(kappa, np.random.default_rng(int(kappa * 100)), edge=1e-12)
    ref, gref = _reference(kappa, u)
    val = np.asarray(jax.jit(log_link_cdf)(jnp.asarray(u), jnp.asarray(kappa)))
    grad = np.asarray(jax.jit(jax.vmap(jax.grad(log_link_cdf), (0, None)))(jnp.asarray(u), jnp.asarray(kappa)))
    assert np.all(np.isfinite(val)) and np.all(np.isfinite(grad))
    ok = ~np.isnan(ref)
    assert ok.sum() > 100
    # For kappa > 1000, terms of size ~kappa cancel down to O(1), costing a
    # few digits; ~1e-10 is still far below anything that matters for MCMC.
    # Gradient errors only affect GHMC efficiency, never exactness.
    tol = 1e-10 if kappa <= 1000 else 5e-10
    assert _rel(val[ok], ref[ok]).max() < tol
    assert _rel(grad[ok], gref[ok]).max() < (1e-10 if kappa <= 1000 else 1e-8)


def test_symmetry():
    """log P(Yea) and log P(Nay) must be log G(u) and log(1 - G(u))."""
    u = jnp.linspace(-PI2 + 1e-3, PI2 - 1e-3, 101)
    for k in [0.3, 4.0, 300.0]:
        total = jnp.exp(log_link_cdf(u, k)) + jnp.exp(log_link_cdf(-u, k))
        np.testing.assert_allclose(np.asarray(total), 1.0, atol=1e-12)


def test_trip_count_is_converged():
    """The fixed trip count must already be converged for kappa <= KAPPA_MAX:
    16 more iterations must not change the result anywhere in u."""
    from circfactor import special

    u = np.concatenate([np.linspace(-PI2, PI2, 2001), np.logspace(-12, 0, 200), -PI2 * (1 - np.logspace(-12, -1, 100))])
    kappa = np.logspace(-3, np.log10(special.KAPPA_MAX), 200)
    U, K = (jnp.asarray(a) for a in np.broadcast_arrays(u[None, :], kappa[:, None]))
    f = jax.jit(special._log_link_cdf_value, static_argnums=2)
    got, more = np.asarray(f(U, K, special.CF_ITERS)), np.asarray(f(U, K, special.CF_ITERS + 16))
    assert np.all(np.isfinite(got))
    assert np.max(np.abs(got - more) / np.maximum(1.0, np.abs(more))) < 1e-14


def test_betaln_helpers():
    from circfactor.special import betaln_half, betaln_sym

    mp.mp.dps = 40
    a = np.concatenate([np.logspace(-3, 6, 200), [9.999999, 10.0, 10.000001]])
    ref_sym = np.array([float(mp.log(mp.beta(x, x))) for x in a])
    ref_half = np.array([float(mp.log(mp.beta(0.5, x))) for x in a])
    np.testing.assert_allclose(np.asarray(betaln_sym(jnp.asarray(a))), ref_sym, rtol=1e-14, atol=1e-14)
    np.testing.assert_allclose(np.asarray(betaln_half(jnp.asarray(a))), ref_half, rtol=1e-14, atol=1e-14)
