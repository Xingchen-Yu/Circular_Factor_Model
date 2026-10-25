"""Numerically stable log-CDF of the circular model's link function.

The link G_kappa is the CDF of a symmetric Beta(kappa, kappa) variable that has
been shifted and scaled to [-pi^2, pi^2] (paper eq. 5): with
s = (pi^2 + u) / (2 pi^2), G_kappa(u) = I_s(kappa, kappa).

Symmetry gives 1 - G_kappa(u) = G_kappa(-u), so the Bernoulli log-likelihood
of a vote is ``log_link_cdf((2y - 1) * e_ij, kappa)``, and we only need
I_lo(kappa, kappa) at the lower-half point lo = (pi^2 - |u|) / (2 pi^2) <= 1/2
(computed from pi^2 - |u| rather than 1 - s to keep the tails accurate).

I_lo is evaluated with the Numerical Recipes continued fraction in one of two
forms, picked per element:

  A. directly:  I_lo(a, a) = lo^a (1-lo)^a / (a B(a, a)) * CF(a, a, lo);
     converges quickly away from the centre (lo << 1/2);
  B. near the centre, through I_x(a, a) = (1 - I_w(1/2, a)) / 2 with
     w = (1 - 2 lo)^2 = (u / pi^2)^2, and
     I_w(1/2, a) = w^(1/2) (1-w)^a / (B(1/2, a) / 2) * CF(1/2, a, w).

Choosing B iff w < min(5 / (a + 2.5), 1/2) bounds the iterations needed for
float64 precision by about 18 for a <= 300, 25 for a <= 1000 and 35 for
a <= 1e4 = KAPPA_MAX.  The loop therefore runs a fixed CF_ITERS = 40 times
(tests/test_special.py checks this is converged).  A static trip count matters
on GPU: a convergence-checked while_loop was ~2-3x slower.
The switch also keeps B away from the tail (B is only used where
I_lo >= 7.8e-4), since (1 - I_w) / 2 cancels badly when I_w is close to 1.

Gradients in ``u`` are supplied analytically through ``custom_jvp``
(d/du log G = g / G), so autodiff never differentiates through the
continued-fraction loop.  ``kappa`` gets no gradient: the sampler updates it by
Metropolis-Hastings.
"""

import jax
import jax.numpy as jnp
import numpy as np
from jax import lax
from jax.scipy.special import gammaln

PI2 = float(np.pi**2)
TWO_PI2 = 2.0 * PI2
LOG_TWO_PI2 = float(np.log(TWO_PI2))
LOG_HALF = float(np.log(0.5))
LOG_FOUR = float(np.log(4.0))

KAPPA_MAX = 1e4  # largest kappa the fixed trip count is validated for
CF_ITERS = 40
SWITCH_C = 5.0
SWITCH_CAP = 0.5


LOG_SQRT_2PI = float(0.5 * np.log(2 * np.pi))
LOG_SQRT_PI = float(0.5 * np.log(np.pi))
LOG_TWO = float(np.log(2.0))
_STIRLING_SWITCH = 10.0


def _stirling_correction(x):
    """lgamma(x) - [(x - 1/2) log x - x + log sqrt(2 pi)], for x >= 10."""
    r = 1.0 / x
    r2 = r * r
    return r * (1 / 12 + r2 * (-1 / 360 + r2 * (1 / 1260 + r2 * (-1 / 1680 + r2 * (1 / 1188 + r2 * (-691 / 360360 + r2 / 156))))))


def betaln_sym(a):
    """log B(a, a), accurate to ~1e-16 relative.

    jax.scipy.special.betaln loses ~5e-10 absolute at a = 1e5 and, for
    asymmetric arguments, ~2e-7 at B(1/2, 10); both matter here.
    """
    small = a < _STIRLING_SWITCH
    a_s = jnp.where(small, a, _STIRLING_SWITCH)
    a_l = jnp.where(small, _STIRLING_SWITCH, a)
    via_gamma = 2.0 * gammaln(a_s) - gammaln(2.0 * a_s)
    via_stirling = (
        -0.5 * jnp.log(a_l)
        - (2.0 * a_l - 0.5) * LOG_TWO
        + LOG_SQRT_2PI
        + 2.0 * _stirling_correction(a_l)
        - _stirling_correction(2.0 * a_l)
    )
    return jnp.where(small, via_gamma, via_stirling)


def betaln_half(a):
    """log B(1/2, a), accurate to ~1e-16 relative (see ``betaln_sym``)."""
    small = a < _STIRLING_SWITCH
    a_s = jnp.where(small, a, _STIRLING_SWITCH)
    a_l = jnp.where(small, _STIRLING_SWITCH, a)
    via_gamma = LOG_SQRT_PI + gammaln(a_s) - gammaln(a_s + 0.5)
    via_stirling = (
        LOG_SQRT_PI
        - 0.5 * jnp.log(a_l)
        + (0.5 - a_l * jnp.log1p(0.5 / a_l))
        + _stirling_correction(a_l)
        - _stirling_correction(a_l + 0.5)
    )
    return jnp.where(small, via_gamma, via_stirling)


def _s_min(dtype):
    """Smallest lower-half argument we evaluate at (keeps logs finite)."""
    return 1e-300 if dtype == jnp.float64 else 1e-35


def _betacf(p, q, x, n_iter):
    """Numerical Recipes continued fraction for I_x(p, q), fixed trip count."""
    dtype = x.dtype
    tiny = jnp.asarray(jnp.finfo(dtype).tiny * 1e4, dtype)

    def _clip(v):
        return jnp.where(jnp.abs(v) < tiny, tiny, v)

    eps = float(jnp.finfo(dtype).eps)
    qab, qap, qam = p + q, p + 1.0, p - 1.0

    def body(m, carry):
        c, d, h, done = carry
        m = m.astype(dtype)
        m2 = 2.0 * m
        aa = m * (q - m) * x / ((qam + m2) * (p + m2))
        d = 1.0 / _clip(1.0 + aa * d)
        c = _clip(1.0 + aa / c)
        h1 = h * d * c
        aa = -(p + m) * (qab + m) * x / ((p + m2) * (qap + m2))
        d = 1.0 / _clip(1.0 + aa * d)
        c = _clip(1.0 + aa / c)
        delta = d * c
        # Stop updating converged entries so round-off doesn't accumulate.
        h = jnp.where(done, h, h1 * delta)
        return c, d, h, done | (jnp.abs(delta - 1.0) <= eps)

    d = 1.0 / _clip(1.0 - qab * x / qap)
    init = (jnp.ones_like(x), d, d, jnp.zeros(x.shape, bool))
    return lax.fori_loop(1, n_iter + 1, body, init)[2]


def _log_link_cdf_value(u, kappa, n_iter=None):
    n_iter = CF_ITERS if n_iter is None else n_iter
    dtype = u.dtype
    a = kappa
    abs_u = jnp.abs(u)
    lo = jnp.maximum((PI2 - abs_u) / TWO_PI2, _s_min(dtype))
    log_lo, log_hi = jnp.log(lo), jnp.log1p(-lo)
    root_w = abs_u / PI2  # = 1 - 2 lo, exactly
    w = root_w * root_w
    use_b = w < jnp.minimum(SWITCH_C / (a + 2.5), SWITCH_CAP)

    # log(lo * hi): lo * hi = (1 - w) / 4 exactly; log1p(-w) avoids the
    # cancellation in log lo + log hi near the centre, which kappa multiplies.
    log_lohi = jnp.where(w < 0.5, jnp.log1p(-w) - LOG_FOUR, log_lo + log_hi)

    # One continued-fraction pass with per-element parameters.
    p = jnp.where(use_b, 0.5, a)
    x = jnp.where(use_b, w, lo)
    cf = _betacf(p, a, x, n_iter)
    log_cf = jnp.log(cf)

    # Form A: log I_lo directly.
    log_a = a * log_lohi - jnp.log(a) - betaln_sym(a) + log_cf
    a_lower, a_upper = log_a, jnp.log1p(-jnp.exp(log_a))

    # Form B: I_w(1/2, a), with log(1 - w) = log 4 + log(lo * hi).
    safe_root = jnp.where(root_w > 0, root_w, 1.0)
    log_iw = jnp.log(safe_root) + a * (LOG_FOUR + log_lohi) - LOG_HALF - betaln_half(a) + log_cf
    iw = jnp.where(root_w > 0, jnp.exp(log_iw), 0.0)
    b_lower, b_upper = LOG_HALF + jnp.log1p(-iw), LOG_HALF + jnp.log1p(iw)

    lower = jnp.where(use_b, b_lower, a_lower)  # log I_lo
    upper = jnp.where(use_b, b_upper, a_upper)  # log (1 - I_lo)
    return jnp.where(u <= 0, lower, upper)


def log_link_pdf(u, kappa):
    """log g_kappa(u): density of the link distribution on [-pi^2, pi^2]."""
    s_min = _s_min(u.dtype)
    log_s = jnp.log(jnp.maximum((PI2 + u) / TWO_PI2, s_min))
    log_1ms = jnp.log(jnp.maximum((PI2 - u) / TWO_PI2, s_min))
    w = (u / PI2) ** 2  # s (1 - s) = (1 - w) / 4
    log_s1s = jnp.where(w < 0.5, jnp.log1p(-w) - LOG_FOUR, log_s + log_1ms)
    return (kappa - 1.0) * log_s1s - betaln_sym(kappa) - LOG_TWO_PI2


@jax.custom_jvp
def log_link_cdf(u, kappa):
    """log G_kappa(u) for u in [-pi^2, pi^2]; broadcasts over inputs."""
    u, kappa = jnp.broadcast_arrays(u, kappa)
    return _log_link_cdf_value(u, kappa)


@log_link_cdf.defjvp
def _log_link_cdf_jvp(primals, tangents):
    u, kappa = primals
    du, _ = tangents
    out = log_link_cdf(u, kappa)
    u_b, kappa_b = jnp.broadcast_arrays(u, kappa)
    dlog = jnp.exp(log_link_pdf(u_b, kappa_b) - out)
    return out, dlog * du


def log_beta_cdf_sym(s, kappa):
    """log I_s(kappa, kappa); convenience wrapper used in tests."""
    return log_link_cdf(TWO_PI2 * s - PI2, kappa)
