"""Geodesic Hamiltonian Monte Carlo on the unit sphere (Byrne & Girolami, 2013).

The kernel advances a *batch* of independent chains at once.  Each row of
``x`` (shape ``(B, d)``) is a point on the unit sphere S^{d-1} embedded in
R^d, and ``log_prob(x)`` must return the ``(B,)`` vector of per-chain log
densities with respect to the Hausdorff (surface) measure, with row ``b``
depending only on ``x[b]``.  That independence is what lets one call to
``jax.vjp`` deliver every chain's gradient at once.

For the circular factor model d = 2, every legislator (or bill) is one chain,
and the batch is the whole block of a Gibbs sweep.

Each chain has its own step size and number of leapfrog steps (the paper's
"jittering").  Chains with fewer steps are frozen once they finish, so the cost
of a call is set by ``max_leap``.

``grad_fn`` may differ from the exact gradient of ``log_prob`` (for example,
the same density evaluated in float32).  The leapfrog map with momentum flip
is reversible and volume preserving for *any* deterministic vector field, and
the Metropolis correction uses ``log_prob``, so the chain targets
``log_prob`` exactly; a poor gradient only lowers the acceptance rate.
"""

from typing import Callable, NamedTuple

import jax
import jax.numpy as jnp
from jax import lax


class GHMCInfo(NamedTuple):
    accepted: jax.Array  # (B,) bool
    accept_prob: jax.Array  # (B,) min(1, exp(dH))
    log_prob: jax.Array  # (B,) log density at the returned state
    aux: object = None  # log_prob's auxiliary output at the returned state


def project(x, v):
    """Project v onto the tangent space of the sphere at x (row-wise)."""
    return v - jnp.sum(x * v, axis=-1, keepdims=True) * x


def _geodesic_flow(x, v, eps):
    """Exact geodesic flow on the sphere for time eps (row-wise)."""
    alpha = jnp.linalg.norm(v, axis=-1, keepdims=True)
    ae = alpha * eps
    cos_ae, sin_ae = jnp.cos(ae), jnp.sin(ae)
    # sin(alpha*eps)/alpha -> eps as alpha -> 0
    safe_alpha = jnp.where(alpha > 0, alpha, 1.0)
    sinc_term = jnp.where(alpha > 0, sin_ae / safe_alpha, eps)
    x_new = x * cos_ae + v * sinc_term
    v_new = v * cos_ae - alpha * x * sin_ae
    # Renormalise to stop round-off drift off the manifold.
    x_new = x_new / jnp.linalg.norm(x_new, axis=-1, keepdims=True)
    return x_new, project(x_new, v_new)


def default_grad_fn(log_prob: Callable) -> Callable:
    def grad_fn(x):
        lp, vjp = jax.vjp(log_prob, x)
        return vjp(jnp.ones_like(lp))[0]

    return grad_fn


def ghmc_step(
    key,
    x,
    log_prob: Callable,
    step_size,
    n_leap,
    max_leap: int,
    grad_fn: Callable | None = None,
    log_prob_x=None,
    has_aux: bool = False,
):
    """One GHMC transition for a batch of independent chains on the sphere.

    Args:
      key: PRNG key.
      x: (B, d) current points, unit norm along the last axis.
      log_prob: maps (B, d) -> (B,) log densities (exact target).
      step_size: (B,) or scalar leapfrog step sizes.
      n_leap: (B,) or scalar integer leapfrog step counts, 1 <= n_leap <= max_leap.
      max_leap: static upper bound on n_leap.
      grad_fn: maps (B, d) -> (B, d) (ambient) gradient used for the momentum
        kicks; defaults to the autodiff gradient of ``log_prob``.
      log_prob_x: optional precomputed ``log_prob(x)`` (with ``has_aux``, the
        pair ``(log_prob, aux)`` at ``x``).
      has_aux: ``log_prob`` returns ``(log_prob, aux)``, where every leaf of
        ``aux`` has the chain axis first; the aux of the returned state is
        reported in ``GHMCInfo.aux`` (accepted rows from the proposal, the
        rest from ``x``).  ``grad_fn`` is then required.

    Returns:
      (x_new, GHMCInfo)
    """
    if grad_fn is None:
        if has_aux:
            raise ValueError("grad_fn is required when has_aux=True")
        grad_fn = default_grad_fn(log_prob)
    lp_fn = log_prob if has_aux else (lambda x: (log_prob(x), None))
    if log_prob_x is not None and not has_aux:
        log_prob_x = (log_prob_x, None)
    B = x.shape[0]
    eps = jnp.broadcast_to(jnp.asarray(step_size, x.dtype), (B,))[:, None]
    n_leap = jnp.broadcast_to(jnp.asarray(n_leap, jnp.int32), (B,))
    half = 0.5 * eps

    key_v, key_u = jax.random.split(key)
    v0 = project(x, jax.random.normal(key_v, x.shape, x.dtype))
    lp0, aux0 = lp_fn(x) if log_prob_x is None else log_prob_x
    h0 = lp0 - 0.5 * jnp.sum(v0 * v0, axis=-1)

    g0 = grad_fn(x).astype(x.dtype)

    def leapfrog(k, carry):
        x, v, g = carry
        active = (k < n_leap)[:, None]
        v_half = project(x, v + half * g)
        x_new, v_new = _geodesic_flow(x, v_half, eps)
        g_new = grad_fn(x_new).astype(x.dtype)
        v_new = project(x_new, v_new + half * g_new)
        return (
            jnp.where(active, x_new, x),
            jnp.where(active, v_new, v),
            jnp.where(active, g_new, g),
        )

    x1, v1, _ = lax.fori_loop(0, max_leap, leapfrog, (x, v0, g0))
    lp1, aux1 = lp_fn(x1)
    h1 = lp1 - 0.5 * jnp.sum(v1 * v1, axis=-1)

    log_ratio = jnp.where(jnp.isnan(h1), -jnp.inf, h1 - h0)
    accept = jnp.log(jax.random.uniform(key_u, (B,), x.dtype)) < log_ratio
    x_out = jnp.where(accept[:, None], x1, x)
    lp_out = jnp.where(accept, lp1, lp0)
    accept_prob = jnp.exp(jnp.minimum(log_ratio, 0.0))
    aux_out = jax.tree.map(lambda a1, a0: jnp.where(accept.reshape((B,) + (1,) * (a1.ndim - 1)), a1, a0), aux1, aux0)
    return x_out, GHMCInfo(accept, accept_prob, lp_out, aux_out)


def angle_to_circle(theta):
    """theta (...,) -> (..., 2) unit vectors (cos, sin)."""
    return jnp.stack([jnp.cos(theta), jnp.sin(theta)], axis=-1)


def circle_to_angle(x):
    """(..., 2) unit vectors -> angles in (-pi, pi]."""
    return jnp.arctan2(x[..., 1], x[..., 0])


def sample_von_mises(key, mu, kappa):
    """One exact draw from vonMises(mu, kappa) (Best & Fisher 1979 rejection).

    Scalar mu and kappa; kappa = 0 gives the uniform distribution.
    """
    kappa = jnp.maximum(jnp.asarray(kappa, jnp.float64), 1e-12)
    s = jnp.sqrt(1.0 + 4.0 * kappa**2)
    tau = 1.0 + s
    # rho = (tau - sqrt(2 tau)) / (2 kappa), rewritten without the 0/0
    # cancellation at small kappa (which would make the loop never accept).
    rho = 2.0 * kappa * tau / ((s + 1.0) * (tau + jnp.sqrt(2.0 * tau)))
    r = (1.0 + rho**2) / (2.0 * rho)

    def cond(state):
        return ~state[0]

    def body(state):
        _, _, key = state
        key, k1, k2 = jax.random.split(key, 3)
        u1, u2 = jax.random.uniform(k1), jax.random.uniform(k2)
        z = jnp.cos(jnp.pi * u1)
        f = (1.0 + r * z) / (r + z)
        c = kappa * (r - f)
        ok = (c * (2.0 - c) - u2 > 0) | (jnp.log(c / u2) + 1.0 - c >= 0)
        return ok, f, key

    key, k_sign = jax.random.split(key)
    _, f, _ = lax.while_loop(cond, body, (jnp.asarray(False), jnp.asarray(0.0), key))
    sign = jnp.where(jax.random.uniform(k_sign) < 0.5, -1.0, 1.0)
    theta = mu + sign * jnp.arccos(jnp.clip(f, -1.0, 1.0))
    return jnp.mod(theta + jnp.pi, 2 * jnp.pi) - jnp.pi
