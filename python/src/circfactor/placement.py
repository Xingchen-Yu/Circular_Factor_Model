"""Place legislators who were left out of a fit (a "cut" posterior).

Members with few votes (partial terms, the Speaker) are excluded from the
fit by the missing-vote rule.  Given the fit's posterior draws
theta^(s) = (psi, zeta, kappa, omega)^(s), each excluded member k is placed by
drawing

    beta_k^(s) ~ p(beta | theta^(s), y_k)
             ∝ prod_j G_kappa_j(y_kj e_kj(beta)) * exp(omega cos(beta - mu)),

one draw per posterior draw s.  The pairs (theta^(s), beta_k^(s)) are draws
from the cut distribution p(theta | y_fit) p(beta_k | theta, y_k) (Plummer
2015): information flows from the roll calls to the member, but the member's
own votes do not update the roll-call parameters.

Because beta is one angle, each conditional is sampled directly on a fine
grid over the circle instead of by an inner MCMC chain.  The naive cut
algorithm, which alternates one MCMC step for beta with a fresh theta, does
not converge to the cut distribution (Plummer 2015); a direct draw for every
theta^(s) has no such problem.  The only approximation is the grid: the density
is treated as constant within each of ``n_grid`` cells (width 2 pi / n_grid,
0.0015 rad for the default 4096).
"""

from __future__ import annotations

from functools import partial

import jax
import jax.numpy as jnp
import numpy as np

from .model import wrap
from .special import log_link_cdf


def _place_one_draw(key, yes, no, kappa, omega, idx, sign, mu, n_grid, dtype):
    """(K,) draws of the placed members' angles for one posterior draw.

    ``idx``/``sign`` are (K, n) compact vote lists: roll-call indices and +-1,
    with sign 0 for padding.
    """
    width = 2 * jnp.pi / n_grid
    grid = (-jnp.pi + width * (jnp.arange(n_grid) + 0.5)).astype(dtype)  # cell centres

    def member_ll(idx_k, sign_k):  # (G,) log-likelihood of one member's votes
        d_yes = wrap(yes.astype(dtype)[idx_k][None, :] - grid[:, None])  # (G, n)
        d_no = wrap(no.astype(dtype)[idx_k][None, :] - grid[:, None])
        u = sign_k.astype(dtype)[None, :] * (d_no * d_no - d_yes * d_yes)
        ll = log_link_cdf(u, kappa.astype(dtype)[idx_k][None, :])
        return jnp.where(sign_k[None, :] != 0, ll, 0.0).sum(1)

    ll = jax.lax.map(lambda a: member_ll(*a), (idx, sign))  # (K, G)
    logits = ll.astype(jnp.float64) + omega * jnp.cos(grid.astype(jnp.float64) - mu)[None, :]
    k_cell, k_off = jax.random.split(key)
    cell = jax.random.categorical(k_cell, logits, axis=-1)  # (K,)
    offset = jax.random.uniform(k_off, cell.shape, minval=-0.5, maxval=0.5)
    return wrap(-jnp.pi + width * (cell + 0.5 + offset))


@partial(jax.jit, static_argnames=("n_grid", "dtype"))
def _place_chunk(keys, yes, no, kappa, omega, idx, sign, mu, n_grid, dtype):
    """_place_one_draw for a chunk of posterior draws (leading axis), in sequence."""
    f = lambda a: _place_one_draw(*a, idx, sign, mu, n_grid, dtype)
    return jax.lax.map(f, (keys, yes, no, kappa, omega))


def place_legislators(
    yes,
    no,
    kappa,
    omega,
    votes,
    mu: float = 0.0,
    n_grid: int = 4096,
    seed: int = 0,
    dtype: str = "float32",
    chunk: int = 64,
) -> np.ndarray:
    """Cut-posterior draws of ideal points for members left out of a fit.

    Args:
      yes, no, kappa: (S, J) posterior draws of psi, zeta, kappa from the fit.
      omega: (S,) posterior draws of omega.
      votes: (K, J) votes of the members to place (1 Yea / 0 Nay / NaN), on
        the same roll calls (columns) as the fit.
      mu: prior mean direction of the ideal points (``Priors.mu``).
      dtype: precision of the log-likelihood on the grid.  float32 is several
        times faster on consumer GPUs; summing a few hundred float32 terms
        changes the log density by ~1e-4, which does not matter for the draw.
      chunk: posterior draws per GPU call.

    Returns:
      (S, K) array of angles, row s paired with posterior draw s.
    """
    votes = np.asarray(votes, dtype=float)
    K = votes.shape[0]
    obs = ~np.isnan(votes)
    n = max(int(obs.sum(1).max()), 1)
    idx = np.zeros((K, n), np.int32)
    sign = np.zeros((K, n))
    for k in range(K):
        j = np.flatnonzero(obs[k])
        idx[k, : len(j)] = j
        sign[k, : len(j)] = 2 * votes[k, j] - 1
    idx_j, sign_j = jnp.asarray(idx), jnp.asarray(sign)
    S = np.shape(yes)[0]
    keys = jax.random.split(jax.random.key(seed), S)
    out = np.empty((S, K))
    for a in range(0, S, chunk):
        b = min(a + chunk, S)
        pad = lambda x: np.concatenate([x[a:b], np.repeat(x[b - 1 : b], chunk - (b - a), 0)])  # fixed shape: no recompiles
        res = _place_chunk(
            jnp.concatenate([keys[a:b], jnp.repeat(keys[b - 1 : b], chunk - (b - a), 0)]),
            jnp.asarray(pad(np.asarray(yes))),
            jnp.asarray(pad(np.asarray(no))),
            jnp.asarray(pad(np.asarray(kappa))),
            jnp.asarray(pad(np.asarray(omega, np.float64))),
            idx_j,
            sign_j,
            mu,
            n_grid,
            dtype,
        )
        out[a:b] = np.asarray(res)[: b - a]
    return out
