"""Grid sampler for placing excluded legislators (``placement.py``).

With the roll-call parameters held at one fixed draw, the placed angles must
follow p(beta | theta, y_k) exactly (up to the grid).  Checked against the
von Mises prior when the member cast no votes, and against the exact CDF from
scipy's incomplete beta function when they did.  (GHMC chains are not a
usable reference here: with few votes the conditional has local modes with
negligible mass, and chains started across the circle get stuck in them.)
"""


import jax.numpy as jnp
import numpy as np
import pytest
from scipy import stats

from scipy.special import betainc

from circfactor.placement import place_legislators
from circfactor.special import log_link_cdf

N = 5000
P_MIN = 1e-3


def wrap(a):
    return (np.asarray(a) + np.pi) % (2 * np.pi) - np.pi


def fixed_theta(J, seed):
    rng = np.random.default_rng(seed)
    return rng.uniform(-np.pi, np.pi, J), rng.uniform(-np.pi, np.pi, J), rng.gamma(2.0, 20.0, J)


def tile(a, S):
    return np.broadcast_to(np.asarray(a), (S,) + np.shape(a))


@pytest.mark.parametrize("mu,omega", [(0.0, 0.0), (0.0, 2.0), (1.0, 50.0)])
def test_no_votes_gives_prior(mu, omega):
    yes, no, kappa = fixed_theta(50, 0)
    votes = np.full((1, 50), np.nan)
    draws = place_legislators(tile(yes, N), tile(no, N), tile(kappa, N), np.full(N, omega), votes, mu=mu, seed=1)
    delta = wrap(draws[:, 0] - mu)
    cdf = (lambda d: (d + np.pi) / (2 * np.pi)) if omega == 0 else (lambda d: stats.vonmises.cdf(d, omega))
    assert stats.kstest(delta, cdf).pvalue > P_MIN


def reference_cdf(yes, no, kappa, y, omega, n=200_000):
    """CDF of p(beta | theta, y) on a fine grid, from scipy.special.betainc.

    Independent of ``log_link_cdf``: G_kappa(u) = I_s(kappa, kappa) with
    s = (pi^2 + u) / (2 pi^2) (paper eq. 5).
    """
    g = -np.pi + 2 * np.pi * (np.arange(n) + 0.5) / n
    o = ~np.isnan(y)
    lp = omega * np.cos(g)
    for j in np.flatnonzero(o):
        u = (2 * y[j] - 1) * (wrap(no[j] - g) ** 2 - wrap(yes[j] - g) ** 2)
        lp = lp + np.log(betainc(kappa[j], kappa[j], (np.pi**2 + u) / (2 * np.pi**2)))
    w = np.exp(lp - lp.max())
    c = np.cumsum(w) / w.sum()
    return lambda x: np.interp(x, g, c, left=0.0, right=1.0)


@pytest.mark.parametrize("n_votes,dtype", [(5, "float32"), (30, "float32"), (200, "float32"), (200, "float64")])
def test_matches_exact_conditional(n_votes, dtype):
    J = 200
    yes, no, kappa = fixed_theta(J, 2)
    rng = np.random.default_rng(3)
    beta_true, omega = 0.8, 3.0
    e = wrap(no - beta_true) ** 2 - wrap(yes - beta_true) ** 2
    p = np.exp(np.asarray(log_link_cdf(jnp.asarray(e), jnp.asarray(kappa))))
    y = (rng.random(J) < p).astype(float)
    y[rng.permutation(J)[n_votes:]] = np.nan
    draws = place_legislators(tile(yes, N), tile(no, N), tile(kappa, N), np.full(N, omega), y[None], seed=4, dtype=dtype)[:, 0]
    assert stats.kstest(draws, reference_cdf(yes, no, kappa, y, omega)).pvalue > P_MIN
