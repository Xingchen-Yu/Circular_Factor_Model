"""Posterior summaries used in the paper: ranks, circular variance, WAIC."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .model import Posterior


def circular_mean(theta, axis=0):
    return np.arctan2(np.sin(theta).mean(axis), np.cos(theta).mean(axis))


def circular_variance(theta, axis=-1):
    """chi = 1 - mean resultant length (paper Section 6.3)."""
    return 1.0 - np.hypot(np.cos(theta).mean(axis), np.sin(theta).mean(axis))


def orient(post: Posterior, party) -> Posterior:
    """Reflect the circle if needed so Democrats sit at negative angles.

    The likelihood is invariant to reflection (beta, psi, zeta) -> -(...);
    the party-based initialisation normally already fixes this.
    """
    party = np.asarray(party)
    dem, rep = party == "D", party == "R"
    if dem.any() and rep.any():
        m = circular_mean(post.beta, axis=0)
        if circular_mean(m[dem]) > circular_mean(m[rep]):
            post.beta, post.yes, post.no = -post.beta, -post.yes, -post.no
    return post


def legislator_summary(post: Posterior, legislators: pd.DataFrame) -> pd.DataFrame:
    """Per-legislator posterior summaries of ideal points and ranks.

    Ranks unwind the circle at the bottom pole (+-pi), as in Section 4.5:
    rank 1 is the most negative angle (most liberal under ``orient``).
    """
    beta = post.beta
    ranks = beta.argsort(axis=1).argsort(axis=1) + 1
    q = lambda a, p: np.quantile(a, p, axis=0)
    df = legislators.reset_index(drop=True).copy()
    df["beta_mean"] = circular_mean(beta, axis=0)
    df["beta_median"] = q(beta, 0.5)
    df["beta_q025"], df["beta_q975"] = q(beta, 0.025), q(beta, 0.975)
    df["beta_circ_var"] = circular_variance(beta.T)
    df["rank_median"] = q(ranks, 0.5)
    df["rank_q025"], df["rank_q975"] = q(ranks, 0.025), q(ranks, 0.975)
    return df


def chamber_summary(post: Posterior, party) -> dict:
    """Circularity chi_0 (all legislators) and within-party chi_D, chi_R per draw."""
    party = np.asarray(party)
    out = {}
    for name, sel in [("chi_0", slice(None)), ("chi_D", party == "D"), ("chi_R", party == "R")]:
        chi = circular_variance(post.beta[:, sel], axis=1)
        out[name] = {"mean": float(chi.mean()), "q025": float(np.quantile(chi, 0.025)), "q975": float(np.quantile(chi, 0.975))}
    out["waic_grouped"] = {k: float(v) for k, v in post.waic_grouped().items()}
    out["waic_pointwise"] = post.waic_pointwise
    out["omega_mean"] = float(post.omega.mean())
    out["accept"] = post.accept
    out["seconds"] = post.seconds
    return out


def bill_summary(post: Posterior, rollcalls: pd.DataFrame) -> pd.DataFrame:
    df = rollcalls.reset_index(drop=True).copy()
    df["psi_mean"] = circular_mean(post.yes, axis=0)
    df["zeta_mean"] = circular_mean(post.no, axis=0)
    df["kappa_mean"] = post.kappa.mean(0)
    df["kappa_q025"], df["kappa_q975"] = np.quantile(post.kappa, [0.025, 0.975], axis=0)
    return df


def ess_bulk(x):
    """Effective sample size of a 1-chain trace (Geyer initial positive sequence)."""
    x = np.asarray(x, float)
    n = len(x)
    x = x - x.mean()
    f = np.fft.rfft(x, 2 * n)
    acf = np.fft.irfft(f * np.conj(f))[:n]
    acf /= acf[0]
    pair = acf[: n - 1 : 2] + acf[1:n:2]
    k = np.argmax(pair < 0) if np.any(pair < 0) else len(pair)
    tau = -1 + 2 * pair[:k].sum()
    return n / max(tau, 1.0)


def split_rhat(chains):
    """Split-R-hat (Gelman et al. 2013) for an (n_chains, n_draws) array."""
    chains = np.asarray(chains, float)
    half = chains.shape[1] // 2
    x = np.concatenate([chains[:, :half], chains[:, half : 2 * half]], axis=0)
    n = x.shape[1]
    within = x.var(axis=1, ddof=1).mean()
    between = n * x.mean(axis=1).var(ddof=1)
    if within == 0:
        return np.nan
    return float(np.sqrt(((n - 1) / n * within + between / n) / within))
