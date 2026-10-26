"""Simulation check of placing low-vote legislators with a cut posterior.

For each (omega, seed) cell: simulate I full-record legislators plus K
low-vote legislators (each keeps only a fraction of their votes, missing at
random) from the circular factor model, then

  * cut:   fit the I full-record members, place the K others on a grid given
           the fit's draws (``circfactor.placement``);
  * joint: fit all I + K members together (what a looser missing-vote cutoff
           would do);
  * naive: (optional) the naive cut algorithm, one or m random-walk steps for
           the K members per randomly chosen posterior draw.

Saves the aligned draws and per-member metrics to <out>/omega{w}_seed{s}/.
``scripts/cut_posterior_report.py`` turns them into tables and figures.

  uv run python scripts/cut_posterior_simulation.py --omega 3 --seed 1 --naive
"""

import argparse
import json
import sys
from pathlib import Path

import jax
import jax.numpy as jnp
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from simulation_study import wrap  # noqa: E402

from circfactor.analysis import circular_mean, circular_variance  # noqa: E402
from circfactor.model import Priors, SamplerConfig, sample  # noqa: E402
from circfactor.placement import place_legislators  # noqa: E402
from circfactor.special import log_link_cdf  # noqa: E402

FRACTIONS = (0.03, 0.10, 0.25, 0.50)  # share of roll calls a low-vote member voted on


def simulate(I, K, J, omega, seed, kappa_shape=2.0, kappa_rate=0.05):
    rng = np.random.default_rng(seed)
    beta = stats.vonmises.rvs(omega, size=I + K, random_state=rng)
    yes, no = rng.uniform(-np.pi, np.pi, J), rng.uniform(-np.pi, np.pi, J)
    kappa = rng.gamma(kappa_shape, 1 / kappa_rate, J)
    e = wrap(no[None] - beta[:, None]) ** 2 - wrap(yes[None] - beta[:, None]) ** 2
    p = np.exp(np.asarray(log_link_cdf(jnp.asarray(e), jnp.asarray(kappa)[None])))
    y = (rng.random((I + K, J)) < p).astype(float)
    y[:I][rng.random((I, J)) < 0.03] = np.nan  # full-record members: ~3% missing, as in the House
    frac = np.array([FRACTIONS[k % len(FRACTIONS)] for k in range(K)])
    for k in range(K):
        y[I + k, rng.permutation(J)[int(round(frac[k] * J)) :]] = np.nan
    return y, dict(beta=beta, yes=yes, no=no, kappa=kappa), frac


def frame(fit_draws, truth_fit):
    """Reflection (global) and per-draw rotation mapping a fit onto the truth.

    The likelihood only sees angle differences, so each draw is rotated to
    best match the true full-record ideal points; the same map is applied to
    every angle of that draw.
    """
    best = None
    for sign in (1, -1):
        shifts = circular_mean(wrap(truth_fit[None] - sign * fit_draws), axis=1)
        err = np.mean(1 - np.cos(truth_fit[None] - (sign * fit_draws + shifts[:, None])))
        if best is None or err < best[0]:
            best = (err, sign, shifts)
    return best[1], best[2]


def apply_frame(draws, sign, shifts):
    return wrap(sign * draws + shifts[:, None])


def arc_interval(draws, level=0.9):
    """Central credible arc: (centre, half-width) per column."""
    c = circular_mean(draws, axis=0)
    return c, np.quantile(np.abs(wrap(draws - c)), level, axis=0)


def member_metrics(draws, truth, fit_aligned, truth_fit):
    """Error, 90% coverage and width for angles; the same for ranks among the fitted members."""
    c, half = arc_interval(draws)
    err = wrap(c - truth)
    rank = (fit_aligned[:, :, None] < draws[:, None, :]).sum(1) + 1  # (S, K) rank among fitted members
    rank_true = (truth_fit[:, None] < truth[None, :]).sum(0) + 1
    lo, hi = np.quantile(rank, [0.05, 0.95], axis=0)
    return dict(
        mean=c,
        err=err,
        covered=np.abs(err) <= half,
        width=2 * half,
        sd=np.sqrt(-2 * np.log(np.hypot(np.cos(draws).mean(0), np.sin(draws).mean(0)))),
        rank_median=np.median(rank, 0),
        rank_true=rank_true,
        rank_covered=(rank_true >= lo) & (rank_true <= hi),
        rank_width=hi - lo,
    )


def naive_cut(post, y_low, init, prop_sd, n_iter, inner, seed):
    """Naive cut MCMC: each iteration picks a random posterior draw theta^(s)
    and takes ``inner`` random-walk Metropolis steps for the low-vote members
    given theta^(s).  Returns (n_iter, K) states."""
    S = post.yes.shape[0]
    K = y_low.shape[0]
    obs = ~np.isnan(y_low)
    n = int(obs.sum(1).max())
    idx = np.zeros((K, n), np.int32)
    sgn = np.zeros((K, n))
    for k in range(K):  # compact vote lists, as in circfactor.placement
        j = np.flatnonzero(obs[k])
        idx[k, : len(j)], sgn[k, : len(j)] = j, 2 * y_low[k, j] - 1
    idx, sgn = jnp.asarray(idx), jnp.asarray(sgn)
    yes, no, kap, om = (jnp.asarray(a, jnp.float64) for a in (post.yes, post.no, post.kappa, post.omega))
    prop_sd = jnp.asarray(prop_sd)

    def logp(b, s):  # (K,)
        u = sgn * (wrap(no[s][idx] - b[:, None]) ** 2 - wrap(yes[s][idx] - b[:, None]) ** 2)
        ll = jnp.where(sgn != 0, log_link_cdf(u, kap[s][idx]), 0.0).sum(1)
        return ll + om[s] * jnp.cos(b)

    def step(b, key):
        k_s, k_in = jax.random.split(key)
        s = jax.random.randint(k_s, (), 0, S)

        def mh(b, k):
            k1, k2 = jax.random.split(k)
            bp = wrap(b + prop_sd * jax.random.normal(k1, b.shape))
            acc = jnp.log(jax.random.uniform(k2, b.shape)) < logp(bp, s) - logp(b, s)
            return jnp.where(acc, bp, b), None

        b, _ = jax.lax.scan(mh, b, jax.random.split(k_in, inner))
        return b, b

    _, trace = jax.jit(lambda b0, keys: jax.lax.scan(step, b0, keys))(jnp.asarray(init), jax.random.split(jax.random.key(seed), n_iter))
    return np.asarray(trace)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--omega", type=float, default=3.0, help="true concentration of the ideal points")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--I", type=int, default=200, help="full-record legislators")
    ap.add_argument("--K", type=int, default=40, help="low-vote legislators")
    ap.add_argument("--J", type=int, default=600)
    ap.add_argument("--burnin", type=int, default=10_000)
    ap.add_argument("--samples", type=int, default=1000)
    ap.add_argument("--thin", type=int, default=5)
    ap.add_argument("--naive", action="store_true", help="also run the naive cut algorithm")
    ap.add_argument("--naive-iter", type=int, default=20_000)
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parents[1] / "results" / "sim_cut")
    args = ap.parse_args()

    out = args.out / f"omega{args.omega:g}_seed{args.seed}"
    out.mkdir(parents=True, exist_ok=True)
    I, K = args.I, args.K
    y, truth, frac = simulate(I, K, args.J, args.omega, args.seed)
    tb = truth["beta"]
    print(f"omega={args.omega}, seed={args.seed}: I={I} + K={K} low-vote, J={args.J}; "
          f"low-vote members cast {np.sort(np.unique((~np.isnan(y[I:])).sum(1)))} votes", flush=True)
    cfg = lambda s: SamplerConfig(burnin=args.burnin, n_samples=args.samples, thin=args.thin, seed=s, print_every=5000)

    # cut: fit the full-record members, then place the low-vote members
    post = sample(y[:I], Priors(), cfg(args.seed))
    placed = place_legislators(post.yes, post.no, post.kappa, post.omega, y[I:], seed=args.seed)
    sign, shifts = frame(post.beta, tb[:I])
    fit_cut = apply_frame(post.beta, sign, shifts)
    low_cut = apply_frame(placed, sign, shifts)

    # joint: fit everyone
    post_j = sample(y, Priors(), cfg(args.seed + 1000))
    sign_j, shifts_j = frame(post_j.beta[:, :I], tb[:I])
    all_j = apply_frame(post_j.beta, sign_j, shifts_j)
    fit_joint, low_joint = all_j[:, :I], all_j[:, I:]

    m_cut = member_metrics(low_cut, tb[I:], fit_cut, tb[:I])
    m_joint = member_metrics(low_joint, tb[I:], fit_joint, tb[:I])
    rows = {"member": np.arange(K), "frac": frac, "n_votes": (~np.isnan(y[I:])).sum(1), "truth": tb[I:]}
    rows |= {f"cut_{k}": v for k, v in m_cut.items()} | {f"joint_{k}": v for k, v in m_joint.items() if k != "rank_true"}

    # effect of including the low-vote members on the full-record members
    c_cut, _ = arc_interval(fit_cut)
    c_joint, _ = arc_interval(fit_joint)
    sd_cut = member_metrics(fit_cut, tb[:I], fit_cut, tb[:I])["sd"]
    summary = {
        "omega": args.omega,
        "seed": args.seed,
        "I": I,
        "K": K,
        "J": args.J,
        "fit_shift_over_sd_median": float(np.median(np.abs(wrap(c_joint - c_cut)) / sd_cut)),
        "fit_shift_over_sd_max": float(np.max(np.abs(wrap(c_joint - c_cut)) / sd_cut)),
        "chi0_cut": float(circular_variance(post.beta, axis=1).mean()),
        "chi0_joint_fitted_only": float(circular_variance(post_j.beta[:, :I], axis=1).mean()),
        "chi0_true": float(circular_variance(tb[:I])),
        "omega_cut": float(post.omega.mean()),
        "omega_joint": float(post_j.omega.mean()),
        "seconds_cut_fit": post.seconds,
        "seconds_joint_fit": post_j.seconds,
    }
    arrays = dict(low_cut=low_cut, low_joint=low_joint, fit_cut=fit_cut, fit_joint=fit_joint, truth=tb, frac=frac,
                  loglik_cut=post.loglik_trace, loglik_joint=post_j.loglik_trace)

    if args.naive:
        # Proposal sd = each member's cut-posterior sd: a well-tuned random walk.
        prop = np.clip(member_metrics(placed, tb[I:], fit_cut, tb[:I])["sd"], 0.01, 1.0)
        init = circular_mean(placed, axis=0)
        arrays["placed_raw"] = placed
        for inner in (1, 10):
            tr = naive_cut(post, y[I:], init, prop, n_iter=args.naive_iter, inner=inner, seed=args.seed)[args.naive_iter // 10 :]
            arrays[f"naive_{inner}"] = tr  # fit frame (not paired with draws), compare with placed_raw
            # Compare in the fit's own frame: the naive draws are not paired with
            # posterior draws, so they cannot be aligned to the truth draw by draw.
            csd = lambda x: np.sqrt(-2 * np.log(np.hypot(np.cos(x).mean(0), np.sin(x).mean(0))))
            sd_grid = csd(placed)
            rows[f"naive{inner}_mean_shift_over_sd"] = np.abs(wrap(circular_mean(tr, 0) - circular_mean(placed, 0))) / sd_grid
            rows[f"naive{inner}_sd_ratio"] = csd(tr) / sd_grid

    pd.DataFrame(rows).to_csv(out / "members.csv", index=False)
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    np.savez_compressed(out / "draws.npz", **arrays)
    df = pd.DataFrame(rows)
    print(df.groupby("frac")[["cut_covered", "joint_covered", "cut_width", "joint_width", "cut_rank_covered"]].mean().round(3))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
