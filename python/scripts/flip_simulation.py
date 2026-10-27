"""Simulation check of the roll-call mirror move (``SamplerConfig.bill_flip``).

Simulates one House-sized data set from the circular factor model (sizes,
omega and kappa calibrated to the 116th House, first session) and fits it
with several chains per sampler variant:

  control:  the paper's scheme (jittering GHMC only)
  flip:     plus the mirror move with kappa rescaled
  flipfix:  plus the mirror move with kappa kept fixed

Each chain starts like the real runs: "Democrats" (true beta < 0) in
(-pi/2, 0), "Republicans" in (0, pi/2), Yea/Nay positions uniform.

  uv run python scripts/flip_simulation.py run --variant control --seed 1
  uv run python scripts/flip_simulation.py report
"""

import argparse
import json
import sys
from pathlib import Path

import jax.numpy as jnp
import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from simulation_study import align, coverage, wrap  # noqa: E402

from circfactor.analysis import circular_mean, circular_variance, ess_bulk, split_rhat  # noqa: E402
from circfactor.model import Data, Posterior, Priors, SamplerConfig, loglik_matrix, sample  # noqa: E402
from circfactor.special import log_link_cdf  # noqa: E402

ROOT = Path(__file__).resolve().parents[1] / "results" / "sim_flip"
VARIANTS = {"control": {}, "flip": dict(bill_flip=True), "flipfix": dict(bill_flip=True, bill_flip_kappa=False)}


def simulate(I=432, J=700, omega=2.0, log_kappa=(4.1, 0.85), seed=1):
    """Ideal points ~ vonMises(0, omega), Yea/Nay uniform, kappa log-normal."""
    rng = np.random.default_rng(seed)
    beta = stats.vonmises.rvs(omega, size=I, random_state=rng)
    yes, no = rng.uniform(-np.pi, np.pi, J), rng.uniform(-np.pi, np.pi, J)
    kappa = np.exp(rng.normal(*log_kappa, J))
    e = wrap(no[None] - beta[:, None]) ** 2 - wrap(yes[None] - beta[:, None]) ** 2
    p = np.exp(np.asarray(log_link_cdf(jnp.asarray(e), jnp.asarray(kappa)[None])))
    y = (rng.random((I, J)) < p).astype(float)
    y[rng.random((I, J)) < 0.03] = np.nan  # ~3% missing, as in the House data
    return y, dict(beta=beta, yes=yes, no=no, kappa=kappa, omega=omega)


def loglik(y, beta, yes, no, kappa):
    mask = ~np.isnan(y)
    data = Data(jnp.asarray(np.where(mask, 2 * np.nan_to_num(y) - 1, 0.0)), jnp.asarray(mask))
    return float(loglik_matrix(*(jnp.asarray(a) for a in (beta, yes, no, kappa)), data).sum())


def run(args):
    y, truth = simulate(seed=args.data_seed)
    out = ROOT / f"data{args.data_seed}" / f"{args.variant}_c{args.seed}"
    out.mkdir(parents=True, exist_ok=True)
    party = np.where(truth["beta"] < 0, "D", "R")
    cfg = SamplerConfig(burnin=args.burnin, n_samples=args.samples, seed=args.seed, print_every=2000, **VARIANTS[args.variant])
    post = sample(y, Priors(), cfg, party=party)
    post.save(out / "posterior.npz")
    print(f"{args.variant} c{args.seed}: accept {post.accept}, {post.seconds / 60:.1f} min")


def report(args):
    base = ROOT / f"data{args.data_seed}"
    y, t = simulate(seed=args.data_seed)
    chi_true = circular_variance(t["beta"])
    ll_true = loglik(y, t["beta"], t["yes"], t["no"], t["kappa"])
    rank_true = t["beta"].argsort().argsort()
    print(f"truth: chi_0 {chi_true:.4f}, omega {t['omega']}, loglik at the true parameters {ll_true:.0f}\n")
    rows = {}
    for variant in VARIANTS:
        dirs = sorted(base.glob(f"{variant}_c*"))
        if not dirs:
            continue
        posts = [Posterior.load(d / "posterior.npz") for d in dirs]
        # Align each chain to the truth (common rotation; reflection is fixed by
        # the start), so that ranks on (-pi, pi] don't wrap near +-pi.
        for p in posts:
            p.beta, sign, shift = align(p.beta, t["beta"])
            p.yes, p.no = wrap(sign * p.yes + shift), wrap(sign * p.no + shift)
        S = posts[0].beta.shape[0]
        chi = np.stack([circular_variance(p.beta, axis=1) for p in posts])
        om = np.stack([p.omega for p in posts])
        ranks = np.stack([p.beta.argsort(1).argsort(1) for p in posts])
        rh = np.array([split_rhat(ranks[:, :, i]) for i in range(ranks.shape[2])])
        rh = rh[~np.isnan(rh)]
        ll_end = [p.loglik_trace[-S // 2 :].mean() for p in posts]
        beta = np.concatenate([p.beta for p in posts])
        yes, no, kappa = (np.concatenate([getattr(p, k) for p in posts]) for k in ("yes", "no", "kappa"))
        # Arrangement of each roll call: is the posterior pair closer to the
        # true (psi, zeta) or to its mirror image (zeta + pi, psi + pi)?
        d_same = lambda p: np.abs(wrap(p.yes - t["yes"])) + np.abs(wrap(p.no - t["no"]))
        d_mirr = lambda p: np.abs(wrap(p.yes - t["no"] - np.pi)) + np.abs(wrap(p.no - t["yes"] - np.pi))
        same = np.stack([(d_same(p) < d_mirr(p)).mean(0) for p in posts])  # (chains, J) share of draws
        lo, hi = np.quantile(kappa, [0.05, 0.95], axis=0)
        r = dict(
            chains=len(posts),
            draws=S,
            minutes=round(float(np.mean([p.seconds for p in posts])) / 60, 1),
            accept_flip=round(float(np.mean([p.accept.get("flip", 0.0) for p in posts])), 4),
            loglik_per_chain=[round(float(v)) for v in ll_end],
            chi0_per_chain=chi.mean(1).round(4).tolist(),
            chi0_mean=round(float(chi.mean()), 4),
            chi0_90ci=np.quantile(chi, [0.05, 0.95]).round(4).tolist(),
            chi0_rhat=round(float(split_rhat(chi)), 3),
            chi0_ess=[round(ess_bulk(c)) for c in chi],
            omega_per_chain=om.mean(1).round(3).tolist(),
            omega_90ci=np.quantile(om, [0.05, 0.95]).round(3).tolist(),
            omega_rhat=round(float(split_rhat(om)), 3),
            rank_rhat_median=round(float(np.median(rh)), 3),
            rank_rhat_max=round(float(rh.max()), 3),
            rank_rhat_gt_1_01=int((rh > 1.01).sum()),
            rank_spearman=round(float(stats.spearmanr(rank_true, np.median(ranks.reshape(-1, ranks.shape[2]), 0)).statistic), 4),
            beta_mean_abs_err=round(float(np.mean(np.abs(wrap(circular_mean(beta, axis=0) - t["beta"])))), 4),
            beta_cov90=round(float(coverage(beta, t["beta"])), 3),
            kappa_cov90=round(float(np.mean((t["kappa"] >= lo) & (t["kappa"] <= hi))), 3),
            kappa_log_corr=round(float(np.corrcoef(np.log(kappa.mean(0)), np.log(t["kappa"]))[0, 1]), 3),
            bills_true_arrangement_per_chain=(same > 0.5).mean(1).round(3).tolist(),
            bills_mixing=int(((same > 0.05) & (same < 0.95)).any(0).sum()),
            bills_chains_disagree=int(((same > 0.5).any(0) & (same <= 0.5).any(0)).sum()),
            psi_cov90=round(float(coverage(yes, t["yes"])), 3),
            zeta_cov90=round(float(coverage(no, t["no"])), 3),
        )
        rows[variant] = r
    keys = list(next(iter(rows.values())))
    w = max(map(len, keys))
    print(f"{'':<{w}}  " + "  ".join(f"{v:<36}" for v in rows))
    for k in keys:
        print(f"{k:<{w}}  " + "  ".join(f"{str(rows[v][k]):<36}" for v in rows))
    (base / "report.json").write_text(json.dumps(dict(truth=dict(chi0=float(chi_true), omega=float(t["omega"]), loglik=ll_true), **rows), indent=1))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--variant", choices=list(VARIANTS), required=True)
    r.add_argument("--seed", type=int, default=1, help="chain seed")
    r.add_argument("--data-seed", type=int, default=1)
    r.add_argument("--burnin", type=int, default=10_000)
    r.add_argument("--samples", type=int, default=10_000)
    p = sub.add_parser("report")
    p.add_argument("--data-seed", type=int, default=1)
    args = ap.parse_args()
    run(args) if args.cmd == "run" else report(args)


if __name__ == "__main__":
    main()
