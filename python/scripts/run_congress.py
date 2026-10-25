"""Fit the circular factor model to one or more Houses of Congress.

Examples:
  # The paper's 116th House analysis (first session only)
  uv run python scripts/run_congress.py 116 --sessions 1
  # 100th to 119th Houses, one after another
  uv run python scripts/run_congress.py 100-119

Outputs go to results/H{congress}[_s{sessions}]/:
  posterior.npz     draws of beta, omega, xi, per-legislator log-lik (and of
                    psi, zeta, kappa with --save-bill-draws)
  legislators.csv   ideal-point and rank summaries
  rollcalls.csv     psi / zeta / kappa summaries
  summary.json      circularity chi_0 / chi_D / chi_R, WAIC, acceptance, timing
"""

import argparse
import json
from pathlib import Path

import numpy as np

from circfactor.analysis import bill_summary, chamber_summary, legislator_summary, orient
from circfactor.data import load_congress
from circfactor.model import Priors, SamplerConfig, sample

ROOT = Path(__file__).resolve().parents[1]


def parse_congresses(spec: str) -> list[int]:
    out = []
    for part in spec.split(","):
        lo, _, hi = part.partition("-")
        out += list(range(int(lo), int(hi or lo) + 1))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("congresses", help="e.g. 116, 112,116 or 100-119")
    ap.add_argument("--chamber", default="House")
    ap.add_argument("--sessions", type=int, nargs="*", help="restrict to these sessions")
    ap.add_argument("--burnin", type=int, default=30_000)
    ap.add_argument("--samples", type=int, default=20_000)
    ap.add_argument("--thin", type=int, default=1)
    ap.add_argument("--seed", type=int, default=2021)
    ap.add_argument("--beta-step", type=float, nargs=2, default=(0.01, 0.04))
    ap.add_argument("--tau-step", type=float, nargs=2, default=(0.01, 0.105))
    ap.add_argument("--no-adapt", action="store_true", help="keep the step-size ranges fixed (as in the paper)")
    ap.add_argument("--rotation", action="store_true", help="enable the common-rotation group move")
    ap.add_argument("--grad-dtype", default="float32", choices=["float32", "float64"])
    ap.add_argument("--omega-mean", type=float, default=10.0, help="mean of the exponential prior on omega")
    ap.add_argument("--xi-rate", type=float, default=25.0, help="rate of the exponential prior on xi = 1/lambda")
    ap.add_argument("--save-bill-draws", action="store_true")
    ap.add_argument("--out", type=Path, default=ROOT / "results")
    ap.add_argument("--print-every", type=int, default=2000)
    args = ap.parse_args()

    priors = Priors(omega_rate=1.0 / args.omega_mean, xi_rate=args.xi_rate)
    for congress in parse_congresses(args.congresses):
        tag = f"{args.chamber[0]}{congress}" + (f"_s{''.join(map(str, args.sessions))}" if args.sessions else "")
        out = args.out / tag
        out.mkdir(parents=True, exist_ok=True)
        data = load_congress(congress, args.chamber, sessions=tuple(args.sessions) if args.sessions else None)
        print(data.summary(), flush=True)
        cfg = SamplerConfig(
            n_samples=args.samples,
            burnin=args.burnin,
            thin=args.thin,
            beta_step=tuple(args.beta_step),
            tau_step=tuple(args.tau_step),
            adapt_steps=not args.no_adapt,
            grad_dtype=args.grad_dtype,
            rotation_move=args.rotation,
            seed=args.seed,
            print_every=args.print_every,
        )
        party = data.legislators["party"].to_numpy()
        post = orient(sample(data.votes, priors, cfg, party=party), party)

        legislator_summary(post, data.legislators).to_csv(out / "legislators.csv", index=False)
        bill_summary(post, data.rollcalls).to_csv(out / "rollcalls.csv", index=False)
        summary = chamber_summary(post, party) | {"congress": congress, "chamber": args.chamber, "I": int(data.shape[0]), "J": int(data.shape[1])}
        (out / "summary.json").write_text(json.dumps(summary, indent=2))
        if not args.save_bill_draws:
            S = post.beta.shape[0]
            post.yes = post.no = post.kappa = np.empty((S, 0), np.float32)
        post.save(out / "posterior.npz")
        print(f"{tag}: chi_0 = {summary['chi_0']['mean']:.3f}, grouped elpd_waic = {summary['waic_grouped']['elpd_waic']:.0f}, "
              f"{post.seconds / 60:.1f} min -> {out}", flush=True)


if __name__ == "__main__":
    main()
