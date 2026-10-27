"""Fit the circular factor model to one or more Houses of Congress.

Examples:
  # The paper's 116th House analysis (first session only)
  uv run python scripts/run_congress.py 116 --sessions 1
  # 100th to 119th Houses, one after another
  uv run python scripts/run_congress.py 100-119
  # Continue those chains for 50k more draws, treating the first run as burn-in
  uv run python scripts/run_congress.py 100-119 --resume results --samples 50000 \
      --discard-old --out results_ext

Outputs go to results/H{congress}[_s{sessions}]/:
  posterior.npz     draws of beta, omega, xi, per-legislator log-lik, a thinned
                    set of psi, zeta, kappa draws (--bill-draws) and the final
                    sampler state (for --resume)
  placed.npz        cut-posterior draws for members excluded by the
                    missing-vote rule (see circfactor/placement.py)
  legislators.csv   ideal-point and rank summaries; placed members have
                    placed=True and are ranked among the fitted members
  rollcalls.csv     psi / zeta / kappa summaries
  summary.json      circularity chi_0 / chi_D / chi_R, WAIC, acceptance, timing
                    (fitted members only)
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from circfactor.analysis import bill_summary, chamber_summary, legislator_summary, orient, placed_summary
from circfactor.data import load_congress
from circfactor.model import Posterior, Priors, SamplerConfig, concatenate, sample
from circfactor.placement import place_legislators

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
    ap.add_argument("--bill-flip", action="store_true", help="enable the per-roll-call Yea/Nay mirror move (kappa rescaled)")
    ap.add_argument("--bill-flip-fixed-kappa", action="store_true", help="with --bill-flip: keep kappa fixed in the mirror move")
    ap.add_argument("--grad-dtype", default="float32", choices=["float32", "float64"])
    ap.add_argument("--omega-mean", type=float, default=10.0, help="mean of the exponential prior on omega")
    ap.add_argument("--xi-rate", type=float, default=25.0, help="rate of the exponential prior on xi = 1/lambda")
    ap.add_argument("--bill-draws", type=int, default=1000, help="evenly spaced psi/zeta/kappa draws to keep (needed for placement)")
    ap.add_argument("--place-min-votes", type=int, default=25,
                    help="place excluded members who cast at least this many votes (0: place nobody)")
    ap.add_argument("--resume", type=Path, help="results directory of an earlier run to continue (same congresses/sessions)")
    ap.add_argument("--discard-old", action="store_true", help="with --resume: treat the earlier draws as burn-in")
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
        party = data.legislators["party"].to_numpy()
        if args.resume:
            prev = Posterior.load(args.resume / tag / "posterior.npz")
            if prev.beta.shape[1] != data.shape[0]:
                raise ValueError(f"{tag}: the data changed since the run being resumed (I {prev.beta.shape[1]} vs {data.shape[0]})")
            # Same sampler settings and priors as the earlier run; only the length changes.
            cfg = SamplerConfig(**prev.config | dict(burnin=0, n_samples=args.samples, n_bill_draws=args.bill_draws, print_every=args.print_every))
            priors = Priors(**prev.priors)
            state = dict(prev.final_state)
            if args.discard_old:
                state = {k: v for k, v in state.items() if not k.startswith("waic_")}
            post = orient(concatenate(prev, sample(data.votes, priors, cfg, resume=state), drop_first=args.discard_old), party)
        else:
            cfg = SamplerConfig(
                n_samples=args.samples,
                burnin=args.burnin,
                thin=args.thin,
                n_bill_draws=args.bill_draws,
                beta_step=tuple(args.beta_step),
                tau_step=tuple(args.tau_step),
                adapt_steps=not args.no_adapt,
                grad_dtype=args.grad_dtype,
                rotation_move=args.rotation,
                bill_flip=args.bill_flip,
                bill_flip_kappa=not args.bill_flip_fixed_kappa,
                seed=args.seed,
                print_every=args.print_every,
            )
            post = orient(sample(data.votes, priors, cfg, party=party), party)

        leg = legislator_summary(post, data.legislators).assign(placed=False)
        leg.insert(leg.columns.get_loc("frac_missing"), "n_votes", (~np.isnan(data.votes)).sum(1))
        n_dropped = (~np.isnan(data.dropped_votes)).sum(1)
        sel = n_dropped >= args.place_min_votes if args.place_min_votes > 0 else np.zeros(len(n_dropped), bool)
        if sel.any() and post.yes.shape[0] > 0:
            beta_pair, omega_pair = post.bill_draw_pairs()
            placed = place_legislators(post.yes, post.no, post.kappa, omega_pair, data.dropped_votes[sel], mu=priors.mu, seed=cfg.seed)
            dropped = data.dropped[sel].assign(n_votes=n_dropped[sel])
            leg = pd.concat([leg, placed_summary(placed, beta_pair, dropped).assign(placed=True)], ignore_index=True)
            np.savez_compressed(out / "placed.npz", beta=placed, person=dropped["person"].to_numpy().astype(str))
            print(f"{tag}: placed {int(sel.sum())} excluded members with >= {args.place_min_votes} votes "
                  f"(skipped {int((~sel).sum())})", flush=True)
        leg.to_csv(out / "legislators.csv", index=False)
        bill_summary(post, data.rollcalls).to_csv(out / "rollcalls.csv", index=False)
        summary = chamber_summary(post, party) | {"congress": congress, "chamber": args.chamber, "I": int(data.shape[0]), "J": int(data.shape[1]),
                                                  "n_placed": int(sel.sum()), "n_skipped": int((~sel).sum())}
        (out / "summary.json").write_text(json.dumps(summary, indent=2))
        post.save(out / "posterior.npz")
        print(f"{tag}: chi_0 = {summary['chi_0']['mean']:.3f}, grouped elpd_waic = {summary['waic_grouped']['elpd_waic']:.0f}, "
              f"{post.seconds / 60:.1f} min -> {out}", flush=True)


if __name__ == "__main__":
    main()
