"""Build each House's reported posterior from the chains chosen in final_choices.json.

Examples:
  uv run python scripts/combine_chains.py            # every House in the choices file
  uv run python scripts/combine_chains.py 102 119    # selected Houses

Outputs go to results/final/H{congress}/:
  legislators.csv   ideal-point and rank summaries over the pooled draws of the
                    chosen chains; placed members (cut posterior) have placed=True
  rollcalls.csv     psi / zeta / kappa summaries
  skipped.csv       members excluded by the missing-vote rule with too few votes to place
  summary.json      chi_0 / chi_D / chi_R, grouped WAIC, convergence diagnostics,
                    the choice and its reason, and provenance (git commit, md5 of
                    every input file, each chain's sampler config)
and results/final/index.csv has one row per House.
"""

import argparse
import hashlib
import json
import subprocess
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from circfactor.analysis import bill_summary, circular_variance, ess_bulk, legislator_summary, placed_summary, split_rhat
from circfactor.data import DEFAULT_DATA_DIR, load_congress
from circfactor.model import Posterior

ROOT = Path(__file__).resolve().parents[1]


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def git_commit():
    run = lambda *a: subprocess.run(["git", *a], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    return {"commit": run("rev-parse", "HEAD"), "dirty": bool(run("status", "--porcelain", "--untracked-files=no"))}


def load_chain(path, drop, I):
    """One chosen chain with its first `drop` draws removed, plus its placed draws on the kept bill draws."""
    post = Posterior.load(path / "posterior.npz")
    if post.beta.shape[1] != I:
        raise ValueError(f"{path}: {post.beta.shape[1]} legislators, data has {I} (data changed since the run?)")
    idx = post.bill_draw_index if post.bill_draw_index is not None else np.arange(post.yes.shape[0])
    keep = idx >= drop
    placed = None
    if (path / "placed.npz").exists():
        z = np.load(path / "placed.npz")
        assert z["beta"].shape[0] == len(idx), f"{path}: placed draws not paired with bill draws"
        placed = {"beta": z["beta"][keep], "person": z["person"]}
    trimmed = replace(post, beta=post.beta[drop:], omega=post.omega[drop:], xi=post.xi[drop:], ll_legislator=post.ll_legislator[drop:],
                      yes=post.yes[keep], no=post.no[keep], kappa=post.kappa[keep], bill_draw_index=idx[keep] - drop)
    return trimmed, placed, post.loglik_trace[-post.beta.shape[0]:][drop:]


# Party used for the within-party measures chi_D / chi_R. Members carry the party they cast most
# of their votes with in that Congress (Voteview's label, also used in the fit). Members who left
# for Independent count with the party they left, or with the party they joined if they joined the
# other party within the same Congress. Sanders (I, never a party member) stays out of both.
PARTY_MEASURE_OVERRIDES = {107: {"G000280": "R"},   # Goode: D -> I (2000) -> R (Aug 2002)
                           116: {"A000367": "R"}}   # Amash: R -> I (July 2019)


def combine(congress, choice, out_dir, place_min_votes=25):
    data = load_congress(congress, "House", download=False)
    I = data.shape[0]
    over = PARTY_MEASURE_OVERRIDES.get(congress, {})
    party = np.array([over.get(p, x) for p, x in zip(data.legislators["person"].astype(str), data.legislators["party"])])
    drop = choice.get("drop_draws", 0)
    chains = [load_chain(ROOT / c, drop, I) for c in choice["chains"]]
    posts = [p for p, _, _ in chains]
    cat = lambda f: np.concatenate([f(p) for p in posts])
    pooled = replace(posts[0], beta=cat(lambda p: p.beta), omega=cat(lambda p: p.omega), xi=cat(lambda p: p.xi),
                     ll_legislator=cat(lambda p: p.ll_legislator), yes=cat(lambda p: p.yes), no=cat(lambda p: p.no), kappa=cat(lambda p: p.kappa),
                     bill_draw_index=None)

    leg = legislator_summary(pooled, data.legislators).assign(placed=False)
    leg.insert(leg.columns.get_loc("frac_missing"), "n_votes", (~np.isnan(data.votes)).sum(1))
    n_dropped = (~np.isnan(data.dropped_votes)).sum(1) if data.dropped_votes is not None else np.zeros(0, int)
    sel = n_dropped >= place_min_votes
    if sel.any():
        dropped = data.dropped[sel].assign(n_votes=n_dropped[sel])
        for _, pl, _ in chains:
            assert pl is not None and list(pl["person"]) == list(dropped["person"].astype(str)), "placed members differ from the data"
        placed = np.concatenate([pl["beta"] for _, pl, _ in chains])
        beta_pair = np.concatenate([p.beta[p.bill_draw_index] for p in posts])
        leg = pd.concat([leg, placed_summary(placed, beta_pair, dropped).assign(placed=True)], ignore_index=True)
    skipped = data.dropped[~sel].assign(n_votes=n_dropped[~sel])

    chi = {name: [circular_variance(p.beta[:, s], axis=1) for p in posts] for name, s in [("chi_0", slice(None)), ("chi_D", party == "D"), ("chi_R", party == "R")]}
    n = min(len(x) for x in chi["chi_0"])
    stack = np.stack([x[-n:] for x in chi["chi_0"]])
    rhat = split_rhat(stack) if len(posts) > 1 else split_rhat(stack[0][: n // 2 * 2].reshape(2, -1))
    summary = {"congress": congress, "chamber": "House", "I": I, "J": int(data.shape[1]),
               "n_fitted": I, "n_placed": int(sel.sum()), "n_skipped": int((~sel).sum())}
    for name, xs in chi.items():
        a = np.concatenate(xs)
        summary[name] = {"mean": float(a.mean()), "q025": float(np.quantile(a, 0.025)), "q975": float(np.quantile(a, 0.975)),
                         "ess": float(sum(ess_bulk(x) for x in xs)), "per_chain_mean": [float(x.mean()) for x in xs]}
    summary["chi_0"]["split_rhat"] = float(rhat)
    summary["party_measure_overrides"] = over
    summary["waic_grouped"] = {k: float(v) for k, v in pooled.waic_grouped().items()}
    summary["omega_mean"] = float(pooled.omega.mean())
    summary["choice"] = choice
    summary["chains"] = [{"path": c, "draws_used": int(p.beta.shape[0]), "loglik_mean": float(t.mean()), "loglik_sd": float(t.std()),
                          "md5_posterior": md5(ROOT / c / "posterior.npz"),
                          "md5_placed": md5(ROOT / c / "placed.npz") if (ROOT / c / "placed.npz").exists() else None,
                          "config": p.config, "priors": p.priors}
                         for c, (p, _, t) in zip(choice["chains"], chains)]
    files = sorted(Path(DEFAULT_DATA_DIR).glob(f"H{congress}_*.csv"))
    summary["provenance"] = {"code": git_commit(), "script": "scripts/combine_chains.py", "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                             "place_min_votes": place_min_votes,
                             "data": [{"file": f.name, "md5": md5(f), "modified": datetime.fromtimestamp(f.stat().st_mtime, timezone.utc).isoformat(timespec="seconds")} for f in files]}

    out = out_dir / f"H{congress}"
    out.mkdir(parents=True, exist_ok=True)
    leg.to_csv(out / "legislators.csv", index=False)
    bill_summary(pooled, data.rollcalls).to_csv(out / "rollcalls.csv", index=False)
    skipped.to_csv(out / "skipped.csv", index=False)
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    return summary


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("congresses", nargs="*", type=int)
    ap.add_argument("--choices", type=Path, default=ROOT / "scripts" / "final_choices.json")
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "final")
    args = ap.parse_args()
    houses = json.loads(args.choices.read_text())["houses"]
    rows = []
    for c in args.congresses or sorted(map(int, houses)):
        s = combine(c, houses[str(c)], args.out)
        rows.append({"congress": c, "chains": len(s["chains"]), "draws": sum(ch["draws_used"] for ch in s["chains"]),
                     **{f"{k}_{q}": s[k][q] for k in ("chi_0", "chi_D", "chi_R") for q in ("mean", "q025", "q975")},
                     "chi_0_ess": s["chi_0"]["ess"], "chi_0_split_rhat": s["chi_0"]["split_rhat"],
                     "n_fitted": s["n_fitted"], "n_placed": s["n_placed"], "n_skipped": s["n_skipped"], "rule": s["choice"]["rule"]})
        print(f"H{c}: chi_0 {s['chi_0']['mean']:.3f} [{s['chi_0']['q025']:.3f}, {s['chi_0']['q975']:.3f}] ess {s['chi_0']['ess']:.0f} "
              f"rhat {s['chi_0']['split_rhat']:.3f} | {rows[-1]['draws']} draws from {len(s['chains'])} chain(s) | "
              f"fitted {s['n_fitted']} placed {s['n_placed']} skipped {s['n_skipped']}", flush=True)
    index = args.out / "index.csv"
    old = pd.read_csv(index) if index.exists() and args.congresses else pd.DataFrame()
    new = pd.DataFrame(rows)
    pd.concat([old[~old.congress.isin(new.congress)] if len(old) else old, new]).sort_values("congress").to_csv(index, index=False)


if __name__ == "__main__":
    main()
