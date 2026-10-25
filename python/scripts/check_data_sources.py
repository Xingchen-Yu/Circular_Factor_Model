"""Check the original local roll-call data against current voteview.com files.

For each House, compares (legislator ICPSR id, roll call) cells of
  * the original Keith Poole .ord files (data/H{n}_votes.ord, used for the paper),
  * the April 2021 bulk dump (data/Hall_votes.csv),
with the current per-congress voteview CSVs, after mapping cast codes to
Yea (1-3) / Nay (4-6) / missing (else).  Also summarises the preprocessed
matrices (I, J, % missing) for every House, including ones without old data.

  uv run python scripts/check_data_sources.py 100-119 --old-dir ../R/data
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from circfactor.data import download_congress, load_congress

OLD_DIR = Path(__file__).resolve().parents[2] / "R" / "data"


def yn(code):
    code = np.asarray(code)
    return np.where((code >= 1) & (code <= 3), 1.0, np.where((code >= 4) & (code <= 6), 0.0, np.nan))


def read_ord(path):
    """{icpsr: vote-code array} from a KH .ord file (votes start at column 37)."""
    rows = {}
    for line in path.read_text().splitlines():
        icpsr = int(line[3:8])
        votes = np.array([int(c) for c in line[36:].replace(" ", "")])
        rows.setdefault(icpsr, []).append(votes)
    return rows


def current_long(paths):
    v = pd.read_csv(paths["votes"], usecols=["chamber", "rollnumber", "icpsr", "cast_code"])
    return v[v["chamber"] == "House"]


def compare(old: pd.DataFrame, new: pd.DataFrame) -> dict:
    """old/new: long frames with icpsr, rollnumber, y (1/0/NaN)."""
    m = old.merge(new, on=["icpsr", "rollnumber"], how="outer", suffixes=("_old", "_new"), indicator=True)
    both = m["_merge"] == "both"
    yo, yn_ = m.loc[both, "y_old"], m.loc[both, "y_new"]
    cast_both = yo.notna() & yn_.notna()
    return {
        "cells_both": int(both.sum()),
        "only_old": int((m["_merge"] == "left_only").sum()),
        "only_new": int((m["_merge"] == "right_only").sum()),
        "vote_flipped": int((yo[cast_both] != yn_[cast_both]).sum()),
        "missing_changed": int((yo.isna() != yn_.isna()).sum()),
        "rollcalls_old": int(old["rollnumber"].nunique()),
        "rollcalls_new": int(new["rollnumber"].nunique()),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("congresses", default="100-119", nargs="?")
    ap.add_argument("--old-dir", type=Path, default=OLD_DIR)
    args = ap.parse_args()
    lo, _, hi = args.congresses.partition("-")
    congresses = range(int(lo), int(hi or lo) + 1)

    hall = None
    if (args.old_dir / "Hall_votes.csv").exists():
        hall = pd.read_csv(args.old_dir / "Hall_votes.csv", usecols=["congress", "chamber", "rollnumber", "icpsr", "cast_code"])
        hall = hall[(hall["chamber"] == "House") & hall["congress"].between(min(congresses), max(congresses))]

    rows = []
    for c in congresses:
        paths = download_congress(c)
        new = current_long(paths).assign(y=lambda d: yn(d["cast_code"]))[["icpsr", "rollnumber", "y"]]
        rec = {"congress": c}
        ord_path = args.old_dir / f"H{c}_votes.ord"
        if ord_path.exists():
            ords = read_ord(ord_path)
            long = [
                (icpsr, j + 1, v)
                for icpsr, vs in ords.items()
                for row in vs
                for j, v in enumerate(yn(row))
                if row[j] != 0  # 0 = not in legislature for that vote
            ]
            old = pd.DataFrame(long, columns=["icpsr", "rollnumber", "y"]).drop_duplicates(["icpsr", "rollnumber"])
            rec |= {f"ord_{k}": v for k, v in compare(old, new).items()}
        if hall is not None and c in set(hall["congress"]):
            old = hall[hall["congress"] == c].assign(y=lambda d: yn(d["cast_code"]))[["icpsr", "rollnumber", "y"]]
            rec |= {f"hall_{k}": v for k, v in compare(old, new).items()}
        d = load_congress(c)
        rec |= {"I": d.shape[0], "J": d.shape[1], "pct_missing": round(100 * float(np.isnan(d.votes).mean()), 2),
                "dropped_absent": len(d.dropped), "merged_switchers": int(d.legislators["merged_ids"].sum()),
                "first_vote": str(d.rollcalls["date"].min()), "last_vote": str(d.rollcalls["date"].max())}
        print(rec, flush=True)
        rows.append(rec)
    out = Path(__file__).resolve().parents[1] / "results" / "data_source_check.csv"
    out.parent.mkdir(exist_ok=True)
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"-> {out}")


if __name__ == "__main__":
    main()
