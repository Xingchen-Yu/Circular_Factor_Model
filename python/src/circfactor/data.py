"""Download and preprocess roll-call data from voteview.com.

Mirrors the preprocessing of the original R code (``get_rollcall_data``):

* cast codes 1-3 are "Yea" (1), 4-6 are "Nay" (0), everything else missing;
* a legislator with several ICPSR ids within one Congress (party switchers such
  as Justin Amash in the 116th) is merged into a single row; here they are
  matched on bioguide id rather than on the formatted name;
* legislators missing more than ``max_missing`` (default 40%) of the votes are
  dropped.
"""

from __future__ import annotations

import urllib.request
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

VOTEVIEW = "https://voteview.com/static/data/out"
DEFAULT_DATA_DIR = Path(__file__).resolve().parents[2] / "data" / "voteview"
PARTY = {100: "D", 200: "R", 328: "I"}
# Non-voting delegates / resident commissioner (they vote only in the Committee
# of the Whole); the original .ord files used for the paper exclude them.
DELEGATE_STATES = {"DC", "PR", "GU", "VI", "AS", "MP"}


@dataclass
class RollCallData:
    """Votes as an (I, J) float array: 1 = Yea, 0 = Nay, NaN = missing."""

    votes: np.ndarray
    legislators: pd.DataFrame  # one row per row of ``votes``
    rollcalls: pd.DataFrame  # one row per column of ``votes``
    chamber: str
    congress: int
    dropped: pd.DataFrame  # legislators removed for excessive absence
    dropped_votes: np.ndarray | None = None  # their votes, one row per row of ``dropped``

    @property
    def shape(self):
        return self.votes.shape

    def summary(self) -> str:
        I, J = self.shape
        miss = np.isnan(self.votes).mean()
        parties = self.legislators["party"].value_counts().to_dict()
        return f"{self.chamber[0]}{self.congress}: I={I} legislators, J={J} roll calls, {miss:.1%} missing, parties={parties}"


def _files(chamber: str, congress: int) -> dict[str, str]:
    c = f"{chamber[0].upper()}{congress:03d}" if congress < 100 else f"{chamber[0].upper()}{congress}"
    return {kind: f"{c}_{kind}.csv" for kind in ("votes", "members", "rollcalls")}


def download_congress(congress: int, chamber: str = "House", data_dir: Path = DEFAULT_DATA_DIR, force: bool = False) -> dict[str, Path]:
    """Fetch the per-congress votes/members/rollcalls CSVs from voteview."""
    data_dir = Path(data_dir)
    data_dir.mkdir(parents=True, exist_ok=True)
    paths = {}
    for kind, name in _files(chamber, congress).items():
        path = data_dir / name
        if force or not path.exists():
            tmp = path.with_suffix(".part")
            urllib.request.urlretrieve(f"{VOTEVIEW}/{kind}/{name}", tmp)
            tmp.rename(path)
        paths[kind] = path
    return paths


def _format_name(bioname: str) -> str:
    """'OCASIO-CORTEZ, Alexandria' -> 'Ocasio-cortez, Alexandria' (as in the R code)."""
    last, _, rest = str(bioname).partition(",")
    last = last[:1] + last[1:].lower()
    return f"{last},{rest}" if rest else last


def load_congress(
    congress: int,
    chamber: str = "House",
    data_dir: Path = DEFAULT_DATA_DIR,
    max_missing: float = 0.4,
    sessions: tuple[int, ...] | None = None,
    download: bool = True,
    include_delegates: bool = False,
) -> RollCallData:
    """Build the legislator x roll-call matrix for one Congress.

    Args:
      sessions: restrict to these sessions (e.g. ``(1,)`` reproduces the
        paper's first-session-only analysis of the 116th House).
    """
    paths = download_congress(congress, chamber, data_dir) if download else {k: Path(data_dir) / v for k, v in _files(chamber, congress).items()}
    votes = pd.read_csv(paths["votes"], usecols=["chamber", "rollnumber", "icpsr", "cast_code"])
    members = pd.read_csv(paths["members"])
    rollcalls = pd.read_csv(paths["rollcalls"], low_memory=False)

    members = members[members["chamber"] == chamber]
    if not include_delegates:
        members = members[~members["state_abbrev"].isin(DELEGATE_STATES)]
    # The votes file also records the President's position; keep members only.
    votes = votes[(votes["chamber"] == chamber) & votes["icpsr"].isin(members["icpsr"])]
    rollcalls = rollcalls[rollcalls["chamber"] == chamber].sort_values("rollnumber")
    if sessions is not None:
        rollcalls = rollcalls[rollcalls["session"].isin(sessions)]
    votes = votes[votes["rollnumber"].isin(rollcalls["rollnumber"])]

    code = votes["cast_code"]
    votes = votes.assign(y=np.where(code.between(1, 3), 1.0, np.where(code.between(4, 6), 0.0, np.nan)))
    wide = votes.pivot_table(index="icpsr", columns="rollnumber", values="y", aggfunc="first")
    wide = wide.reindex(columns=rollcalls["rollnumber"].to_numpy())

    # One person per row: merge multiple ICPSR ids of the same bioguide id.
    members = members.drop_duplicates("icpsr").set_index("icpsr").loc[wide.index]
    person = members["bioguide_id"].fillna(members["bioname"]).rename("person")
    n_cast = wide.notna().sum(axis=1)
    merged = wide.groupby(person).first()  # first non-missing vote per roll call
    clash = (wide.notna().groupby(person).sum() > 1).to_numpy().sum()
    if clash:
        raise ValueError(f"{clash} roll calls with conflicting votes for merged legislators")

    # Metadata from the ICPSR id with the most cast votes for each person.
    main = members.assign(n_cast=n_cast).sort_values("n_cast", ascending=False)
    main = main.assign(person=person).drop_duplicates("person").set_index("person").loc[merged.index]
    icpsrs = person.reset_index().groupby("person")["icpsr"].apply(lambda s: ",".join(map(str, sorted(s))))
    legislators = pd.DataFrame(
        {
            "person": merged.index,
            "icpsr": icpsrs.loc[merged.index].to_numpy(),
            "bioname": main["bioname"].to_numpy(),
            "party": main["party_code"].map(PARTY).fillna("O").to_numpy(),
            "state": main["state_abbrev"].to_numpy(),
            "district": main["district_code"].astype("Int64").to_numpy(),
            "nominate_dim1": main["nominate_dim1"].to_numpy(),
            "nominate_dim2": main["nominate_dim2"].to_numpy(),
        }
    )
    legislators["label"] = [
        f"{_format_name(b)} ({p} {s}-{d})" for b, p, s, d in zip(legislators["bioname"], legislators["party"], legislators["state"], legislators["district"])
    ]
    party_switch = person.reset_index().groupby("person")["icpsr"].nunique() > 1
    legislators["merged_ids"] = party_switch.loc[merged.index].to_numpy()

    y = merged.to_numpy(dtype=float)
    frac_missing = np.isnan(y).mean(axis=1)
    keep = frac_missing <= max_missing
    dropped = legislators[~keep].assign(frac_missing=frac_missing[~keep]).reset_index(drop=True)
    dropped_votes = y[~keep]
    legislators = legislators[keep].assign(frac_missing=frac_missing[keep]).reset_index(drop=True)
    y = y[keep]

    # Order rows like voteview (by state, then district) for readability.
    order = np.lexsort((legislators["district"].fillna(0).to_numpy(), legislators["state"].to_numpy()))
    legislators = legislators.iloc[order].reset_index(drop=True)
    y = y[order]

    rc_cols = [c for c in ["rollnumber", "date", "session", "bill_number", "vote_result", "vote_question", "vote_desc", "yea_count", "nay_count"] if c in rollcalls]
    return RollCallData(
        votes=y,
        legislators=legislators,
        rollcalls=rollcalls[rc_cols].reset_index(drop=True),
        chamber=chamber,
        congress=congress,
        dropped=dropped,
        dropped_votes=dropped_votes,
    )
