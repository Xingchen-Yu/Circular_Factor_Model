"""Export results/final/ to the JSON files the Full Circle website reads, and check them.

    uv run python scripts/export_website.py            # writes ../website/data/

Inputs: results/final/ (combine_chains.py), data/voteview/ (votes, roll calls),
data/congress_legislators/ (display names, public domain, github.com/unitedstates/congress-legislators).
Optional: scripts/featured_votes.json, hand-picked featured votes per House
({"116": [{"roll": 797, "title": "...", "text": "..."}]}); Houses without an entry get the
top-scoring "ends against the middle" votes.

Outputs (website/data/):
  index.json      one entry per House: years, majority, chi_0/chi_D/chi_R, counts, featured votes
  members.json    one entry per person: name and their record in every House (search, member pages)
  H{c}.json       members of House c (display angle, rank, interval) and its roll calls
  V{c}.json       how each member voted on each roll call (loaded when a vote is picked)
  csv/            the per-House result files, for download
Also results/website/featured_candidates.csv: the top-scoring votes per House, for review.
"""

import json
import re
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FINAL = ROOT / "results" / "final"
VOTEVIEW = ROOT / "data" / "voteview"
NAMES = ROOT / "data" / "congress_legislators"
OUT = ROOT.parent / "website" / "data"
CONGRESSES = range(100, 120)
# House majority party by Congress
MAJORITY = {c: "D" for c in [100, 101, 102, 103, 110, 111, 116, 117]} | {c: "R" for c in [104, 105, 106, 107, 108, 109, 112, 113, 114, 115, 118, 119]}
PARTY_CODE = {100: "D", 200: "R", 328: "I"}
PARTY_NAME = {"D": "Democrat", "R": "Republican", "I": "Independent"}
MONTHS = "January February March April May June July August September October November December".split()


def names():
    cl = pd.concat([pd.read_csv(NAMES / f"legislators-{k}.csv", dtype=str) for k in ("current", "historical")])
    cl = cl.drop_duplicates("bioguide_id").set_index("bioguide_id").fillna("")
    out = {}
    for b, r in cl.iterrows():
        n = f"{r.nickname or r.first_name} {r.last_name}"
        out[b] = n + (f" {r.suffix}" if r.suffix else "")
    return out


def ids(icpsr):
    return [int(float(i)) for i in re.split(r"[,; ]+", str(icpsr)) if i]


def circ_mean(a):
    return float(np.arctan2(np.sin(a).mean(), np.cos(a).mean()))


def wrap(a):
    return np.arctan2(np.sin(a), np.cos(a))


def switch_note(c, leg, members, votes, rollcalls):
    """Party switch note for a person merged from several ICPSR ids, None if the party did not change."""
    spans = []
    for i in ids(leg.icpsr):
        v = votes[(votes.icpsr == i) & votes.cast_code.between(1, 6)]
        d = rollcalls.set_index("rollnumber").loc[v.rollnumber, "date"]
        code = int(members.loc[members.icpsr == i, "party_code"].iloc[0])
        if len(d):
            spans.append((d.min(), PARTY_CODE.get(code, "I")))
    spans.sort()
    if len(spans) < 2 or spans[0][1] == spans[1][1]:
        return None
    y, m = spans[1][0][:4], int(spans[1][0][5:7])
    return f"Switched from {PARTY_NAME[spans[0][1]]} to {PARTY_NAME[spans[1][1]]} during this Congress (first vote as {PARTY_NAME[spans[1][1]]}: {MONTHS[m - 1]} {y})"


def vote_strings(leg, votes, rolls):
    """One string per roll call, one character per member: Y yes, N no, P present / not voting, _ not in office."""
    col = {}
    for k, r in enumerate(leg.itertuples()):
        for i in ids(r.icpsr):
            col[i] = k
    v = votes[votes.icpsr.isin(col) & votes.rollnumber.isin(rolls)]
    v = v.assign(k=v.icpsr.map(col), ch=np.select([v.cast_code.between(1, 3), v.cast_code.between(4, 6), v.cast_code.between(7, 9)], ["Y", "N", "P"], "_"))
    # merged switchers: keep the cast vote if one of their ids has one
    v = v.assign(pri=v.ch.map({"Y": 0, "N": 0, "P": 1, "_": 2})).sort_values("pri").drop_duplicates(["rollnumber", "k"])
    row = {r: j for j, r in enumerate(rolls)}
    grid = np.full((len(rolls), len(leg)), "_")
    grid[v.rollnumber.map(row).to_numpy(), v.k.to_numpy()] = v.ch.to_numpy()
    return ["".join(row) for row in grid]


def ends_scores(leg, vs):
    """How strongly each roll call puts both ends of the House against both parties' middles (0-1)."""
    n = len(leg)
    order = leg.rank_median.rank(method="first").to_numpy()
    q = max(1, round(n / 10))
    group = np.where(order <= q, "L", np.where(order > n - q, "R", np.where(leg.party == "D", "Dm", np.where(leg.party == "R", "Rm", ""))))
    masks = {g: group == g for g in ("L", "Dm", "Rm", "R")}
    out = []
    for s in vs:
        a = np.frombuffer(s.encode(), "S1")
        share = {}
        for g, m in masks.items():
            y, nn = (a[m] == b"Y").sum(), (a[m] == b"N").sum()
            share[g] = y / (y + nn) if y + nn >= 5 else np.nan
        if any(np.isnan(x) for x in share.values()):
            out.append(0.0)
            continue
        ends_yes = min(share["L"], share["R"], 1 - share["Dm"], 1 - share["Rm"])
        ends_no = min(1 - share["L"], 1 - share["R"], share["Dm"], share["Rm"])
        out.append(round(float(max(ends_yes, ends_no)), 3))
    return out


def describe(r):
    for k in ("vote_desc", "dtl_desc"):
        if isinstance(r[k], str) and r[k].strip():
            return r[k].strip()[:400]
    return ""


def main():
    nm = names()
    index = pd.read_csv(FINAL / "index.csv").set_index("congress")
    featured_manual = json.loads((ROOT / "scripts" / "featured_votes.json").read_text()) if (ROOT / "scripts" / "featured_votes.json").exists() else {}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "csv").mkdir(exist_ok=True)
    people, houses, candidates = {}, [], []
    for c in CONGRESSES:
        d = FINAL / f"H{c}"
        leg = pd.read_csv(d / "legislators.csv")
        skipped = pd.read_csv(d / "skipped.csv")
        summ = json.loads((d / "summary.json").read_text())
        rc = pd.read_csv(d / "rollcalls.csv")
        vv_rc = pd.read_csv(VOTEVIEW / f"H{c}_rollcalls.csv")
        rc = rc.merge(vv_rc[["rollnumber", "dtl_desc"]], on="rollnumber", how="left")
        members = pd.read_csv(VOTEVIEW / f"H{c}_members.csv")
        votes = pd.read_csv(VOTEVIEW / f"H{c}_votes.csv")
        n_rank = int(summ["n_fitted"])

        # checks on the inputs
        assert leg.person.is_unique, f"H{c}: duplicated person"
        assert len(leg) == summ["n_fitted"] + summ["n_placed"] and len(skipped) == summ["n_skipped"], f"H{c}: member counts differ from summary.json"
        assert set(leg.person).isdisjoint(skipped.person), f"H{c}: member both fitted and skipped"
        assert leg.person.isin(nm.keys()).all() and skipped.person.isin(nm.keys()).all(), f"H{c}: member without a display name"
        assert leg.rank_median.between(1, n_rank).all(), f"H{c}: rank outside 1..{n_rank}"
        dmed, rmed = leg.loc[leg.party == "D", "rank_median"].median(), leg.loc[leg.party == "R", "rank_median"].median()
        assert dmed < rmed, f"H{c}: rank direction wrong (Democrats should rank more liberal)"

        # display angle: the midpoint between the two parties' typical members at the top of the circle
        beta = leg.beta_median.to_numpy()
        mid = circ_mean(np.array([circ_mean(beta[leg.party == "D"]), circ_mean(beta[leg.party == "R"])]))
        ang = wrap(beta - mid)
        assert np.median(ang[leg.party == "D"]) < 0 < np.median(ang[leg.party == "R"]), f"H{c}: Democrats should sit on the left"

        rolls = rc.rollnumber.tolist()
        vs = vote_strings(leg, votes, rolls)
        known = {i for x in pd.concat([leg.icpsr, skipped.icpsr]) for i in ids(x)}
        # the President's positions and non-voting delegates (Committee of the Whole votes, 103rd) are not in the fit
        outside = set(members.loc[(members.chamber == "President") | members.state_abbrev.isin(["DC", "PR", "GU", "VI", "AS", "MP"]), "icpsr"])
        cast = set(votes.loc[votes.cast_code.between(1, 6) & votes.rollnumber.isin(rolls), "icpsr"]) - outside
        assert cast <= known, f"H{c}: votes from members missing in the results: {sorted(cast - known)}"
        for s, r in zip(vs, rc.itertuples()):  # Voteview's yea_count counts only plain yes votes (code 1), we also count paired/announced
            y = s.count("Y")
            assert r.yea_count - len(skipped) - 6 <= y <= r.yea_count + 25, f"H{c} roll {r.rollnumber}: {y} yes in export, {r.yea_count} in Voteview"
        score = ends_scores(leg, vs)
        rc["score"] = score

        notes = {}
        for r in leg[leg.merged_ids == True].itertuples():  # noqa: E712
            notes[r.person] = switch_note(c, r, members, votes, vv_rc)

        mrows = []
        for k, r in enumerate(leg.itertuples()):
            mrows.append([r.person, r.party, r.state, int(r.district), round(float(ang[k]), 4), int(r.rank_median), int(r.rank_q025), int(r.rank_q975), int(r.n_votes), int(bool(r.placed))])
            people.setdefault(r.person, {"id": r.person, "name": nm[r.person], "h": []})["h"].append(
                [c, r.party, r.state, int(r.district), int(r.rank_median), int(r.rank_q025), int(r.rank_q975), n_rank, int(r.n_votes), int(bool(r.placed)), notes.get(r.person)])
        for r in skipped.itertuples():
            people.setdefault(r.person, {"id": r.person, "name": nm[r.person], "h": []})["h"].append(
                [c, r.party, r.state, int(r.district), None, None, None, n_rank, int(r.n_votes), 0, None])

        rcrows = [[int(r.rollnumber), r.date, r.bill_number if isinstance(r.bill_number, str) else "", r.vote_question if isinstance(r.vote_question, str) else "",
                   describe(r._asdict()), r.vote_result if isinstance(r.vote_result, str) else "", int(r.yea_count), int(r.nay_count), r.score] for r in rc.itertuples()]
        (OUT / f"H{c}.json").write_text(json.dumps({"congress": c, "nRank": n_rank, "cols": ["id", "party", "state", "district", "angle", "rank", "lo", "hi", "votes", "placed"],
                                                     "members": mrows, "rollcols": ["roll", "date", "bill", "question", "desc", "result", "yea", "nay", "score"], "rolls": rcrows},
                                                    separators=(",", ":"), ensure_ascii=False))
        (OUT / f"V{c}.json").write_text(json.dumps({"rolls": rolls, "votes": vs}, separators=(",", ":")))

        # featured votes: hand-picked if present, else the best-scoring roll call per bill
        top = rc.sort_values("score", ascending=False)
        top = top[~top.bill_number.duplicated() | top.bill_number.isna()].head(8)
        for r in top.itertuples():
            candidates.append({"congress": c, "roll": r.rollnumber, "date": r.date, "bill": r.bill_number, "score": r.score, "yea": r.yea_count, "nay": r.nay_count,
                               "question": r.vote_question, "desc": describe(r._asdict())})
        if str(c) in featured_manual:
            featured = featured_manual[str(c)]
            assert all(f["roll"] in rolls for f in featured), f"H{c}: featured roll call not in the data"
        else:
            featured = [{"roll": int(r.rollnumber)} for r in top.itertuples() if r.score >= 0.4][:4]

        ix = index.loc[c]
        for k in ("chi_0", "chi_D", "chi_R"):
            assert abs(ix[f"{k}_mean"] - summ[k]["mean"]) < 1e-12, f"H{c}: index.csv and summary.json disagree on {k}"
        houses.append({"congress": c, "years": [1787 + 2 * c, 1789 + 2 * c], "majority": MAJORITY[c],
                       **{k: [round(summ[k]["mean"], 4), round(summ[k]["q025"], 4), round(summ[k]["q975"], 4)] for k in ("chi_0", "chi_D", "chi_R")},
                       "nFitted": summ["n_fitted"], "nPlaced": summ["n_placed"], "nSkipped": summ["n_skipped"],
                       "nD": int((leg.party == "D").sum()), "nR": int((leg.party == "R").sum()), "nI": int((leg.party == "I").sum()),
                       "nRolls": len(rolls), "first": rc.date.min(), "last": rc.date.max(), "featured": featured})
        for f in ("legislators.csv", "rollcalls.csv", "skipped.csv", "summary.json"):
            shutil.copy(d / f, OUT / "csv" / f"H{c}_{f}")
        print(f"H{c}: {len(leg)} members, {len(skipped)} skipped, {len(rolls)} roll calls, top ends score {max(score):.2f}", flush=True)

    # every person once, names checked
    for p in people.values():
        p["h"].sort(key=lambda h: h[0])
    for c in CONGRESSES:
        rows = [(p["name"], h[2]) for p in people.values() for h in p["h"] if h[0] == c]
        dup = pd.Series(rows).duplicated()
        assert not dup.any(), f"H{c}: two members share a name and state: {pd.Series(rows)[dup].tolist()}"
    (OUT / "members.json").write_text(json.dumps(sorted(people.values(), key=lambda p: p["name"]), separators=(",", ":"), ensure_ascii=False))
    (OUT / "index.json").write_text(json.dumps({"houses": houses, "lastVote": houses[-1]["last"], "source": "Voteview (voteview.com), downloaded 2026-10"}, separators=(",", ":")))
    cand = ROOT / "results" / "website"
    cand.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(candidates).to_csv(cand / "featured_candidates.csv", index=False)
    print(f"{len(people)} people; wrote {OUT}")


if __name__ == "__main__":
    main()
