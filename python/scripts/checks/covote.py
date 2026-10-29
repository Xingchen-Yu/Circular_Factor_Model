"""Voting profile of one member: party unity, defections, who defects with them.
uv run python covote.py 102 Schroeder"""
import sys, numpy as np, pandas as pd
from circfactor.data import load_congress
c, name = int(sys.argv[1]), sys.argv[2]
d = load_congress(c, "House"); V = d.votes; leg = d.legislators.reset_index(drop=True); rc = d.rollcalls.reset_index(drop=True)
rk = {s: pd.read_csv(f"results/full/c{s}/H{c}/legislators.csv").set_index("icpsr").rank_median for s in (1, 2)}
leg["r1"] = leg.icpsr.astype(int).map(rk[1]); leg["r2"] = leg.icpsr.astype(int).map(rk[2])
i = leg.index[leg.label.str.startswith(name)][0]; me = V[i]
D = (leg.party == "D").values; R = (leg.party == "R").values
dy = np.nanmean(V[D], 0); ry = np.nanmean(V[R], 0)
party_vote = ((dy > .5) != (ry > .5))
dpos = (dy > .5).astype(float)
voted = ~np.isnan(me)
print(f"{leg.label[i]}: ranks c1 {leg.r1[i]:.0f}, c2 {leg.r2[i]:.0f}; voted {voted.sum()}/{len(me)}")
pv = party_vote & voted
unity = (me[pv] == dpos[pv]).mean()
# party unity of every Democrat for comparison
du = np.array([(V[k][party_vote & ~np.isnan(V[k])] == dpos[party_vote & ~np.isnan(V[k])]).mean() for k in np.where(D)[0]])
print(f"party-line votes: {pv.sum()}; votes with D majority {unity:.1%}  (Dem median {np.median(du):.1%}, rank {(du > unity).sum()+1} of {D.sum()} in loyalty)")
dis = voted & (me != dpos)  # any vote against the D majority
print(f"votes against the D majority (all roll calls): {dis.sum()}")
# on those, the side she joins: how liberal are the Democrats and how conservative the Republicans who join her?
for lab, cols in [("all votes against D majority", dis), ("party-line defections", dis & party_vote)]:
    if not cols.any(): continue
    same = (V[:, cols] == me[cols]); n = (~np.isnan(V[:, cols])).sum(1)
    agree = same.sum(1) / np.maximum(n, 1)
    leg[f"ag"] = agree
    print(f"\n== {lab} ({cols.sum()} roll calls): who votes with her")
    print(f"   Dems joining per vote: median {np.median((V[D][:, cols] == me[cols]).sum(0)):.0f} of {D.sum()};  Reps joining: median {np.median((V[R][:, cols] == me[cols]).sum(0)):.0f} of {R.sum()}")
    for p in "DR":
        t = leg[(leg.party == p) & (leg.index != i)].sort_values("ag", ascending=False).head(10)
        print(f"   top {p} co-voters:  " + "; ".join(f"{r.label.split(' (')[0].split(',')[0]} {r.ag:.0%} (c1 {r.r1:.0f}/c2 {r.r2:.0f})" for r in t.itertuples()))
    # correlation of agreement with rank among Republicans: far right or moderate Rs?
    rr = leg[R]
    print(f"   Republicans: corr(agreement, c1 rank) = {np.corrcoef(rr.ag, rr.r1)[0,1]:+.2f}  (+ = far-right Rs join her, - = moderate Rs)")
    dd = leg[D & (leg.index != i)]
    print(f"   Democrats:   corr(agreement, c1 rank) = {np.corrcoef(dd.ag, dd.r1)[0,1]:+.2f}  (- = liberal Ds join her)")
    if lab.startswith("party"):
        q = rc[cols].assign(me=me[cols], dems=(V[D][:, cols] == me[cols]).sum(0), reps=(V[R][:, cols] == me[cols]).sum(0))
        print("   examples:")
        for r in q.sample(min(12, len(q)), random_state=0).itertuples():
            print(f"     {'YEA' if r.me else 'NAY'} with {r.dems} D + {r.reps} R | {str(r.bill_number):10s} {str(r.vote_question)[:30]:30s} | {str(r.vote_desc)[:90]}")
