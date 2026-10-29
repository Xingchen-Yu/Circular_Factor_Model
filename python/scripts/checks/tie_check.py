"""Domain check for members whose placement differs between equally good chains.
For each member: agreement with the 20 members ranked around them in chain A vs chain B (neighbours by that chain's ranks),
and on the member's party-line defections, the mean rank of the other-party members who join them."""
import numpy as np, pandas as pd
from circfactor.data import load_congress
cases = [c for c in [(102, "brute/c6", "brute/c7", ["Jacobs, Andrew"]),
         (103, "brute/c8", "full/c2", ["Schroeder", "Archer"]),
         (104, "brute/c6", "full/c2", ["Neumann", "Roth, Toby"])] if c[0] == 104]
for c, A, B, names in cases:
    d = load_congress(c, "House"); V = d.votes; leg = d.legislators.reset_index(drop=True)
    rk = {k: leg.person.map(pd.read_csv(f"results/{k}/H{c}/legislators.csv").set_index("person").rank_median).values for k in (A, B)}
    D = (leg.party == "D").values; R = (leg.party == "R").values
    dpos = (np.nanmean(V[D], 0) > .5).astype(float); rpos = (np.nanmean(V[R], 0) > .5).astype(float); pv = dpos != rpos
    for n in names:
        i = leg.index[leg.label.str.startswith(n)][0]; me = V[i]; own = dpos if D[i] else rpos; other = R if D[i] else D
        both = ~np.isnan(V) & ~np.isnan(me); ag = np.where(both, V == me, False).sum(1) / both.sum(1); ag[i] = np.nan
        defect = pv & ~np.isnan(me) & (me != own)
        join = (V[:, defect] == me[defect])  # who votes with them on their defections
        print(f"H{c} {leg.label[i]}: votes with own party on {1 - defect.sum() / (pv & ~np.isnan(me)).sum():.0%} of party-line votes ({defect.sum()} defections)")
        for k in (A, B):
            nb = np.abs(rk[k] - rk[k][i]) <= 10; nb[i] = False
            oth = other & (join.sum(1) > 0)
            w = join[other].mean(1)  # share of their defections each other-party member joins
            corr = np.corrcoef(w, rk[k][other])[0, 1]
            print(f"   {k:9s} rank {rk[k][i]:4.0f}: agreement with rank neighbours {np.nanmean(ag[nb]):.1%}")
        top = np.argsort(-np.nan_to_num(ag, nan=-1))[:8]
        print("   most-agreeing members (rank A/B):", "; ".join(f"{leg.label[t].split(',')[0]} ({leg.party[t]}) {ag[t]:.0%} {rk[A][t]:.0f}/{rk[B][t]:.0f}" for t in top))
        w = join[other].mean(1); o = np.where(other)[0]; topo = o[np.argsort(-w)[:6]]
        print(f"   other-party members joining their defections most (rank A/B):", "; ".join(f"{leg.label[t].split(',')[0]} {join[t].mean():.0%} {rk[A][t]:.0f}/{rk[B][t]:.0f}" for t in topo),
              f"| corr(join share, rank) {np.corrcoef(w, rk[A][other])[0,1]:+.2f}")
