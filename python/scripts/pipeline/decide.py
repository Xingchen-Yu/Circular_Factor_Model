"""Per-House decision after the main chains (results/full/c*/H{c}).
Prints one line: 'extend <k>' | 'none' | 'brute', then the evidence.
  extend: best chain's mean loglik beats the worst by > 3 within-chain sd -> continue the best chain.
  brute:  fits tie but some member's rank differs by > 50 (non-overlapping 95% intervals) -> more seeds needed.
  none:   fits tie and ranks agree.
"""
import sys, glob, numpy as np, pandas as pd
from circfactor.model import Posterior
c = int(sys.argv[1])
runs = sorted(glob.glob(f"results/full/c*/H{c}/summary.json"))
ll, sd, rk = {}, {}, {}
for r in runs:
    k = r.split("/")[2][1:]
    p = Posterior.load(r.replace("summary.json", "posterior.npz"))
    t = p.loglik_trace[-p.beta.shape[0]:]
    ll[k], sd[k] = t.mean(), t.std()
    rk[k] = pd.read_csv(r.replace("summary.json", "legislators.csv")).set_index("icpsr")[["rank_median", "rank_q025", "rank_q975"]]
best = max(ll, key=ll.get); worst = min(ll, key=ll.get)
gap, s = ll[best] - ll[worst], max(sd.values())
# a member counts as moved only if the chains' 95% rank intervals don't overlap (few-vote members have wide ones)
def moved(a, b):
    a, b = a.align(b, join="inner")[0], a.align(b, join="inner")[1]
    sep = (a.rank_q975 < b.rank_q025) | (b.rank_q975 < a.rank_q025)
    return (a.rank_median - b.rank_median).abs()[sep].max() if sep.any() else 0.0
moves = max(moved(rk[a], rk[b]) for a in rk for b in rk if a < b)
ev = f"H{c}: " + " ".join(f"c{k} {ll[k]:.0f}" for k in sorted(ll)) + f" | gap {gap:.0f} sd {s:.0f} | max rank move {moves:.0f}"
print("extend " + best if gap > 3 * s else "brute" if moves > 50 else "none")
print(ev)
