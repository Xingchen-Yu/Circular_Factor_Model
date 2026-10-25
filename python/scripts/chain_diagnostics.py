"""Between-chain diagnostics for several runs of the same data set.

  uv run python scripts/chain_diagnostics.py results/chains/c*/H116_s1
"""

import sys

import numpy as np
import pandas as pd

from circfactor.analysis import ess_bulk, split_rhat
from circfactor.model import Posterior

dirs = sys.argv[1:]
posts = [Posterior.load(f"{d}/posterior.npz") for d in dirs]
leg = pd.read_csv(f"{dirs[0]}/legislators.csv")
ranks = np.stack([p.beta.argsort(1).argsort(1) + 1 for p in posts])  # (chains, S, I)
chi = np.stack([1 - np.hypot(np.cos(p.beta).mean(1), np.sin(p.beta).mean(1)) for p in posts])

rh = np.array([split_rhat(ranks[:, :, i]) for i in range(ranks.shape[2])])
ok = ~np.isnan(rh)
print(f"{len(posts)} chains x {ranks.shape[1]} draws")
print(f"rank split-Rhat: median {np.median(rh[ok]):.3f}, 95th pct {np.quantile(rh[ok], 0.95):.3f}, max {rh[ok].max():.3f}; "
      f"{(rh[ok] > 1.01).sum()} legislators > 1.01, {(rh[ok] > 1.05).sum()} > 1.05")
print(f"chi_0 split-Rhat {split_rhat(chi):.3f}; per-chain means {chi.mean(1).round(4)}; ESS per chain {[round(ess_bulk(c)) for c in chi]}")
print(f"omega split-Rhat {split_rhat(np.stack([p.omega for p in posts])):.3f}")
print("worst legislators:")
for i in np.argsort(-np.nan_to_num(rh))[:10]:
    print(f"  {leg.label[i]:<45} Rhat {rh[i]:.3f}  per-chain median rank {[int(np.median(r)) for r in ranks[:, :, i]]}")
for name in ["Gaetz", "Rooney", "Lofgren", "Defazio", "Pressley", "Omar", "Tlaib", "Ocasio", "Amash", "Massie"]:
    i = int(leg.index[leg.label.str.startswith(name)][0])
    print(f"  {leg.label[i]:<45} Rhat {rh[i]:.3f}  per-chain median rank {[int(np.median(r)) for r in ranks[:, :, i]]}")
