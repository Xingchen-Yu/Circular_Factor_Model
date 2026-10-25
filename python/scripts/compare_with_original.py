"""Compare a new H116 (session 1) fit with the original R/C++ posterior.

Usage: uv run python scripts/compare_with_original.py results/H116_s1 /path/to/H116_posterior_beta_estimates.csv
"""

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from circfactor.model import Posterior


def key(label: str) -> str:
    """'Ocasio-cortez, Alexandria (D NY-14)' / 'Ocasio-cort (D NY-14)' -> 'ocasio-cort|NY-14'."""
    last = re.split(r"[,(]", label)[0].strip().lower()[:11]
    district = re.search(r"\(\w+ ([A-Z]{2}-\d+)\)", label).group(1)
    return f"{last}|{district}"


def main(new_dir, old_csv):
    new_dir = Path(new_dir)
    leg = pd.read_csv(new_dir / "legislators.csv")
    post = Posterior.load(new_dir / "posterior.npz")
    old = pd.read_csv(old_csv)
    old_keys = [key(c) for c in old.columns]
    new_keys = [key(c) for c in leg["label"]]
    common = [k for k in new_keys if k in set(old_keys)]
    print(f"matched {len(common)} of {len(new_keys)} legislators (old run has {len(old_keys)})")
    i_new = np.array([new_keys.index(k) for k in common])
    i_old = np.array([old_keys.index(k) for k in common])

    b_new, b_old = post.beta[:, i_new], old.to_numpy()[:, i_old]
    rank = lambda b: b.argsort(1).argsort(1) + 1
    r_new, r_old = rank(b_new), rank(b_old)
    m_new, m_old = np.median(b_new, 0), np.median(b_old, 0)
    print(f"posterior median beta: Pearson r = {np.corrcoef(m_new, m_old)[0, 1]:.4f}, max |diff| = {np.abs(m_new - m_old).max():.3f} rad")
    print(f"posterior median rank: Spearman = {stats.spearmanr(np.median(r_new, 0), np.median(r_old, 0)).statistic:.4f}, "
          f"mean |diff| = {np.abs(np.median(r_new, 0) - np.median(r_old, 0)).mean():.2f}")
    sd_new, sd_old = b_new.std(0), b_old.std(0)
    print(f"posterior sd of beta: new/old ratio median {np.median(sd_new / sd_old):.3f} (IQR {np.quantile(sd_new / sd_old, 0.25):.3f}-{np.quantile(sd_new / sd_old, 0.75):.3f})")
    print(f"circular variance chi_0: new {np.mean(1 - np.hypot(np.cos(b_new).mean(1), np.sin(b_new).mean(1))):.4f}, "
          f"old {np.mean(1 - np.hypot(np.cos(b_old).mean(1), np.sin(b_old).mean(1))):.4f}")
    print("\nPaper Tables 2-3 (median rank, 95% interval):")
    for name in ["pressley", "omar", "tlaib", "ocasio-cort", "amash", "massie", "gaetz"]:
        j = next(k for k, c in enumerate(common) if c.startswith(name + "|"))
        f = lambda r: f"{np.median(r[:, j]):.0f} ({np.quantile(r[:, j], 0.025):.0f},{np.quantile(r[:, j], 0.975):.0f})"
        print(f"  {common[j]:<18} new {f(r_new):<16} old {f(r_old)}")


if __name__ == "__main__":
    main(*sys.argv[1:3])
