"""Source-result checks per House on the chains listed in scripts/final_choices.json.
uv run python source_checks.py [congress ...]   -> results/checks/source_checks.csv (+ printed table)
Orientation, membership accounting, party labels, classification vs NOMINATE, rank agreement with NOMINATE dim1."""
import sys, os, numpy as np, pandas as pd, jax.numpy as jnp
from scipy.stats import spearmanr
from circfactor.data import load_congress
from circfactor.model import Posterior, Data, loglik_matrix

import json
CHOICES = json.load(open("scripts/final_choices.json"))["houses"]
def final_dirs(c):
    return CHOICES[str(c)]["chains"]

rows = []
for c in [int(a) for a in sys.argv[1:]] or range(100, 120):
    dirs = final_dirs(c)
    d = load_congress(c, "House"); V = d.votes; mask = ~np.isnan(V)
    data = Data(jnp.asarray(np.where(mask, 2 * np.nan_to_num(V) - 1, 0.0)), jnp.asarray(mask))
    # classification: posterior-mean probability of the observed vote, over thinned bill draws of every final chain
    p_obs = np.zeros(V.shape); n = 0
    for r in dirs:
        p = Posterior.load(f"{r}/posterior.npz"); beta, _ = p.bill_draw_pairs()
        for k in np.linspace(0, len(beta) - 1, 100).astype(int):
            p_obs += np.exp(np.asarray(loglik_matrix(jnp.asarray(beta[k]), jnp.asarray(p.yes[k]), jnp.asarray(p.no[k]), jnp.asarray(p.kappa[k]), data))); n += 1
    p_obs /= n
    acc = (p_obs[mask] > 0.5).mean()
    minority = np.minimum(np.nansum(V, 0), mask.sum(0) - np.nansum(V, 0)).sum()
    apre = (minority - (p_obs[mask] <= 0.5).sum()) / minority
    # NOMINATE classification on the same votes (Voteview prob = NOMINATE prob of the cast vote, in %)
    vv = pd.read_csv(f"data/voteview/H{c}_votes.csv", usecols=["rollnumber", "icpsr", "prob"])
    nom = vv.dropna(subset=["prob"]).merge(d.rollcalls[["rollnumber"]], on="rollnumber")
    nom = nom[nom.icpsr.astype(str).isin(d.legislators.icpsr.astype(str))]
    acc_nom = (nom.prob > 50).mean() if len(nom) else float("nan")
    leg = pd.read_csv(f"{dirs[0]}/legislators.csv")
    fit = leg[~leg.placed]
    # orientation: Democrats on the low-rank side
    dmed, rmed = fit[fit.party == "D"].rank_median.median(), fit[fit.party == "R"].rank_median.median()
    # membership accounting vs the Voteview roster (non-delegate House members)
    mem = pd.read_csv(f"data/voteview/H{c}_members.csv"); mem = mem[(mem.chamber == "House") & ~mem.state_abbrev.isin(["DC", "PR", "GU", "VI", "AS", "MP"])]
    n_roster = mem.icpsr.nunique()
    n_fit, n_placed = int((~leg.placed).sum()), int(leg.placed.sum())
    n_dropped = len(d.dropped); n_skipped = n_dropped - n_placed
    merged = int(fit.merged_ids.sum()) if "merged_ids" in fit else 0
    dup = int(leg.icpsr.duplicated().sum())
    parties = ",".join(f"{k}:{v}" for k, v in leg.party.value_counts().items())
    rho = spearmanr(fit.rank_median, fit.nominate_dim1, nan_policy="omit")[0]
    # per-member fit: members whose votes the model predicts worst
    mem_acc = np.array([(p_obs[i][mask[i]] > 0.5).mean() for i in range(V.shape[0])])
    worst = d.legislators.assign(acc=mem_acc).nsmallest(3, "acc")
    rows.append(dict(H=c, chains=len(dirs), acc=acc, apre=apre, acc_nominate=acc_nom, D_med_rank=dmed, R_med_rank=rmed, oriented=dmed < rmed,
                     roster=n_roster, fitted=n_fit, placed=n_placed, skipped=n_skipped, accounted=n_fit + n_dropped,
                     merged_ids=merged, dup_icpsr=dup, parties=parties, spearman_nom=rho,
                     worst=" / ".join(f"{r.label.split(' (')[0]} {r.acc:.0%}" for r in worst.itertuples())))
    r = rows[-1]
    print(f"H{c} acc {acc:.3f} (NOMINATE {acc_nom:.3f}) APRE {apre:.3f} | oriented {r['oriented']} | roster {n_roster} fitted {n_fit} placed {n_placed} skipped {n_skipped} -> {n_fit+n_dropped} | merged {merged} dup {dup} | {parties} | rho_nom {rho:.3f} | worst: {r['worst']}", flush=True)
pd.DataFrame(rows).to_csv("results/checks/source_checks.csv", index=False)
