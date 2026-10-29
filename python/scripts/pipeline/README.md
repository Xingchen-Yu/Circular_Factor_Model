# How the 100th-119th House results were produced (2026-10-27 to 10-29)

These are the scripts exactly as run (paths are those of the PC,
`~/circular_factor_model_remake/python`, and of the RunPod pods,
`/workspace/Circular_Factor_Model/python`). Every run used
`scripts/run_congress.py` with `--bill-flip --bill-shift` (jittering only, no
rotation; tempering and the legislator shift were tested and not used).
Each chain's full sampler config and seed is saved in its `posterior.npz` and
copied into `results/final/H*/summary.json`.

| Step | Script | Output |
|---|---|---|
| 1. Main runs: 2 chains per House, seeds 2021/2022, 20k burn-in + 40k draws, one House at a time | `run_lane.sh` (PC), `pod_lane.sh` (pods; `pod_setup.sh` sets a pod up) | `results/full/c{1,2}/H*` |
| 1b. H100 only: 2 more chains (seeds 2023/2024) | `run_congress.py` by hand | `results/full/c{3,4}/H100` |
| 2. Per House decision: best-vs-worst mean loglik gap > 3 within-chain sd -> extend the best chain +40k; fits tie but a member's 95% rank intervals don't overlap (and medians differ by > 50) -> 4 extra seeds; else nothing | `decide.py`, `decide_lane.sh` | `results/full_ext/c*/H*`, `results/brute/c{5..8}/H*` (seeds 2025-2028, 20k + 10k) |
| 3. Houses with extra seeds: extend the best of all chains +40k | `pick.py`, `pick_lane.sh` | `results/full_ext/c*/H{102,103,104,109}` |
| 4. Low chi_0 ESS: +80k for H102, H113, H115 (best chain), +40k per chain for H116, H118 | `run_congress.py --resume` by hand | `results/full_ext2/c*/H{102,113,115}`, `results/full_ext/c*/H{116,118}` |
| 5. H119: both chains +40k (ongoing Congress, chains pooled by user decision) | `run_congress.py --resume` by hand | `results/full_ext/c{1,2}/H119` |
| 6. Final posterior per House from the chosen chains | `../combine_chains.py` with `../final_choices.json` (choice + reason per House) | `results/final/H*`, `results/final/index.csv` |

`brute_lane.sh` was queued once and replaced by `decide_lane.sh` before it ran.

Checks (in `../checks/`):
* `source_checks.py`: orientation, membership accounting against the Voteview
  roster, classification vs NOMINATE, rank agreement with NOMINATE dim1, worst-fit
  members -> `results/checks/source_checks.csv`.
* `covote.py`, `tie_check.py`: voting-record check used to confirm the chosen
  chain where members' placements differ between chains (party-line defections,
  who joins them, agreement with rank neighbours under each chain).

`results/` and `logs/` are not in git; `summary.json` records the md5 of every
input posterior and data file so the outputs can be matched to them.
