# circfactor: circular factor model in JAX

A Python/JAX port of the circular (spherical) latent factor model from

> Yu, X. & Rodríguez, A. (2021). *Spatial voting models in circular spaces: A
> case study of the U.S. House of Representatives.* Annals of Applied
> Statistics 15(4), 1897–1922.

It replaces the R + Rcpp/RcppArmadillo + snowfall implementation
(`SLFM()` in `../R/main_script/Circular_Factor_Model_tidyverse.R`). The sampler is
the same Metropolis-within-Gibbs scheme with geodesic HMC for the angles, but
each block is now a single vectorised GPU call, gradients come from autodiff,
and the numerics are tested against high-precision references.

## Install

All commands are run from this `python/` directory.

```bash
uv sync            # Python 3.12, JAX with CUDA 12
uv run pytest      # full validation suite (~3 min on a GPU)
```

On WSL, XLA's default 75% GPU preallocation fails, so the package sets
`XLA_PYTHON_CLIENT_PREALLOCATE=false` unless you set it yourself.

## Usage

```bash
# Paper's 116th House analysis (first session: I=432, J=700, as in Table 4)
uv run python scripts/run_congress.py 116 --sessions 1

# 100th through 119th Houses (data is downloaded from voteview.com)
uv run python scripts/run_congress.py 100-119

# Compare with the original R results
uv run python scripts/compare_with_original.py results/H116_s1 \
    path/to/H116_posterior_beta_estimates.csv   # from the original R run (not in the repo)

# Parameter recovery on simulated data
uv run python scripts/simulation_study.py --I 200 --J 600 --omega 3
```

From Python:

```python
from circfactor.data import load_congress
from circfactor.model import Priors, SamplerConfig, sample
from circfactor.analysis import orient, legislator_summary, chamber_summary

d = load_congress(116, sessions=(1,))
party = d.legislators["party"].to_numpy()
post = orient(sample(d.votes, Priors(), SamplerConfig(burnin=30_000, n_samples=20_000), party=party), party)
legislator_summary(post, d.legislators)   # ideal points + rank medians / 95% intervals
chamber_summary(post, party)              # chi_0, chi_D, chi_R, grouped WAIC
```

## Layout

| Module | Replaces | Contents |
|---|---|---|
| `special.py` | `betaInc_log*`, `betaContFrac`, Bessel code | log-CDF of the link G_κ, with an analytic `custom_jvp` |
| `ghmc.py` | the leapfrog loops in `update_beta` / `update_tau_*` | batched geodesic HMC on S^{d-1} (any d) |
| `model.py` | `SLFM()`, `update_kappa`, `update_omega`, `waic_cpp` | priors, likelihood, Gibbs sweep, sampler, WAIC |
| `data.py` | `get_rollcall_data`, `readKH2`, `ymat_spit` | voteview download and preprocessing |
| `analysis.py` | parts of `tables_figures.R` | ranks, circular variance, summaries |

## What the sampler does

Each iteration, in the original order:

1. κ_j: log-normal random-walk MH (proposal sd tuned toward 30–60% acceptance during burn-in, as in `update_tsig`)
2. ξ = 1/λ: conjugate Gamma draw
3. β_i: geodesic HMC, all legislators at once
4. ψ_j, then ζ_j: geodesic HMC, all roll calls at once
5. ω: log-normal random-walk MH

Step sizes and leapfrog counts are "jittered" every 50 iterations, with
the paper's ranges (β: U(0.01, 0.04); ψ, ζ: U(0.01, 0.105); L ~ U{1..10}).

## Validation

`uv run pytest` runs:

* **`test_ghmc_vonmises.py`: the GHMC kernel recovers von Mises
  distributions.** Thousands of independent chains target
  vonMises(μ, ω) for μ ∈ {0, 2, π−0.05, −π/2} and ω ∈ {0, 0.1, 1, 4, 16, 64,
  256}; final states pass KS tests against the exact CDF, plus checks on the
  mean resultant length and mean direction. Further tests cover the
  paper's jitter ranges, deliberately coarse steps (so the Metropolis
  correction matters), a deliberately *wrong* gradient (the chain must stay
  exact), von Mises–Fisher on S², and **posterior recovery**: for data
  x₁..xₙ ~ vonMises(μ, κ), the posterior of μ is exactly vonMises(x̄, κR),
  and GHMC must reproduce it and give ~90% coverage of the true μ over 40
  simulated datasets. Injecting bugs (always-accept, wrong kinetic-energy sign)
  makes these tests fail.
* **`test_special.py`**: log G_κ and its derivative against mpmath/scipy for
  κ ∈ [10⁻³, 10⁴] including tails down to 10⁻¹²: relative error ≤ 1e-10;
  convergence of the fixed-length continued fraction; `betaln` helpers.
* **`test_model.py`**: the likelihood matrix against a direct scipy evaluation of
  paper eq. (4)–(5); a **Geweke (2004) joint-distribution test** of the whole Gibbs sweep
  (2000 chains started from the prior, alternating sweeps with fresh
  y | θ, must stay at the prior); and a consistency check of the carried
  log-likelihood matrix. Injected bugs (a dropped Jacobian, an off-by-one in
  the ξ conjugate, stale likelihood bookkeeping) make these tests fail.

## Differences from the original code

* **Missing votes** are dropped from the likelihood (missing completely at
  random, as in the paper) instead of being imputed every 50 iterations. Both
  target the same posterior for the parameters. WAIC is computed over observed
  votes only, so it is not directly comparable with Table 4 if that included
  imputed entries.
* **Step-size adaptation** (on by default, burn-in only): each coordinate's
  jitter range is rescaled so its acceptance lands in [0.6, 0.9], the window
  the paper tuned toward by hand. `--no-adapt` reproduces fixed ranges.
* **Mixed precision**: leapfrog gradients are computed in float32, while the
  Hamiltonian/acceptance step uses float64. The leapfrog map is reversible and
  volume-preserving for any deterministic gradient, so the chain stays exact
  (the wrong-gradient test checks this). This makes it about 4x faster than
  all-float64. Use `--grad-dtype float64` to disable it.
* **κ_j ≤ 10⁴**: the fixed-length continued fraction is validated up to this
  value, so proposals above it are rejected. That truncates the κ prior where
  it has mass ~exp(−10⁴ξ), which is negligible for any plausible ξ.
* **Optional rotation move** (`--rotation`, off by default): an exact Gibbs
  draw of a common rotation of all angles. It is exact (Geweke test; identical
  stationary distribution on simulated data), but on the 116th House every
  chain run with it climbed more slowly in log-likelihood, so it is not used.
* Party switchers are merged by bioguide id rather than by formatted name.
* Non-voting delegates (DC, PR, GU, VI, AS, MP) are excluded, as in the
  original `.ord` files (`include_delegates=True` keeps them).

## Run length and mixing

The sampler mixes slowly on the House data, both here and in the original R
code; this is a property of the model and the data, not of the port.

* **Paper:** 20,000 draws kept after a burn-in of 20,000 to 80,000
  iterations depending on the House (median ~30,000). Convergence was judged
  from the log-likelihood trace and Gelman-Rubin across several chains.
* **Original R code:** `R/README.md` recommends 100,000 iterations with
  80,000 burn-in for the 112th and 116th Houses (1-3 days on 12 CPU cores).
* **This port, 116th House, first session** (one seed): the log-likelihood
  by 10k block is -29318, -28769, -28625, -28607, -28594, so it levels off
  only after ~30-40k iterations. Chains with different seeds can still
  disagree on individual legislators (e.g. Gaetz's median rank 400 vs 411
  after 30k + 20k).

The defaults (`--burnin 30000 --samples 20000`) match the paper's median.
For final results, run several chains per House with
`--burnin 80000` or more and check them with `scripts/chain_diagnostics.py`.
A single 116th House chain runs at about 22 it/s on an RTX 5070 Ti, so
100,000 iterations take about 75 minutes.

## Data source check (2026-10-24)

`scripts/check_data_sources.py` compares the original data in `../R/data` with the
current voteview.com files (last updated 2026-10-01); the output is in
`results/data_source_check.csv`.

* **100th–116th Houses: every cast vote is identical** to both the original
  `.ord` files and the April 2021 `Hall_votes.csv` dump: 0 flipped votes and
  0 changes in missingness over ~9.4M legislator×roll-call cells. The only extra
  cells in the current files are non-voting delegates (Committee of the Whole
  votes in the 103rd, 110th, 111th and 116th) and presidential positions;
  both are excluded. The old 116th `.ord` stops at roll call 808 (June 2020);
  the current file has all 952.
* After preprocessing, (I, J) matches paper Table 4 exactly for 100–103, 105,
  107, 109, 110, 112–114 and the 116th first session (432, 700). Differences:
  * 104th (+3), 106th (+1), 108th (+1), 111th (+1): party switchers. The old
    pipeline kept each ICPSR id as its own row, so each half-record failed the
    40% absence rule; they are now merged and kept.
  * 115th: 427 vs 450 in Table 4. The old `H115_votes.ord` has 450 rows, and
    applying the 40% rule to it also gives 427, so 450 appears to be the count
    before filtering.
* New Houses: 117th I=427, J=996; 118th I=428, J=1235; 119th (to
  2026-09-16, in progress) I=428, J=674.
