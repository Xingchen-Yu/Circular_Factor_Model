"""Parameter recovery on data simulated from the circular factor model.

Simulates I legislators and J roll calls with known (beta, psi, zeta, kappa),
fits the model, aligns the fit to the truth (common rotation / reflection,
which the likelihood cannot see) and reports recovery and interval coverage.

  uv run python scripts/simulation_study.py --I 200 --J 600 --omega 3
"""

import argparse

import jax
import jax.numpy as jnp
import numpy as np
from scipy import stats

from circfactor.analysis import circular_mean
from circfactor.model import Priors, SamplerConfig, sample
from circfactor.special import log_link_cdf


def wrap(a):
    return (a + np.pi) % (2 * np.pi) - np.pi


def simulate(I, J, omega, kappa_shape, kappa_rate, seed):
    rng = np.random.default_rng(seed)
    beta = stats.vonmises.rvs(omega, size=I, random_state=rng)
    yes, no = rng.uniform(-np.pi, np.pi, J), rng.uniform(-np.pi, np.pi, J)
    kappa = rng.gamma(kappa_shape, 1 / kappa_rate, J)
    e = wrap(no[None] - beta[:, None]) ** 2 - wrap(yes[None] - beta[:, None]) ** 2
    p = np.exp(np.asarray(log_link_cdf(jnp.asarray(e), jnp.asarray(kappa)[None])))
    y = (rng.random((I, J)) < p).astype(float)
    y[rng.random((I, J)) < 0.03] = np.nan  # ~3% missing, as in the House data
    return y, dict(beta=beta, yes=yes, no=no, kappa=kappa)


def align(draws, truth):
    """Best common rotation/reflection of posterior draws onto the truth."""
    m = circular_mean(draws, axis=0)
    best = None
    for sign in (1, -1):
        shift = circular_mean(wrap(truth - sign * m)[None], axis=1)[0]
        err = np.mean(1 - np.cos(truth - (sign * m + shift)))
        if best is None or err < best[0]:
            best = (err, sign, shift)
    _, sign, shift = best
    return wrap(sign * draws + shift), sign, shift


def circular_corr(a, b):
    """Jammalamadaka-SenGupta circular correlation coefficient."""
    sa = np.sin(a - circular_mean(a[None], axis=1)[0])
    sb = np.sin(b - circular_mean(b[None], axis=1)[0])
    return np.sum(sa * sb) / np.sqrt(np.sum(sa**2) * np.sum(sb**2))


def coverage(draws, truth, level=0.9):
    """Coverage of central credible arcs around the posterior circular mean."""
    center = circular_mean(draws, axis=0)
    half = np.quantile(np.abs(wrap(draws - center)), level, axis=0)
    return np.mean(np.abs(wrap(truth - center)) <= half)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--I", type=int, default=200)
    ap.add_argument("--J", type=int, default=600)
    ap.add_argument("--omega", type=float, default=3.0, help="true concentration of the ideal points")
    ap.add_argument("--kappa-shape", type=float, default=2.0)
    ap.add_argument("--kappa-rate", type=float, default=0.05)
    ap.add_argument("--burnin", type=int, default=5000)
    ap.add_argument("--samples", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    y, truth = simulate(args.I, args.J, args.omega, args.kappa_shape, args.kappa_rate, args.seed)
    print(f"simulated I={args.I}, J={args.J}, omega={args.omega}; Yea rate {np.nanmean(y):.2f}")
    post = sample(y, Priors(), SamplerConfig(burnin=args.burnin, n_samples=args.samples, seed=args.seed, print_every=2500))

    beta, sign, shift = align(post.beta, truth["beta"])
    yes = wrap(sign * post.yes + shift)
    no = wrap(sign * post.no + shift)
    print(f"alignment: reflection={'yes' if sign < 0 else 'no'}, rotation={shift:+.3f} rad")
    for name, d, t in [("beta", beta, truth["beta"]), ("psi", yes, truth["yes"]), ("zeta", no, truth["no"])]:
        m = circular_mean(d, axis=0)
        print(f"{name:>5}: circular corr(mean, truth) = {circular_corr(m, t):.3f}, "
              f"mean |error| = {np.mean(np.abs(wrap(m - t))):.3f} rad, 90% coverage = {coverage(d, t):.2f}")
    ranks_true = truth["beta"].argsort().argsort()
    ranks_est = np.median(beta.argsort(1).argsort(1), 0)
    print(f" rank: Spearman(true, posterior median rank) = {stats.spearmanr(ranks_true, ranks_est).statistic:.3f}")
    k_mean = post.kappa.mean(0)
    lo, hi = np.quantile(post.kappa, [0.05, 0.95], axis=0)
    print(f"kappa: corr(log mean, log truth) = {np.corrcoef(np.log(k_mean), np.log(truth['kappa']))[0, 1]:.3f}, "
          f"90% coverage = {np.mean((truth['kappa'] >= lo) & (truth['kappa'] <= hi)):.2f}")
    print(f"omega: posterior mean {post.omega.mean():.2f} (truth {args.omega}), 90% CI {np.quantile(post.omega, [0.05, 0.95]).round(2)}")
    print(f"accept {post.accept}, {post.seconds / 60:.1f} min")


if __name__ == "__main__":
    main()
