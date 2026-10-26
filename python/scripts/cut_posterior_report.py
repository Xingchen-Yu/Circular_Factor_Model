"""Tables and figures for docs/cut_posterior_simulation.md.

Reads every results/sim_cut/omega*_seed*/ written by cut_posterior_simulation.py.

  uv run python scripts/cut_posterior_report.py
"""

import json
from pathlib import Path

import matplotlib
import matplotlib.ticker

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
SIM = ROOT / "results" / "sim_cut"
FIG = ROOT / "docs" / "figures"

# Reference palette (dataviz skill, light mode)
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE = "#2a78d6", "#eb6834"
RAMP = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]  # ordinal blue: fewer -> more votes


def style(ax):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.grid(color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def wrap(a):
    return (a + np.pi) % (2 * np.pi) - np.pi


def md(df):
    """DataFrame as a GitHub Markdown table (avoids the tabulate dependency)."""
    fmt = lambda v: f"{v:.3f}" if isinstance(v, float) else str(v)
    rows = [" | ".join(map(str, df.columns))] + [" | ".join("---" for _ in df.columns)]
    rows += [" | ".join(fmt(v) for v in r) for r in df.itertuples(index=False)]
    return "\n".join(f"| {r} |" for r in rows)


def load():
    members, summaries = [], []
    for d in sorted(SIM.glob("omega*_seed*")):
        if not (d / "summary.json").exists():
            continue
        s = json.loads((d / "summary.json").read_text())
        summaries.append(s)
        members.append(pd.read_csv(d / "members.csv").assign(omega=s["omega"], seed=s["seed"]))
    return pd.concat(members, ignore_index=True), pd.DataFrame(summaries)


def tables(m, s):
    m = m.assign(
        cut_abs_err=m.cut_err.abs(),
        joint_abs_err=m.joint_err.abs(),
        diff_over_sd=np.abs(wrap(m.cut_mean - m.joint_mean)) / m.joint_sd,
        width_ratio=m.cut_width / m.joint_width,
    )
    g = m.groupby(["omega", "n_votes"])
    t1 = pd.DataFrame(
        {
            "members": g.size(),
            "cover_cut": g.cut_covered.mean(),
            "cover_joint": g.joint_covered.mean(),
            "abs_err_cut": g.cut_abs_err.mean(),
            "abs_err_joint": g.joint_abs_err.mean(),
            "width_cut": g.cut_width.mean(),
            "width_joint": g.joint_width.mean(),
            "rank_cover_cut": g.cut_rank_covered.mean(),
            "rank_cover_joint": g.joint_rank_covered.mean(),
            "rank_width_cut": g.cut_rank_width.mean(),
        }
    ).reset_index()
    t2 = pd.DataFrame(
        {"median_diff_over_sd": g.diff_over_sd.median(), "max_diff_over_sd": g.diff_over_sd.max(), "median_width_ratio": g.width_ratio.median()}
    ).reset_index()
    tot = m.groupby("omega")[["cut_covered", "joint_covered", "cut_rank_covered", "joint_rank_covered"]].agg(["mean", "size"])
    return t1, t2, tot, naive_table()


def naive_table():
    """Naive cut vs grid, both in the fit's own frame (recomputed from draws.npz)."""
    csd = lambda x: np.sqrt(-2 * np.log(np.hypot(np.cos(x).mean(0), np.sin(x).mean(0))))
    cm = lambda x: np.arctan2(np.sin(x).mean(0), np.cos(x).mean(0))
    out = []
    for d in sorted(SIM.glob("omega*_seed*")):
        z = np.load(d / "draws.npz")
        if "naive_1" not in z:
            continue
        mem = pd.read_csv(d / "members.csv")
        g = z["placed_raw"]
        df = pd.DataFrame({"omega": json.loads((d / "summary.json").read_text())["omega"], "n_votes": mem.n_votes})
        for k in (1, 10):
            n = z[f"naive_{k}"]
            df[f"naive{k}_shift_over_sd"] = np.abs(wrap(cm(n) - cm(g))) / csd(g)
            df[f"naive{k}_sd_ratio"] = csd(n) / csd(g)
        out.append(df)
    if not out:
        return None
    df = pd.concat(out)
    g = df.groupby(["omega", "n_votes"])
    t = g.median().reset_index()
    for k in (1, 10):
        t[f"naive{k}_sd_ratio_min"] = g[f"naive{k}_sd_ratio"].min().to_numpy()
    return t


def fig_agreement(m):
    """Cut vs joint posterior means for the low-vote members."""
    omegas = sorted(m.omega.unique())
    fig, axes = plt.subplots(1, len(omegas), figsize=(4.6 * len(omegas), 4.4), facecolor=SURFACE, squeeze=False)
    votes = sorted(m.n_votes.unique())
    for ax, w in zip(axes[0], omegas):
        style(ax)
        d = m[m.omega == w]
        lim = np.pi
        ax.plot([-lim, lim], [-lim, lim], color=INK2, linewidth=1, zorder=1)
        for v, c in zip(votes, RAMP):
            e = d[d.n_votes == v]
            ax.scatter(e.joint_mean, e.cut_mean, s=36, color=c, edgecolor=SURFACE, linewidth=1.5, zorder=3, label=f"{v} votes")
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_xlabel("joint fit: posterior mean (rad)", color=INK2, fontsize=9)
        ax.set_ylabel("cut placement: posterior mean (rad)", color=INK2, fontsize=9)
        ax.set_title(f"true ω = {w:g}", color=INK, fontsize=11, loc="left")
    axes[0][-1].legend(frameon=False, fontsize=9, labelcolor=INK2, loc="upper left", title="votes cast", title_fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG / "cut_vs_joint.png", dpi=160, facecolor=SURFACE)
    plt.close(fig)


def fig_truth(m):
    """Truth vs cut posterior mean with 90% arcs."""
    omegas = sorted(m.omega.unique())
    votes = sorted(m.n_votes.unique())
    fig, axes = plt.subplots(1, len(omegas), figsize=(4.6 * len(omegas), 4.4), facecolor=SURFACE, squeeze=False)
    for ax, w in zip(axes[0], omegas):
        style(ax)
        d = m[m.omega == w]
        ax.plot([-np.pi, np.pi], [-np.pi, np.pi], color=INK2, linewidth=1, zorder=1)
        for v, c in zip(votes, RAMP):
            e = d[d.n_votes == v]
            ax.errorbar(e.truth, e.cut_mean, yerr=e.cut_width / 2, fmt="o", ms=5, color=c, mec=SURFACE, mew=1.2, elinewidth=1.2, capsize=0, zorder=3, label=f"{v} votes")
        ax.set_xlim(-np.pi, np.pi)
        ax.set_ylim(-np.pi, np.pi)
        ax.set_aspect("equal")
        ax.set_xlabel("true ideal point (rad)", color=INK2, fontsize=9)
        ax.set_ylabel("cut placement: mean ± 90% arc (rad)", color=INK2, fontsize=9)
        ax.set_title(f"true ω = {w:g}", color=INK, fontsize=11, loc="left")
    axes[0][-1].legend(frameon=False, fontsize=9, labelcolor=INK2, loc="upper left", title="votes cast", title_fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG / "cut_vs_truth.png", dpi=160, facecolor=SURFACE)
    plt.close(fig)


def fig_naive():
    """Grid (exact cut) vs naive cut draws for a few members of one cell."""
    for d in sorted(SIM.glob("omega*_seed1")):
        z = np.load(d / "draws.npz")
        if "naive_1" not in z:
            continue
        mem = pd.read_csv(d / "members.csv")
        # the member of each vote count where the naive chain is furthest off
        picks = mem.loc[mem.groupby("n_votes").naive1_sd_ratio.idxmin(), ["member", "n_votes"]].to_numpy()
        fig, axes = plt.subplots(1, len(picks), figsize=(3.3 * len(picks), 2.9), facecolor=SURFACE, squeeze=False)
        for ax, (k, v) in zip(axes[0], picks):
            style(ax)
            g, n1 = z["placed_raw"][:, k], z["naive_1"][:, k]
            c = np.arctan2(np.sin(g).mean(), np.cos(g).mean())
            g, n1 = wrap(g - c), wrap(n1 - c)
            lo, hi = np.quantile(np.concatenate([g, n1]), [0.002, 0.998])
            bins = np.linspace(lo, hi, 40)
            ax.hist(g, bins=bins, density=True, color=BLUE, alpha=0.55, label="grid (exact cut)")
            ax.hist(n1, bins=bins, density=True, histtype="step", color=ORANGE, linewidth=2, label="naive, 1 step")
            ax.set_title(f"{v} votes", color=INK, fontsize=10, loc="left")
            ax.set_xlabel("angle − grid mean (rad)", color=INK2, fontsize=9)
            ax.set_yticks([])
            ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
        axes[0][0].legend(frameon=False, fontsize=8, labelcolor=INK2, loc="upper left")
        fig.tight_layout()
        fig.savefig(FIG / f"naive_vs_grid_{d.name}.png", dpi=160, facecolor=SURFACE)
        plt.close(fig)


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    m, s = load()
    t1, t2, tot, t3 = tables(m, s)
    pd.set_option("display.width", 200)
    print("== per (omega, votes cast): 90% coverage, mean |error| (rad), 90% width (rad), rank coverage ==")
    print(md(t1))
    print("\n== cut vs joint: |difference of means| / joint sd, width ratio cut/joint ==")
    print(md(t2))
    print("\n== overall coverage ==")
    print(tot.round(3))
    print("\n== effect on fitted members and chamber summaries ==")
    print(md(s))
    if t3 is not None:
        print("\n== naive cut vs grid (median over members): |mean shift| / sd, sd ratio ==")
        print(md(t3))
    fig_agreement(m)
    fig_truth(m)
    fig_naive()
    print(f"\nfigures in {FIG}")


if __name__ == "__main__":
    main()
