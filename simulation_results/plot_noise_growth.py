#!/usr/bin/env python3
"""Plot noise growth from the CSV written by writeOpChainCSV().

Usage:
    python plot_noise_growth.py noise_growth.csv
    python plot_noise_growth.py noise_growth.csv --x numMults --show
    python plot_noise_growth.py noise_growth.csv --experiments mult add

Creates two figures in --out-dir (default: plots/):
    noise_growth_per_experiment.png  one panel per experiment: every chain (thin lines),
                                     median and min-max band of log2|e|max, RMS noise,
                                     decryption threshold and failures
    noise_growth_overlay.png         all experiments in one plot, plus the remaining
                                     noise budget (threshold - noise) in bits
"""

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def read_metadata(path: Path) -> dict:
    meta = {}
    with open(path) as f:
        for line in f:
            if not line.startswith("#"):
                break
            key, _, value = line[1:].strip().partition("=")
            meta[key.strip()] = value.strip()
    return meta


def summarize(df: pd.DataFrame) -> None:
    print(f"{'experiment':<20}{'chains':>8}{'failed':>8}{'mean steps ok':>15}{'max steps ok':>14}{'first wrap':>12}")
    for name, g in df.groupby("experiment", sort=False):
        ok_steps = g[g.correct == 1].groupby("chain").step.max()
        failed_chains = g[g.correct == 0].chain.nunique()
        wrapped_steps = g[g.wrapped > 0].step
        first_wrap = str(int(wrapped_steps.min())) if not wrapped_steps.empty else "-"
        print(f"{name:<20}{g.chain.nunique():>8}{failed_chains:>8}"
              f"{ok_steps.mean():>15.1f}{ok_steps.max():>14}{first_wrap:>12}")


def plot_experiment(ax, g: pd.DataFrame, x: str, name: str) -> None:
    for _, chain in g.groupby("chain"):
        ax.plot(chain[x], chain.log2MaxNoise, color="tab:blue", alpha=0.15, linewidth=0.8)

    stats = g.groupby(x).agg(
        med=("log2MaxNoise", "median"),
        lo=("log2MaxNoise", "min"),
        hi=("log2MaxNoise", "max"),
        rms=("log2RmsNoise", "median"),
        thr=("log2Threshold", "median"),
    ).reset_index()

    ax.fill_between(stats[x], stats.lo, stats.hi, color="tab:blue", alpha=0.2, label="min-max over chains")
    ax.plot(stats[x], stats.med, color="tab:blue", marker="o", markersize=3, label="median log2 max|e|")
    ax.plot(stats[x], stats.rms, color="tab:green", linestyle="--", label="median log2 RMS(e)")
    ax.step(stats[x], stats.thr, where="post", color="black", linestyle=":", label="threshold log2(Q/2t)")

    failures = g[g.correct == 0]
    if not failures.empty:
        ax.scatter(failures[x], failures.log2MaxNoise, color="tab:red", marker="x", s=50, zorder=5,
                   label=f"decryption failure ({failures.chain.nunique()} chains)")

    # If the threshold is far above the noise (e.g. pure additions), zoom in on the noise
    # and print the threshold instead of drawing it, so the growth stays visible.
    noise_top = stats.hi.max()
    if stats.thr.min() > 3 * max(noise_top, 1.0):
        ax.set_ylim(0, noise_top * 1.3 + 2)
        ax.text(0.98, 0.04, f"threshold ~ {stats.thr.max():.0f} bits (off-scale)", transform=ax.transAxes,
                ha="right", fontsize=8, color="black")

    ax.set_title(name)
    ax.set_xlabel(x)
    ax.set_ylabel("bits (log2)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, loc="best")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("csv", type=Path, nargs="?", default=Path("noise_growth.csv"))
    parser.add_argument("--x", choices=["step", "numMults", "numAdds"], default="step",
                        help="quantity on the x axis (default: step = total number of operations)")
    parser.add_argument("--experiments", nargs="*", help="only plot these experiments")
    parser.add_argument("--out-dir", type=Path, default=Path("plots"))
    parser.add_argument("--show", action="store_true", help="also open the figures in a window")
    args = parser.parse_args()

    meta = read_metadata(args.csv)
    df = pd.read_csv(args.csv, comment="#")
    if args.experiments:
        df = df[df.experiment.isin(args.experiments)]
    if df.empty:
        raise SystemExit("No rows to plot.")

    summarize(df)
    subtitle = ", ".join(f"{k}={v}" for k, v in meta.items())
    args.out_dir.mkdir(parents=True, exist_ok=True)
    experiments = list(dict.fromkeys(df.experiment))  # keep CSV order

    # Figure 1: one panel per experiment
    cols = min(2, len(experiments))
    rows = math.ceil(len(experiments) / cols)
    fig1, axes = plt.subplots(rows, cols, figsize=(6.5 * cols, 4.5 * rows), squeeze=False)
    for ax, name in zip(axes.flat, experiments):
        plot_experiment(ax, df[df.experiment == name], args.x, name)
    for ax in list(axes.flat)[len(experiments):]:
        ax.set_visible(False)
    fig1.suptitle(f"Noise growth per experiment\n{subtitle}", fontsize=10)
    fig1.tight_layout()
    out1 = args.out_dir / "noise_growth_per_experiment.png"
    fig1.savefig(out1, dpi=150)

    # Figure 2: overlay of medians and remaining noise budget
    fig2, (ax_noise, ax_budget) = plt.subplots(1, 2, figsize=(13, 4.5))
    for name in experiments:
        g = df[df.experiment == name].copy()
        g["budget"] = g.log2Threshold - g.log2MaxNoise
        stats = g.groupby(args.x)[["log2MaxNoise", "budget"]].median().reset_index()
        line, = ax_noise.plot(stats[args.x], stats.log2MaxNoise, marker="o", markersize=3, label=name)
        ax_budget.plot(stats[args.x], stats.budget, marker="o", markersize=3, color=line.get_color(), label=name)

    ax_noise.set_title("Median log2 max|e|")
    ax_budget.set_title("Remaining noise budget (median)")
    ax_budget.axhline(0, color="tab:red", linestyle="--", linewidth=1, label="decryption fails below 0")
    for ax, ylabel in ((ax_noise, "bits (log2)"), (ax_budget, "bits left")):
        ax.set_xlabel(args.x)
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig2.suptitle(subtitle, fontsize=10)
    fig2.tight_layout()
    out2 = args.out_dir / "noise_growth_overlay.png"
    fig2.savefig(out2, dpi=150)

    print(f"\nSaved {out1} and {out2}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
