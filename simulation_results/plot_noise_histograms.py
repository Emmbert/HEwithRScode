#!/usr/bin/env python3
"""Histograms of the error coefficients e_i after k operations, from writeNoiseHistogramCSV().

Usage:
    python plot_noise_histograms.py noise_histograms.csv
    python plot_noise_histograms.py noise_histograms.csv --experiments mult --steps 10 11 12 13 14
    python plot_noise_histograms.py noise_histograms.csv --rs-capacity 100 --show

For every experiment, three figures are written to --out-dir (default: plots/):

  hist_grid_<exp>.png     one histogram of log2|e_i| per step (all N coefficients of all chains),
                          with the decryption bound log2(Q/2t). Coefficients whose noise crossed the
                          bound wrap around mod Q; they are drawn in red with their decrypted value
                          log2(|v_i mod Q|/t), which always lies at or below the bound.
  hist_overlay_<exp>.png  all selected steps in one plot, to see the distribution move toward the threshold.
  wrapped_<exp>.png       left:  safety margin in bits, i.e. threshold minus the largest and the RMS noise.
                                 Under a Gaussian model, the first coefficient of N wraps once
                                 T/sigma drops below about sqrt(2 ln N) (dashed line).
                          right: wrapped coefficients per ciphertext (observed, mean over chains) and the
                                 Gaussian prediction N * P(|e| > T). --rs-capacity draws the number of
                                 symbol errors your RS code can correct.

x axis: log2|e_i| in bits, y axis: number of coefficients (log scale; use --linear to change).
"""

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ----------------------------------------------------------------------------- data loading

def read_metadata(path: Path) -> dict:
    meta = {}
    with open(path) as f:
        for line in f:
            if not line.startswith("#"):
                break
            key, sep, value = line[1:].strip().partition("=")
            if sep:
                meta[key.strip()] = value.strip()
    return meta


def load(path: Path):
    """Returns (hist, special).
    hist:    numeric bins with column kind = 'ok' (log2|e_i| below the bound) or
             'wrapped' (log2(|v_i mod q|/t) of the coefficients above the bound)
    special: totals per measurement with bin = 'zero' or 'wrapped'
    Older CSVs without the kind column are read as well."""
    df = pd.read_csv(path, comment="#", dtype={"bin": str})
    if "kind" not in df:
        df["kind"] = np.where(pd.to_numeric(df["bin"], errors="coerce").notna(), "ok", "special")
    numeric = df.kind.isin(["ok", "wrapped"])
    hist = df[numeric].assign(bin=pd.to_numeric(df.loc[numeric, "bin"]).astype(float))
    special = df[~numeric].copy()
    special.loc[special.kind == "zero", "bin"] = "zero"
    special.loc[special.kind == "wrappedCount", "bin"] = "wrapped"
    return hist, special


def per_measurement(hist: pd.DataFrame, special: pd.DataFrame, N: int) -> pd.DataFrame:
    """One row per (experiment, chain, step): wrapped count, log2 max / RMS of the non-wrapped noise."""
    keys = ["experiment", "chain", "step"]
    meta_cols = ["numMults", "numAdds", "log2Threshold"]

    def stats(g):
        w = g.binWidth.iloc[0]
        mids = g.bin.to_numpy() + w / 2
        counts = g["count"].to_numpy(dtype=float)
        top = mids.max()
        s = (counts * np.exp2(2 * (mids - top))).sum()
        return pd.Series({"log2Max": g.bin.max() + w, "sumSq": s, "top": top, "measured": counts.sum()})

    ok = hist[hist.kind == "ok"]
    m = ok.groupby(keys).apply(stats, include_groups=False).reset_index()
    first = pd.concat([hist, special]).groupby(keys)[meta_cols].first().reset_index()
    m = first.merge(m, on=keys, how="left")

    sp = special.pivot_table(index=keys, columns="bin", values="count", aggfunc="sum").reset_index()
    m = m.merge(sp, on=keys, how="left")
    for col in ("wrapped", "zero"):
        if col not in m:
            m[col] = 0
        m[col] = m[col].fillna(0)
    m["measured"] = m["measured"].fillna(0)

    nonwrapped = m.measured + m.zero
    with np.errstate(divide="ignore", invalid="ignore"):
        m["log2Rms"] = m.top + 0.5 * np.log2(m.sumSq / nonwrapped)
    # Without any measured coefficient (everything wrapped) the noise is unknown.
    m.loc[m.measured == 0, ["log2Max", "log2Rms"]] = np.nan
    m["N"] = N
    return m


def aggregate_hist(hist: pd.DataFrame, experiment: str, step: int, kind: str = "ok") -> pd.DataFrame:
    g = hist[(hist.experiment == experiment) & (hist.step == step) & (hist.kind == kind)]
    return g.groupby("bin", as_index=False)["count"].sum().sort_values("bin")


def choose_steps(steps, wanted, max_panels, first_wrap):
    steps = sorted(steps)
    if wanted:
        return [s for s in steps if s in set(wanted)]
    if len(steps) <= max_panels:
        return steps
    idx = set(np.linspace(0, len(steps) - 1, max_panels).round().astype(int))
    chosen = {steps[i] for i in idx}
    if first_wrap is not None:  # always show the transition into wrapping
        pos = steps.index(first_wrap)
        chosen.update(steps[max(0, pos - 1):pos + 1])
    return sorted(chosen)


def fmt_pct(x: float) -> str:
    if x == 0 or x >= 1:
        return f"{x:.1f}%" if x < 99.95 or x == 100 else f"{x:.2f}%"
    return f"{x:.2g}%"


def dense(h: pd.DataFrame, width: float):
    """Bin edges and counts on a contiguous grid (missing bins = 0), for ax.stairs."""
    k = np.round(h.bin.to_numpy() / width).astype(int)
    grid = np.arange(k.min(), k.max() + 1)
    counts = np.zeros(len(grid))
    counts[k - k.min()] = h["count"].to_numpy()
    edges = np.append(grid, grid[-1] + 1) * width
    return edges, counts


def gaussian_tail(r: float) -> float:
    """P(|X| > r*sigma) for X ~ N(0, sigma^2)."""
    return math.erfc(r / math.sqrt(2.0))


# ----------------------------------------------------------------------------- plots

def draw_hist_panel(ax, hist, meas, experiment, step, N, linear=False, title=None):
    """Histogram of log2|e_i| for (experiment, step), summed over all chains.
    Blue: coefficients below the bound (their noise e_i).
    Red:  coefficients above the bound, shown with their decrypted value log2(|v_i mod q| / t)
          -- after the reduction mod q they always land at or below the bound.
    The two are stacked (red on top of blue), so every bin shows its total count."""
    h = aggregate_hist(hist, experiment, step, "ok")
    hw = aggregate_hist(hist, experiment, step, "wrapped")
    ms = meas[(meas.experiment == experiment) & (meas.step == step)]
    if ms.empty:
        ax.set_visible(False)
        return
    thr = ms.log2Threshold.median()
    wrapped = int(ms.wrapped.sum())
    total = int(len(ms) * N)
    both = pd.concat([h, hw])
    width = both.bin.drop_duplicates().sort_values().diff().min() if len(both) > 1 else 0.25
    width = width if np.isfinite(width) and width > 0 else 0.25
    baseline = 0 if linear else 0.5

    # Stacked: blue = below the bound, red stacked on top = above the bound (same bins).
    if not both.empty:
        k_all = np.round(both.bin.to_numpy() / width).astype(int)
        grid = np.arange(k_all.min(), k_all.max() + 1)
        edges = np.append(grid, grid[-1] + 1) * width

        def on_grid(frame):
            c = np.zeros(len(grid))
            if not frame.empty:
                k = np.round(frame.bin.to_numpy() / width).astype(int) - grid[0]
                np.add.at(c, k, frame["count"].to_numpy())
            return c

        blue, red = on_grid(h), on_grid(hw)
        if red.any():
            ax.stairs(blue + red, edges, fill=True, color="tab:red", alpha=0.75, baseline=baseline,
                      label="above bound: |v_i mod q|/t")
        if blue.any():
            ax.stairs(blue, edges, fill=True, color="tab:blue", baseline=baseline,
                      label="below bound: |e_i|")
    if hw.empty and wrapped:  # older CSV without values of wrapped coefficients: show the count only
        ax.bar(thr, wrapped, width=width * 4, align="edge", color="tab:red", label="above bound (count)")

    if not both.empty:
        lo, hi = both.bin.min(), both.bin.max() + width
    else:
        lo, hi = thr - 4, thr

    spread = max(hi - lo, 1.0)
    far = (thr - hi) > 3 * spread and wrapped == 0
    if far:  # bound far away: zoom to the data and print the margin instead
        ax.set_xlim(lo - 0.1 * spread, hi + 0.1 * spread)
    else:
        ax.axvline(thr, color="black", linestyle="--", linewidth=1, label="bound log2(q/2p)")
        ax.set_xlim(min(lo, thr) - 0.05 * spread, max(hi, thr) + 0.05 * spread)

    if not linear:
        ax.set_yscale("log")
    nm, na = int(ms.numMults.iloc[0]), int(ms.numAdds.iloc[0])
    margin = thr - ms.log2Max.max() if ms.log2Max.notna().any() else float("nan")
    info = f"above bound: {wrapped}/{total} ({fmt_pct(100 * wrapped / total)})"
    if wrapped == 0 and np.isfinite(margin):
        info += f", margin {margin:.1f} bits"
    if far:
        info += f"\nbound at {thr:.0f} bits (off-scale →)"
    ax.set_title(f"{title or f'step {step} ({nm} mult, {na} add)'}\n{info}", fontsize=8,
                 color="tab:red" if wrapped else "black")
    if wrapped:
        ax.legend(fontsize=6, loc="upper left")
    ax.set_xlabel("log2|e_i|  (bits)", fontsize=8)
    ax.set_ylabel("coefficients", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.grid(alpha=0.3)


def plot_grid(hist, meas, experiment, steps, args, subtitle, N):
    cols = min(4, len(steps))
    rows = math.ceil(len(steps) / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(4.2 * cols, 3.6 * rows), squeeze=False)
    for ax, step in zip(axes.flat, steps):
        draw_hist_panel(ax, hist, meas, experiment, step, N, args.linear)
    for ax in list(axes.flat)[len(steps):]:
        ax.set_visible(False)
    fig.suptitle(f"{experiment}: distribution of the error coefficients per step "
                 f"(all coefficients of all chains)\n{subtitle}", fontsize=10)
    fig.tight_layout()
    return fig


def plot_overlay(hist, meas, experiment, steps, args, subtitle):
    fig, ax = plt.subplots(figsize=(11, 4.8))
    cmap = plt.get_cmap("viridis")
    for k, step in enumerate(steps):
        h = aggregate_hist(hist, experiment, step)
        color = cmap(k / max(1, len(steps) - 1))
        ms = meas[(meas.experiment == experiment) & (meas.step == step)]
        thr = ms.log2Threshold.median()
        if not h.empty:
            width = h.bin.diff().min() if len(h) > 1 else 0.25
            edges, counts = dense(h, width)
            ax.stairs(counts, edges, color=color, linewidth=1.2, label=f"step {step}")
        ax.axvline(thr, color=color, linestyle=":", linewidth=1)
        hw = aggregate_hist(hist, experiment, step, "wrapped")
        if not hw.empty:
            width = hw.bin.diff().min() if len(hw) > 1 else 0.25
            edges, counts = dense(hw, width)
            ax.stairs(counts, edges, color=color, linewidth=1.2, linestyle="--")
    ax.plot([], [], "--", color="gray", label="above bound: |v_i mod q|/t")
    ax.plot([], [], ":", color="gray", label="threshold log2(Q/2t)")
    if not args.linear:
        ax.set_yscale("log")
    ax.set_xlabel("log2|e_i|  (bits)")
    ax.set_ylabel("coefficients")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, ncol=2 if len(steps) > 8 else 1, loc="upper left", bbox_to_anchor=(1.01, 1))
    ax.set_title(f"{experiment}: error distribution moving toward the threshold\n{subtitle}", fontsize=10)
    fig.tight_layout()
    return fig


def plot_wrapped(meas, experiment, args, subtitle, N):
    m = meas[meas.experiment == experiment]
    per_step = m.groupby("step").agg(
        thr=("log2Threshold", "median"),
        log2Max=("log2Max", "median"),
        log2Rms=("log2Rms", "median"),
        wrappedMean=("wrapped", "mean"),
        anyWrapped=("wrapped", lambda w: (w > 0).mean()),
        chains=("chain", "nunique"),
    ).reset_index()
    # Once coefficients wrap, the measured noise only covers the survivors (biased low).
    # For those steps, extrapolate log2 RMS linearly from the last two clean steps.
    per_step = per_step.sort_values("step").reset_index(drop=True)
    clean = per_step.wrappedMean == 0
    rms_used = per_step.log2Rms.copy()
    extrapolated = pd.Series(False, index=per_step.index)
    clean_idx = list(per_step.index[clean])
    for i in per_step.index[~clean]:
        prev = [j for j in clean_idx if j < i]
        if len(prev) >= 2:
            a, b = prev[-2], prev[-1]
            slope = (per_step.log2Rms[b] - per_step.log2Rms[a]) / (per_step.step[b] - per_step.step[a])
            rms_used[i] = per_step.log2Rms[b] + slope * (per_step.step[i] - per_step.step[b])
            extrapolated[i] = True
        else:
            rms_used[i] = np.nan
    per_step["marginMax"] = (per_step.thr - per_step.log2Max).where(clean)
    per_step["marginRms"] = per_step.thr - rms_used
    per_step["extrapolated"] = extrapolated

    # Gaussian model: e_i ~ N(0, sigma^2), sigma = RMS  ->  expected wrapped = N * P(|e| > T).
    per_step["predicted"] = [
        N * gaussian_tail(2.0 ** mr) if np.isfinite(mr) and mr < 10 else 0.0
        for mr in per_step.marginRms
    ]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4.6))

    c = per_step[~per_step.extrapolated]
    x = per_step[per_step.extrapolated]
    ax1.plot(c.step, c.marginMax, marker="o", label="threshold − log2 max|e_i|")
    ax1.plot(c.step, c.marginRms, marker="s", linestyle="--", label="threshold − log2 RMS(e)")
    if not x.empty:
        ax1.plot(x.step, x.marginRms, marker="s", linestyle="none", markerfacecolor="white",
                 color="tab:orange", label="RMS margin, extrapolated (after wrapping)")
    first_wrap_level = math.log2(math.sqrt(2 * math.log(N)))
    ax1.axhline(first_wrap_level, color="tab:red", linestyle="--", linewidth=1,
                label=f"RMS margin where the first of N wraps (≈{first_wrap_level:.1f} bits)")
    ax1.axhline(0, color="black", linewidth=0.8)
    ax1.set_yscale("symlog", linthresh=4)
    ax1.set_xlabel("step")
    ax1.set_ylabel("safety margin (bits, symlog)")
    ax1.set_title("Distance of the noise to the decryption threshold")
    ax1.grid(alpha=0.3)
    ax1.legend(fontsize=7)

    floor = 0.5
    observed = per_step.wrappedMean.clip(lower=floor)
    ax2.plot(per_step.step, observed, marker="o", color="tab:red", label="observed (mean per ciphertext)")
    zero = per_step.wrappedMean == 0
    ax2.scatter(per_step.step[zero], [floor] * zero.sum(), facecolors="white", edgecolors="tab:red",
                zorder=5, label="observed 0")
    pred = per_step.predicted.where(per_step.predicted > 1e-3)
    ax2.plot(per_step.step, pred, marker="x", linestyle=":", color="tab:purple",
             label="Gaussian prediction N·P(|e|>T)")
    if args.rs_capacity:
        ax2.axhline(args.rs_capacity, color="tab:green", linestyle="--",
                    label=f"RS capacity ({args.rs_capacity} errors)")
    ax2.axhline(N, color="gray", linestyle=":", linewidth=1, label=f"N = {N}")
    ax2.set_yscale("log")
    ax2.set_ylim(bottom=floor * 0.6)
    ax2.set_xlabel("step")
    ax2.set_ylabel("wrapped coefficients per ciphertext")
    ax2.set_title("Coefficients whose noise exceeded the threshold")
    ax2.grid(alpha=0.3)
    ax2.legend(fontsize=7)

    fig.suptitle(f"{experiment}\n{subtitle}", fontsize=10)
    fig.tight_layout()
    return fig, per_step


# ----------------------------------------------------------------------------- main

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("csv", type=Path, nargs="?", default=Path("noise_histograms.csv"))
    parser.add_argument("--experiments", nargs="*", help="only these experiments (default: all)")
    parser.add_argument("--steps", nargs="*", type=int, help="only these steps (default: up to --max-panels)")
    parser.add_argument("--max-panels", type=int, default=16, help="max histograms per grid figure")
    parser.add_argument("--rs-capacity", type=int, help="number of coefficient errors your RS code corrects")
    parser.add_argument("--linear", action="store_true", help="linear instead of log y axis for counts")
    parser.add_argument("--out-dir", type=Path, default=Path("plots"))
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()

    meta = read_metadata(args.csv)
    N = int(meta.get("ringDimension", 0))
    hist, special = load(args.csv)
    if N == 0:
        raise SystemExit("ringDimension missing in the CSV header.")
    meas = per_measurement(hist, special, N)

    experiments = list(dict.fromkeys(pd.concat([hist, special]).experiment))
    if args.experiments:
        experiments = [e for e in experiments if e in args.experiments]
    if not experiments:
        raise SystemExit("No matching experiments.")

    subtitle = ", ".join(f"{k}={v}" for k, v in meta.items() if k != "bin")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    for exp in experiments:
        m = meas[meas.experiment == exp]
        wrapped_steps = sorted(m[m.wrapped > 0].step.unique())
        first_wrap = wrapped_steps[0] if wrapped_steps else None
        steps = choose_steps(m.step.unique(), args.steps, args.max_panels, first_wrap)
        if not steps:
            print(f"\n=== {exp} === no measurements at the requested steps, skipped.")
            continue

        safe = exp.replace("/", "_").replace("+", "_")
        fig1 = plot_grid(hist, meas, exp, steps, args, subtitle, N)
        fig1.savefig(args.out_dir / f"hist_grid_{safe}.png", dpi=150)
        fig2 = plot_overlay(hist, meas, exp, steps, args, subtitle)
        fig2.savefig(args.out_dir / f"hist_overlay_{safe}.png", dpi=150)
        fig3, per_step = plot_wrapped(meas, exp, args, subtitle, N)
        fig3.savefig(args.out_dir / f"wrapped_{safe}.png", dpi=150)

        print(f"\n=== {exp} ===")
        print(f"{'step':>5}{'chains':>8}{'margin max':>12}{'margin RMS':>12}{'wrapped/ct':>12}{'P(ct has wrap)':>16}")
        for _, r in per_step.iterrows():
            print(f"{int(r.step):>5}{int(r.chains):>8}{r.marginMax:>12.1f}{r.marginRms:>12.1f}"
                  f"{r.wrappedMean:>12.1f}{r.anyWrapped:>16.2f}")
        if first_wrap is not None:
            print(f"First wrapped coefficient at step {first_wrap}.")
        else:
            print("No coefficient exceeded the threshold.")

    print(f"\nSaved figures to {args.out_dir}/")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
