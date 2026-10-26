#!/usr/bin/env python3
"""Table and plots for a parameter sweep written by runParameterSweep().

Usage:
    python plot_param_sweep.py sweep                    # sweep = SweepConfig::outDir
    python plot_param_sweep.py sweep --latex --full     # also a LaTeX table and all per-set plots
    python plot_param_sweep.py sweep --rs-capacity 100 --show

Definitions (per parameter set, over all chains):
  error                 a decryption failure OR at least one error coefficient above the decryption
                        bound log2(q/2p) (a "wrapped" coefficient; with coefficient packing this can
                        happen before the n message coefficients are affected)
  #mult fully correct   largest k such that EVERY chain is still error-free after k multiplications
  #mult first error     smallest k at which ANY chain has an error, plus how many chains fail there
  % coeffs > bound      share of the N error coefficients above the bound at the first-error step
                        (mean over chains, [min-max])
  margin                bits between the largest error coefficient and the bound at the last fully
                        correct step (median over chains) -- how close the last correct step came

Outputs in --out-dir (default: <sweep>/plots):
  sweep_table.csv / sweep_table.md (/ sweep_table.tex with --latex)
  sweep_margin.png       remaining noise budget vs. number of multiplications, all parameter sets
  sweep_histograms.png   one row per parameter set: error distribution after 0 multiplications, at the
                         last fully correct step, at the first error, and one step later
  <id>/...               with --full: the per-set plots of plot_noise_histograms.py
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import plot_noise_histograms as pnh  # noqa: E402  (same folder)


# ----------------------------------------------------------------------------- analysis

def analyse(growth: pd.DataFrame, N: int) -> dict:
    g = growth.copy()
    g["error"] = (g.wrapped > 0) | (g.correct == 0)
    chains = sorted(g.chain.unique())
    max_step = int(g.step.max())

    first_err = {}
    for c in chains:
        e = g[(g.chain == c) & g.error]
        first_err[c] = int(e.step.min()) if not e.empty else None

    errs = [v for v in first_err.values() if v is not None]
    res = {"chains": len(chains), "maxStep": max_step, "chainsWithError": len(errs)}

    if not errs:
        res.update(fullyCorrect=max_step, firstError=None, chainsFailingAtFirst=0,
                   wrappedPct=np.nan, wrappedPctMin=np.nan, wrappedPctMax=np.nan)
        last_ok = max_step
    else:
        k = min(errs)
        at_k = g[g.step == k]
        pct = 100.0 * at_k.wrapped / N
        res.update(fullyCorrect=k - 1, firstError=k,
                   chainsFailingAtFirst=int(at_k.error.sum()),
                   wrappedPct=pct.mean(), wrappedPctMin=pct.min(), wrappedPctMax=pct.max())
        res["firstErrorRange"] = (min(errs), max(errs))
        last_ok = k - 1

    ok = g[(g.step == last_ok) & ~g.error]
    res["marginBits"] = float((ok.log2Threshold - ok.log2MaxNoise).median()) if not ok.empty else np.nan
    res["lastCorrectStep"] = last_ok
    return res


def build_table(index: pd.DataFrame, results: dict) -> pd.DataFrame:
    rows = []
    for _, r in index.iterrows():
        base = {
            "id": r.id, "label": r.label,
            "n": int(r.ringDim) if pd.notna(r.ringDim) else "—",
            "log2 q": r.log2Q if pd.notna(r.log2Q) else "—", "sigma": r.sigma,
            "p": r.plaintextModulus, "depth": r.multDepth, "secure": "yes" if r.secure == 1 else "no",
        }
        if r.status != "ok":
            base.update({"#mult fully correct": "—", "#mult first error": "—",
                         "chains failing": "—", "% coeffs > bound": "—", "margin (bits)": "—",
                         "note": str(r.status).split("(): ")[-1][:100]})
            rows.append(base)
            continue
        a = results[r.id]
        if a["firstError"] is None:
            first, failing, pct = f"> {a['maxStep']}", "0/" + str(a["chains"]), "—"
        else:
            first = str(a["firstError"])
            lo, hi = a["firstErrorRange"]
            if hi != lo:
                first += f" (chains: {lo}–{hi})"
            failing = f"{a['chainsFailingAtFirst']}/{a['chains']}"
            pct = f"{a['wrappedPct']:.3g}"
            if a["wrappedPctMax"] - a["wrappedPctMin"] > 0.05:
                pct += f" [{a['wrappedPctMin']:.3g}–{a['wrappedPctMax']:.3g}]"
        base.update({
            "#mult fully correct": a["fullyCorrect"],
            "#mult first error": first,
            "chains failing": failing,
            "% coeffs > bound": pct,
            "margin (bits)": f"{a['marginBits']:.1f}" if np.isfinite(a["marginBits"]) else "—",
            "note": "",
        })
        rows.append(base)
    return pd.DataFrame(rows)


def to_markdown(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
    return "\n".join(lines)


def to_latex(df: pd.DataFrame) -> str:
    cols = [c for c in df.columns if c not in ("label", "note")]
    esc = lambda s: str(s).replace("%", r"\%").replace("#", r"\#").replace("_", r"\_").replace("–", "--")
    out = [r"\begin{tabular}{" + "l" * len(cols) + "}", r"\toprule",
           " & ".join(esc(c) for c in cols) + r" \\", r"\midrule"]
    for _, r in df.iterrows():
        out.append(" & ".join(esc(r[c]) for c in cols) + r" \\")
    out += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(out)


# ----------------------------------------------------------------------------- plots

def plot_margins(index, growths, results):
    fig, ax = plt.subplots(figsize=(11, 5))
    ok = index[index.status == "ok"].reset_index(drop=True)
    for i, r in ok.iterrows():
        g = growths[r.id]
        clean = g[(g.wrapped == 0) & (g.correct == 1)]
        med = (clean.log2Threshold - clean.log2MaxNoise).groupby(clean.step).median()
        line, = ax.plot(med.index, med.values, marker="o", markersize=3, label=f"{r.id}: {r.label}")
        k = results[r.id]["firstError"]
        if k is not None:
            # small horizontal offset so markers of sets failing at the same step stay visible
            dx = (i - (len(ok) - 1) / 2) * 0.06
            last = results[r.id]["lastCorrectStep"]
            if last in med.index:
                ax.plot([last, k + dx], [med[last], 0], linestyle=":", color=line.get_color())
            ax.plot(k + dx, 0, marker="X", markersize=10, color=line.get_color(), linestyle="none",
                    markeredgecolor="white")
    ax.plot([], [], "X", color="gray", label="first error (q/2p exceeded)")
    ax.axhline(0, color="tab:red", linestyle="--", linewidth=1)
    ax.set_yscale("symlog", linthresh=8)
    ax.set_ylim(bottom=-0.6)
    ax.set_xlabel("number of multiplications")
    ax.set_ylabel("remaining noise budget (bits, symlog)\n= log2(q/2p) − log2 max|e_i|")
    ax.set_title("Noise budget per parameter set (median over chains; error-free steps only)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1))
    fig.tight_layout()
    return fig


def plot_histogram_rows(index, hists, results, args):
    ok = index[index.status == "ok"]
    fig, axes = plt.subplots(len(ok), 4, figsize=(17, 3.4 * len(ok)), squeeze=False)
    for row, (_, r) in zip(axes, ok.iterrows()):
        hist, meas, N = hists[r.id]
        a = results[r.id]
        available = set(meas.step.unique())
        k = a["firstError"]
        columns = [
            (0, "fresh (0 mult)"),
            (a["lastCorrectStep"], f"last fully correct ({a['lastCorrectStep']} mult)"),
            (k, f"first error ({k} mult)" if k is not None else None),
            (k + 1 if k is not None else None, f"one more ({k + 1} mult)" if k is not None else None),
        ]
        for ax, (step, title) in zip(row, columns):
            if step is None or step not in available:
                ax.set_visible(False)
                continue
            pnh.draw_hist_panel(ax, hist, meas, r.label, step, N, args.linear, title=title)
        row[0].set_ylabel(f"{r.id}\ncoefficients", fontsize=8)
        row[0].annotate(f"{r.id}: {r.label}", xy=(0, 1.32), xycoords="axes fraction", fontsize=9,
                        fontweight="bold")
    fig.suptitle("Distribution of the error coefficients log2|e_i| (all chains) — "
                 "red: coefficients above the decryption bound", fontsize=11)
    fig.tight_layout()
    return fig


# ----------------------------------------------------------------------------- main

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("sweep_dir", type=Path, nargs="?", default=Path("sweep"))
    parser.add_argument("--out-dir", type=Path)
    parser.add_argument("--latex", action="store_true", help="also write sweep_table.tex")
    parser.add_argument("--full", action="store_true", help="also all per-set histogram plots")
    parser.add_argument("--rs-capacity", type=int, help="for --full: errors your RS code corrects")
    parser.add_argument("--linear", action="store_true", help="linear y axis for histogram counts")
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()

    out = args.out_dir or (args.sweep_dir / "plots")
    out.mkdir(parents=True, exist_ok=True)
    index = pd.read_csv(args.sweep_dir / "sweep_index.csv")

    growths, hists, results = {}, {}, {}
    for _, r in index[index.status == "ok"].iterrows():
        growth = pd.read_csv(args.sweep_dir / r.growthFile, comment="#")
        meta = pnh.read_metadata(args.sweep_dir / r.histFile)
        N = int(meta["ringDimension"])
        hist, special = pnh.load(args.sweep_dir / r.histFile)
        meas = pnh.per_measurement(hist, special, N)
        growths[r.id], hists[r.id] = growth, (hist, meas, N)
        results[r.id] = analyse(growth, N)

    table = build_table(index, results)
    table.to_csv(out / "sweep_table.csv", index=False)
    md = to_markdown(table)
    (out / "sweep_table.md").write_text(md + "\n")
    print(md)
    if args.latex:
        (out / "sweep_table.tex").write_text(to_latex(table) + "\n")

    if results:
        plot_margins(index, growths, results).savefig(out / "sweep_margin.png", dpi=150)
        plot_histogram_rows(index, hists, results, args).savefig(out / "sweep_histograms.png", dpi=130)

    if args.full:
        for _, r in index[index.status == "ok"].iterrows():
            hist, meas, N = hists[r.id]
            sub = out / r.id
            sub.mkdir(exist_ok=True)
            steps = pnh.choose_steps(meas.step.unique(), None, 16, results[r.id]["firstError"])
            subtitle = f"{r.id}: {r.label}"
            pnh.plot_grid(hist, meas, r.label, steps, args, subtitle, N).savefig(sub / "hist_grid.png", dpi=130)
            pnh.plot_overlay(hist, meas, r.label, steps, args, subtitle).savefig(sub / "hist_overlay.png", dpi=130)
            fig, _ = pnh.plot_wrapped(meas, r.label, args, subtitle, N)
            fig.savefig(sub / "wrapped.png", dpi=130)
            plt.close("all")

    print(f"\nSaved table and figures to {out}/")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
