#!/usr/bin/env python3
"""Compare candidate parameter sets after the exact sweep: ciphertext size vs. errors at step K.

Usage:
    python compare_candidates.py final --K 10 --candidates candidates.csv
    python compare_candidates.py final --K 10 --latex --show

Errors are always counted in the coefficients of the decrypted polynomial (before the conversion
to slots with SIMD packing). Per parameter set (all chains of the sweep):
  ct size        theoretical size of a fresh ciphertext, 2 * n * log2(q) bits (KiB)
  ct RNS         size in memory, 2 * n * towers * 8 bytes (KiB)
  rel. size      ct size relative to the smallest set with ALL K multiplications error-free
  fully correct  multiplications without any coefficient above the bound (min over chains)
  % > bound @K   coefficients above the bound after K multiplications (mean [min-max] over chains)
  % correct @K   100 - the above: coefficients still below the bound
  category       'K ok'            every chain error-free through K multiplications
                 'K-1 ok, partial' every chain error-free through K-1, some coefficients
                                   above the bound at K, but not all
                 'K-1 ok, all'     error-free through K-1, but (almost) every coefficient wrong at K
                 'fails < K'       an error before multiplication K

Outputs in --out-dir (default: <sweep>/plots): compare_table.csv/.md(/.tex), compare_tradeoff.png
(size vs. % correct at K), compare_histograms.png (error distribution at step K per set).
"""

import argparse
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import plot_noise_histograms as pnh  # noqa: E402


def analyse(growth: pd.DataFrame, N: int, K: int) -> dict:
    g = growth.copy()
    g["error"] = (g.wrapped > 0) | (g.correct == 0)
    fully, at_k, at_km1 = [], [], []
    for _, c in g.groupby("chain"):
        errs = c[c.error].step
        fully.append(int(errs.min()) - 1 if not errs.empty else int(c.step.max()))
        k_row = c[c.step == K]
        at_k.append(100.0 * k_row.wrapped.iloc[0] / N if not k_row.empty else np.nan)
        km1 = c[c.step == K - 1]
        at_km1.append(bool(km1.error.iloc[0]) if not km1.empty else True)
    at_k = np.array(at_k, dtype=float)
    res = {"chains": len(fully), "fullyMin": min(fully), "fullyMax": max(fully),
           "pctK": np.nanmean(at_k) if np.isfinite(at_k).any() else np.nan,
           "pctKmin": np.nanmin(at_k) if np.isfinite(at_k).any() else np.nan,
           "pctKmax": np.nanmax(at_k) if np.isfinite(at_k).any() else np.nan}
    if res["fullyMin"] >= K:
        res["category"] = "K ok"
    elif res["fullyMin"] >= K - 1:
        res["category"] = "K-1 ok, partial" if res["pctKmax"] < 99.9 else "K-1 ok, all"
    else:
        res["category"] = "fails < K"
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sweep_dir", type=Path)
    ap.add_argument("--K", type=int, required=True, help="target number of multiplications")
    ap.add_argument("--candidates", type=Path, help="candidates.csv from find_candidates.py (adds targets)")
    ap.add_argument("--out-dir", type=Path)
    ap.add_argument("--latex", action="store_true")
    ap.add_argument("--linear", action="store_true", help="linear y axis in the histograms")
    ap.add_argument("--show", action="store_true")
    args = ap.parse_args()

    out = args.out_dir or (args.sweep_dir / "plots")
    out.mkdir(parents=True, exist_ok=True)
    index = pd.read_csv(args.sweep_dir / "sweep_index.csv")
    cand = pd.read_csv(args.candidates) if args.candidates else None

    rows, hists = [], {}
    packing = "?"
    for _, r in index[index.status == "ok"].iterrows():
        packing = pnh.read_metadata(args.sweep_dir / r.histFile).get("packing", "?")
        break
    print(f"Packing: {packing} (errors are counted in the coefficients of the decrypted polynomial)")
    for _, r in index.iterrows():
        base = {"id": r.id, "label": r.label}
        if cand is not None:
            m = cand[cand.label == r.label]
            base["target"] = m.target.iloc[0] if not m.empty else ""
            base["pred. % > bound @K"] = m.predictedPctAboveBoundAtK.iloc[0] if not m.empty else np.nan
        if r.status != "ok":
            base.update({"category": "rejected by OpenFHE", "note": str(r.status)[:90]})
            rows.append(base)
            continue
        N = int(r.ringDim)
        growth = pd.read_csv(args.sweep_dir / r.growthFile, comment="#")
        a = analyse(growth, N, args.K)
        hist, special = pnh.load(args.sweep_dir / r.histFile)
        hists[r.id] = (hist, pnh.per_measurement(hist, special, N), N, r.label)
        base.update({
            "n": N, "log2 q": float(r.log2Q), "towers": int(r.numTowers),
            "digits": int(r.numLargeDigits) if pd.notna(r.get("numLargeDigits")) else "",
            "log2 P": r.get("log2P", ""), "sigma": float(r.sigma),
            "ct size KiB": 2 * N * float(r.log2Q) / 8 / 1024,
            "ct RNS KiB": 2 * N * int(r.numTowers) * 8 / 1024,
            "fully correct": a["fullyMin"] if a["fullyMin"] == a["fullyMax"]
            else f"{a['fullyMin']}–{a['fullyMax']}",
            "pctK": a["pctK"], "pctKmin": a["pctKmin"], "pctKmax": a["pctKmax"],
            "chains": a["chains"], "category": a["category"], "note": "",
        })
        rows.append(base)

    df = pd.DataFrame(rows)
    ok = df[df.category == "K ok"]
    ref = ok["ct size KiB"].min() if not ok.empty else np.nan
    df["rel. size"] = df["ct size KiB"] / ref * 100 if np.isfinite(ref) else np.nan
    df = df.sort_values("ct size KiB", na_position="last")

    def fmt_pct(r):
        if not np.isfinite(r.get("pctK", np.nan)):
            return "—"
        s = f"{r.pctK:.3g}"
        if r.pctKmax - r.pctKmin > 0.05:
            s += f" [{r.pctKmin:.3g}–{r.pctKmax:.3g}]"
        return s

    table = pd.DataFrame({
        "label": df.label,
        **({"target": df.target} if "target" in df else {}),
        "n": df.get("n"), "log2 q": df.get("log2 q").round(1) if "log2 q" in df else None,
        "towers": df.get("towers"), "digits": df.get("digits"), "log2 P": df.get("log2 P"),
        "sigma": df.get("sigma"),
        "ct size KiB": df.get("ct size KiB").round(0) if "ct size KiB" in df else None,
        "ct RNS KiB": df.get("ct RNS KiB").round(0) if "ct RNS KiB" in df else None,
        "fully correct": df.get("fully correct"),
        "rel. size %": df["rel. size"].round(1),
        f"% > bound @{args.K}": [fmt_pct(r) for _, r in df.iterrows()],
        f"errors @{args.K} (coeffs)": [f"{r.pctK * r.n / 100:.4g}" if np.isfinite(r.get("pctK", np.nan)) else "—"
                                     for _, r in df.iterrows()],
        f"% correct @{args.K}": [f"{100 - r.pctK:.3g}" if np.isfinite(r.get("pctK", np.nan)) else "—"
                                 for _, r in df.iterrows()],
        **({f"pred. % > bound @{args.K}": df["pred. % > bound @K"]} if "pred. % > bound @K" in df else {}),
        "category": df.category,
        "note": df.get("note"),
    }).fillna("—")

    table.to_csv(out / "compare_table.csv", index=False)
    cols = list(table.columns)
    md = ["| " + " | ".join(cols) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    md += ["| " + " | ".join(str(v) for v in row) + " |" for row in table.itertuples(index=False)]
    (out / "compare_table.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))
    if args.latex:
        esc = lambda s: str(s).replace("%", r"\%").replace("_", r"\_").replace("–", "--").replace("#", r"\#")
        tex_cols = [c for c in cols if c not in ("label", "note")]
        tex = [r"\begin{tabular}{" + "l" * len(tex_cols) + "}", r"\toprule",
               " & ".join(esc(c) for c in tex_cols) + r" \\", r"\midrule"]
        tex += [" & ".join(esc(row[c]) for c in tex_cols) + r" \\" for _, row in table.iterrows()]
        tex += [r"\bottomrule", r"\end{tabular}"]
        (out / "compare_table.tex").write_text("\n".join(tex) + "\n")

    good = df[df.category != "rejected by OpenFHE"].dropna(subset=["pctK"])
    if not good.empty:
        # trade-off: size vs. correct coefficients at K
        fig, ax = plt.subplots(figsize=(9, 5))
        for n, g in good.groupby("n"):
            ax.scatter(g["ct size KiB"], 100 - g.pctK, s=60, label=f"n = {int(n)}")
            for _, r in g.iterrows():
                ax.annotate(r.get("target", r.id) or r.id, (r["ct size KiB"], 100 - r.pctK),
                            textcoords="offset points", xytext=(5, 4), fontsize=7)
        ax.set_xscale("log")
        ax.set_xlabel("ciphertext size (KiB, 2·n·log2 q)")
        ax.set_ylabel(f"polynomial coefficients below the bound after {args.K} mults (%)")
        ax.set_title(f"Size vs. correctness after {args.K} multiplications")
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(out / "compare_tradeoff.png", dpi=150)

        # error distribution at step K for every set
        ids = [i for i in good.id if i in hists]
        cols_n = min(4, len(ids))
        rows_n = math.ceil(len(ids) / cols_n)
        fig, axes = plt.subplots(rows_n, cols_n, figsize=(4.4 * cols_n, 3.6 * rows_n), squeeze=False)
        for ax, i in zip(axes.flat, ids):
            hist, meas, N, label = hists[i]
            r = good[good.id == i].iloc[0]
            title = f"{label}\nct {r['ct size KiB']:.0f} KiB"
            pnh.draw_hist_panel(ax, hist, meas, label, args.K, N, args.linear, title=title)
        for ax in list(axes.flat)[len(ids):]:
            ax.set_visible(False)
        fig.suptitle(f"Error distribution after {args.K} multiplications (blue: below the bound, "
                     f"red: above, shown as |v_i mod q|/t)", fontsize=10)
        fig.tight_layout()
        fig.savefig(out / "compare_histograms.png", dpi=130)

    print(f"\nSaved table and figures to {out}/")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
