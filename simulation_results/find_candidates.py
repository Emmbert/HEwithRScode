#!/usr/bin/env python3
"""Find secure BGV parameter sets for a fixed plaintext modulus p and K multiplications.

For each (n, sigma) the script computes the ciphertext modulus q that puts the noise after
K multiplications at a chosen position relative to the decryption bound q/(2p):

  base        all K multiplications error-free (safety margin, like a normal parameter choice)
  f=25%       K-1 multiplications error-free, the K-th pushes ~25% of the polynomial coefficients
              over the bound (75% still correct) -- one row per value of --targets
  errors=10   optional (--target-errors): on average 10 coefficients over the bound at step K
  maxq        optional (--max-q): the largest q that is still secure for this n
  For error targets, q is rounded up to whole bits and sigma is raised slightly to hit the target
  (sigma only grows: more secure). Errors are always counted in the coefficients of the decrypted
  polynomial -- with SIMD packing that is before the conversion to slots.

q is realized as RNS towers (firstModSize + depth * scalingModSize bits), and the security of
(n, q, P) is checked exactly like OpenFHE does (HE standard, ternary secret, 128-bit classic,
including the key-switching modulus P). The number of key-switching digits is chosen so that P is
as small as possible -- this often decides whether a smaller n is allowed.

Typical workflow (K = 10, p = 65537):

  1) python find_candidates.py calibrate --p 65537 --K 10 --n 16384 32768 --out calib_params.csv
     C++:  runParameterSweep(loadHEParamsCSV("calib_params.csv"), cfg)   with cfg.outDir = "calib",
           cfg.maxMults = 10, cfg.numChains = 10
     (generous q, so every chain survives K multiplications; this measures the exact noise level --
      without ModReduce the noise does not depend on q)

  2) python find_candidates.py select --p 65537 --K 10 --calib calib --targets 0.1 0.25 0.5 \\
         --out candidates.csv
     C++:  runParameterSweep(loadHEParamsCSV("candidates.csv"), cfg)     with cfg.outDir = "final",
           cfg.maxMults = 11 (K + 1), cfg.stepsAfterFailure = 1

  3) python compare_candidates.py final --K 10 --candidates candidates.csv

Step 2 also works without --calib; it then uses the analytical model (less precise, ~0.5 bit per
multiplication, i.e. a few bits at K = 10 -- fine for 'base', too coarse to hit 25% exactly).
"""

import argparse
import csv
import math
from statistics import NormalDist
from pathlib import Path

# HE standard, ternary secret, 128-bit classic: max log2(Q*P) per ring dimension
# (OpenFHE src/core/lib/lattice/stdlatticeparms.cpp)
HE_STD_TERNARY_128 = {1024: 27, 2048: 54, 4096: 109, 8192: 218, 16384: 438, 32768: 881,
                      65536: 1747, 131072: 3523}
AUX_BITS = 60        # size of the special primes in P (DCRT_MODULUS::MAX_SIZE)
MAX_MOD_BITS = 60


# ----------------------------------------------------------------------------- noise model

def model_rms(n: int, sigma: float, p: int, k: int) -> float:
    """log2 RMS of the error after k multiplications (no ModReduce), fitted to measurements:
    fresh:   log2 sigma + 1/2 log2 n + 0.23
    growth:  log2 p + log2 rms_fresh + 1/2 log2 n + 0.6   per multiplication"""
    rms0 = math.log2(sigma) + 0.5 * math.log2(n) + 0.23
    growth = math.log2(p) + rms0 + 0.5 * math.log2(n) + 0.6
    return rms0 + k * growth


def read_calibration(calib_dir: Path, K: int):
    """{(n, sigma, p): (median log2 RMS at step K, std over chains)} from a calibration sweep."""
    import pandas as pd
    index = pd.read_csv(calib_dir / "sweep_index.csv")
    out = {}
    for _, r in index[index.status == "ok"].iterrows():
        g = pd.read_csv(calib_dir / r.growthFile, comment="#")
        g = g[(g.step == K) & (g.wrapped == 0) & (g.correct == 1)]
        if g.empty:
            print(f"  warning: calibration {r.id} has no error-free measurement at step {K}; skipped")
            continue
        key = (int(r.ringDim), round(float(r.sigma), 4), int(r.plaintextModulus))
        out[key] = (float(g.log2RmsNoise.median()), float(g.log2RmsNoise.std(ddof=0)))
    return out


# ----------------------------------------------------------------------------- OpenFHE replicas

def default_digits(depth: int) -> int:
    """ComputeNumLargeDigits() for numLargeDigits = 0 (src/pke/include/scheme/scheme-utils.h)."""
    if depth > 3:
        return 3
    if depth > 0:
        return 2
    return 1


def estimate_log_p(tower_bits, digits):
    """CryptoParametersRNS::EstimateLogP() for FIXEDMANUAL (no extra modulus, addOne = False).
    Returns log2 P, or None if the towers cannot be split into this many digits."""
    size_q = len(tower_bits)
    digits = min(digits, size_q)
    per = math.ceil(size_q / digits)
    if size_q <= per * (digits - 1):
        return None
    max_bits = 0
    for j in range(digits):
        part = tower_bits[j * per: min((j + 1) * per, size_q)]
        max_bits = max(max_bits, int(sum(part)))
    return math.ceil(max_bits / AUX_BITS) * AUX_BITS


def min_prime_bits(n: int, p: int) -> int:
    """Moduli must be primes = 1 mod lcm(2n, p) (FIXEDMANUAL); keep a few bits of room above it."""
    odd, pow2 = p, 1
    while odd % 2 == 0:
        odd //= 2
        pow2 *= 2
    order = max(pow2, 2 * n) * odd
    return math.ceil(math.log2(order)) + 4


def realize_towers(q_bits: int, n: int, p: int):
    """Split q_bits into firstModSize + depth * scalingModSize with as few towers as possible."""
    lo = min_prime_bits(n, p)
    for towers in range(max(1, math.ceil(q_bits / MAX_MOD_BITS)), 64):
        s = math.ceil(q_bits / towers)
        f = q_bits - (towers - 1) * s
        if s <= MAX_MOD_BITS and lo <= f <= MAX_MOD_BITS and s >= lo:
            return f, s, towers - 1
    return None


def best_digits(f, s, depth):
    """Smallest digit count that gives the smallest P (more digits = larger keys, not ciphertexts)."""
    towers = [f] + [s] * depth
    best = None
    for d in range(1, len(towers) + 1):
        lp = estimate_log_p(towers, d)
        if lp is not None and (best is None or lp < best[1]):
            best = (d, lp)
    return best


def min_secure_n(q_bound: int):
    for n in sorted(HE_STD_TERNARY_128):
        if q_bound <= HE_STD_TERNARY_128[n]:
            return n
    return None


# ----------------------------------------------------------------------------- commands

PARAM_COLUMNS = ["label", "ringDim", "multDepth", "firstModSize", "scalingModSize", "sigma",
                 "plaintextModulus", "secure", "numLargeDigits"]


def cmd_calibrate(args):
    args.n = args.n or [16384, 32768]
    rows = []
    for n in args.n:
        for sigma in args.sigma:
            rms = model_rms(n, sigma, args.p, args.K)
            q_bits = math.ceil(rms + 40 + math.log2(args.p) + 1)  # 40 bits of slack: nobody fails
            towers = realize_towers(q_bits, n, args.p)
            if towers is None:
                print(f"  n={n} sigma={sigma}: cannot realize {q_bits} bits of q; skipped")
                continue
            f, s, depth = towers
            rows.append({"label": f"calib n={n} sigma={sigma} p={args.p}", "ringDim": n, "multDepth": depth,
                         "firstModSize": f, "scalingModSize": s, "sigma": sigma,
                         "plaintextModulus": args.p, "secure": 0, "numLargeDigits": 0})
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=PARAM_COLUMNS)
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {len(rows)} calibration sets to {args.out} (insecure on purpose: q is oversized).")
    print(f"Run them with cfg.maxMults = {args.K}, then: find_candidates.py select --calib <outDir> ...")


def is_prime(x: int) -> bool:
    if x < 2:
        return False
    i = 2
    while i * i <= x:
        if x % i == 0:
            return False
        i += 1
    return True


def expected_errors(n, q_bits, p_bits, rms, nd):
    """Expected number of polynomial coefficients above the bound q/(2p) at step K."""
    return n * 2 * (1 - nd.cdf(2 ** ((q_bits - 1 - p_bits) - rms)))


def max_secure_q(n, p):
    """Largest q (whole bits) that passes the security check for n with the best digit count."""
    bound = HE_STD_TERNARY_128.get(n)
    if bound is None:
        return None
    for q_bits in range(bound - AUX_BITS, 20, -1):
        towers = realize_towers(q_bits, n, p)
        if towers is None:
            continue
        f, s, depth = towers
        d, lp = best_digits(f, s, depth)
        if q_bits + lp <= bound:
            return q_bits
    return None


def cmd_select(args):
    calib = read_calibration(Path(args.calib), args.K) if args.calib else {}
    p_bits = math.log2(args.p)
    nd = NormalDist()
    simd = args.packing == "simd"

    # Targets: 'base', shares of erroneous coefficients, optional absolute counts and 'maxq'.
    targets = [("base", None, None)]
    targets += [(f"f={100 * f:g}%", f, None) for f in args.targets]
    targets += [(f"errors={c:g}", None, c) for c in (args.target_errors or [])]
    if args.max_q:
        targets.append(("maxq", None, None))

    sources = []
    if calib:
        for (n, sigma, p), (rms, std) in sorted(calib.items()):
            if p == args.p and (not args.n or n in args.n):
                sources.append((n, sigma, rms, max(std, 0.1), "measured"))
    else:
        for n in (args.n or [16384, 32768]):
            for sigma in args.sigma:
                sources.append((n, sigma, model_rms(n, sigma, args.p, args.K), 0.25, "model"))

    rows = []
    for n, sigma, rms, std, origin in sources:
        if simd and not (is_prime(args.p) and (args.p - 1) % (2 * n) == 0):
            print(f"  n={n}: SIMD packing needs a prime p = 1 mod 2n; p={args.p} is not -- skipped")
            continue
        for name, f, count in targets:
            if name == "maxq":
                q_bits = max_secure_q(n, args.p)
                if q_bits is None:
                    continue
                log2_t = q_bits - 1 - p_bits
            else:
                if f is None and count is None:  # base: (almost) no coefficient crosses at step K
                    z = nd.inv_cdf(1 - args.base_failure / (2 * n))
                    log2_t = rms + math.log2(z) + args.margin + 2 * std
                else:
                    share = f if f is not None else count / n
                    z = nd.inv_cdf(1 - share / 2)  # P(|e| > T) = share  for e ~ N(0, rms^2)
                    log2_t = rms + math.log2(z)
                q_bits = math.ceil(log2_t + 1 + p_bits)  # T = q / (2p), rounded up to whole bits
            towers = realize_towers(q_bits, n, args.p)
            if towers is None:
                continue

            # Fine tuning: q has whole bits only. Raise sigma so that the noise after K mults grows
            # by the missing fraction of a bit: shift = (K+1) * log2(sigma'/sigma). More secure.
            sigma_used, rms_used = sigma, rms
            if (f is not None or count is not None) and not args.no_sigma_tuning:
                shift = (q_bits - 1 - p_bits) - log2_t  # in [0, 1)
                sigma_used = sigma * 2 ** (shift / (args.K + 1))
                rms_used = rms + shift

            fbits, s, depth = towers
            d_def = default_digits(depth)
            lp_def = estimate_log_p([fbits] + [s] * depth, d_def)
            d_best, lp_best = best_digits(fbits, s, depth)
            bound_def = q_bits + lp_def if lp_def is not None else math.inf
            bound_best = q_bits + lp_best
            lam = expected_errors(n, q_bits, p_bits, rms_used, nd)
            rows.append({
                "label": f"K{args.K} {name} n={n} sigma={sigma_used:.3f}",
                "target": name, "ringDim": n, "sigma": round(sigma_used, 3), "plaintextModulus": args.p,
                "noiseSource": origin, "log2RmsAtK": round(rms_used, 2), "rmsStd": round(std, 2),
                "log2q": q_bits, "multDepth": depth, "towers": depth + 1,
                "firstModSize": fbits, "scalingModSize": s,
                "numLargeDigits": d_best, "log2P": lp_best, "log2QP": bound_best,
                "boundForN": HE_STD_TERNARY_128.get(n, ""),
                "secureDefaultDigits": int(bound_def <= HE_STD_TERNARY_128.get(n, 0)),
                "secure": int(bound_best <= HE_STD_TERNARY_128.get(n, 0)),
                "minSecureN": min_secure_n(bound_best),
                "predictedPctAboveBoundAtK": round(100 * lam / n, 4),
                "predictedErrorsAtK": round(lam, 3),
                "ctSizeKiB": round(2 * n * q_bits / 8 / 1024, 1),
                "ctRnsKiB": round(2 * n * (depth + 1) * 8 / 1024, 1),
            })

    rows.sort(key=lambda r: (not r["secure"], r["ctSizeKiB"]))
    kept = [r for r in rows if r["secure"] or args.include_insecure]
    # drop duplicates (e.g. maxq == base), keep the first
    seen, unique = set(), []
    for r in kept:
        k = (r["ringDim"], r["log2q"], r["sigma"])
        if k not in seen:
            seen.add(k)
            unique.append(r)

    cols = ["label", "ringDim", "multDepth", "firstModSize", "scalingModSize", "sigma",
            "plaintextModulus", "secure", "numLargeDigits", "target", "log2q", "towers", "log2P",
            "log2QP", "boundForN", "secureDefaultDigits", "minSecureN", "noiseSource", "log2RmsAtK",
            "rmsStd", "predictedPctAboveBoundAtK", "predictedErrorsAtK", "ctSizeKiB", "ctRnsKiB"]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(unique)

    hdr = ["label", "secure", "log2q", "towers", "digits", "log2(QP)/bound", "default digits ok",
           "pred. % > bound @K", "pred. errors @K", "ct size KiB"]
    print("| " + " | ".join(hdr) + " |")
    print("|" + "|".join("---" for _ in hdr) + "|")
    for r in rows:
        sec = "yes" if r["secure"] else f"no (needs n={r['minSecureN']})"
        print(f"| {r['label']} | {sec} | {r['log2q']} | {r['towers']} | {r['numLargeDigits']} | "
              f"{r['log2QP']}/{r['boundForN']} | {'yes' if r['secureDefaultDigits'] else 'no'} | "
              f"{r['predictedPctAboveBoundAtK']:g} | {r['predictedErrorsAtK']:.4g} | {r['ctSizeKiB']} |")
    print(f"\nWrote {len(unique)} candidates to {args.out}"
          + ("" if args.include_insecure else " (secure only; --include-insecure to keep the rest)"))
    if not calib:
        print("Noise from the analytical model -- run 'calibrate' first for bit-precise targets.")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--p", type=int, default=65537, help="plaintext modulus")
    common.add_argument("--K", type=int, default=10, help="number of multiplications")
    common.add_argument("--n", type=int, nargs="*",
                        help="ring dimensions (calibrate: default 16384 32768; select with --calib: "
                             "default all calibrated ones)")
    common.add_argument("--sigma", type=float, nargs="*", default=[3.19], help="error std deviations")

    c = sub.add_parser("calibrate", parents=[common], help="write calibration parameter sets")
    c.add_argument("--out", default="calib_params.csv")
    c.set_defaults(func=cmd_calibrate)

    s = sub.add_parser("select", parents=[common], help="compute candidate parameter sets")
    s.add_argument("--calib", help="output directory of the calibration sweep")
    s.add_argument("--packing", choices=["simd", "coef"], default="simd",
                   help="plaintext packing (must match PACKING_MODE in crypto_context.h); simd only adds "
                        "the check that p is a prime = 1 mod 2n")
    s.add_argument("--targets", type=float, nargs="*", default=[0.1, 0.25, 0.5],
                   help="shares of polynomial coefficients above the bound at step K")
    s.add_argument("--target-errors", type=float, nargs="*",
                   help="optional: expected NUMBER of polynomial coefficients above the bound at step K")
    s.add_argument("--max-q", action="store_true", help="also the largest secure q per n")
    s.add_argument("--base-failure", type=float, default=0.01,
                   help="'base': allowed probability that any coefficient crosses at step K")
    s.add_argument("--margin", type=float, default=1.0, help="'base': extra safety bits")
    s.add_argument("--include-insecure", action="store_true")
    s.add_argument("--no-sigma-tuning", action="store_true",
                   help="keep sigma fixed (targets are then only hit to within one bit of q)")
    s.add_argument("--out", default="candidates.csv")
    s.set_defaults(func=cmd_select)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
