#include "rs/decoders/welch_berlekamp.hpp"

#include <stdexcept>

namespace rs {
namespace {

// One element Q0(x) + Q1(x)*y of the interpolation module, with cached degrees.
struct Candidate {
    Poly q0, q1;
    int deg0 = -1, deg1 = -1;

    // (1, w)-weighted degree, w = k - 1: max(deg Q0, deg Q1 + w).
    long wdeg(long w) const {
        const long a = deg0;
        const long b = deg1 < 0 ? -1 : deg1 + w;
        return a > b ? a : b;
    }
    void refreshDegrees() {
        deg0 = degree(q0);
        deg1 = degree(q1);
    }
};

inline Elem evalUpTo(const PrimeField& F, const Poly& p, int deg, Elem x) {
    Elem acc = 0;
    for (int i = deg; i >= 0; --i) acc = F.add(F.mul(acc, x), p[static_cast<std::size_t>(i)]);
    return acc;
}

// p <- (x - a) * p
inline void mulLinear(const PrimeField& F, Poly& p, int& deg, Elem a) {
    if (deg < 0) return;
    const std::size_t need = static_cast<std::size_t>(deg) + 2;
    if (p.size() < need) p.resize(need, 0);
    for (int i = deg + 1; i >= 1; --i) {
        const std::size_t u = static_cast<std::size_t>(i);
        p[u] = F.sub(p[u - 1], F.mul(p[u], a));
    }
    p[0] = F.neg(F.mul(p[0], a));
    ++deg;
}

// dst <- s1 * dst - s2 * src
inline void combine(const PrimeField& F, Poly& dst, int dstDeg, Elem s1, const Poly& src, int srcDeg, Elem s2) {
    const int top = dstDeg > srcDeg ? dstDeg : srcDeg;
    if (top < 0) return;
    if (dst.size() < static_cast<std::size_t>(top) + 1) dst.resize(static_cast<std::size_t>(top) + 1, 0);
    for (int i = 0; i <= top; ++i) {
        const std::size_t u = static_cast<std::size_t>(i);
        const Elem a = i <= dstDeg ? F.mul(s1, dst[u]) : Elem{0};
        const Elem b = i <= srcDeg ? F.mul(s2, src[u]) : Elem{0};
        dst[u] = F.sub(a, b);
    }
}

}  // namespace

DecodeResult WelchBerlekampDecoder::decode(const ReceivedWord& rx) const {
    const PrimeField& F = code_.field();
    const std::size_t n = code_.n(), k = code_.k();
    const auto& xs = code_.evalPoints();
    if (rx.hard.size() != n) throw std::invalid_argument("WelchBerlekamp: received word has wrong length");

    const long w = static_cast<long>(k) - 1;

    // Basis of the module before any interpolation constraint: {1, y}.
    Candidate g[2];
    g[0].q0 = {1};
    g[0].deg0 = 0;
    g[1].q1 = {1};
    g[1].deg1 = 0;
    for (auto& c : g) {
        c.q0.reserve(n + 2);
        c.q1.reserve(n + 2);
    }

    // Process the points one by one: O(n) per point -> O(n^2) in total.
    for (std::size_t i = 0; i < n; ++i) {
        const Elem xi = xs[i], yi = rx.hard[i] % F.p();
        Elem delta[2];
        for (int j = 0; j < 2; ++j)
            delta[j] = F.add(evalUpTo(F, g[j].q0, g[j].deg0, xi),
                             F.mul(evalUpTo(F, g[j].q1, g[j].deg1, xi), yi));

        int js = -1;  // pivot: nonzero discrepancy with minimal weighted degree
        for (int j = 0; j < 2; ++j)
            if (delta[j] != 0 && (js < 0 || g[j].wdeg(w) < g[js].wdeg(w))) js = j;
        if (js < 0) continue;  // both already vanish at (xi, yi)

        for (int j = 0; j < 2; ++j) {
            if (j == js || delta[j] == 0) continue;
            // g_j <- delta_js * g_j - delta_j * g_js   (vanishes at (xi, yi), keeps wdeg)
            combine(F, g[j].q0, g[j].deg0, delta[js], g[js].q0, g[js].deg0, delta[j]);
            combine(F, g[j].q1, g[j].deg1, delta[js], g[js].q1, g[js].deg1, delta[j]);
            g[j].refreshDegrees();
        }
        // g_js <- (x - xi) * g_js   (raises wdeg by exactly one)
        mulLinear(F, g[js].q0, g[js].deg0, xi);
        mulLinear(F, g[js].q1, g[js].deg1, xi);
    }

    // Minimal weighted degree element: Q0 + Q1*y.
    const Candidate& best = (g[0].wdeg(w) <= g[1].wdeg(w)) ? g[0] : g[1];
    DecodeResult res;
    if (best.deg1 < 0) return res;  // Q1 == 0: no valid solution

    Poly q0 = best.q0, q1 = best.q1;
    trim(q0);
    trim(q1);
    Poly rem;
    Poly m = divmod(F, q0, q1, rem);  // Q0 + Q1*m = 0  ->  m = -Q0 / Q1
    if (!rem.empty()) return res;     // not divisible: more than t errors
    if (degree(m) >= static_cast<int>(k)) return res;
    for (auto& c : m) c = F.neg(c);

    m.resize(k, 0);
    Codeword c = code_.encodePoly(m);  // negacyclic NTT, O(n log n)
    const std::size_t dist = hammingDistance(c, rx.hard);
    if (dist > code_.t()) return res;  // stay a bounded-distance decoder

    res.status = DecodeStatus::Success;
    res.codeword = std::move(c);
    res.message = std::move(m);
    res.numCorrected = dist;
    return res;
}

}  // namespace rs
