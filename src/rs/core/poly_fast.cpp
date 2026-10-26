#include "rs/core/poly_fast.hpp"

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace rs {

namespace {

std::size_t nextPow2(std::size_t x) {
    std::size_t p = 1;
    while (p < x) p <<= 1;
    return p;
}

// a div x^h
Poly shiftDown(const Poly& a, std::size_t h) {
    if (a.size() <= h) return {};
    return Poly(a.begin() + static_cast<std::ptrdiff_t>(h), a.end());
}

void truncateTo(Poly& a, std::size_t m) {
    if (a.size() > m) a.resize(m);
    trim(a);
}

}  // namespace

// ---------------------------------------------------------------------------
// PolyMatrix
// ---------------------------------------------------------------------------
PolyMatrix PolyMatrix::identity() {
    PolyMatrix M;
    M.m00 = {1};
    M.m11 = {1};
    return M;
}

bool PolyMatrix::isIdentity() const {
    return m00 == Poly{1} && m11 == Poly{1} && m01.empty() && m10.empty();
}

// ---------------------------------------------------------------------------
// Setup and NTT
// ---------------------------------------------------------------------------
FastPolyArith::FastPolyArith(const PrimeField& F, std::size_t schoolbookThreshold,
                             std::size_t halfGcdThreshold)
    : F_(F),
      schoolbook_(std::max<std::size_t>(1, schoolbookThreshold)),
      hgcdThreshold_(std::max<std::size_t>(2, halfGcdThreshold)) {
    const std::uint64_t pm1 = F_.p() - 1;
    while (pm1 % (2 * maxLen_) == 0 && maxLen_ < (std::size_t{1} << 20)) maxLen_ *= 2;

    const Elem w = F_.rootOfUnity(maxLen_), wi = F_.inv(w);
    roots_.resize(maxLen_ / 2);
    iroots_.resize(maxLen_ / 2);
    Elem a = 1, b = 1;
    for (std::size_t j = 0; j < maxLen_ / 2; ++j) {
        roots_[j] = a;
        iroots_[j] = b;
        a = F_.mul(a, w);
        b = F_.mul(b, wi);
    }
}

void FastPolyArith::ntt(std::vector<Elem>& a, bool inverse) const {
    const std::size_t n = a.size();
    for (std::size_t i = 1, j = 0; i < n; ++i) {  // bit-reversal permutation
        std::size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }
    const std::vector<Elem>& tw = inverse ? iroots_ : roots_;
    for (std::size_t len = 2; len <= n; len <<= 1) {
        const std::size_t half = len / 2, step = maxLen_ / len;
        for (std::size_t i = 0; i < n; i += len) {
            for (std::size_t j = 0; j < half; ++j) {
                const Elem u = a[i + j];
                const Elem v = F_.mul(a[i + j + half], tw[j * step]);
                a[i + j] = F_.add(u, v);
                a[i + j + half] = F_.sub(u, v);
            }
        }
    }
    if (inverse) {
        const Elem nInv = F_.inv(static_cast<Elem>(n % F_.p()));
        for (auto& x : a) x = F_.mul(x, nInv);
    }
}

// ---------------------------------------------------------------------------
// Basic arithmetic
// ---------------------------------------------------------------------------
Poly FastPolyArith::add(const Poly& a, const Poly& b) const {
    Poly c(std::max(a.size(), b.size()), 0);
    for (std::size_t i = 0; i < c.size(); ++i)
        c[i] = F_.add(i < a.size() ? a[i] : 0, i < b.size() ? b[i] : 0);
    trim(c);
    return c;
}

Poly FastPolyArith::sub(const Poly& a, const Poly& b) const {
    Poly c(std::max(a.size(), b.size()), 0);
    for (std::size_t i = 0; i < c.size(); ++i)
        c[i] = F_.sub(i < a.size() ? a[i] : 0, i < b.size() ? b[i] : 0);
    trim(c);
    return c;
}

Poly FastPolyArith::mulSchoolbook(const Poly& a, const Poly& b) const {
    const int da = degree(a), db = degree(b);
    if (da < 0 || db < 0) return {};
    Poly c(static_cast<std::size_t>(da + db + 1), 0);
    for (int i = 0; i <= da; ++i) {
        const Elem ai = a[static_cast<std::size_t>(i)];
        if (ai == 0) continue;
        for (int j = 0; j <= db; ++j) {
            Elem& t = c[static_cast<std::size_t>(i + j)];
            t = F_.add(t, F_.mul(ai, b[static_cast<std::size_t>(j)]));
        }
    }
    trim(c);
    return c;
}

Poly FastPolyArith::mul(const Poly& a, const Poly& b) const {
    const int da = degree(a), db = degree(b);
    if (da < 0 || db < 0) return {};
    const std::size_t na = static_cast<std::size_t>(da) + 1, nb = static_cast<std::size_t>(db) + 1;
    if (std::min(na, nb) <= schoolbook_) return mulSchoolbook(a, b);

    const std::size_t need = na + nb - 1, L = nextPow2(need);
    if (L <= maxLen_) {  // single NTT product
        std::vector<Elem> fa(a.begin(), a.begin() + static_cast<std::ptrdiff_t>(na));
        std::vector<Elem> fb(b.begin(), b.begin() + static_cast<std::ptrdiff_t>(nb));
        fa.resize(L, 0);
        fb.resize(L, 0);
        ntt(fa, false);
        ntt(fb, false);
        for (std::size_t i = 0; i < L; ++i) fa[i] = F_.mul(fa[i], fb[i]);
        ntt(fa, true);
        fa.resize(need);
        trim(fa);
        return fa;
    }

    // Product longer than the largest NTT: split both operands into blocks of
    // h = maxLen/2 coefficients; each block product fits into one NTT of size
    // maxLen. Block products with equal i + j are summed in the transform domain.
    const std::size_t h = maxLen_ / 2;
    auto toBlocks = [&](const Poly& x, std::size_t nx) {
        std::vector<std::vector<Elem>> B((nx + h - 1) / h);
        for (std::size_t i = 0; i < B.size(); ++i) {
            B[i].assign(maxLen_, 0);
            const std::size_t end = std::min(nx, (i + 1) * h);
            std::copy(x.begin() + static_cast<std::ptrdiff_t>(i * h),
                      x.begin() + static_cast<std::ptrdiff_t>(end), B[i].begin());
            ntt(B[i], false);
        }
        return B;
    };
    const auto A = toBlocks(a, na), B = toBlocks(b, nb);
    Poly c(need + maxLen_, 0);
    std::vector<Elem> acc(maxLen_);
    for (std::size_t s = 0; s + 1 < A.size() + B.size(); ++s) {
        std::fill(acc.begin(), acc.end(), 0);
        for (std::size_t i = 0; i < A.size(); ++i) {
            if (s < i || s - i >= B.size()) continue;
            const auto& Ai = A[i];
            const auto& Bj = B[s - i];
            for (std::size_t t = 0; t < maxLen_; ++t) acc[t] = F_.add(acc[t], F_.mul(Ai[t], Bj[t]));
        }
        ntt(acc, true);
        for (std::size_t t = 0; t < maxLen_; ++t) c[s * h + t] = F_.add(c[s * h + t], acc[t]);
    }
    c.resize(need);
    trim(c);
    return c;
}

// ---------------------------------------------------------------------------
// Division
// ---------------------------------------------------------------------------
Poly FastPolyArith::invSeries(const Poly& b, std::size_t m) const {
    if (b.empty() || b[0] == 0) throw std::domain_error("invSeries: constant term is zero");
    Poly g{F_.inv(b[0])};
    for (std::size_t len = 1; len < m; len *= 2) {
        // Newton step: g <- g * (2 - b*g) mod x^(2 len)
        Poly bl(b.begin(), b.begin() + static_cast<std::ptrdiff_t>(std::min(b.size(), 2 * len)));
        Poly t = mul(bl, g);
        truncateTo(t, 2 * len);
        t.resize(std::max<std::size_t>(t.size(), 1), 0);
        for (auto& x : t) x = F_.neg(x);
        t[0] = F_.add(t[0], 2 % F_.p());
        trim(t);
        g = mul(g, t);
        truncateTo(g, 2 * len);
    }
    truncateTo(g, m);
    return g;
}

Poly FastPolyArith::divmod(const Poly& a, const Poly& b, Poly& rem) const {
    const int da = degree(a), db = degree(b);
    if (db < 0) throw std::domain_error("divmod: division by zero polynomial");
    if (da < db) {
        rem = a;
        trim(rem);
        return {};
    }
    const std::size_t mq = static_cast<std::size_t>(da - db) + 1;  // number of quotient coefficients
    if (mq <= schoolbook_ || static_cast<std::size_t>(db) <= schoolbook_)
        return rs::divmod(F_, a, b, rem);  // classical, O(mq * db)

    // rev(q) = rev(a) * rev(b)^-1 mod x^mq
    Poly ar(mq), br(std::min<std::size_t>(mq, static_cast<std::size_t>(db) + 1));
    for (std::size_t i = 0; i < mq; ++i) ar[i] = a[static_cast<std::size_t>(da) - i];
    for (std::size_t i = 0; i < br.size(); ++i) br[i] = b[static_cast<std::size_t>(db) - i];
    Poly qr = mul(ar, invSeries(br, mq));
    qr.resize(mq, 0);
    Poly q(mq);
    for (std::size_t i = 0; i < mq; ++i) q[i] = qr[mq - 1 - i];
    trim(q);

    rem = sub(a, mul(q, b));
    truncateTo(rem, static_cast<std::size_t>(db));  // deg rem < deg b (exact up to rounding-free arithmetic)
    return q;
}

// ---------------------------------------------------------------------------
// Matrices
// ---------------------------------------------------------------------------
PolyMatrix FastPolyArith::mul(const PolyMatrix& A, const PolyMatrix& B) const {
    PolyMatrix C;
    C.m00 = add(mul(A.m00, B.m00), mul(A.m01, B.m10));
    C.m01 = add(mul(A.m00, B.m01), mul(A.m01, B.m11));
    C.m10 = add(mul(A.m10, B.m00), mul(A.m11, B.m10));
    C.m11 = add(mul(A.m10, B.m01), mul(A.m11, B.m11));
    return C;
}

void FastPolyArith::apply(const PolyMatrix& M, Poly& a, Poly& b) const {
    Poly a2 = add(mul(M.m00, a), mul(M.m01, b));
    Poly b2 = add(mul(M.m10, a), mul(M.m11, b));
    a = std::move(a2);
    b = std::move(b2);
}

PolyMatrix FastPolyArith::euclidStep(const PolyMatrix& M, const Poly& q) const {
    // [[0, 1], [1, -q]] * M
    PolyMatrix R;
    R.m00 = M.m10;
    R.m01 = M.m11;
    R.m10 = sub(M.m00, mul(q, M.m10));
    R.m11 = sub(M.m01, mul(q, M.m11));
    return R;
}

// ---------------------------------------------------------------------------
// Half-GCD and partial GCD
// ---------------------------------------------------------------------------
PolyMatrix FastPolyArith::halfGcd(Poly a, Poly b) const {
    trim(a);
    trim(b);
    const int da = degree(a);
    const int m = (da + 1) / 2;  // ceil(deg a / 2)
    if (degree(b) < m) return PolyMatrix::identity();

    if (static_cast<std::size_t>(da) <= hgcdThreshold_) {  // small: classical Euclid
        PolyMatrix M = PolyMatrix::identity();
        while (degree(b) >= m) {
            Poly r;
            Poly q = divmod(a, b, r);
            M = euclidStep(M, q);
            a = std::move(b);
            b = std::move(r);
        }
        return M;
    }

    // 1) Reduce using only the upper half of the coefficients.
    PolyMatrix M1 = halfGcd(shiftDown(a, static_cast<std::size_t>(m)), shiftDown(b, static_cast<std::size_t>(m)));
    apply(M1, a, b);
    if (degree(b) < m) return M1;

    // 2) One explicit Euclid step.
    Poly r;
    Poly q = divmod(a, b, r);
    M1 = euclidStep(M1, q);
    a = std::move(b);
    b = std::move(r);
    if (degree(b) < m) return M1;

    // 3) Second recursive call, again on the relevant upper coefficients.
    const int k = 2 * m - degree(a);  // >= 0
    PolyMatrix M2 = halfGcd(shiftDown(a, static_cast<std::size_t>(k)), shiftDown(b, static_cast<std::size_t>(k)));
    return mul(M2, M1);
}

PartialGcdResult FastPolyArith::partialGcd(Poly a, Poly b, int s) const {
    trim(a);
    trim(b);
    PartialGcdResult res;
    res.M = PolyMatrix::identity();
    if (degree(b) > degree(a)) {  // ensure deg a >= deg b
        std::swap(a, b);
        res.M = PolyMatrix{{}, {1}, {1}, {}};
    }

    while (degree(b) >= s) {
        const int da = degree(a);
        const int h = std::max(0, 2 * s - da);  // truncation so that halfGcd stops near degree s
        bool advanced = false;

        if (da - h > static_cast<int>(hgcdThreshold_)) {
            const PolyMatrix Mh = halfGcd(shiftDown(a, static_cast<std::size_t>(h)),
                                          shiftDown(b, static_cast<std::size_t>(h)));
            if (!Mh.isIdentity()) {
                Poly a2 = a, b2 = b;
                apply(Mh, a2, b2);
                // Accept only if we did not jump past the first remainder below s.
                if (degree(a2) >= s) {
                    a = std::move(a2);
                    b = std::move(b2);
                    res.M = mul(Mh, res.M);
                    advanced = true;
                }
            }
        }
        if (!advanced) {  // single classical step (never overshoots)
            Poly r;
            Poly q = divmod(a, b, r);
            res.M = euclidStep(res.M, q);
            a = std::move(b);
            b = std::move(r);
        }
    }
    res.r0 = std::move(a);
    res.r1 = std::move(b);
    return res;
}

}  // namespace rs
