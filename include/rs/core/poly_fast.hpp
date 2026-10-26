#pragma once
// Fast polynomial arithmetic over GF(p) for decoders that need O(n log n)
// multiplication and O(n log^2 n) (partial) GCDs.
//
//   mul          NTT-based multiplication, O(n log n). Products longer than the
//                largest NTT length of the field (2-adic part of p-1) are split into
//                blocks, so any p works (e.g. p = 40961 allows length 8192).
//   divmod       Newton iteration on reversed polynomials, O(n log n);
//                classical division for small quotients or divisors.
//   halfGcd      Half-GCD (Knuth / Schoenhage / Thull-Yap), O(n log^2 n).
//   partialGcd   Extended Euclid on (a, b) stopped at the FIRST remainder with
//                degree < s, including the cofactor matrix. This is the key step of
//                the Gao decoder and of other Euclid-based RS decoders.
//
// Conventions: Poly p[i] is the coefficient of x^i; all results are trimmed
// (no leading zeros, zero polynomial = empty vector).
#include <cstddef>
#include <vector>

#include "rs/core/field.hpp"
#include "rs/core/poly.hpp"

namespace rs {

// 2x2 polynomial matrix acting on pairs: (a', b') = M * (a, b).
struct PolyMatrix {
    Poly m00, m01, m10, m11;
    static PolyMatrix identity();
    bool isIdentity() const;
};

struct PartialGcdResult {
    Poly r0, r1;   // consecutive remainders with deg r0 >= s > deg r1
    PolyMatrix M;  // M * (a, b) = (r0, r1); in particular r1 = M.m10 * a + M.m11 * b
};

class FastPolyArith {
public:
    explicit FastPolyArith(const PrimeField& F, std::size_t schoolbookThreshold = 32,
                           std::size_t halfGcdThreshold = 64);

    const PrimeField& field() const { return F_; }
    std::size_t maxNttLength() const { return maxLen_; }

    Poly add(const Poly& a, const Poly& b) const;
    Poly sub(const Poly& a, const Poly& b) const;
    Poly mul(const Poly& a, const Poly& b) const;
    Poly mulSchoolbook(const Poly& a, const Poly& b) const;

    // Power series inverse: b * result = 1 (mod x^m); requires b[0] != 0.
    Poly invSeries(const Poly& b, std::size_t m) const;
    // a = q * b + rem with deg rem < deg b. Returns q.
    Poly divmod(const Poly& a, const Poly& b, Poly& rem) const;

    PolyMatrix mul(const PolyMatrix& A, const PolyMatrix& B) const;  // A * B
    void apply(const PolyMatrix& M, Poly& a, Poly& b) const;          // (a, b) <- M (a, b)

    // Requires deg a >= deg b. Returns M with M (a, b) = (r_j, r_{j+1}),
    // deg r_j >= ceil(deg a / 2) > deg r_{j+1}.
    PolyMatrix halfGcd(Poly a, Poly b) const;
    // Extended Euclid on (a, b) until the first remainder of degree < s.
    PartialGcdResult partialGcd(Poly a, Poly b, int s) const;

private:
    void ntt(std::vector<Elem>& a, bool inverse) const;  // size: power of two <= maxLen_
    PolyMatrix euclidStep(const PolyMatrix& M, const Poly& q) const;  // [[0,1],[1,-q]] * M

    PrimeField F_;
    std::size_t schoolbook_, hgcdThreshold_;
    std::size_t maxLen_ = 1;            // largest power of two dividing p-1 (capped at 2^20)
    std::vector<Elem> roots_, iroots_;  // w^j, w^-j for j < maxLen_/2, w of order maxLen_
};

}  // namespace rs
