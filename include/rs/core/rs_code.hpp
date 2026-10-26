#pragma once
// Negacyclic Reed-Solomon code over GF(p).
//
//   Ring:      R_p = Z_p[x] / (x^n + 1),  n = 2^d,  p prime with 2n | p - 1
//   psi:       primitive 2n-th root of unity in Z_p  (psi^n = -1)
//   Locators:  x_i = psi^(2i+1), i = 0..n-1  -> exactly the n roots of x^n + 1,
//              i.e. the primitive 2n-th roots of unity
//   Message:   m(x) in R_p with deg m < k  (k coefficients)
//   Codeword:  c_i = m(x_i)  ->  c = NTT_psi(m), the negacyclic NTT
//
// Consequences of this choice:
//   * encode() is a negacyclic NTT, O(n log n).
//   * evaluate(a * b mod x^n+1) == evaluate(a) .* evaluate(b) (pointwise), so ring
//     products map to pointwise products of evaluation vectors.
//   * The code is an evaluation (GRS) code, so all evaluation-form decoders
//     (Welch-Berlekamp, Gao, Guruswami-Sudan, Koetter-Vardy) apply directly.
#include <cstdint>
#include <random>
#include <string>
#include <vector>
#include "rs/core/field.hpp"
#include "rs/core/poly.hpp"

namespace rs {

using Message = std::vector<Elem>;   // k coefficients m_0..m_{k-1}
using Codeword = std::vector<Elem>;  // n evaluations

class RSCode {
public:
    // psi = 0 selects a primitive 2n-th root of unity automatically.
    RSCode(std::uint32_t p, std::size_t n, std::size_t k, Elem psi = 0);

    // Smallest prime p = 1 (mod 2n), e.g. 40961 for n = 4096.
    static std::uint32_t smallestNegacyclicPrime(std::size_t n);

    const PrimeField& field() const { return field_; }
    std::size_t n() const { return n_; }
    std::size_t k() const { return k_; }
    std::size_t redundancy() const { return n_ - k_; }
    std::size_t d() const { return n_ - k_ + 1; }
    std::size_t t() const { return (n_ - k_) / 2; }
    double rate() const { return static_cast<double>(k_) / static_cast<double>(n_); }
    Elem psi() const { return psi_; }
    const std::vector<Elem>& evalPoints() const { return points_; }  // x_i = psi^(2i+1)

    Codeword encode(const Message& msg) const;  // msg.size() == k
    Codeword encodePoly(const Poly& f) const;   // deg f < k
    // Negacyclic NTT of an arbitrary polynomial; it is first reduced mod x^n + 1.
    Codeword evaluate(const Poly& f) const;
    // Inverse negacyclic NTT: the ring element a (size n) with evaluate(a) == c.
    Poly interpolate(const Codeword& c) const;
    bool isCodeword(const Codeword& c) const { return degree(interpolate(c)) < static_cast<int>(k_); }

    Message randomMessage(std::mt19937_64& rng) const;
    std::string describe() const;

private:
    void cyclicNtt(std::vector<Elem>& a, const std::vector<Elem>& twiddles) const;

    PrimeField field_;
    std::size_t n_, k_;
    Elem psi_;
    Elem nInv_;
    std::vector<Elem> points_;
    std::vector<Elem> psiPow_, psiInvPow_;      // psi^j, psi^-j, j < n
    std::vector<Elem> omegaPow_, omegaInvPow_;  // omega = psi^2, powers j < n/2
};

std::size_t hammingDistance(const std::vector<Elem>& a, const std::vector<Elem>& b);

}  // namespace rs
