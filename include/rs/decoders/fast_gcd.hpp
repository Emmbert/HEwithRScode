#pragma once
// Gao decoder with fast partial GCD (Half-GCD), O(n log^2 n), hard decision,
// bounded distance up to t.
//
// Algorithm (S. Gao, "A new algorithm for decoding Reed-Solomon codes", 2003):
//   1. g0(x) = prod_i (x - x_i). For the negacyclic code the locators x_i = psi^(2i+1)
//      are exactly the roots of x^n + 1, so g0 = x^n + 1 (no product tree needed).
//   2. g1(x) = interpolation polynomial of the received word, g1(x_i) = y_i:
//      the inverse negacyclic NTT, O(n log n).
//   3. Partial extended Euclid on (g0, g1), stopped at the first remainder
//      g = u*g0 + v*g1 with deg g < (n + k) / 2. Computed with Half-GCD, O(n log^2 n).
//   4. m = g / v. Success if the division is exact, deg m < k and the re-encoded
//      codeword is within distance t of the received word.
// Here v is (up to a scalar) the error locator polynomial.
#include "rs/core/poly_fast.hpp"
#include "rs/decoder.hpp"

namespace rs {

class FastGcdDecoder final : public Decoder {
public:
    explicit FastGcdDecoder(const RSCode& code);
    std::string name() const override { return "Gao (fast Half-GCD)"; }
    DecodeResult decode(const ReceivedWord& rx) const override;

private:
    FastPolyArith arith_;
    Poly g0_;  // x^n + 1
};

}  // namespace rs
