#pragma once
// Welch-Berlekamp decoder (hard decision, bounded distance up to t) for the
// negacyclic RS code over GF(p).
//
// Key equation (no syndromes needed, works directly on the evaluation form):
//   find E(x) (error locator, deg E <= t) and N(x) (deg N < k + deg E) with
//   N(x_i) = y_i * E(x_i) for all locators x_i = psi^(2i+1), then m = N / E.
//
// Implementation: O(n^2) point-by-point rational interpolation (Welch-Berlekamp /
// Koetter, i.e. Guruswami-Sudan with multiplicity 1 and list size 1). We look for
// Q(x, y) = Q0(x) + Q1(x)*y with Q(x_i, y_i) = 0 and minimal (1, k-1)-weighted
// degree; then m = -Q0 / Q1.
// Nothing here relies on the locators being roots of x^n + 1, so the decoder works
// for any evaluation points; the negacyclic structure is only exploited for
// O(n log n) re-encoding via the NTT.
#include "rs/decoder.hpp"

namespace rs {

class WelchBerlekampDecoder final : public Decoder {
public:
    using Decoder::Decoder;
    std::string name() const override { return "Welch-Berlekamp"; }
    DecodeResult decode(const ReceivedWord& rx) const override;
};

}  // namespace rs
