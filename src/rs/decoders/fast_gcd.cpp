#include "rs/decoders/fast_gcd.hpp"

#include <stdexcept>

namespace rs {

FastGcdDecoder::FastGcdDecoder(const RSCode& code) : Decoder(code), arith_(code.field()) {
    g0_.assign(code.n() + 1, 0);
    g0_[0] = 1;          // x^n + 1 = prod (x - psi^(2i+1))
    g0_[code.n()] = 1;
}

DecodeResult FastGcdDecoder::decode(const ReceivedWord& rx) const {
    const std::size_t n = code_.n(), k = code_.k();
    if (rx.hard.size() != n) throw std::invalid_argument("FastGcdDecoder: received word has wrong length");

    // Step 2: interpolation polynomial via inverse NTT
    Poly g1 = code_.interpolate(rx.hard);
    trim(g1);

    // Step 3: partial GCD, stop at the first remainder with deg < (n+k)/2
    const int s = static_cast<int>((n + k + 1) / 2);  // deg g < (n+k)/2  <=>  deg g < ceil((n+k)/2)
    const PartialGcdResult pg = arith_.partialGcd(g0_, g1, s);
    const Poly& g = pg.r1;    // g = u*g0 + v*g1
    const Poly& v = pg.M.m11;

    DecodeResult res;
    if (degree(v) < 0) return res;

    // Step 4: m = g / v
    Poly rem;
    Poly m = arith_.divmod(g, v, rem);
    if (!rem.empty()) return res;                      // more than t errors
    if (degree(m) >= static_cast<int>(k)) return res;

    m.resize(k, 0);
    Codeword c = code_.encodePoly(m);  // negacyclic NTT
    const std::size_t dist = hammingDistance(c, rx.hard);
    if (dist > code_.t()) return res;  // bounded-distance decoder

    res.status = DecodeStatus::Success;
    res.codeword = std::move(c);
    res.message = std::move(m);
    res.numCorrected = dist;
    return res;
}

}  // namespace rs
