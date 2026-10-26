#include "rs/decoders/registry.hpp"

#include "rs/decoders/fast_gcd.hpp"
#include "rs/decoders/welch_berlekamp.hpp"
// #include "rs/decoders/berlekamp_massey.hpp"

namespace rs {

std::vector<DecoderPtr> makeAllDecoders(const RSCode& code) {
    std::vector<DecoderPtr> v;
    v.push_back(std::make_unique<WelchBerlekampDecoder>(code));
    v.push_back(std::make_unique<FastGcdDecoder>(code));
    // v.push_back(std::make_unique<BerlekampMasseyDecoder>(code));
    return v;
}

}  // namespace rs
