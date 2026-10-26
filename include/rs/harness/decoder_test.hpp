#pragma once
// Generic correctness test usable by every decoder.
//   e <= t : decoder must return exactly the transmitted codeword and message.
//   e >  t : if the decoder claims success, the result must be consistent
//            (re-encoding the message gives the codeword) and lie within
//            decodingRadius() of the received word. Miscorrections are allowed.
#include <cstdint>
#include <ostream>
#include <string>
#include <vector>
#include "rs/decoder.hpp"

namespace rs {

struct DecoderTestReport {
    std::size_t trials = 0;
    std::size_t failures = 0;
    std::string firstError;
    bool passed() const { return failures == 0; }
};

std::vector<std::size_t> defaultTestErrorCounts(const RSCode& code);

DecoderTestReport testDecoder(const Decoder& dec, std::size_t trialsPerErrorCount,
                              std::uint64_t seed, std::ostream* log = nullptr);

}  // namespace rs
