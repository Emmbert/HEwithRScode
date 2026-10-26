#pragma once
// Central list of all decoders. Adding a decoder = one line in makeAllDecoders().
#include <vector>
#include "rs/decoder.hpp"

namespace rs {

std::vector<DecoderPtr> makeAllDecoders(const RSCode& code);

}  // namespace rs
