#pragma once
// Channel models over Z_p. Every decoder receives a ReceivedWord: hard decisions
// are always present; soft decoders may use the fractional part of the observation.
#include <random>
#include <vector>
#include "rs/core/rs_code.hpp"

namespace rs {

struct ReceivedWord {
    std::vector<Elem> hard;   // n hard-decision symbols in [0, p)
    // Optional soft information (Gaussian channel): the real observation is
    // y_i = hard_i + frac_i (mod p), frac_i in [-0.5, 0.5], noise std. deviation sigma.
    std::vector<float> frac;
    double sigma = 0.0;       // 0 = no soft information
};

// Adds exactly numErrors errors at random positions with uniformly random nonzero values.
ReceivedWord symbolErrorChannel(const RSCode& code, const Codeword& cw,
                                std::size_t numErrors, std::mt19937_64& rng);

// Additive Gaussian noise on every symbol, y_i = c_i + N(0, sigma^2), hard = round(y) mod p.
// Models small coefficient noise as in RLWE-type rings; a symbol is wrong once |noise| >= 0.5.
ReceivedWord gaussianChannel(const RSCode& code, const Codeword& cw,
                             double sigma, std::mt19937_64& rng);

}  // namespace rs
