// Correctness tests for the Gao decoder with fast Half-GCD.
// Usage: test_fast_gcd [--large]   (--large adds n = 4096 codes)
#include <cstring>
#include <random>

#include "rs/channel.hpp"
#include "rs/decoders/fast_gcd.hpp"
#include "rs/decoders/welch_berlekamp.hpp"
#include "test_util.hpp"

int main(int argc, char** argv) {
    using namespace rs;
    const bool large = argc > 1 && std::strcmp(argv[1], "--large") == 0;

    struct Case { std::uint32_t p; std::size_t n, k, trials; };
    std::vector<Case> cases = {
        {0, 8, 4, 100},    {0, 8, 1, 50},    {0, 8, 8, 10},   // k = 1 and t = 0 edge cases
        {0, 16, 8, 50},    {0, 16, 15, 50},  {0, 256, 224, 20},
        {0, 256, 128, 10}, {0, 256, 2, 10},  {0, 1024, 512, 3},
    };
    if (large) {
        cases.push_back({0, 4096, 2048, 2});
        cases.push_back({0, 4096, 3584, 2});
        cases.push_back({65537, 4096, 3072, 2});
    }

    for (const Case& c : cases) {
        const std::uint32_t p = c.p ? c.p : RSCode::smallestNegacyclicPrime(c.n);
        RSCode code(p, c.n, c.k);
        FastGcdDecoder dec(code);
        rs::test::checkDecoder(dec, c.trials, 4321 + c.n + c.k);
    }

    // Cross-check: Gao and Welch-Berlekamp must agree on every received word,
    // also beyond the radius (both are bounded-distance decoders).
    {
        RSCode code(RSCode::smallestNegacyclicPrime(256), 256, 128);
        FastGcdDecoder gao(code);
        WelchBerlekampDecoder wb(code);
        std::mt19937_64 rng(5);
        for (std::size_t e = 0; e <= 2 * code.t(); e += 7) {
            const Codeword cw = code.encode(code.randomMessage(rng));
            const ReceivedWord rx = symbolErrorChannel(code, cw, e, rng);
            const DecodeResult a = gao.decode(rx), b = wb.decode(rx);
            RS_CHECK(a.status == b.status);
            if (a.status == DecodeStatus::Success && b.status == DecodeStatus::Success)
                RS_CHECK(a.codeword == b.codeword);
        }
    }
    return rs::test::finish();
}
