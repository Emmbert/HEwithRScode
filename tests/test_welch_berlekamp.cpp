// Correctness tests for the Welch-Berlekamp decoder.
// Usage: test_welch_berlekamp [--large]   (--large adds n = 4096 codes)
#include <cstring>

#include "rs/decoders/welch_berlekamp.hpp"
#include "test_util.hpp"

int main(int argc, char** argv) {
    using namespace rs;
    const bool large = argc > 1 && std::strcmp(argv[1], "--large") == 0;

    struct Case { std::uint32_t p; std::size_t n, k, trials; };
    std::vector<Case> cases = {
        {0, 8, 4, 100},    {0, 8, 1, 50},    {0, 8, 8, 10},   // k = 1 and t = 0 edge cases
        {0, 16, 8, 50},    {0, 16, 15, 50},  {0, 256, 224, 20},
        {0, 256, 128, 10}, {0, 1024, 512, 2},
    };
    if (large) {
        cases.push_back({0, 4096, 2048, 1});
        cases.push_back({0, 4096, 3584, 1});
        cases.push_back({65537, 4096, 3072, 1});
    }

    for (const Case& c : cases) {
        const std::uint32_t p = c.p ? c.p : RSCode::smallestNegacyclicPrime(c.n);
        RSCode code(p, c.n, c.k);
        WelchBerlekampDecoder dec(code);
        rs::test::checkDecoder(dec, c.trials, 1234 + c.n + c.k);
    }
    return rs::test::finish();
}
