#pragma once
// Minimal test helpers (no external framework). A test executable returns 0 on
// success, so CTest can run it directly.
#include <cstdlib>
#include <iostream>
#include <string>

#include "rs/harness/decoder_test.hpp"

namespace rs::test {

inline int g_failures = 0;

#define RS_CHECK(cond)                                                                 \
    do {                                                                               \
        if (!(cond)) {                                                                 \
            ++rs::test::g_failures;                                                    \
            std::cerr << __FILE__ << ":" << __LINE__ << ": CHECK failed: " #cond "\n"; \
        }                                                                              \
    } while (0)

// Runs the generic decoder test (rs::testDecoder) and records the result.
inline void checkDecoder(const Decoder& dec, std::size_t trials, std::uint64_t seed) {
    const DecoderTestReport rep = testDecoder(dec, trials, seed, &std::cerr);
    std::cout << "[" << (rep.passed() ? "PASS" : "FAIL") << "] " << dec.name() << " on "
              << dec.code().describe() << " (" << rep.trials << " trials)\n";
    if (!rep.passed()) ++g_failures;
}

inline int finish() {
    if (g_failures) std::cerr << g_failures << " check(s) failed\n";
    return g_failures ? EXIT_FAILURE : EXIT_SUCCESS;
}

}  // namespace rs::test
