#include "rs/harness/decoder_test.hpp"

#include <set>
#include <sstream>

#include "rs/channel.hpp"

namespace rs {

std::vector<std::size_t> defaultTestErrorCounts(const RSCode& code) {
    const std::size_t t = code.t(), n = code.n();
    std::set<std::size_t> s;
    if (t <= 16) {
        for (std::size_t e = 0; e <= t; ++e) s.insert(e);
    } else {
        for (std::size_t e : {std::size_t{0}, std::size_t{1}, std::size_t{2}, t / 4, t / 2, 3 * t / 4, t - 1, t})
            s.insert(e);
    }
    for (std::size_t e : {t + 1, t + 1 + t / 2, 2 * t + 1, n})  // beyond the radius
        if (e <= n) s.insert(e);
    return {s.begin(), s.end()};
}

DecoderTestReport testDecoder(const Decoder& dec, std::size_t trialsPerErrorCount,
                              std::uint64_t seed, std::ostream* log) {
    const RSCode& code = dec.code();
    const std::size_t t = code.t();
    std::mt19937_64 rng(seed);
    DecoderTestReport rep;

    auto fail = [&](std::size_t e, const std::string& what) {
        ++rep.failures;
        std::ostringstream os;
        os << dec.name() << " on " << code.describe() << ", e = " << e << ": " << what;
        if (rep.firstError.empty()) rep.firstError = os.str();
        if (log) *log << "  FAIL: " << os.str() << "\n";
    };

    for (std::size_t e : defaultTestErrorCounts(code)) {
        for (std::size_t trial = 0; trial < trialsPerErrorCount; ++trial) {
            ++rep.trials;
            const Message msg = code.randomMessage(rng);
            const Codeword cw = code.encode(msg);
            const ReceivedWord rx = symbolErrorChannel(code, cw, e, rng);

            DecodeResult r;
            try {
                r = dec.decode(rx);
            } catch (const std::exception& ex) {
                fail(e, std::string("exception: ") + ex.what());
                continue;
            }

            if (e <= t) {
                if (r.status != DecodeStatus::Success) fail(e, "reported failure within radius t");
                else if (r.codeword != cw) fail(e, "wrong codeword within radius t");
                else if (r.message != msg) fail(e, "wrong message within radius t");
                else if (r.numCorrected != e) fail(e, "numCorrected mismatch");
            } else if (r.status == DecodeStatus::Success) {
                if (r.codeword.size() != code.n() || r.message.size() != code.k() ||
                    code.encode(r.message) != r.codeword)
                    fail(e, "returned codeword/message inconsistent");
                else if (hammingDistance(r.codeword, rx.hard) > dec.decodingRadius())
                    fail(e, "returned codeword outside decoding radius");
            }
        }
    }
    return rep;
}

}  // namespace rs
