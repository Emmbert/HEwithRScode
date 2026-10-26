#include "noise_simulations/noise_growth_experiment.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>

using namespace lbcrypto;

namespace {

int64_t reduceBalanced(int64_t x, int64_t t) {
    int64_t r = x % t;
    if (r < 0) r += t;
    if (r > t / 2) r -= t;  // same convention as encryption.cpp
    return r;
}

// log2 of an arbitrarily large non-negative integer, without overflowing a double.
template <typename BigInt>
double log2Big(const BigInt& x) {
    const auto msb = x.GetMSB();  // 1-based position of the highest set bit, 0 for x == 0
    if (msb == 0) return -std::numeric_limits<double>::infinity();
    if (msb <= 63) return std::log2(static_cast<double>(x.template ConvertToInt<uint64_t>()));
    const uint32_t shift = msb - 63;
    return std::log2(static_cast<double>((x >> shift).template ConvertToInt<uint64_t>())) + shift;
}

// The plaintext polynomial p (coefficient domain, balanced mod t) that the ciphertext
// encrypts, i.e. v = c0 + c1*s = sf*p + t*e (mod Q).
//   COEFFICIENT packing: p is the message itself, zero-padded to N.
//   SIMD packing:        p is the inverse NTT (mod t) of the slot values, so every
//                        coefficient of p depends on all slots -- and vice versa, one
//                        wrapped coefficient corrupts essentially all slots.
std::vector<int64_t> expectedPlaintextPoly(const Context& ctx, const std::vector<int64_t>& message) {
    const int64_t t = ctx.cc->GetCryptoParameters()->GetPlaintextModulus();
    const size_t N = ctx.cc->GetRingDimension();
    std::vector<int64_t> p(N, 0);

    if (PACKING_MODE == PackingMode::COEFFICIENT) {
        for (size_t i = 0; i < std::min(N, message.size()); ++i) p[i] = reduceBalanced(message[i], t);
    } else {
        Plaintext pt = ctx.cc->MakePackedPlaintext(message);
        DCRTPoly encoded = pt->GetElement<DCRTPoly>();  // copy; encoded mod each q_i, values < t
        encoded.SetFormat(Format::COEFFICIENT);
        const auto& limb0 = encoded.GetElementAtIndex(0).GetValues();
        for (size_t i = 0; i < N; ++i) p[i] = reduceBalanced(static_cast<int64_t>(limb0[i].ConvertToInt()), t);
    }
    return p;
}

bool decryptsCorrectly(const Context& ctx, const EncryptedSample& s) {
    try {
        return isCorrect(s.original, decryptMessage(ctx, s.ciphertext, s.original.size()));
    } catch (const std::exception&) {
        return false;  // a badly broken ciphertext may make decoding throw
    }
}

}  // namespace

NoiseMeasurement measureNoise(const Context& ctx, const EncryptedSample& sample, double histBinWidth) {
    const auto& cv = sample.ciphertext->GetElements();
    const int64_t t = ctx.cc->GetCryptoParameters()->GetPlaintextModulus();
    const size_t N = ctx.cc->GetRingDimension();
    const double log2t = std::log2(static_cast<double>(t));

    // The secret key lives at the full modulus Q; after ModReduce the ciphertext
    // only has the first few RNS towers, so drop the matching towers of s.
    DCRTPoly s = ctx.keyPair.secretKey->GetPrivateElement();
    const size_t towers = cv[0].GetNumOfElements();
    if (s.GetNumOfElements() > towers) s.DropLastElements(s.GetNumOfElements() - towers);

    // v = c0 + c1*s + c2*s^2 + ...  (more than 2 elements only without relinearization)
    DCRTPoly v = cv[0];
    DCRTPoly sPow = s;
    for (size_t j = 1; j < cv.size(); ++j) {
        v += cv[j] * sPow;
        if (j + 1 < cv.size()) sPow = sPow * s;
    }
    v.SetFormat(Format::COEFFICIENT);

    const auto big = v.CRTInterpolate();
    const auto Q = big.GetModulus();
    const auto halfQ = Q >> 1;
    const decltype(Q) tBig(static_cast<uint64_t>(t));
    const auto& coeffs = big.GetValues();

    // sf stays 1 under FIXEDMANUAL (moduli are 1 mod t); FLEXIBLEAUTO* change it.
    int64_t sf = static_cast<int64_t>(sample.ciphertext->GetScalingFactorInt().ConvertToInt());
    if (sf == 0) sf = 1;
    std::vector<int64_t> p = expectedPlaintextPoly(ctx, sample.original);
    for (auto& pi : p) pi = reduceBalanced(pi * sf, t);

    NoiseMeasurement out;
    out.numTowers = towers;
    out.numCoeffs = N;
    out.histBinWidth = histBinWidth;
    out.log2Threshold = log2Big(Q) - 1.0 - log2t;

    auto countWrapped = [&](size_t i) {
        out.numWrapped++;
        if (i < sample.original.size()) out.numWrappedInMessage++;
    };

    const double negInf = -std::numeric_limits<double>::infinity();
    double maxLog = negInf;
    std::vector<double> log2e;
    log2e.reserve(N);

    for (size_t i = 0; i < N; ++i) {
        const auto& c = coeffs[i];
        const bool negative = c > halfQ;
        const auto mag = negative ? (Q - c) : c;  // |v_i| in the centered representation

        double l;
        if (mag.GetMSB() < 62) {
            // Small values: exact integer arithmetic.
            int64_t vi = static_cast<int64_t>(mag.template ConvertToInt<uint64_t>());
            if (negative) vi = -vi;
            const int64_t diff = vi - p[i];
            if (diff % t != 0) { countWrapped(i); continue; }
            const int64_t e = diff / t;
            if (e == 0) { out.numZero++; continue; }
            l = std::log2(std::abs(static_cast<double>(e)));
        } else {
            // Large values: check v_i == p_i (mod t) exactly, then |e| ~ |v| / t
            // (|p_i| <= t/2 is negligible next to |v_i| >= 2^61).
            int64_t r = static_cast<int64_t>(mag.Mod(tBig).template ConvertToInt<uint64_t>());
            if (negative) r = -r;
            if (reduceBalanced(r - p[i], t) != 0) { countWrapped(i); continue; }
            l = log2Big(mag) - log2t;
        }

        log2e.push_back(l);
        maxLog = std::max(maxLog, l);
        if (histBinWidth > 0.0) out.log2Histogram[static_cast<int>(std::floor(l / histBinWidth))]++;
    }

    if (!log2e.empty() && !std::isinf(maxLog)) {
        // RMS in the log domain without overflow (zeros contribute 0 to the sum).
        double sum = 0.0;
        for (double l : log2e) sum += std::exp2(2.0 * (l - maxLog));
        const double nonWrapped = static_cast<double>(log2e.size() + out.numZero);
        out.log2MaxNoise = std::max(0.0, maxLog);
        out.log2RmsNoise = std::max(0.0, maxLog + 0.5 * std::log2(sum / nonWrapped));
    }
    // A wrapped coefficient's true noise is at least the threshold.
    if (out.numWrapped > 0) out.log2MaxNoise = std::max(out.log2MaxNoise, out.log2Threshold);
    return out;
}

EncryptedSample addSamples(const Context& ctx, const EncryptedSample& a, const EncryptedSample& b) {
    const int64_t t = ctx.cc->GetCryptoParameters()->GetPlaintextModulus();
    const size_t len = std::max(a.original.size(), b.original.size());

    std::vector<int64_t> expected(len, 0);
    for (size_t i = 0; i < len; ++i) {
        int64_t x = (i < a.original.size()) ? a.original[i] : 0;
        int64_t y = (i < b.original.size()) ? b.original[i] : 0;
        expected[i] = reduceBalanced(x + y, t);
    }

    auto ciphertext = ctx.cc->EvalAdd(a.ciphertext, b.ciphertext);
    return EncryptedSample{expected, ciphertext};
}

std::vector<HeOp> repeatOp(HeOp op, size_t count) {
    return std::vector<HeOp>(count, op);
}

std::vector<HeOp> mixedOps(size_t numMults, size_t addsPerMult) {
    std::vector<HeOp> ops;
    for (size_t k = 0; k < numMults; ++k) {
        ops.push_back(HeOp::Mult);
        for (size_t a = 0; a < addsPerMult; ++a) ops.push_back(HeOp::Add);
    }
    return ops;
}

std::vector<OpChainRow> runOpChainExperiment(const Context& ctx, const OpChainConfig& cfg) {
    std::vector<OpChainRow> rows;
    const size_t stride = std::max<size_t>(1, cfg.measureStride);

    for (size_t chain = 0; chain < cfg.numChains; ++chain) {
        EncryptedSample acc = encryptRandomMessage(ctx, cfg.n);
        size_t numAdds = 0, numMults = 0;
        bool failed = false;
        size_t failedAtStep = 0;

        rows.push_back({cfg.name, chain, 0, "fresh", 0, 0, measureNoise(ctx, acc, cfg.histBinWidth),
                        decryptsCorrectly(ctx, acc)});

        for (size_t step = 1; step <= cfg.ops.size(); ++step) {
            const HeOp op = cfg.ops[step - 1];
            // Fresh operand at the top level; OpenFHE aligns levels automatically if acc was ModReduced.
            EncryptedSample fresh = encryptRandomMessage(ctx, cfg.n);

            if (op == HeOp::Mult) {
                acc = multiplySamples(ctx, acc, fresh);
                numMults++;
                if (cfg.modReduceAfterMult && acc.ciphertext->GetElements()[0].GetNumOfElements() > 1) {
                    ctx.cc->ModReduceInPlace(acc.ciphertext);
                }
            } else {
                acc = addSamples(ctx, acc, fresh);
                numAdds++;
            }

            const bool correct = decryptsCorrectly(ctx, acc);
            const bool measure = !correct || step % stride == 0 || step == cfg.ops.size();
            bool error = !correct;
            if (measure) {
                rows.push_back({cfg.name, chain, step, op == HeOp::Mult ? "mult" : "add", numAdds, numMults,
                                measureNoise(ctx, acc, cfg.histBinWidth), correct});
                // A wrapped coefficient is an error even if the n message coefficients
                // happen to decrypt correctly (COEFFICIENT packing).
                error = error || rows.back().noise.numWrapped > 0;
            }
            if (error && !failed) {
                failed = true;
                failedAtStep = step;
            }
            if (failed && (cfg.stopAtFailure || step - failedAtStep >= cfg.stepsAfterFailure)) break;
        }

        std::cout << "[" << cfg.name << "] chain " << (chain + 1) << "/" << cfg.numChains << " done" << std::endl;
    }
    return rows;
}

namespace {

void writeMetadata(std::ofstream& file, const Context& ctx) {
    const auto params = ctx.cc->GetCryptoParameters();
    file << "# ringDimension=" << ctx.cc->GetRingDimension() << "\n";
    file << "# plaintextModulus=" << params->GetPlaintextModulus() << "\n";
    file << "# log2Q=" << params->GetElementParams()->GetModulus().GetMSB() << "\n";
    if (auto rlwe = std::dynamic_pointer_cast<CryptoParametersRLWE<DCRTPoly>>(params))
        file << "# sigma=" << rlwe->GetDistributionParameter() << "\n";
    file << "# packing=" << (PACKING_MODE == PackingMode::SIMD ? "SIMD" : "COEFFICIENT") << "\n";
}

}  // namespace

void writeOpChainCSV(const std::string& path, const Context& ctx, const std::vector<OpChainRow>& rows) {
    std::ofstream file(path);
    if (!file.is_open()) {
        std::cerr << "Error: could not open " << path << " for writing." << std::endl;
        return;
    }

    writeMetadata(file, ctx);
    file << "experiment,chain,step,op,numAdds,numMults,towers,"
            "log2MaxNoise,log2RmsNoise,log2Threshold,wrapped,correct\n";

    size_t suspicious = 0;
    for (const auto& r : rows) {
        file << r.experiment << "," << r.chain << "," << r.step << "," << r.op << "," << r.numAdds << ","
             << r.numMults << "," << r.noise.numTowers << "," << r.noise.log2MaxNoise << ","
             << r.noise.log2RmsNoise << "," << r.noise.log2Threshold << "," << r.noise.numWrapped << ","
             << (r.correct ? 1 : 0) << "\n";
        // A wrapped coefficient corrupts essentially all slots (SIMD packing), or exactly
        // message coefficient i (COEFFICIENT packing; positions >= message length are not
        // part of the message). In those cases decryption must have failed.
        const size_t relevant = (PACKING_MODE == PackingMode::SIMD) ? r.noise.numWrapped : r.noise.numWrappedInMessage;
        if (r.correct && relevant > 0) suspicious++;
    }

    std::cout << "Wrote " << rows.size() << " rows to " << path << std::endl;
    if (suspicious > 0) {
        std::cerr << "Warning: " << suspicious << " rows decrypted correctly although coefficients were "
                  << "flagged as wrapped -- please check the noise reconstruction." << std::endl;
    }
}

void writeNoiseHistogramCSV(const std::string& path, const Context& ctx, const std::vector<OpChainRow>& rows) {
    std::ofstream file(path);
    if (!file.is_open()) {
        std::cerr << "Error: could not open " << path << " for writing." << std::endl;
        return;
    }

    writeMetadata(file, ctx);
    file << "# bin = lower edge of a log2|e_i| bin of width binWidth; 'zero' = e_i == 0; "
            "'wrapped' = noise exceeded the threshold\n";
    file << "experiment,chain,step,op,numAdds,numMults,log2Threshold,binWidth,bin,count\n";

    size_t lines = 0;
    for (const auto& r : rows) {
        const auto& m = r.noise;
        if (m.histBinWidth <= 0.0) continue;
        const std::string prefix = r.experiment + "," + std::to_string(r.chain) + "," + std::to_string(r.step) +
                                   "," + r.op + "," + std::to_string(r.numAdds) + "," +
                                   std::to_string(r.numMults) + "," + std::to_string(m.log2Threshold) + "," +
                                   std::to_string(m.histBinWidth) + ",";
        for (const auto& [bin, count] : m.log2Histogram) {
            file << prefix << (bin * m.histBinWidth) << "," << count << "\n";
            lines++;
        }
        if (m.numZero > 0) { file << prefix << "zero," << m.numZero << "\n"; lines++; }
        if (m.numWrapped > 0) { file << prefix << "wrapped," << m.numWrapped << "\n"; lines++; }
    }
    std::cout << "Wrote " << lines << " histogram rows to " << path << std::endl;
}