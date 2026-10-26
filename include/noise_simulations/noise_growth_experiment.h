#pragma once

// Noise growth under chains of homomorphic additions and multiplications.
//
// Noise is measured as log2 values (doubles), so it stays correct when |e| grows
// far beyond 2^63, which happens after one or two multiplications.
//
// Wrapped coefficients: after decryption every coefficient is reduced mod Q, so a
// *measured* |e_i| can never exceed the threshold Q/(2t). When the true noise of a
// coefficient crosses it, the coefficient wraps around and decrypts to a wrong value
// mod t. measureNoise() detects exactly these coefficients (v_i mod t != expected
// plaintext coefficient) and counts them in numWrapped instead of reporting a bogus
// noise value. Works for COEFFICIENT and SIMD packing.

#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "noise_simulations/crypto_context.h"
#include "noise_simulations/encryption.h"

enum class HeOp { Add, Mult };

struct NoiseMeasurement {
    double log2MaxNoise  = 0.0;  // log2(max_i |e_i|) over non-wrapped coefficients; >= log2Threshold if any wrapped
    double log2RmsNoise  = 0.0;  // log2(sqrt(mean_i e_i^2)) over non-wrapped coefficients
    double log2Threshold = 0.0;  // log2(Q_l / (2t)) at the ciphertext's current modulus
    size_t numTowers     = 0;    // RNS towers (shrinks after ModReduce)
    size_t numCoeffs     = 0;    // N
    size_t numWrapped    = 0;    // coefficients whose noise exceeded the threshold
    size_t numWrappedInMessage = 0;  // of those, positions i < message length (COEFFICIENT packing)
    size_t numZero       = 0;    // coefficients with e_i == 0 (only possible for tiny noise)

    // Optional histogram of log2|e_i| over the non-wrapped, non-zero coefficients:
    // key k counts coefficients with log2|e_i| in [k*histBinWidth, (k+1)*histBinWidth).
    double histBinWidth = 0.0;
    std::map<int, size_t> log2Histogram;
};

// histBinWidth > 0 also fills the histogram (in bits, e.g. 0.25).
NoiseMeasurement measureNoise(const Context& ctx, const EncryptedSample& sample, double histBinWidth = 0.0);

// Homomorphic addition with the matching expected plaintext (counterpart of multiplySamples).
EncryptedSample addSamples(const Context& ctx, const EncryptedSample& a, const EncryptedSample& b);

struct OpChainConfig {
    std::string name;                // label in the CSVs, e.g. "mult", "add"
    std::vector<HeOp> ops;           // step k combines acc with a fresh ciphertext
    size_t numChains = 10;           // independent repetitions
    size_t n = 16;                   // length of each random message
    bool modReduceAfterMult = false; // call ModReduce after every multiplication
    size_t measureStride = 1;        // measure every k steps (failures and the last step always)
    bool stopAtFailure = true;       // stop at the first error (decryption failure or any wrapped
                                     // coefficient); false: keep going (see stepsAfterFailure)
    double histBinWidth = 0.25;      // histogram resolution in bits; 0 disables histograms
    size_t stepsAfterFailure = SIZE_MAX;  // with stopAtFailure == false: how many steps to continue
                                          // after the first error (default: until the end)
};

struct OpChainRow {
    std::string experiment;
    size_t chain = 0;
    size_t step = 0;                 // 0 = fresh ciphertext, k = after k operations
    std::string op;                  // "fresh", "add" or "mult"
    size_t numAdds = 0;
    size_t numMults = 0;
    NoiseMeasurement noise;
    bool correct = true;             // OpenFHE's Decrypt returned the expected message
};

std::vector<HeOp> repeatOp(HeOp op, size_t count);

// numMults blocks of: one multiplication followed by addsPerMult additions.
std::vector<HeOp> mixedOps(size_t numMults, size_t addsPerMult);

std::vector<OpChainRow> runOpChainExperiment(const Context& ctx, const OpChainConfig& cfg);

// One row per measurement (summary statistics) -> plot_noise_growth.py
void writeOpChainCSV(const std::string& path, const Context& ctx, const std::vector<OpChainRow>& rows);

// One row per (measurement, histogram bin), plus "zero" and "wrapped" pseudo-bins
// -> plot_noise_histograms.py
void writeNoiseHistogramCSV(const std::string& path, const Context& ctx, const std::vector<OpChainRow>& rows);