#pragma once

// Parameter sweep: for each HE parameter set, run multiplication chains until the
// noise exceeds the decryption bound, and write per-parameter-set CSVs plus an index
// for plot_param_sweep.py.
//
// Output in SweepConfig::outDir:
//   sweep_index.csv          one row per parameter set (requested + actual parameters, files)
//   pXX_growth.csv           per-step summary (same format as writeOpChainCSV)
//   pXX_hist.csv             per-step histograms of log2|e_i| (same format as writeNoiseHistogramCSV)

#include <string>
#include <vector>

#include "noise_simulations/crypto_context.h"
#include "noise_simulations/noise_growth_experiment.h"

struct HEParams {
    std::string label;               // free text for tables and plots; auto-generated if empty
    uint32_t ringDim = 0;            // n; 0 = let OpenFHE choose (needs enforce128BitSecurity)
    uint32_t multDepth = 10;         // number of levels -> together with the mod sizes determines q
    uint32_t firstModSize = 0;       // bits of the first RNS modulus; 0 = OpenFHE default
    uint32_t scalingModSize = 0;     // bits of the other RNS moduli; 0 = OpenFHE default
    double sigma = 3.19;             // standard deviation of the error distribution
    uint32_t plaintextModulus = 65537;  // p (= t); SIMD packing needs a prime with p = 1 mod 2n
    bool enforce128BitSecurity = true;  // false: HEStd_NotSet, any n allowed (simulation only!)
    uint32_t numLargeDigits = 0;     // key-switching digits (HYBRID); 0 = OpenFHE default.
                                     // More digits -> smaller auxiliary modulus P -> a larger q
                                     // passes the security check for the same n.
};

// Reads parameter sets from a CSV with a header line, e.g. written by find_candidates.py.
// Known columns (any order, missing ones keep their defaults, unknown ones are ignored):
//   label, ringDim, multDepth, firstModSize, scalingModSize, sigma, plaintextModulus,
//   secure (1/0), numLargeDigits
std::vector<HEParams> loadHEParamsCSV(const std::string& path);

// Builds a context like generateContext(), but from explicit parameters.
// Throws (OpenFHE exception) for invalid or insecure combinations.
Context generateContext(const HEParams& p);

struct SweepConfig {
    size_t numChains = 10;           // independent multiplication chains per parameter set
    size_t n = 16;                   // message length
    size_t maxMults = 40;            // upper bound on the chain length
    bool modReduceAfterMult = false; // ModReduce after every multiplication
    size_t stepsAfterFailure = 1;    // continue this many multiplications after the first error
    double histBinWidth = 0.25;      // histogram resolution in bits
    std::string outDir = "sweep";
};

void runParameterSweep(const std::vector<HEParams>& paramSets, const SweepConfig& cfg);