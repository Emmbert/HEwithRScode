#pragma once

#include "encryption.h"
#include <string>

// ---------------------------------------------------------------------
// Aggregated result of running many independent encrypt/decrypt round
// trips and collecting every reconstructed noise coefficient.
// ---------------------------------------------------------------------
struct NoiseExperimentResult {
    std::vector<int64_t> allNoiseValues;   // flattened noise samples, all trials & coefficients
    size_t numTrials = 0;
    size_t numDecryptionFailures = 0;      // trials where isCorrect() was false
    size_t numOverThreshold = 0;           // trials where some |noise| > maxNoiseThreshold
    size_t numInconsistentTrials = 0;      // trials where the two checks above disagreed
};

bool exceedsThreshold(const std::vector<int64_t>& noise, double maxNoiseThreshold);

// Runs `numTrials` independent encrypt/decrypt round trips of `n` packed
// values each. Set useSymmetric=true to encrypt with the secret key
// (isolates a single fresh noise term, directly comparable to sigma);
// false (default) uses public-key encryption, matching real-world usage
// but with larger, convolution-amplified noise (see selfTestNoiseReconstruction
// for the derivation). Collects every reconstructed noise coefficient, and
// cross-checks two independent failure indicators against each other:
//   (a) isCorrect() on the decrypted vs. original message
//   (b) whether any reconstructed |noise| exceeded computeMaxNoiseThreshold()
// These should always agree; numInconsistentTrials counts cases where they
// don't, which would indicate a bug in decryptWithNoise() or in the
// threshold formula itself.
NoiseExperimentResult runNoiseExperiment(const Context& ctx, size_t numTrials, size_t n, bool useSymmetric = false);

// Writes the collected noise samples to a CSV file for external plotting.
// Sigma and the noise threshold are written as leading "#" comment lines
// (Python's csv/pandas readers can skip these) so the plot script has the
// theoretical parameters without needing to recompute them.
void writeNoiseCSV(const std::string& path, const NoiseExperimentResult& result, double sigma, double maxNoiseThreshold);

// Theoretical stddev of the *total* noise in a fresh public-key ciphertext,
// i.e. the noise selfTestNoiseReconstruction(ctx, /*symmetric=*/false)
// should match. Unlike symmetric encryption (a single fresh e ~ N(0,sigma)),
// public-key encryption combines an ephemeral secret u, the key-generation
// noise e0, and the secret key s via ring multiplication:
//   (v - m)/t = u*e0 + e1 + e2*s
// Using the standard coefficient-wise noise-growth heuristic
// Var(a*b) = N * Var(a) * Var(b) for independent, zero-mean length-N
// polynomials, and assuming u and s share the same secret key distribution:
//   Var(total) = sigma^2 * (1 + 2*N*Var(s))
// Returns -1.0 if the secret key distribution isn't one we have a variance
// formula for (rather than silently guessing).
double computeExpectedPublicKeyNoiseStdDev(const Context& ctx);

// Lightweight, dependency-free self-check: runs a moderate-size experiment
// and verifies the empirical standard deviation of the reconstructed noise
// is close to the theoretical value -- sigma directly if symmetric=true,
// or computeExpectedPublicKeyNoiseStdDev() if symmetric=false (the default,
// matching real-world public-key usage). Prints PASS/FAIL with details and
// returns the result.
//
// This lives in its own function specifically so it can be called from
// main() today, and from a real unit-test runner (GoogleTest/Catch2) later,
// without any changes.
bool selfTestNoiseReconstruction(const Context& ctx, bool symmetric = false);