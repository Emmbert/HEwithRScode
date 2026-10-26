#pragma once

#include "encryption.h"
#include "noise_experiment.h"
#include <string>

// ---------------------------------------------------------------------
// Runs `numTrials` independent single-multiplication trials: encrypts two
// fresh random messages of length n each, multiplies them homomorphically
// (multiplySamples), and reconstructs the noise on the resulting ciphertext.
//
// Reuses NoiseExperimentResult / writeNoiseCSV from noise_experiment.h, so
// plot_noise.py works unchanged on this data too -- just point it at
// whichever CSV file you wrote this result to.
//
// Note: under COEFFICIENT packing the product is longer than the inputs
// (2n-1), which this function accounts for internally.
// ---------------------------------------------------------------------
NoiseExperimentResult runNoiseAfterMultExperiment(const Context& ctx, size_t numTrials, size_t n);

// ---------------------------------------------------------------------
// Repeatedly multiplies an accumulator ciphertext by fresh random messages
// until decryption fails, or maxMultiplications is reached first. Doing
// this for `numChains` independent chains gives a distribution of "depth
// until failure" you can plot as a histogram.
//
// CAVEAT: under COEFFICIENT packing, the accumulator's message length
// grows by (n-1) with every multiplication. This is only a valid model as
// long as the true polynomial degree stays below the ring dimension N.
// ---------------------------------------------------------------------

// One row of noise-vs-depth data: the noise magnitude observed at a given
// multiplication depth, within a given chain. maxAbsNoise is the largest
// |noise coefficient| across all N ring positions at that depth -- a
// single summary number per step, since storing the full N-length noise
// vector at every depth of every chain would be enormous (numChains *
// maxMultiplications * N values).
struct NoiseAtDepth {
    size_t chain;
    size_t depth;
    int64_t maxAbsNoise;
};

struct MultiplicationDepthResult {
    std::vector<size_t> depthsUntilFailure;  // one entry per chain
    std::vector<NoiseAtDepth> noiseHistory;  // one entry per (chain, depth) sampled
    size_t maxMult = 0;
    size_t numChains = 0;
    size_t numCorrect = 0;  // chains that reached maxMultiplications without failing
    size_t numFailed = 0;   // chains that hit a decryption failure before maxMultiplications
};

// noiseSampleStride: only reconstruct+record noise every `noiseSampleStride`
// multiplications (default 1 = every step). Noise reconstruction is the
// expensive part (CRTInterpolate over the full ciphertext modulus), so for
// large maxMultiplications you may want e.g. 10 or 50 to keep runtime
// reasonable while still getting a representative growth curve. The
// correctness check (isCorrect) still runs every single step regardless --
// only the noise *recording* is subsampled.
MultiplicationDepthResult runMultiplicationDepthExperiment(const Context& ctx, size_t numChains, size_t n,
                                                             size_t maxMultiplications, size_t noiseSampleStride = 1);

void writeDepthsCSV(const std::string& path, const MultiplicationDepthResult& result);

void writeNoiseHistoryCSV(const std::string& path, const MultiplicationDepthResult& result);