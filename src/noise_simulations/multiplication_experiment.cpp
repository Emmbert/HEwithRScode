#include "../../include/noise_simulations/multiplication_experiment.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cstdlib>

NoiseExperimentResult runNoiseAfterMultExperiment(const Context& ctx, size_t numTrials, size_t n) {
    NoiseExperimentResult result;
    result.numTrials = numTrials;

    double maxNoiseThreshold = computeMaxNoiseThreshold(ctx);
    result.allNoiseValues.reserve(numTrials * (2 * n));  // rough capacity hint

    for (size_t trial = 0; trial < numTrials; ++trial) {
        EncryptedSample a = encryptRandomMessage(ctx, n);
        EncryptedSample b = encryptRandomMessage(ctx, n);
        EncryptedSample product = multiplySamples(ctx, a, b);

        size_t productLength = product.original.size();
        DecryptedSampleWithNoise decryptedSample = decryptWithNoise(ctx, product, productLength);

        bool correct = isCorrect(product.original, decryptedSample.decrypted);
        bool overThreshold = exceedsThreshold(decryptedSample.noise, maxNoiseThreshold);

        if (!correct) {
            result.numDecryptionFailures++;
        }
        if (overThreshold) {
            result.numOverThreshold++;
        }
        if (correct == overThreshold) {
            result.numInconsistentTrials++;
        }

        result.allNoiseValues.insert(result.allNoiseValues.end(), decryptedSample.noise.begin(),
                                      decryptedSample.noise.end());
    }

    return result;
}

MultiplicationDepthResult runMultiplicationDepthExperiment(const Context& ctx, size_t numChains, size_t n,
                                                             size_t maxMultiplications, size_t noiseSampleStride) {
    MultiplicationDepthResult result;
    result.maxMult = maxMultiplications;
    result.numChains = numChains;

    for (size_t chain = 0; chain < numChains; ++chain) {
        EncryptedSample acc = encryptRandomMessage(ctx, n);
        size_t successfulMultiplications = 0;
        bool failed = false;

        for (size_t k = 0; k < maxMultiplications; ++k) {
            EncryptedSample next = encryptRandomMessage(ctx, n);
            EncryptedSample product = multiplySamples(ctx, acc, next);

            bool shouldSampleNoise = (k % noiseSampleStride == 0);
            bool correct;

            if (shouldSampleNoise) {
                // Heavier path: also reconstructs the noise polynomial.
                DecryptedSampleWithNoise decryptedSample = decryptWithNoise(ctx, product, product.original.size());
                correct = isCorrect(product.original, decryptedSample.decrypted);

                int64_t maxAbsNoise = 0;
                for (auto v : decryptedSample.noise) {
                    maxAbsNoise = std::max(maxAbsNoise, std::abs(v));
                }
                result.noiseHistory.push_back(NoiseAtDepth{chain, k, maxAbsNoise});
            } else {
                // Cheaper path: only checks correctness, no noise reconstruction.
                std::vector<int64_t> decrypted = decryptMessage(ctx, product.ciphertext, product.original.size());
                correct = isCorrect(product.original, decrypted);
            }

            if (!correct) {
                failed = true;
                break;
            }

            acc = product;
            successfulMultiplications++;
        }

        result.depthsUntilFailure.push_back(successfulMultiplications);
        if (!failed) {
            result.numCorrect++;
        } else {
            result.numFailed++;
        }
    }

    return result;
}

void writeDepthsCSV(const std::string& path, const MultiplicationDepthResult& result) {
    std::ofstream file(path);
    if (!file.is_open()) {
        std::cerr << "Error: could not open " << path << " for writing." << std::endl;
        return;
    }
    file << "# maxMult=" << result.maxMult << "\n";
    file << "# numChains=" << result.numChains << "\n";
    file << "# numCorrect=" << result.numCorrect << "\n";
    file << "# numFailed=" << result.numFailed << "\n";
    file << "depth\n";
    for (auto d : result.depthsUntilFailure) {
        file << d << "\n";
    }

    file.close();
    std::cout << "Wrote " << result.depthsUntilFailure.size() << " chain depths to " << path << std::endl;
}

void writeNoiseHistoryCSV(const std::string& path, const MultiplicationDepthResult& result) {
    std::ofstream file(path);
    if (!file.is_open()) {
        std::cerr << "Error: could not open " << path << " for writing." << std::endl;
        return;
    }

    file << "# maxMult=" << result.maxMult << "\n";
    file << "# numChains=" << result.numChains << "\n";
    file << "chain,depth,maxAbsNoise\n";
    for (const auto& row : result.noiseHistory) {
        file << row.chain << "," << row.depth << "," << row.maxAbsNoise << "\n";
    }

    file.close();
    std::cout << "Wrote " << result.noiseHistory.size() << " noise-history rows to " << path << std::endl;
}