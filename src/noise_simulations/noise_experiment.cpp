#include "../../include/noise_simulations/noise_experiment.h"
#include <cmath>
#include <fstream>
#include <iostream>

bool exceedsThreshold(const std::vector<int64_t>& noise, double maxNoiseThreshold) {
    for (auto e : noise) {
        if (std::abs(static_cast<double>(e)) > maxNoiseThreshold) {
            return true;
        }
    }
    return false;
}

namespace {

    struct MeanStdDev {
        double mean;
        double stdDev;
    };

    MeanStdDev computeMeanStdDev(const std::vector<int64_t>& values) {
        if (values.empty()) {
            return {0.0, 0.0};
        }

        double sum = 0.0;
        for (auto v : values) sum += static_cast<double>(v);
        double mean = sum / static_cast<double>(values.size());

        double sqDiffSum = 0.0;
        for (auto v : values) {
            double diff = static_cast<double>(v) - mean;
            sqDiffSum += diff * diff;
        }
        double stdDev = std::sqrt(sqDiffSum / static_cast<double>(values.size()));

        return {mean, stdDev};
    }

}  // namespace


double computeExpectedPublicKeyNoiseStdDev(const Context& ctx) {
    auto cryptoParamsRLWE = std::dynamic_pointer_cast<CryptoParametersRLWE<DCRTPoly>>(ctx.cc->GetCryptoParameters());
    if (!cryptoParamsRLWE) {
        std::cerr << "Warning: could not access CryptoParametersRLWE." << std::endl;
        return -1.0;
    }

    double sigma = cryptoParamsRLWE->GetDistributionParameter();
    double N = static_cast<double>(ctx.cc->GetRingDimension());
    SecretKeyDist skDist = cryptoParamsRLWE->GetSecretKeyDist();

    double varSecret;
    switch (skDist) {
        case UNIFORM_TERNARY:
            varSecret = 2.0 / 3.0;  // coefficients uniform in {-1,0,1}: E[s^2] = 2/3
            break;
        case GAUSSIAN:
            // Best-effort assumption: the Gaussian secret key uses the same
            // sigma as the noise distribution. OpenFHE doesn't expose a
            // separate getter for the key's own stddev, so treat this case
            // as approximate rather than exact.
            varSecret = sigma * sigma;
            break;
        default:
            std::cerr << "Warning: no variance formula implemented for this "
                         "SecretKeyDist (e.g. SPARSE_TERNARY) -- returning -1."
                      << std::endl;
            return -1.0;
    }

    // Var(total) = sigma^2 * (1 + N*Var(u) + N*Var(s)), with u and s sharing
    // the same distribution => Var(u) == Var(s) == varSecret.
    double variance = sigma * sigma * (1.0 + 2.0 * N * varSecret);
    return std::sqrt(variance);
}

NoiseExperimentResult runNoiseExperiment(const Context& ctx, size_t numTrials, size_t n, bool useSymmetric) {
    NoiseExperimentResult result;
    result.numTrials = numTrials;

    double maxNoiseThreshold = computeMaxNoiseThreshold(ctx);

    result.allNoiseValues.reserve(numTrials * n);

    for (size_t trial = 0; trial < numTrials; ++trial) {
        EncryptedSample sample = encryptRandomMessage(ctx, n, useSymmetric);
        DecryptedSampleWithNoise decryptedSample = decryptWithNoise(ctx, sample, n);

        bool correct = isCorrect(sample.original, decryptedSample.decrypted);
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

        result.allNoiseValues.insert(result.allNoiseValues.end(),
                                      decryptedSample.noise.begin(),
                                      decryptedSample.noise.end());
    }

    return result;
}

void writeNoiseCSV(const std::string& path, const NoiseExperimentResult& result, double sigma, double maxNoiseThreshold) {
    std::ofstream file(path);
    if (!file.is_open()) {
        std::cerr << "Error: could not open " << path << " for writing." << std::endl;
        return;
    }

    file << "# numTrials=" << result.numTrials << "\n";
    file << "# sigma=" << sigma << "\n";
    file << "# maxNoiseThreshold=" << maxNoiseThreshold << "\n";
    file << "# numDecryptionFailures=" << result.numDecryptionFailures << "\n";
    file << "# numOverThreshold=" << result.numOverThreshold << "\n";
    file << "# numInconsistentTrials=" << result.numInconsistentTrials << "\n";
    file << "noise\n";
    for (auto v : result.allNoiseValues) {
        file << v << "\n";
    }

    file.close();
    std::cout << "Wrote " << result.allNoiseValues.size() << " noise samples to " << path << std::endl;
}

bool selfTestNoiseReconstruction(const Context& ctx, bool symmetric) {
    constexpr size_t kTrials = 200;
    constexpr size_t kSlotsPerTrial = 32;
    constexpr double kRelativeTolerance = 0.20;

    auto cryptoParamsRLWE = std::dynamic_pointer_cast<CryptoParametersRLWE<DCRTPoly>>(ctx.cc->GetCryptoParameters());
    double sigma = cryptoParamsRLWE ? cryptoParamsRLWE->GetDistributionParameter() : -1.0;

    // Symmetric encryption isolates a single fresh e ~ N(0, sigma); public-key
    // encryption's noise is amplified by convolution with the secret key.
    double expectedStdDev = symmetric ? sigma : computeExpectedPublicKeyNoiseStdDev(ctx);

    NoiseExperimentResult result = runNoiseExperiment(ctx, kTrials, kSlotsPerTrial, symmetric);
    MeanStdDev stats = computeMeanStdDev(result.allNoiseValues);

    double relativeError = (expectedStdDev > 0.0) ? std::abs(stats.stdDev - expectedStdDev) / expectedStdDev : -1.0;
    bool passed = (result.numInconsistentTrials == 0) && (relativeError >= 0.0) && (relativeError <= kRelativeTolerance);

    std::cout << "----- selfTestNoiseReconstruction (" << (symmetric ? "symmetric" : "public-key") << ") -----" << std::endl;
    std::cout << "Trials:                  " << kTrials << std::endl;
    std::cout << "Theoretical sigma:       " << sigma << std::endl;
    std::cout << "Expected stddev(e):      " << expectedStdDev << std::endl;
    std::cout << "Empirical mean(e):       " << stats.mean << std::endl;
    std::cout << "Empirical stddev(e):     " << stats.stdDev << std::endl;
    std::cout << "Relative error:          " << (relativeError * 100.0) << "%" << std::endl;
    std::cout << "Inconsistent trials:     " << result.numInconsistentTrials << std::endl;
    std::cout << "Result:                  " << (passed ? "PASS" : "FAIL") << std::endl;
    std::cout << "----------------------------------------" << std::endl;

    return passed;
}