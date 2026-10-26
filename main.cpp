#include "include/noise_simulations/crypto_context.h"
#include "include/noise_simulations/encryption.h"
#include "include/noise_simulations/noise_experiment.h"
#include "include/noise_simulations/multiplication_experiment.h"
#include "include/noise_simulations/noise_growth_experiment.h"
#include "include/noise_simulations/param_sweep.h"

#include <iostream>

int performSingleMultiplication() {
    Context ctx = generateContext();
    printContextDetails(ctx);

    constexpr size_t n = 2;

    // --- Step 1: encrypt two fresh random values ---
    EncryptedSample a = encryptRandomMessage(ctx, n);
    EncryptedSample b = encryptRandomMessage(ctx, n);

    DecryptedSampleWithNoise decryptedA = decryptWithNoise(ctx, a, n);
    DecryptedSampleWithNoise decryptedB = decryptWithNoise(ctx, b, n);

    std::cout << "===== Value A =====" << std::endl;
    std::cout << "Original:  " << a.original[0] << std::endl;
    std::cout << "Decrypted: " << decryptedA.decrypted[0] << std::endl;
    std::cout << "Noise (e): " << decryptedA.noise[0] << std::endl;

    std::cout << "===== Value B =====" << std::endl;
    std::cout << "Original:  " << b.original[0] << std::endl;
    std::cout << "Decrypted: " << decryptedB.decrypted[0] << std::endl;
    std::cout << "Noise (e): " << decryptedB.noise[0] << std::endl;

    // --- Step 2: multiply homomorphically ---
    EncryptedSample product = multiplySamples(ctx, a, b);

    // Under COEFFICIENT packing, multiplying two length-1 vectors still
    // gives a length-1 result (1+1-1=1), so this stays a single number too.
    size_t productLength = product.original.size();
    DecryptedSampleWithNoise decryptedProduct = decryptWithNoise(ctx, product, productLength);

    std::cout << "===== A * B =====" << std::endl;
    std::cout << "Expected (A*B mod t, balanced): ";
    for (auto v : product.original) std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "Decrypted: ";
    for (auto v : decryptedProduct.decrypted) std::cout << v << " ";
    std::cout << std::endl;

    std::cout << "Noise (e): ";
    for (auto v : decryptedProduct.noise) std::cout << v << " ";
    std::cout << std::endl;

    bool correct = isCorrect(product.original, decryptedProduct.decrypted);
    std::cout << "Correct decryption: " << (correct ? "YES" : "NO") << std::endl;

    return 0;
}

int runEncryptionNoiseExperiment() {

    Context ctx = generateContext();
    printContextDetails(ctx);

    constexpr size_t n = 16384;
    constexpr size_t numTrials = 100;

    NoiseExperimentResult result = runNoiseExperiment(ctx, numTrials, n, true);

    double sigma = 0.0;
    if (auto p = std::dynamic_pointer_cast<CryptoParametersRLWE<DCRTPoly>>(ctx.cc->GetCryptoParameters())) {
        sigma = p->GetDistributionParameter();
    }
    double maxNoiseThreshold = computeMaxNoiseThreshold(ctx);

    std::cout << "----- Noise Experiment Summary -----" << std::endl;
    std::cout << "Trials:                " << result.numTrials << std::endl;
    std::cout << "Decryption failures:   " << result.numDecryptionFailures << std::endl;
    std::cout << "Over-threshold trials: " << result.numOverThreshold << std::endl;
    std::cout << "Inconsistent trials:   " << result.numInconsistentTrials << std::endl;
    std::cout << "Total noise samples:   " << result.allNoiseValues.size() << std::endl;
    std::cout << "-------------------------------------" << std::endl;

    writeNoiseCSV("noise_samples.csv", result, sigma, maxNoiseThreshold);
    return 0;
}

int runMultiplicationNoiseExperiment() {
    Context ctx = generateContext();
    printContextDetails(ctx);

    constexpr size_t n = 16384;
    constexpr size_t numTrials = 100;

    double maxNoiseThreshold = computeMaxNoiseThreshold(ctx);

    // --- Experiment 2: noise after one multiplication ---
    NoiseExperimentResult multResult = runNoiseAfterMultExperiment(ctx, numTrials, n);

    std::cout << "----- Noise After One Multiplication Summary -----" << std::endl;
    std::cout << "Trials:                " << multResult.numTrials << std::endl;
    std::cout << "Decryption failures:   " << multResult.numDecryptionFailures << std::endl;
    std::cout << "Over-threshold trials: " << multResult.numOverThreshold << std::endl;
    std::cout << "Inconsistent trials:   " << multResult.numInconsistentTrials << std::endl;
    std::cout << "Total noise samples:   " << multResult.allNoiseValues.size() << std::endl;
    std::cout << "----------------------------------------------------" << std::endl;

    // We don't have a closed-form theoretical sigma for post-multiplication
    // noise (unlike the fresh case), so we pass -1.0 as a placeholder --
    // plot_noise.py will simply skip drawing a theoretical curve.
    writeNoiseCSV("noise_after_mult.csv", multResult, -1.0, maxNoiseThreshold);
    return 0;
}

// Your current experiment: noise growth with the default parameters
int runNoiseGrowthExperiment() {
    Context ctx = generateContext();
    printContextDetails(ctx);

    std::vector<OpChainRow> rows;
    auto run = [&](const OpChainConfig& cfg) {
        auto r = runOpChainExperiment(ctx, cfg);
        rows.insert(rows.end(), r.begin(), r.end());
    };
    //   name               operations                  chains  n   modReduce  measure every
    run({"mult",           repeatOp(HeOp::Mult, 20),   10,     16, false,     1});
    run({"mult+modreduce", repeatOp(HeOp::Mult, 20),   10,     16, true,      1});
    run({"add",            repeatOp(HeOp::Add, 200),   10,     16, false,     10});
    run({"mult+3add",      mixedOps(8, 3),             10,     16, false,     1});
    writeOpChainCSV("simulation_results/noise_growth.csv", ctx, rows);
    writeNoiseHistogramCSV("simulation_results/noise_histograms.csv", ctx, rows);
    return 0;
}

// New: the parameter sweep -- edit this list to choose the parameters
int runParameterSweepExperiment() {
    std::vector<HEParams> sweep = {
        // label  n      depth  firstMod  scalingMod  sigma  p       secure
        {"",      0,     10,    0,        0,          3.19,  65537,  true },  // your current default
        {"",      0,     10,    0,        0,          1.19,  65537,  true },  // your current default
        {"",      0,     10,    0,        0,          6.19,  65537,  true },  // your current default
        {"",      16384, 3,     0,        0,          3.19,  65537,  true },
        {"",      4096,  2,     0,        0,          3.19,  65537,  false},
        {"",      4096,  2,     0,        0,          50.0,  65537,  false},
        {"",      4096,  2,     0,        0,          3.19,  257,    false},
    };

    SweepConfig cfg;
    cfg.numChains = 10;          // chains per parameter set
    cfg.maxMults = 40;           // upper bound on the chain length
    cfg.stepsAfterFailure = 1;   // one more multiplication after the first error
    cfg.outDir = "simulation_results/sweep";
    runParameterSweep(sweep, cfg);
    return 0;
}

// Runs the sweep simulation for all parameter sets listed in a CSV file.
//   csv      : file written by find_candidates.py (calib_params.csv or candidates.csv)
//   outDir   : folder for the results (sweep_index.csv, pXX_growth.csv, pXX_hist.csv)
//   maxMults : maximum number of multiplications per chain
int runSweepFromCSV(const std::string& csv, const std::string& outDir, size_t maxMults) {
    SweepConfig cfg;
    cfg.numChains = 10;          // chains per parameter set
    cfg.maxMults = maxMults;
    cfg.stepsAfterFailure = 1;   // one more multiplication after the first error
    cfg.outDir = outDir;
    runParameterSweep(loadHEParamsCSV(csv), cfg);
    return 0;
}

int main(int argc, char** argv) {
    if (argc == 4) {
        // ./HEwithRScodes <params.csv> <outDir> <maxMults>
        return runSweepFromCSV(argv[1], argv[2], std::stoul(argv[3]));
    }
    return runParameterSweepExperiment();   // or whatever you run by default

    //return runParameterSweepExperiment();
    // return runNoiseGrowthExperiment();
    // return runEncryptionNoiseExperiment();
    // return runMultiplicationNoiseExperiment();
    // return performSingleMultiplication();
}