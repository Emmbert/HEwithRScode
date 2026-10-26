#include "../../include/noise_simulations/crypto_context.h"
#include <iostream>
#include <string>

Context generateContext() {
    CCParams<CryptoContextBGVRNS> parameters;

    uint32_t multDepth = 10;
    parameters.SetPlaintextModulus(65537);
    parameters.SetMultiplicativeDepth(multDepth);
    parameters.SetSecurityLevel(HEStd_128_classic);
    // Use the original ("unscaled") BGV variant, not the default GHS-scaled
    // FLEXIBLEAUTOEXT. This makes v = c0 + c1*s = m + t*e (mod Q) hold exactly,
    // which decryptWithNoise() in encryption.cpp relies on.
    parameters.SetScalingTechnique(FIXEDMANUAL);

    /*CCParams<CryptoContextBGVRNS> parameters;
    uint32_t multDepth = 0;
    parameters.SetPlaintextModulus(65537);
    parameters.SetScalingTechnique(FIXEDMANUAL);

    parameters.SetSecurityLevel(HEStd_NotSet);   // <-- unlocks manual n
    parameters.SetRingDim(1024);                 // n
    parameters.SetStandardDeviation(50.0f);      // sigma (exaggerated, for testing)

    parameters.SetMultiplicativeDepth(multDepth);
    parameters.SetFirstModSize(28);              // controls size of q
    parameters.SetScalingModSize(28);*/

    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

    //CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);

    KeyPair<DCRTPoly> keyPair = cc->KeyGen();
    cc->EvalMultKeyGen(keyPair.secretKey);

    return Context{cc, keyPair, multDepth};
}

double computeMaxNoiseThreshold(const Context& ctx) {
    auto cryptoParams = ctx.cc->GetCryptoParameters();
    auto ciphertextModulus = cryptoParams->GetElementParams()->GetModulus();
    auto plaintextModulus = cryptoParams->GetPlaintextModulus();

    double approxQ = ciphertextModulus.ConvertToDouble();
    return approxQ / (2.0 * static_cast<double>(plaintextModulus));
}

void printContextDetails(const Context& ctx) {
    auto cryptoParams = ctx.cc->GetCryptoParameters();

    uint32_t ringDimension = ctx.cc->GetRingDimension();
    auto ciphertextModulus = cryptoParams->GetElementParams()->GetModulus();
    auto plaintextModulus = cryptoParams->GetPlaintextModulus();

    auto cryptoParamsRLWE = std::dynamic_pointer_cast<CryptoParametersRLWE<DCRTPoly>>(cryptoParams);
    double stdDev = cryptoParamsRLWE ? cryptoParamsRLWE->GetDistributionParameter() : -1.0;

    std::string packingDescription = (PACKING_MODE == PackingMode::SIMD)
        ? "SIMD-style (Packed) encoding via MakePackedPlaintext"
        : "Coefficient-wise encoding via MakeCoefPackedPlaintext";

    double approxMaxNoise = computeMaxNoiseThreshold(ctx);

    std::cout << "----- Crypto Context Details -----" << std::endl;
    std::cout << "Encryption scheme:        BGV" << std::endl;
    std::cout << "Ring dimension (N):       " << ringDimension << std::endl;
    std::cout << "Ciphertext modulus (Q):   " << ciphertextModulus << std::endl;
    std::cout << "Plaintext modulus (t):    " << plaintextModulus << std::endl;
    std::cout << "Standard deviation (sigma): " << stdDev << std::endl;
    std::cout << "Plaintext packing method: " << packingDescription << std::endl;
    std::cout << "Multiplicative depth:     " << ctx.multiplicativeDepth << std::endl;
    std::cout << "Max decryptable noise |e| (approx, ~Q/(2t)): " << approxMaxNoise << std::endl;
    std::cout << "-----------------------------------" << std::endl;
}