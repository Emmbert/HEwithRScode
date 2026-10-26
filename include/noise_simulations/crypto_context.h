
//#ifndef HEWITHRSCODE_CRYPTOCONTEXT_H
//#define HEWITHRSCODE_CRYPTOCONTEXT_H


#pragma once

#include "openfhe.h"

using namespace lbcrypto;

enum class PackingMode { SIMD, COEFFICIENT };
constexpr PackingMode PACKING_MODE = PackingMode::COEFFICIENT;

struct Context {
    CryptoContext<DCRTPoly> cc;
    KeyPair<DCRTPoly> keyPair;
    uint32_t multiplicativeDepth;
};

Context generateContext();

double computeMaxNoiseThreshold(const Context& ctx);

void printContextDetails(const Context& ctx);

//#endif //HEWITHRSCODE_CRYPTOCONTEXT_H