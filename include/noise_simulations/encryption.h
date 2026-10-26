#pragma once

#include "crypto_context.h"
#include <vector>
#include <cstdint>

struct EncryptedSample {
    std::vector<int64_t> original;
    Ciphertext<DCRTPoly> ciphertext;
};

EncryptedSample encryptMessageSymmetric(const Context& ctx, std::vector<int64_t> message_vec);

EncryptedSample encryptMessage(const Context& ctx, std::vector<int64_t> message_vec);

EncryptedSample encryptRandomMessage(const Context& ctx, size_t n, bool useSymmetric = true);

std::vector<int64_t> decryptMessage(const Context& ctx, const Ciphertext<DCRTPoly>& ciphertext, size_t n);

// noise.size() == ring dimension N (not n!) -- every ring coefficient is
// sampled, since the noise polynomial is independent of the message and
// COEFFICIENT packing zero-pads every ring position past the message.
struct DecryptedSampleWithNoise {
    std::vector<int64_t> decrypted;
    std::vector<int64_t> noise;
};

// Only valid 1:1 with COEFFICIENT packing (ring-position <-> message-index
// correspondence breaks under SIMD packing -- see earlier discussion).
DecryptedSampleWithNoise decryptWithNoise(const Context& ctx, const EncryptedSample& encrypted_sample, size_t n);

EncryptedSample multiplySamples(const Context& ctx, const EncryptedSample& a, const EncryptedSample& b);

bool isCorrect(const std::vector<int64_t>& original, const std::vector<int64_t>& decrypted);