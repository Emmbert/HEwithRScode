#include "../../include/noise_simulations/encryption.h"
#include <random>
#include <iostream>
#include <algorithm>

namespace {

// Reduces x into the balanced representation [-t/2, t/2), matching how
// OpenFHE's decoders represent plaintext values.
int64_t reduceBalanced(int64_t x, int64_t t) {
    int64_t r = x % t;
    if (r < 0) r += t;
    if (r > t / 2) r -= t;  // strict '>', matching OpenFHE's own centering convention
    return r;
}

// What EvalMult produces at the plaintext level, given how PACKING_MODE
// interprets the ring:
//   - SIMD: elementwise (slot-wise) product, same length as the inputs.
//   - COEFFICIENT: the two messages multiplied as polynomials in the ring
//     Z[X]/(X^N+1) -- a genuine NEGACYCLIC convolution. Terms landing at
//     degree >= N wrap around with a sign flip (since X^N = -1), not just
//     truncate. Result length = min(a.size()+b.size()-1, N).
std::vector<int64_t> computeExpectedProduct(const std::vector<int64_t>& a, const std::vector<int64_t>& b, int64_t t,
                                             size_t ringDimension) {
    if (PACKING_MODE == PackingMode::SIMD) {
        std::vector<int64_t> result(a.size());
        for (size_t i = 0; i < a.size(); ++i) {
            result[i] = reduceBalanced(a[i] * b[i], t);
        }
        return result;
    } else {
        size_t rawLength = a.size() + b.size() - 1;
        size_t resultLength = std::min(rawLength, ringDimension);
        std::vector<int64_t> result(resultLength, 0);

        for (size_t i = 0; i < a.size(); ++i) {
            for (size_t j = 0; j < b.size(); ++j) {
                size_t k = i + j;
                int64_t term = a[i] * b[j];
                if (k < ringDimension) {
                    result[k] += term;
                } else {
                    // X^k = X^(k - N) * X^N = -X^(k - N)  (negacyclic wraparound)
                    result[k - ringDimension] -= term;
                }
            }
        }

        for (auto& v : result) v = reduceBalanced(v, t);
        return result;
    }
}

}  // namespace

// ---------------------------------------------------------------------
// Symmetric-key encryption: c1 = a (uniform), c0 = -a*s + t*e + m.
// No ephemeral u, no public-key generation error e0 involved -- the
// reconstructed noise here is a single fresh discrete Gaussian sample,
// directly comparable to sigma. Use this to validate the noise
// reconstruction math itself; use encryptMessage() (public-key) for
// noise that reflects how ciphertexts are actually used in practice.
// ---------------------------------------------------------------------
EncryptedSample encryptMessageSymmetric(const Context& ctx, std::vector<int64_t> message_vec) {
    Plaintext plaintext;
    if (PACKING_MODE == PackingMode::SIMD) {
        plaintext = ctx.cc->MakePackedPlaintext(message_vec);
    } else {
        plaintext = ctx.cc->MakeCoefPackedPlaintext(message_vec);
    }

    auto ciphertext = ctx.cc->Encrypt(ctx.keyPair.secretKey, plaintext);

    return EncryptedSample{message_vec, ciphertext};
}

EncryptedSample encryptMessage(const Context& ctx, std::vector<int64_t> message_vec) {

    Plaintext plaintext;
    if (PACKING_MODE == PackingMode::SIMD) {
        plaintext = ctx.cc->MakePackedPlaintext(message_vec);
    } else {
        plaintext = ctx.cc->MakeCoefPackedPlaintext(message_vec);
    }

    auto ciphertext = ctx.cc->Encrypt(ctx.keyPair.publicKey, plaintext);

    return EncryptedSample{message_vec, ciphertext};
}

EncryptedSample encryptRandomMessage(const Context& ctx, size_t n, bool useSymmetric) {
    static std::mt19937 rng(std::random_device{}());

    int64_t t = ctx.cc->GetCryptoParameters()->GetPlaintextModulus();
    std::vector<int64_t> vec(n);
    std::uniform_int_distribution<int64_t> dist(-t / 2, t / 2 - 1);
    for (auto& v : vec) v = dist(rng);

    return useSymmetric ? encryptMessageSymmetric(ctx, vec) : encryptMessage(ctx, vec);
}

std::vector<int64_t> decryptMessage(const Context& ctx, const Ciphertext<DCRTPoly>& ciphertext, size_t n) {
    Plaintext result;
    ctx.cc->Decrypt(ctx.keyPair.secretKey, ciphertext, &result);
    result->SetLength(n);

    if (PACKING_MODE == PackingMode::SIMD) {
        return result->GetPackedValue();
    } else {
        return result->GetCoefPackedValue();
    }
}

// ---------------------------------------------------------------------
// See crypto_context.h/.cpp for background: v = c0 + c1*s = m + t*e (mod Q)
// under FIXEDMANUAL scaling. We redo that linear step ourselves and solve
// for e. Only valid 1:1 with COEFFICIENT packing (see earlier discussion).
//
// IMPORTANT: the noise polynomial has N (ring dimension) coefficients,
// generated independently of the message -- not just n of them. Since
// COEFFICIENT packing zero-pads every position i >= n, we know the exact
// expected plaintext value at every one of the N positions (original[i]
// for i < n, else 0), so we sample noise at ALL N positions rather than
// only the n "meaningful" ones. This means every single encryption gives
// you N noise samples "for free" instead of just n.
// ---------------------------------------------------------------------
DecryptedSampleWithNoise decryptWithNoise(const Context& ctx, const EncryptedSample& encrypted_sample, size_t n){
    Ciphertext<DCRTPoly> ciphertext = encrypted_sample.ciphertext;
    const std::vector<int64_t>& original = encrypted_sample.original;
    std::vector<int64_t> decrypted = decryptMessage(ctx, ciphertext, n);

    const std::vector<DCRTPoly>& cv = ciphertext->GetElements();
    DCRTPoly s = ctx.keyPair.secretKey->GetPrivateElement();

    DCRTPoly raw = cv[0] + cv[1] * s;
    raw.SetFormat(COEFFICIENT);

    auto rawBig = raw.CRTInterpolate();
    auto Qmod = rawBig.GetModulus();
    const auto& coeffs = rawBig.GetValues();

    int64_t t = ctx.cc->GetCryptoParameters()->GetPlaintextModulus();
    size_t ringDimension = ctx.cc->GetRingDimension();

    std::vector<int64_t> noise(ringDimension);
    for (size_t i = 0; i < ringDimension; ++i) {
        auto val = coeffs[i];

        int64_t signedVal;
        if (val > (Qmod - val)) {
            signedVal = -static_cast<int64_t>((Qmod - val).ConvertToInt<uint64_t>());
        } else {
            signedVal = static_cast<int64_t>(val.ConvertToInt<uint64_t>());
        }

        // Expected plaintext value at ring position i: the original message
        // value if i falls within it, otherwise 0 (zero-padding).
        int64_t expectedMessageAtI = (i < original.size()) ? original[i] : 0;

        int64_t diff = signedVal - expectedMessageAtI;

        if (diff % t != 0) {
            std::cerr << "Warning: noise reconstruction inconsistent at index " << i
                       << " (diff=" << diff << ", t=" << t << ")" << std::endl;
        }

        noise[i] = diff / t;
    }

    return DecryptedSampleWithNoise{decrypted, noise};
}

bool isCorrect(const std::vector<int64_t>& original, const std::vector<int64_t>& decrypted) {
    if (original.size() != decrypted.size()) {
        return false;
    }
    for (size_t i = 0; i < original.size(); ++i) {
        if (original[i] != decrypted[i]) {
            return false;
        }
    }
    return true;
}

EncryptedSample multiplySamples(const Context& ctx, const EncryptedSample& a, const EncryptedSample& b) {
    int64_t t = ctx.cc->GetCryptoParameters()->GetPlaintextModulus();
    size_t ringDimension = ctx.cc->GetRingDimension();

    auto ciphertext = ctx.cc->EvalMult(a.ciphertext, b.ciphertext);
    std::vector<int64_t> expected = computeExpectedProduct(a.original, b.original, t, ringDimension);

    return EncryptedSample{expected, ciphertext};
}