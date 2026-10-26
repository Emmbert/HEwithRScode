#pragma once
// Common interface for all decoders. Add a new decoder by:
//   1. include/rs/decoders/<name>.hpp + src/rs/decoders/<name>.cpp deriving from Decoder,
//   2. one line in src/rs/decoders/registry.cpp,
//   3. tests/test_<name>.cpp + one line in tests/CMakeLists.txt.
#include <memory>
#include <string>
#include "rs/channel.hpp"
#include "rs/core/rs_code.hpp"

namespace rs {

enum class DecodeStatus { Success, Failure };

struct DecodeResult {
    DecodeStatus status = DecodeStatus::Failure;
    Codeword codeword;             // valid when status == Success
    Message message;               // k coefficients, valid when status == Success
    std::size_t numCorrected = 0;  // symbols changed relative to rx.hard
};

class Decoder {
public:
    explicit Decoder(const RSCode& code) : code_(code) {}  // code must outlive the decoder
    virtual ~Decoder() = default;

    virtual std::string name() const = 0;
    virtual bool isSoftDecision() const { return false; }
    // Max. distance between returned codeword and hard decisions (list/soft decoders may exceed t).
    virtual std::size_t decodingRadius() const { return code_.t(); }

    // Must be thread-safe and free of side effects (all scratch memory local).
    virtual DecodeResult decode(const ReceivedWord& rx) const = 0;

    const RSCode& code() const { return code_; }

protected:
    const RSCode& code_;
};

using DecoderPtr = std::unique_ptr<Decoder>;

}  // namespace rs
