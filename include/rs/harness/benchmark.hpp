#pragma once
// Runs all decoders on identical received words and compares runtime and
// decoding performance.
#include <cstdint>
#include <ostream>
#include <string>
#include <vector>
#include "rs/decoder.hpp"

namespace rs {

struct ChannelPoint {
    enum class Kind { SymbolErrors, Gaussian };
    Kind kind = Kind::SymbolErrors;
    double value = 0;  // number of symbol errors, or noise sigma
    std::string label() const;
};

struct BenchmarkConfig {
    std::vector<ChannelPoint> points;
    std::size_t trials = 100;
    std::uint64_t seed = 1;
};

struct BenchmarkRow {
    std::string decoder, channel;
    std::size_t trials = 0, success = 0, failure = 0, miscorrection = 0;
    double meanUs = 0, medianUs = 0, maxUs = 0;
    double pct(std::size_t x) const { return trials ? 100.0 * double(x) / double(trials) : 0.0; }
};

std::vector<BenchmarkRow> runAllDecoders(const RSCode& code,
                                         const std::vector<DecoderPtr>& decoders,
                                         const BenchmarkConfig& cfg);

void printComparison(std::ostream& os, const RSCode& code, const std::vector<BenchmarkRow>& rows);
void writeCsv(const std::string& path, const RSCode& code, const std::vector<BenchmarkRow>& rows);

}  // namespace rs
