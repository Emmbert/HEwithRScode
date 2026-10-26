#include "rs/harness/benchmark.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>

#include "rs/channel.hpp"

namespace rs {

std::string ChannelPoint::label() const {
    std::ostringstream os;
    if (kind == Kind::SymbolErrors) os << "e=" << static_cast<std::size_t>(value);
    else os << "sigma=" << std::fixed << std::setprecision(3) << value;
    return os.str();
}

std::vector<BenchmarkRow> runAllDecoders(const RSCode& code,
                                         const std::vector<DecoderPtr>& decoders,
                                         const BenchmarkConfig& cfg) {
    using Clock = std::chrono::steady_clock;
    std::mt19937_64 rng(cfg.seed);
    std::vector<BenchmarkRow> rows;

    for (const ChannelPoint& pt : cfg.points) {
        std::vector<BenchmarkRow> local(decoders.size());
        std::vector<std::vector<double>> times(decoders.size());
        for (std::size_t d = 0; d < decoders.size(); ++d) {
            local[d].decoder = decoders[d]->name();
            local[d].channel = pt.label();
            times[d].reserve(cfg.trials);
        }

        for (std::size_t trial = 0; trial < cfg.trials; ++trial) {
            // Every decoder sees the same received word -> fair comparison.
            const Message msg = code.randomMessage(rng);
            const Codeword cw = code.encode(msg);
            const ReceivedWord rx = (pt.kind == ChannelPoint::Kind::SymbolErrors)
                ? symbolErrorChannel(code, cw, static_cast<std::size_t>(pt.value), rng)
                : gaussianChannel(code, cw, pt.value, rng);

            for (std::size_t d = 0; d < decoders.size(); ++d) {
                BenchmarkRow& row = local[d];
                DecodeResult r;
                const auto t0 = Clock::now();
                try {
                    r = decoders[d]->decode(rx);
                } catch (const std::exception&) {
                    r.status = DecodeStatus::Failure;
                }
                const auto t1 = Clock::now();
                times[d].push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());

                ++row.trials;
                if (r.status != DecodeStatus::Success) ++row.failure;
                else if (r.codeword == cw) ++row.success;
                else ++row.miscorrection;
            }
        }

        for (std::size_t d = 0; d < decoders.size(); ++d) {
            auto& ts = times[d];
            if (!ts.empty()) {
                local[d].meanUs = std::accumulate(ts.begin(), ts.end(), 0.0) / double(ts.size());
                local[d].maxUs = *std::max_element(ts.begin(), ts.end());
                std::nth_element(ts.begin(), ts.begin() + long(ts.size() / 2), ts.end());
                local[d].medianUs = ts[ts.size() / 2];
            }
            rows.push_back(local[d]);
        }
    }
    return rows;
}

void printComparison(std::ostream& os, const RSCode& code, const std::vector<BenchmarkRow>& rows) {
    const double infoBits = double(code.k()) * std::log2(double(code.field().p()));
    std::size_t w = 8;
    for (const auto& r : rows) w = std::max(w, r.decoder.size());

    os << "\n" << code.describe() << "\n\n" << std::left << std::setw(int(w)) << "Decoder"
       << "  " << std::setw(16) << "Channel" << std::right << std::setw(7) << "Trials"
       << std::setw(10) << "Correct%" << std::setw(8) << "Fail%" << std::setw(10) << "Miscorr%"
       << std::setw(12) << "mean[us]" << std::setw(12) << "median[us]" << std::setw(12) << "max[us]"
       << std::setw(12) << "Mbit/s" << "\n"
       << std::string(w + 2 + 16 + 7 + 10 + 8 + 10 + 12 * 4, '-') << "\n";

    os << std::fixed;
    for (const auto& r : rows) {
        os << std::left << std::setw(int(w)) << r.decoder << "  " << std::setw(16) << r.channel
           << std::right << std::setw(7) << r.trials << std::setprecision(1)
           << std::setw(10) << r.pct(r.success) << std::setw(8) << r.pct(r.failure)
           << std::setw(10) << r.pct(r.miscorrection) << std::setprecision(1)
           << std::setw(12) << r.meanUs << std::setw(12) << r.medianUs << std::setw(12) << r.maxUs
           << std::setprecision(2) << std::setw(12) << (r.meanUs > 0 ? infoBits / r.meanUs : 0.0)
           << "\n";
    }
    os << "\nMbit/s = information bits (k*log2 p) per mean decoding time.\n";
}

void writeCsv(const std::string& path, const RSCode& code, const std::vector<BenchmarkRow>& rows) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("writeCsv: cannot open " + path);
    f << "p,n,k,decoder,channel,trials,success,failure,miscorrection,mean_us,median_us,max_us\n";
    for (const auto& r : rows)
        f << code.field().p() << ',' << code.n() << ',' << code.k() << ",\"" << r.decoder << "\","
          << r.channel << ',' << r.trials << ',' << r.success << ',' << r.failure << ','
          << r.miscorrection << ',' << r.meanUs << ',' << r.medianUs << ',' << r.maxUs << '\n';
}

}  // namespace rs
