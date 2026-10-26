// rs_bench: performance comparison of all registered decoders
// for negacyclic RS codes over GF(p), ring Z_p[x]/(x^n + 1).
//
// Usage: rs_bench [n k trials p] [--psi v] [--errors e]... [--sigma s]... [--csv file]
// Correctness tests are run separately with ctest.
//   p = 0 or omitted -> smallest prime p = 1 (mod 2n) (40961 for n = 4096)
// Example: rs_bench 4096 3072 50 65537 --errors 256 --errors 512 --errors 513 --sigma 0.2
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "rs/decoders/registry.hpp"
#include "rs/harness/benchmark.hpp"

int main(int argc, char** argv) {
    std::size_t n = 4096, k = 2048, trials = 20;
    std::uint32_t p = 0;
    rs::Elem psi = 0;
    std::string csv;
    std::vector<std::size_t> errors;
    std::vector<double> sigmas;

    std::vector<std::string> pos;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--csv" && i + 1 < argc) csv = argv[++i];
        else if (a == "--errors" && i + 1 < argc) errors.push_back(std::stoul(argv[++i]));
        else if (a == "--sigma" && i + 1 < argc) sigmas.push_back(std::stod(argv[++i]));
        else if (a == "--psi" && i + 1 < argc) psi = static_cast<rs::Elem>(std::stoul(argv[++i]));
        else pos.push_back(a);
    }
    if (pos.size() > 0) n = std::stoul(pos[0]);
    if (pos.size() > 1) k = std::stoul(pos[1]);
    if (pos.size() > 2) trials = std::stoul(pos[2]);
    if (pos.size() > 3) p = static_cast<std::uint32_t>(std::stoul(pos[3]));

    try {
        if (p == 0) p = rs::RSCode::smallestNegacyclicPrime(n);
        rs::RSCode code(p, n, k, psi);
        auto decoders = rs::makeAllDecoders(code);

        rs::BenchmarkConfig cfg;
        cfg.trials = trials;
        if (errors.empty() && sigmas.empty()) {
            const std::size_t t = code.t();
            std::set<std::size_t> def{0, t / 2, t};
            if (t + 1 <= n) def.insert(t + 1);
            errors.assign(def.begin(), def.end());
        }
        for (std::size_t e : errors)
            cfg.points.push_back({rs::ChannelPoint::Kind::SymbolErrors, double(e)});
        for (double s : sigmas)
            cfg.points.push_back({rs::ChannelPoint::Kind::Gaussian, s});

        const auto rows = rs::runAllDecoders(code, decoders, cfg);
        rs::printComparison(std::cout, code, rows);
        if (!csv.empty()) rs::writeCsv(csv, code, rows);
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
        return 1;
    }
    return 0;
}
