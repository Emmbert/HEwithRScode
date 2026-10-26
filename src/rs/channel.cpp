#include "rs/channel.hpp"

#include <cmath>
#include <numeric>
#include <stdexcept>

namespace rs {

ReceivedWord symbolErrorChannel(const RSCode& code, const Codeword& cw,
                                std::size_t numErrors, std::mt19937_64& rng) {
    const std::size_t n = code.n();
    const PrimeField& F = code.field();
    if (cw.size() != n) throw std::invalid_argument("symbolErrorChannel: wrong codeword length");
    if (numErrors > n) throw std::invalid_argument("symbolErrorChannel: numErrors > n");

    ReceivedWord rx;
    rx.hard = cw;
    std::vector<std::size_t> idx(n);
    std::iota(idx.begin(), idx.end(), std::size_t{0});
    std::uniform_int_distribution<Elem> nonzero(1, F.p() - 1);
    for (std::size_t i = 0; i < numErrors; ++i) {  // partial Fisher-Yates
        std::uniform_int_distribution<std::size_t> pick(i, n - 1);
        std::swap(idx[i], idx[pick(rng)]);
        rx.hard[idx[i]] = F.add(rx.hard[idx[i]], nonzero(rng));
    }
    rx.frac.assign(n, 0.0f);
    return rx;
}

ReceivedWord gaussianChannel(const RSCode& code, const Codeword& cw,
                             double sigma, std::mt19937_64& rng) {
    const PrimeField& F = code.field();
    const std::int64_t p = F.p();
    std::normal_distribution<double> noise(0.0, sigma);

    ReceivedWord rx;
    rx.sigma = sigma;
    rx.hard.resize(cw.size());
    rx.frac.resize(cw.size());
    for (std::size_t i = 0; i < cw.size(); ++i) {
        const double e = noise(rng);
        const double r = std::round(e);
        std::int64_t v = (static_cast<std::int64_t>(cw[i]) + static_cast<std::int64_t>(r)) % p;
        if (v < 0) v += p;
        rx.hard[i] = static_cast<Elem>(v);
        rx.frac[i] = static_cast<float>(e - r);
    }
    return rx;
}

}  // namespace rs
