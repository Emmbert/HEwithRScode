#include "rs/core/rs_code.hpp"

#include <sstream>
#include <stdexcept>

namespace rs {

namespace {
bool isPowerOfTwo(std::size_t x) { return x && !(x & (x - 1)); }
}  // namespace

std::uint32_t RSCode::smallestNegacyclicPrime(std::size_t n) {
    const std::uint64_t step = 2 * static_cast<std::uint64_t>(n);
    for (std::uint64_t c = step + 1; c < (1ull << 31); c += step)
        if (PrimeField::isPrime(c)) return static_cast<std::uint32_t>(c);
    throw std::invalid_argument("smallestNegacyclicPrime: no prime < 2^31 found");
}

RSCode::RSCode(std::uint32_t p, std::size_t n, std::size_t k, Elem psi)
    : field_(p), n_(n), k_(k), psi_(psi) {
    const PrimeField& F = field_;
    if (!isPowerOfTwo(n) || n < 2) throw std::invalid_argument("RSCode: n must be a power of two (ring x^n + 1)");
    if (k < 1 || k > n) throw std::invalid_argument("RSCode: need 1 <= k <= n");
    if ((static_cast<std::uint64_t>(p) - 1) % (2 * static_cast<std::uint64_t>(n)) != 0)
        throw std::invalid_argument("RSCode: need p = 1 (mod 2n) for a negacyclic code");

    if (psi_ == 0) psi_ = F.rootOfUnity(2 * n);
    // For n a power of two, psi^n == -1 is equivalent to psi having order exactly 2n.
    if (psi_ >= p || F.pow(psi_, n) != p - 1)
        throw std::invalid_argument("RSCode: psi is not a primitive 2n-th root of unity");

    const Elem omega = F.mul(psi_, psi_), psiInv = F.inv(psi_), omegaInv = F.inv(omega);
    nInv_ = F.inv(static_cast<Elem>(n % p));

    psiPow_.resize(n);
    psiInvPow_.resize(n);
    points_.resize(n);
    Elem a = 1, b = 1;
    for (std::size_t j = 0; j < n; ++j) {
        psiPow_[j] = a;
        psiInvPow_[j] = b;
        points_[j] = F.mul(F.mul(a, a), psi_);  // psi^(2j+1)
        a = F.mul(a, psi_);
        b = F.mul(b, psiInv);
    }
    omegaPow_.resize(n / 2);
    omegaInvPow_.resize(n / 2);
    a = b = 1;
    for (std::size_t j = 0; j < n / 2; ++j) {
        omegaPow_[j] = a;
        omegaInvPow_[j] = b;
        a = F.mul(a, omega);
        b = F.mul(b, omegaInv);
    }
}

// In-place cyclic NTT of length n: A_i = sum_j a_j w^(ij), twiddles[j] = w^j (j < n/2).
void RSCode::cyclicNtt(std::vector<Elem>& a, const std::vector<Elem>& tw) const {
    const PrimeField& F = field_;
    const std::size_t n = n_;
    for (std::size_t i = 1, j = 0; i < n; ++i) {  // bit-reversal permutation
        std::size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }
    for (std::size_t len = 2; len <= n; len <<= 1) {
        const std::size_t half = len / 2, step = n / len;
        for (std::size_t i = 0; i < n; i += len) {
            for (std::size_t j = 0; j < half; ++j) {
                const Elem u = a[i + j];
                const Elem v = F.mul(a[i + j + half], tw[j * step]);
                a[i + j] = F.add(u, v);
                a[i + j + half] = F.sub(u, v);
            }
        }
    }
}

Codeword RSCode::evaluate(const Poly& f) const {
    const PrimeField& F = field_;
    std::vector<Elem> a(n_, 0);
    for (std::size_t j = 0; j < f.size(); ++j) {  // reduce mod x^n + 1 (x^n = -1)
        const std::size_t idx = j % n_;
        a[idx] = ((j / n_) & 1) ? F.sub(a[idx], f[j] % F.p()) : F.add(a[idx], f[j] % F.p());
    }
    for (std::size_t j = 0; j < n_; ++j) a[j] = F.mul(a[j], psiPow_[j]);  // twist
    cyclicNtt(a, omegaPow_);  // A_i = sum_j a_j psi^j omega^(ij) = a(psi^(2i+1))
    return a;
}

Poly RSCode::interpolate(const Codeword& c) const {
    if (c.size() != n_) throw std::invalid_argument("RSCode::interpolate: wrong length");
    const PrimeField& F = field_;
    std::vector<Elem> a = c;
    cyclicNtt(a, omegaInvPow_);
    for (std::size_t j = 0; j < n_; ++j) a[j] = F.mul(F.mul(a[j], nInv_), psiInvPow_[j]);
    return a;
}

Codeword RSCode::encodePoly(const Poly& f) const {
    if (degree(f) >= static_cast<int>(k_)) throw std::invalid_argument("RSCode::encodePoly: deg f >= k");
    return evaluate(f);
}

Codeword RSCode::encode(const Message& msg) const {
    if (msg.size() != k_) throw std::invalid_argument("RSCode::encode: message length != k");
    return evaluate(msg);
}

Message RSCode::randomMessage(std::mt19937_64& rng) const {
    std::uniform_int_distribution<Elem> sym(0, field_.p() - 1);
    Message msg(k_);
    for (auto& s : msg) s = sym(rng);
    return msg;
}

std::string RSCode::describe() const {
    std::ostringstream os;
    os << "Negacyclic RS(" << n_ << ", " << k_ << ") over GF(" << field_.p() << "), psi = " << psi_
       << ", d = " << d() << ", t = " << t() << ", R = " << rate();
    return os.str();
}

std::size_t hammingDistance(const std::vector<Elem>& a, const std::vector<Elem>& b) {
    if (a.size() != b.size()) throw std::invalid_argument("hammingDistance: length mismatch");
    std::size_t d = 0;
    for (std::size_t i = 0; i < a.size(); ++i) d += (a[i] != b[i]);
    return d;
}

}  // namespace rs
