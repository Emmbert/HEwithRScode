#include "rs/core/field.hpp"

namespace rs {

bool PrimeField::isPrime(std::uint64_t x) {
    if (x < 2) return false;
    for (std::uint64_t d = 2; d * d <= x; ++d)
        if (x % d == 0) return false;
    return true;
}

PrimeField::PrimeField(std::uint32_t p) : p_(p) {
    if (p < 3 || p >= (1u << 31) || !isPrime(p))
        throw std::invalid_argument("PrimeField: p must be an odd prime < 2^31");

    std::uint64_t r = p - 1;  // factor p-1
    for (std::uint64_t d = 2; d * d <= r; ++d) {
        if (r % d == 0) {
            factors_.push_back(d);
            while (r % d == 0) r /= d;
        }
    }
    if (r > 1) factors_.push_back(r);

    for (Elem g = 2; g < p; ++g) {  // smallest generator of Z_p^*
        bool ok = true;
        for (std::uint64_t f : factors_)
            if (pow(g, (p - 1) / f) == 1) { ok = false; break; }
        if (ok) { gen_ = g; break; }
    }
}

Elem PrimeField::pow(Elem a, std::uint64_t e) const {
    Elem result = 1 % p_, base = a % p_;
    while (e) {
        if (e & 1) result = mul(result, base);
        base = mul(base, base);
        e >>= 1;
    }
    return result;
}

Elem PrimeField::inv(Elem a) const {
    if (a % p_ == 0) throw std::domain_error("PrimeField::inv(0)");
    return pow(a, p_ - 2);
}

Elem PrimeField::rootOfUnity(std::uint64_t order) const {
    if (order == 0 || (p_ - 1) % order != 0)
        throw std::invalid_argument("PrimeField::rootOfUnity: order does not divide p-1");
    return pow(gen_, (p_ - 1) / order);
}

}  // namespace rs
