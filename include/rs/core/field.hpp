#pragma once
// Prime field GF(p) = Z_p, p prime, p < 2^31.
//
// Why not GF(2^m)? Negacyclic codes over Z_p[x]/(x^n + 1) with n = 2^d need a
// primitive 2n-th root of unity, i.e. 2n | (field size - 1). In characteristic 2
// the multiplicative group has odd order 2^m - 1, so the only root of unity of
// power-of-two order is 1 (x^(2^d) + 1 = (x + 1)^(2^d)). Hence a prime field with
// p = 1 (mod 2n) is required, e.g. p = 40961 or 65537 for n = 4096.
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace rs {

using Elem = std::uint32_t;  // field element in [0, p)

class PrimeField {
public:
    explicit PrimeField(std::uint32_t p);

    std::uint32_t p() const { return p_; }

    Elem add(Elem a, Elem b) const {  // p < 2^31: a + b cannot overflow
        const Elem s = a + b;
        return s >= p_ ? s - p_ : s;
    }
    Elem sub(Elem a, Elem b) const { return a >= b ? a - b : a + p_ - b; }
    Elem neg(Elem a) const { return a ? p_ - a : 0; }
    Elem mul(Elem a, Elem b) const {
        return static_cast<Elem>(static_cast<std::uint64_t>(a) * b % p_);
    }
    Elem pow(Elem a, std::uint64_t e) const;
    Elem inv(Elem a) const;
    Elem div(Elem a, Elem b) const { return mul(a, inv(b)); }

    Elem primitiveRoot() const { return gen_; }
    // Element of exact multiplicative order `order` (throws if order does not divide p-1).
    Elem rootOfUnity(std::uint64_t order) const;

    static bool isPrime(std::uint64_t x);

private:
    std::uint32_t p_;
    Elem gen_ = 0;
    std::vector<std::uint64_t> factors_;  // distinct prime factors of p-1
};

}  // namespace rs
