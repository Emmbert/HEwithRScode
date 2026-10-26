#pragma once
// Univariate polynomials over GF(p). p[i] is the coefficient of x^i.
#include <utility>
#include <vector>
#include "rs/core/field.hpp"

namespace rs {

using Poly = std::vector<Elem>;

// Degree of p, -1 for the zero polynomial.
inline int degree(const Poly& p) {
    for (int i = static_cast<int>(p.size()) - 1; i >= 0; --i)
        if (p[static_cast<std::size_t>(i)] != 0) return i;
    return -1;
}

inline void trim(Poly& p) {
    while (!p.empty() && p.back() == 0) p.pop_back();
}

// Horner evaluation.
inline Elem eval(const PrimeField& F, const Poly& p, Elem x) {
    Elem acc = 0;
    for (std::size_t i = p.size(); i-- > 0;) acc = F.add(F.mul(acc, x), p[i]);
    return acc;
}

// Polynomial division: returns quotient, writes remainder to rem.
inline Poly divmod(const PrimeField& F, Poly num, const Poly& den, Poly& rem) {
    const int dd = degree(den);
    if (dd < 0) throw std::domain_error("poly divmod: division by zero polynomial");
    const int nd = degree(num);
    Poly quo;
    if (nd >= dd) {
        quo.assign(static_cast<std::size_t>(nd - dd + 1), 0);
        const Elem invLead = F.inv(den[static_cast<std::size_t>(dd)]);
        for (int i = nd; i >= dd; --i) {
            const Elem c = num[static_cast<std::size_t>(i)];
            if (c == 0) continue;
            const Elem f = F.mul(c, invLead);
            quo[static_cast<std::size_t>(i - dd)] = f;
            for (int j = 0; j <= dd; ++j) {
                Elem& t = num[static_cast<std::size_t>(i - dd + j)];
                t = F.sub(t, F.mul(f, den[static_cast<std::size_t>(j)]));
            }
        }
    }
    trim(num);
    rem = std::move(num);
    return quo;
}

}  // namespace rs
