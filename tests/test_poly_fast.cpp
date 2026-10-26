// Tests for FastPolyArith: every fast routine is compared with its classical
// counterpart on random inputs, including products longer than the largest NTT.
#include <random>

#include "rs/core/poly_fast.hpp"
#include "test_util.hpp"

namespace {

rs::Poly randomPoly(std::size_t len, const rs::PrimeField& F, std::mt19937_64& rng) {
    std::uniform_int_distribution<rs::Elem> sym(0, F.p() - 1);
    rs::Poly a(len);
    for (auto& x : a) x = sym(rng);
    if (len) a.back() = 1 + sym(rng) % (F.p() - 1);  // exact degree len-1
    return a;
}

// Classical partial extended Euclid (reference implementation).
rs::PartialGcdResult classicalPartialGcd(const rs::FastPolyArith& A, rs::Poly a, rs::Poly b, int s) {
    rs::PartialGcdResult res;
    res.M = rs::PolyMatrix::identity();
    while (rs::degree(b) >= s) {
        rs::Poly r;
        rs::Poly q = rs::divmod(A.field(), a, b, r);
        rs::PolyMatrix N;
        N.m00 = res.M.m10;
        N.m01 = res.M.m11;
        N.m10 = A.sub(res.M.m00, A.mulSchoolbook(q, res.M.m10));
        N.m11 = A.sub(res.M.m01, A.mulSchoolbook(q, res.M.m11));
        res.M = N;
        a = std::move(b);
        b = std::move(r);
    }
    res.r0 = a;
    res.r1 = b;
    return res;
}

}  // namespace

int main() {
    using namespace rs;
    std::mt19937_64 rng(11);

    // p = 17: maxNtt = 16 (block multiplication kicks in early),
    // p = 40961: maxNtt = 8192 (exactly 2n for n = 4096)
    for (std::uint32_t p : {17u, 7681u, 40961u, 65537u}) {
        PrimeField F(p);
        FastPolyArith A(F);
        std::cout << "p = " << p << ", max NTT length = " << A.maxNttLength() << "\n";

        // Multiplication (incl. lengths beyond maxNttLength)
        for (std::size_t na : {1u, 5u, 40u, 300u, 3000u, 6000u}) {
            for (std::size_t nb : {1u, 33u, 500u, 5000u}) {
                if (na * nb > 4000000u) continue;
                const Poly a = randomPoly(na, F, rng), b = randomPoly(nb, F, rng);
                RS_CHECK(A.mul(a, b) == A.mulSchoolbook(a, b));
            }
        }

        // Power series inverse
        {
            const Poly b = randomPoly(700, F, rng);
            if (b[0] != 0) {
                Poly prod = A.mul(b, A.invSeries(b, 500));
                prod.resize(500, 0);
                Poly one(500, 0);
                one[0] = 1;
                RS_CHECK(prod == one);
            }
        }

        // Division
        for (auto [na, nb] : {std::pair{2000u, 700u}, std::pair{2000u, 1990u}, std::pair{900u, 40u},
                              std::pair{100u, 300u}, std::pair{5000u, 2500u}}) {
            const Poly a = randomPoly(na, F, rng), b = randomPoly(nb, F, rng);
            Poly r1, r2;
            const Poly q1 = A.divmod(a, b, r1);
            const Poly q2 = rs::divmod(F, a, b, r2);
            RS_CHECK(q1 == q2);
            RS_CHECK(r1 == r2);
        }

        // Partial GCD: fast == classical (remainders and cofactor matrix)
        for (std::size_t n : {64u, 300u, 1024u, 2048u}) {
            if (n > 1024 && p < 7681) continue;
            Poly g0(n + 1, 0);
            g0[0] = 1;
            g0[n] = 1;  // x^n + 1, as in the decoder
            const Poly g1 = randomPoly(n - 1, F, rng);
            for (int s : {static_cast<int>(n) / 2, static_cast<int>(3 * n / 4), static_cast<int>(n) - 3, 5}) {
                const PartialGcdResult f = A.partialGcd(g0, g1, s);
                const PartialGcdResult c = classicalPartialGcd(A, g0, g1, s);
                RS_CHECK(f.r0 == c.r0);
                RS_CHECK(f.r1 == c.r1);
                RS_CHECK(f.M.m10 == c.M.m10);
                RS_CHECK(f.M.m11 == c.M.m11);
                RS_CHECK(degree(f.r0) >= s && degree(f.r1) < s);
            }
        }
    }
    return rs::test::finish();
}
