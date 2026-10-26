// Tests for the negacyclic RS code: locators, NTT, inverse NTT, ring products.
#include <random>

#include "rs/core/rs_code.hpp"
#include "test_util.hpp"

namespace {

// Naive product in Z_p[x]/(x^n + 1).
rs::Poly negacyclicProduct(const rs::PrimeField& F, const rs::Poly& a, const rs::Poly& b) {
    const std::size_t n = a.size();
    rs::Poly c(n, 0);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            const rs::Elem t = F.mul(a[i], b[j]);
            const std::size_t idx = (i + j) % n;
            c[idx] = (i + j >= n) ? F.sub(c[idx], t) : F.add(c[idx], t);
        }
    return c;
}

}  // namespace

int main() {
    using namespace rs;
    RS_CHECK(RSCode::smallestNegacyclicPrime(4096) == 40961);
    RS_CHECK(RSCode::smallestNegacyclicPrime(1024) == 12289);

    std::mt19937_64 rng(7);
    for (std::size_t n : {8u, 16u, 256u, 4096u}) {
        RSCode code(RSCode::smallestNegacyclicPrime(n), n, n / 2);
        const PrimeField& F = code.field();
        std::uniform_int_distribution<Elem> sym(0, F.p() - 1);
        Poly a(n), b(n);
        for (auto& x : a) x = sym(rng);
        for (auto& x : b) x = sym(rng);

        // Locators are the n roots of x^n + 1 and pairwise distinct
        std::vector<bool> seen(F.p(), false);
        for (Elem x : code.evalPoints()) {
            RS_CHECK(F.pow(x, n) == F.p() - 1);
            RS_CHECK(!seen[x]);
            seen[x] = true;
        }

        // NTT == direct evaluation, inverse NTT round-trips
        const Codeword A = code.evaluate(a);
        for (std::size_t i = 0; i < n && i < 64; ++i) RS_CHECK(A[i] == eval(F, a, code.evalPoints()[i]));
        RS_CHECK(code.interpolate(A) == a);

        // Encoding produces codewords, random words are (almost surely) not
        RS_CHECK(code.isCodeword(code.encode(code.randomMessage(rng))));
        RS_CHECK(!code.isCodeword(A));

        // Reduction mod x^n + 1: x^n evaluates like -1
        Poly xn(n + 1, 0);
        xn[n] = 1;
        RS_CHECK(code.evaluate(xn) == code.evaluate(Poly{F.p() - 1}));

        // Ring product <-> pointwise product
        if (n <= 256) {
            const Codeword B = code.evaluate(b), C = code.evaluate(negacyclicProduct(F, a, b));
            for (std::size_t i = 0; i < n; ++i) RS_CHECK(C[i] == F.mul(A[i], B[i]));
        }
        std::cout << "[DONE] " << code.describe() << "\n";
    }

    // Invalid parameters must throw
    auto throws = [](auto&& f) { try { f(); } catch (const std::invalid_argument&) { return true; } return false; };
    RS_CHECK(throws([] { RSCode c(40961, 3000, 100); }));  // n not a power of two
    RS_CHECK(throws([] { RSCode c(12289, 4096, 100); }));  // 2n does not divide p-1
    RS_CHECK(throws([] { RSCode c(40961, 4096, 0); }));    // k = 0

    return rs::test::finish();
}
