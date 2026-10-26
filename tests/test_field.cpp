// Unit tests for PrimeField and polynomial helpers.
#include "rs/core/field.hpp"
#include "rs/core/poly.hpp"
#include "test_util.hpp"

int main() {
    using namespace rs;

    RS_CHECK(PrimeField::isPrime(40961));
    RS_CHECK(PrimeField::isPrime(65537));
    RS_CHECK(!PrimeField::isPrime(8193));

    for (std::uint32_t p : {17u, 97u, 7681u, 12289u, 40961u, 65537u}) {
        PrimeField F(p);
        // Field axioms on a few elements
        for (Elem a = 1; a < 50 && a < p; ++a) {
            RS_CHECK(F.mul(a, F.inv(a)) == 1);
            RS_CHECK(F.add(a, F.neg(a)) == 0);
            RS_CHECK(F.sub(F.add(a, 7 % p), 7 % p) == a);
        }
        // Generator has order exactly p-1
        RS_CHECK(F.pow(F.primitiveRoot(), p - 1) == 1);
        RS_CHECK(F.pow(F.primitiveRoot(), (p - 1) / 2) == p - 1);
    }

    // Invalid inputs must throw
    bool threw = false;
    try { PrimeField F(15); } catch (const std::invalid_argument&) { threw = true; }
    RS_CHECK(threw);

    // Polynomial division: (x^2 - 1) / (x - 1) = x + 1 over GF(17)
    PrimeField F(17);
    Poly num{16, 0, 1}, den{16, 1}, rem;
    Poly quo = divmod(F, num, den, rem);
    RS_CHECK(quo == (Poly{1, 1}));
    RS_CHECK(rem.empty());
    RS_CHECK(eval(F, num, 1) == 0);

    return rs::test::finish();
}
