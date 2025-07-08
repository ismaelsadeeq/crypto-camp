#include "fast_exp.hpp"
#include "multiplicative_inverse.hpp"
#include "elgamal.hpp"

#include <NTL/ZZ.h>
#include "cassert"
#include "iostream"

void fast_exp_tests() {
    // ================================
    // Test without modulo (small range)
    assert(fast_exponent(NTL::ZZ(3), NTL::ZZ(218), NTL::ZZ(1000)) == NTL::ZZ(489));
    assert(fast_exponent(NTL::ZZ(2), NTL::ZZ(10)) == NTL::ZZ(1024));
    assert(fast_exponent(NTL::ZZ(3), NTL::ZZ(3)) == NTL::ZZ(27));
    assert(fast_exponent(NTL::ZZ(5), NTL::ZZ(0)) == NTL::ZZ(1));
    assert(fast_exponent(NTL::ZZ(0), NTL::ZZ(5)) == NTL::ZZ(0));
    assert(fast_exponent(NTL::ZZ(0), NTL::ZZ(0)) == NTL::ZZ(1));  // defined as 1 in your implementation

    // ================================
    // With modulo
    assert(fast_exponent(NTL::ZZ(3), NTL::ZZ(2), NTL::ZZ(4)) == NTL::ZZ(1));
    assert(fast_exponent(NTL::ZZ(2), NTL::ZZ(10), NTL::ZZ(1000)) == NTL::ZZ(24)); // 1024 % 1000 = 24
    assert(fast_exponent(NTL::ZZ(7), NTL::ZZ(4), NTL::ZZ(5)) == NTL::ZZ(1));      // 2401 % 5 = 1
    assert(fast_exponent(NTL::ZZ(10), NTL::ZZ(0), NTL::ZZ(6)) == NTL::ZZ(1));     // 1 % 6 = 1
    assert(fast_exponent(NTL::ZZ(0), NTL::ZZ(5), NTL::ZZ(7)) == NTL::ZZ(0));      // Base zero, result 0

    // Negative base, modulo
    assert(fast_exponent(NTL::ZZ(-3), NTL::ZZ(7), NTL::ZZ(7)) == NTL::ZZ(4));     // (-3)^7 = -2187 % 7 = 4

    // ================================
    // Test with huge numbers
    NTL::ZZ huge_base = NTL::conv<NTL::ZZ>("84265675725482892459719348378630146162719620409152809167814480007059199482163");
    NTL::ZZ huge_exp  = NTL::conv<NTL::ZZ>("123456789123456789");
    NTL::ZZ huge_mod  = NTL::conv<NTL::ZZ>("115792089237316195423570985008687907853269984665640564039457584007913129639937"); // 2^256 - 1

    NTL::ZZ huge_result = NTL::conv<NTL::ZZ>("27255166557789855811328668216188941042231292868742423663488474786644472725517");
    assert(fast_exponent(huge_base, huge_exp, huge_mod) == huge_result);

    std::cout << "All fast exponent test cases passed!" << std::endl;
}

void test_multiplicative_inverse() {
    // ============================================
    // Well-known small primes
    assert(get_multiplicative_inverse(NTL::ZZ(2), NTL::ZZ(7)) == NTL::ZZ(4)); // 2 * 4 ≡ 1 mod 7
    assert(get_multiplicative_inverse(NTL::ZZ(3), NTL::ZZ(7)) == NTL::ZZ(5)); // 3 * 5 ≡ 1 mod 7
    assert(get_multiplicative_inverse(NTL::ZZ(4), NTL::ZZ(7)) == NTL::ZZ(2));
    assert(get_multiplicative_inverse(NTL::ZZ(5), NTL::ZZ(7)) == NTL::ZZ(3));
    assert(get_multiplicative_inverse(NTL::ZZ(6), NTL::ZZ(7)) == NTL::ZZ(6)); // self-inverse

    // ============================================
    // Edge cases
    assert(get_multiplicative_inverse(NTL::ZZ(1), NTL::ZZ(7)) == NTL::ZZ(1));
    assert(get_multiplicative_inverse(NTL::ZZ(6), NTL::ZZ(7)) == NTL::ZZ(6));

    // ============================================
    // Medium primes
    assert(get_multiplicative_inverse(NTL::ZZ(10), NTL::ZZ(17)) == NTL::ZZ(12));
    assert(get_multiplicative_inverse(NTL::ZZ(4), NTL::ZZ(17)) == NTL::ZZ(13));

    // ============================================
    // Large primes
    assert(get_multiplicative_inverse(NTL::ZZ(123), NTL::ZZ(1009)) == NTL::ZZ(484));
    assert(get_multiplicative_inverse(NTL::ZZ(456), NTL::ZZ(10007)) == NTL::ZZ(1602));

    // ============================================
    // Negative base
    assert(get_multiplicative_inverse(NTL::ZZ(-3), NTL::ZZ(7)) == NTL::ZZ(2)); // -3 ≡ 4 mod 7; inverse of 4 is 2

    std::cout << "All multiplicative inverse test cases passed!" << std::endl;
}

void test_elgamal() {
    struct Vector {
        NTL::ZZ p, g, x, k, m, a, b;
    };

    std::vector<Vector> test_vectors = {
        {NTL::ZZ(71), NTL::ZZ(33), NTL::ZZ(62), NTL::ZZ(31), NTL::ZZ(15), NTL::ZZ(62), NTL::ZZ(18)},
        {NTL::ZZ(23), NTL::ZZ(11), NTL::ZZ(6), NTL::ZZ(3), NTL::ZZ(10), NTL::ZZ(20), NTL::ZZ(22)},
        {NTL::ZZ(809), NTL::ZZ(3), NTL::ZZ(68), NTL::ZZ(89), NTL::ZZ(100), NTL::ZZ(345), NTL::ZZ(517)},
        {NTL::ZZ(17), NTL::ZZ(6), NTL::ZZ(5), NTL::ZZ(10), NTL::ZZ(13), NTL::ZZ(15), NTL::ZZ(9)},
        {
            NTL::conv<NTL::ZZ>("84265675725482892459719348378630146162719620409152809167814480007059199482163"),
            NTL::ZZ(5),
            NTL::conv<NTL::ZZ>("2799014790424892046701478888900891009403869701173893426"),
            NTL::conv<NTL::ZZ>("23517683968368899022119256606644551548285683288848885921"),
            NTL::conv<NTL::ZZ>("87521618088882658227876453"),
            NTL::conv<NTL::ZZ>("22954586883013884818653063688294540134886732496160582262267014428782771199687"),
            NTL::conv<NTL::ZZ>("56046128113101346099694619669629128017849277484825379502821514323706183544424")
        },
        {
            NTL::conv<NTL::ZZ>("12658517083168187407924345155971956101250996576825115113297001855799796437288935576230034157578333666497170430505565580165565829633685607504706642034926119"),
            NTL::conv<NTL::ZZ>(7),
            NTL::conv<NTL::ZZ>("2001688878140630728014209681954697141876038523595247208"),
            NTL::conv<NTL::ZZ>("5446024688717452254835115775456957961297236108858862823"),
            NTL::conv<NTL::ZZ>("87521618088882658227876453"),
            NTL::conv<NTL::ZZ>("2150519483988769855483983776372336742288374425191291528256965705108393490638750082340115568718132372731853110762124400441550538499580316268601341087676203"),
            NTL::conv<NTL::ZZ>("1540471266850557563382406324432354117072109094950140952195099581066490559252112349492583688225692526496193879919152401794896907007394565292272866724291488")
        }
    };

    for (size_t i = 0; i < test_vectors.size(); ++i) {
        const auto& v = test_vectors[i];

        NTL::ZZ public_key = fast_exponent(v.g, v.x, v.p);
        auto [c1, c2] = basic_elgamal_encrypt(public_key, v.m, v.k, v.g, v.p);

        // Check against known ciphertexts
        assert(c1 == v.a);
        assert(c2 == v.b);

        // Now decrypt and check the original message
        NTL::ZZ decrypted = basic_elgamal_decrypt(c1, c2, v.x, v.p);
        assert(decrypted == v.m);

        std::cout << "Test vector " << i + 1 << " passed" << std::endl;
    }

    std::cout << "All ElGamal test vectors passed successfully!" << std::endl;
}

int main() {
    fast_exp_tests();
    test_multiplicative_inverse();
    test_elgamal();
    return 0;
}