#include "iostream"
#include "fast_exp.hpp"
#include "multiplicative_inverse.hpp"
#include "elgamal.hpp"


void fast_exp_tests() {
    // ================================
    // Test without modulo (results within int range)
    assert(fast_exponent(3, 218, 1000) == 489);
    assert(fast_exponent(2, 10) == 1024);
    assert(fast_exponent(3, 3) == 27);
    assert(fast_exponent(5, 0) == 1);
    assert(fast_exponent(0, 5) == 0);
    assert(fast_exponent(0, 0) == 1); 

    // Test with modulo
    assert(fast_exponent(3, 2, 4) == 1);
    assert(fast_exponent(2, 10, 1000) == 24); // 1024 % 1000 = 24
    assert(fast_exponent(7, 4, 5) == 1);      // 2401 % 5 = 1
    assert(fast_exponent(10, 0, 6) == 1 % 6); // Should be 1 % 6 = 1
    assert(fast_exponent(0, 5, 7) == 0);      // Base zero, should return 0
    assert(fast_exponent(-3, 7, 7) == 4);
    std::cout << "All fast exponent test cases passed!" << std::endl;
}

void test_multiplicative_inverse() {
    // ============================================
    // Well-known small primes
    assert(get_multiplicative_inverse(2, 7) == 4); // 2 * 4 ≡ 1 mod 7
    assert(get_multiplicative_inverse(3, 7) == 5); // 3 * 5 ≡ 15 ≡ 1 mod 7
    assert(get_multiplicative_inverse(4, 7) == 2); // 4 * 2 ≡ 8 ≡ 1 mod 7
    assert(get_multiplicative_inverse(5, 7) == 3); // 5 * 3 ≡ 15 ≡ 1 mod 7
    assert(get_multiplicative_inverse(6, 7) == 6); // Self-inverse: 6 * 6 ≡ 36 ≡ 1 mod 7

    // ============================================
    // Edge cases
    assert(get_multiplicative_inverse(1, 7) == 1); // 1 is its own inverse
    assert(get_multiplicative_inverse(6, 7) == 6); // p - 1 is always self-inverse in prime mod

    // ============================================
    // Medium primes
    assert(get_multiplicative_inverse(10, 17) == 12); // Verified: 10 * 12 ≡ 120 ≡ 1 mod 17
    assert(get_multiplicative_inverse(4, 17) == 13);  // Verified: 4 * 13 ≡ 52 ≡ 1 mod 17

    // // ============================================
    // TODO Large primes
    // assert(get_multiplicative_inverse(123, 1009) == 824);
    // assert(get_multiplicative_inverse(456, 10007) == 9847);

    // ============================================
    // TODO: Negative base (if your fast_exponent handles negatives correctly)
    assert(get_multiplicative_inverse(-3, 7) == 2); // -3 ≡ 4 mod 7, inverse of 4 is 2

    // ============================================
    std::cout << "All multiplicative inverse test cases passed!" << std::endl;
}

void test_elgamal_vectors() {
    // Test Vectors: p, g, x, k, m, expected c1 (a), expected c2 (b)
    struct TestVector {
        int p;
        int g;
        int x;
        int k;
        int m;
        int expected_c1;
        int expected_c2;
    };

    // TODO: Handle large primes
    TestVector vectors[] = {
        {71, 33, 62, 31, 15, 62, 18},
        {23, 11, 6, 3, 10, 20, 22},
        {809, 3, 68, 89, 100, 345, 517},
        {17, 6, 5, 10, 13, 15, 9}
    };

    for (const auto& test : vectors) {
        int public_key = fast_exponent(test.g, test.x, test.p);

        auto cipher = basic_elgamal_encrypt(public_key, test.m, test.k, test.g, test.p);

        assert(cipher.first == test.expected_c1);
        assert(cipher.second == test.expected_c2);

        int decrypted_message = basic_elgamal_decrypt(cipher.first, cipher.second, test.x, test.p);

        assert(decrypted_message == test.m);
    }

    // ============================================
    std::cout << "All elgamal encryption and decryption test cases passed!" << std::endl;
}

int main() {
    fast_exp_tests();
    test_multiplicative_inverse();
    test_elgamal_vectors();
    return 0;
}