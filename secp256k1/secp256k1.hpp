#include "../common/big_int.hpp"
#include "../common/fast_exp.hpp"
#include "../common/multiplicative_inverse.hpp"

#ifndef SECP256K1
#define SECP256K1

static const big_int p = conv<big_int>("115792089237316195423570985008687907853269984665640564039457584007908834671663");

typedef std::pair<big_int, big_int> point;

bool point_at_infinity(point& a)
{
    return a.first == 0 && a.second == 0;
}

bool point_are_equal(point& a, point& b)
{
    return a.first == b.first && a.second == b.second;
}

/*
    The curve has points whose coordinates are elements of a finite field F=(Z/pZ,+,×),
    where p = 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f
    The curve equation follows y^2 = x^3 + 7
*/

std::pair<big_int, big_int> affine_point_addition(point P, point Q) {
    // Handle point at infinity (identity element in elliptic curve group)
    if (point_at_infinity(P)) return Q;
    if (point_at_infinity(Q)) return P;

    // Check if P and Q are vertical reflections of each other
    if (P.first == Q.first && (mod(P.second + Q.second, p) == 0)) {
        return std::make_pair(big_int(0), big_int(0));  // Point at infinity
    }

    big_int lambda;

    if (point_are_equal(P, Q)) {
        // Point doubling: P == Q
        // lambda = (3 * x^2 + a) / (2 * y) mod p
        big_int numerator = mod(3 * fast_exponent(P.first, conv<big_int>(2), p), p);
        big_int denominator = mod(2 * P.second, p);
        lambda = mod(numerator * get_multiplicative_inverse(denominator, p), p);
    } else {
        // Point addition: P != Q
        // lambda = (y2 - y1) / (x2 - x1) mod p
        big_int numerator = mod(Q.second - P.second, p);
        big_int denominator = mod(Q.first - P.first, p);
        lambda = mod(numerator * get_multiplicative_inverse(denominator, p), p);
    }

    // Compute resulting x and y coordinates
    // x3 = lambda^2 - (x1 + x2) mod p
    big_int x3 = mod(fast_exponent(lambda, conv<big_int>(2), p) - (P.first + Q.first), p);

    // y3 = lambda * (x1 - x3) - y1 mod p
    big_int y3 = mod(lambda * (P.first - x3) - P.second, p);

    return std::make_pair(x3, y3);
}

#endif
