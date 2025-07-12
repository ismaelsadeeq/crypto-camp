#include "../common/big_int.hpp"
#include "../common/fast_exp.hpp"
#include "../common/multiplicative_inverse.hpp"

#include <iostream>

#ifndef SECP256K1
#define SECP256K1

static const big_int p = conv<big_int>("115792089237316195423570985008687907853269984665640564039457584007908834671663");

typedef std::pair<big_int, big_int> point;

static const point G = {conv<big_int>("55066263022277343669578718895168534326250603453777594175500187360389116729240"), conv<big_int>("32670510020758816978083085130507043184471273380659243275938904335757337482424")};
point POINT_AT_INFINITY{conv<big_int>(0), conv<big_int>(0)};

bool point_are_equal(point& a, point& b)
{
    return a.first == b.first && a.second == b.second;
}

bool point_at_infinity(point& a)
{
    return point_are_equal(a, POINT_AT_INFINITY);
}

/*
    The curve has points whose coordinates are elements of a finite field F=(Z/pZ,+,×),
    where p = 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f
    The curve equation follows y^2 = x^3 + 7
*/

point affine_point_addition(point P, point Q) {
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


/*
 * Multiply a point P by a scalar using the double and add method.
 * Efficiently computes Q = scalar * P using point doubling and conditional addition.
*/
point affine_scalar_multiplication(big_int scalar, point P)
{
    big_int base = conv<big_int>(2);
    big_int current_base = base;
    point med_res = P;
    point result = POINT_AT_INFINITY;
    while (scalar > 0) {
        if (bit(scalar, 0)) {
            result = affine_point_addition(result, med_res);
        }
        scalar >>= 1;
        med_res = affine_point_addition(med_res, med_res);
    }
    return result;
}

#endif
