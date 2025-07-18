#include "../common/big_int.hpp"
#include "../common/fast_exp.hpp"
#include "../common/multiplicative_inverse.hpp"

#include <iostream>

#ifndef SECP256K1
#define SECP256K1

using std::make_tuple;

static const big_int p = conv<big_int>("115792089237316195423570985008687907853269984665640564039457584007908834671663");

typedef std::pair<big_int, big_int> point;

typedef std::tuple<big_int, big_int, big_int> jacobian_point;

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
 * Multiply a point P by a scalar using the double and add method using affine point addition algorithm.
 * Efficiently computes Q = scalar * P using point doubling and conditional addition.
*/
point affine_scalar_multiplication(big_int scalar, point P)
{
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

// Jacobian point doubling: https://en.wikibooks.org/wiki/Cryptography/Prime_Curve/Jacobian_Coordinates
jacobian_point jacobian_point_doubling(jacobian_point P) {
    big_int x = std::get<0>(P);
    big_int y = std::get<1>(P);
    big_int z = std::get<2>(P);

    if (y == 0) {
        return make_tuple(conv<big_int>(0), conv<big_int>(0), conv<big_int>(0)); // Point at infinity
    }

    big_int s = mod(4 * x * (y * y), p);
    big_int m = mod(3 * (x * x), p); // a = 0 for secp256k1
    big_int x2 = mod((m * m) - 2 * s, p);
    big_int y2 = mod(m * (s - x2) - 8 * (y * y * y * y), p);
    big_int z2 = mod(2 * y * z, p);

    return make_tuple(x2, y2, z2);
}

// Jacobian point addition: https://en.wikibooks.org/wiki/Cryptography/Prime_Curve/Jacobian_Coordinates
jacobian_point jacobian_point_addition(jacobian_point P, jacobian_point Q) {
    
    big_int x1 = std::get<0>(P), y1 = std::get<1>(P), z1 = std::get<2>(P);
    big_int x2 = std::get<0>(Q), y2 = std::get<1>(Q), z2 = std::get<2>(Q);

    if (z1 == 0) return Q;
    if (z2 == 0) return P;

    big_int z1_sq = mod(z1 * z1, p);
    big_int z2_sq = mod(z2 * z2, p);

    big_int u1 = mod(x1 * z2_sq, p);
    big_int u2 = mod(x2 * z1_sq, p);

    big_int z1_cu = mod(z1_sq * z1, p);
    big_int z2_cu = mod(z2_sq * z2, p);

    big_int s1 = mod(y1 * z2_cu, p);
    big_int s2 = mod(y2 * z1_cu, p);

    big_int h = mod(u2 - u1, p);
    big_int r = mod(s2 - s1, p);

    if (h == 0) {
        if (r == 0) {
            return jacobian_point_doubling(P);
        } else {
            return make_tuple(conv<big_int>(0), conv<big_int>(0), conv<big_int>(0)); // Point at infinity
        }
    }

    big_int h_sq = mod(h * h, p);
    big_int h_cu = mod(h_sq * h, p);
    big_int u1_h_sq = mod(u1 * h_sq, p);

    big_int x3 = mod(r * r - h_cu - 2 * u1_h_sq, p);
    big_int y3 = mod(r * (u1_h_sq - x3) - s1 * h_cu, p);
    big_int z3 = mod(h * z1 * z2, p);

    return make_tuple(x3, y3, z3);
}

// Converts Jacobian (X, Y, Z) to affine (x, y): https://en.wikibooks.org/wiki/Cryptography/Prime_Curve/Jacobian_Coordinates
point convert_jacobian_to_affine(jacobian_point P) {
    big_int x = std::get<0>(P);
    big_int y = std::get<1>(P);
    big_int z = std::get<2>(P);

    if (z == 0) {
        return POINT_AT_INFINITY;
    }

    big_int z_sq = mod(z * z, p);
    big_int z_cu = mod(z_sq * z, p);

    big_int inv_z_sq = get_multiplicative_inverse(z_sq, p);
    big_int inv_z_cu = get_multiplicative_inverse(z_cu, p);

    big_int x_affine = mod(x * inv_z_sq, p);
    big_int y_affine = mod(y * inv_z_cu, p);

    return std::make_pair(x_affine, y_affine);
}

jacobian_point convert_affine_to_jacobian(const point& P) {
    return std::make_tuple(P.first, P.second, big_int(1));
}

/*
 * Multiply a point P by a scalar using the double and add method using jacobian point addition algorithm.
 * Efficiently computes Q = scalar * P using point doubling and conditional addition.
*/
jacobian_point jacobian_scalar_multiplication(big_int scalar, jacobian_point P)
{
    jacobian_point med_res = P;
    jacobian_point result = make_tuple(conv<big_int>(0), conv<big_int>(0), conv<big_int>(0));
    while (scalar > 0) {
        if (bit(scalar, 0)) {
            result = jacobian_point_addition(result, med_res);
        }
        scalar >>= 1;
        med_res = jacobian_point_doubling(med_res);
    }
    return result;
}

#endif
