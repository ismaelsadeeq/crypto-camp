#include "fast_exp.hpp"

#ifndef MULTIPLICATIVE_INVERSE
#define MULTIPLICATIVE_INVERSE

/* This function returns the multiplicative inverse of 
   an integer a in ZP^* 

   It takes the integer a and the prime number p
   and returns the multiplicative inverse by applying
   Fermat's Little Theorem which states that
   a^(p-1) ≡ 1 mod p
   which means that a * a^(p-2) ≡ 1 mod p
   and since the group ZP* has 1 as the identity,
   the multiplicative inverse is a^(p-2) mod p.
 */
int get_multiplicative_inverse(int a, int p)
{
    return fast_exponent(a, p - 2, p);
}

#endif
