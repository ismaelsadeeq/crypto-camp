#include <NTL/ZZ.h>
#include <optional>

#ifndef FAST_EXPONENTIATION
#define FAST_EXPONENTIATION

NTL::ZZ mod(NTL::ZZ a, NTL::ZZ n) {
    return ( (a % n) + n ) % n;
}

/*
 *  Fast exponentiation using iterative binary exponentiation.
 *  This efficiently computes (base ^ exponent) % modulo (if provided) in O(log exponent) time.
 *  When modulo is provided, intermediate values stay small and avoid overflow.
 *  Without modulo, results can quickly exceed the storage limits of the int type.
*/
NTL::ZZ fast_exponent(NTL::ZZ base, NTL::ZZ exponent, std::optional<NTL::ZZ> modulo = std::nullopt){
    // Edge case: Any number to the power of 0 is 1.
    // Even with modulo, we return 1 % modulo (to respect modulo when provided).
    if (exponent == 0) return NTL::conv<NTL::ZZ>(1) % modulo.value_or(NTL::conv<NTL::ZZ>(2));

    // Edge case: 0 raised to any positive exponent is 0.
    if (base == 0) return NTL::conv<NTL::ZZ>(0);

    // Result accumulator: starts at 1 because it's the multiplicative identity.
    NTL::ZZ result = NTL::conv<NTL::ZZ>(1);

    // Current base: we will repeatedly square this value.
    NTL::ZZ current_base = base;

    // Loop while exponent is not zero.
    // Each iteration processes one bit of the exponent.
    while (exponent > 0) {
        // If the current least significant bit (LSB) is 1, the exponent is odd.
        // We multiply the result by the current base.
        if (bit(exponent, 0)) {
            result *= current_base;

            // Apply modulo if provided, to prevent integer overflow.
            if (modulo) result = mod(result, modulo.value());
        }

        // Square the base for the next bit position.
        current_base = current_base * current_base;

        // Apply modulo if provided, to keep numbers small and avoid overflow.
        if (modulo) current_base = mod(current_base, modulo.value());

        // Move to the next higher bit by right-shifting the exponent.
        exponent >>= 1;
    }
    return result;
}
#endif