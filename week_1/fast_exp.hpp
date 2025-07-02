#include <optional>

#ifndef FAST_EXPONENTIATION
#define FAST_EXPONENTIATION

int mod(int a, int n) {
    return ( (a % n) + n ) % n;
}

/*
 *  Fast exponentiation using iterative binary exponentiation.
 *  This efficiently computes (base ^ exponent) % modulo (if provided) in O(log exponent) time.
 *  When modulo is provided, intermediate values stay small and avoid overflow.
 *  Without modulo, results can quickly exceed the storage limits of the int type.
    INT_MAX is typically 2,147,483,647 for 32-bit int.
    For example: fast_exponent(10, 10) = 10,000,000,000 > INT_MAX
    This will result in integer overflow and undefined behavior in C++.
    NOTE: This does not handle negatives
*/
int fast_exponent(int base, int exponent, std::optional<int> modulo = std::nullopt)
{
    // Edge case: Any number to the power of 0 is 1.
    // Even with modulo, we return 1 % modulo (to respect modulo when provided).
    if (exponent == 0) return 1 % modulo.value_or(2);

    // Edge case: 0 raised to any positive exponent is 0.
    if (base == 0) return 0;

    // Result accumulator: starts at 1 because it's the multiplicative identity.
    int result = 1;

    // Current base: we will repeatedly square this value.
    int current_base = base;

    // Loop while exponent is not zero.
    // Each iteration processes one bit of the exponent.
    while (exponent > 0) {
        // If the current least significant bit (LSB) is 1, the exponent is odd.
        // We multiply the result by the current base.
        if (exponent & 1) {
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