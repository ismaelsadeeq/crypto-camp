#include "../common/big_int.hpp"
#include "../common/fast_exp.hpp"
#include "../common/multiplicative_inverse.hpp"

#include "utility"

#ifndef ELGAMAL_ENCRYPTION
#define ELGAMAL_ENCRYPTION

/*
 * Basic ElGamal Encryption Algorithm
 *
 * This implementation provides a simple version of the ElGamal encryption 
 * scheme, which is based on the Diffie-Hellman key exchange concept.
 *
 * Encryption:
 * - Computes the first part  (c1): generator^ephemeral_key mod prime_field
 * - Computes the second part (c2): message * (public_key^ephemeral_key) mod prime_field
 *
 * The exponentiation uses the binary (fast) exponentiation algorithm.
 *
 * Returns a pair consisting of c1 and c2.
 */
std::pair<big_int, big_int> basic_elgamal_encrypt(big_int public_key, big_int message, big_int ephemeral_key, big_int generator, big_int prime_field)
{
    big_int c1 = fast_exponent(generator, ephemeral_key, prime_field);
    big_int public_key_ephemeral = fast_exponent(public_key, ephemeral_key, prime_field);
    big_int cipher_component = message * public_key_ephemeral;
    return std::make_pair(c1, mod(cipher_component, prime_field));
}

/*
 * Basic ElGamal Decryption Algorithm
 *
 * This function decrypts the ciphertext using the private key and prime field.
 *
 * Decryption:
 * - Computes x = c1^private_key mod prime_field
 * - Computes the multiplicative inverse of x modulo prime_field
 * - Recovers the original message: (multiplicative_inverse * c2) mod prime_field
 * 
 * Returns:
 * - The decrypted message.
 */
big_int basic_elgamal_decrypt(big_int c1, big_int c2, big_int private_key, big_int prime_field)
{
    big_int shared_secret = fast_exponent(c1, private_key, prime_field);
    big_int shared_secret_inverse = get_multiplicative_inverse(shared_secret, prime_field);
    big_int message = shared_secret_inverse * c2;
    return mod(message, prime_field);
}

#endif
