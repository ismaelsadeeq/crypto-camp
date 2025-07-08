#include "fast_exp.hpp"
#include "multiplicative_inverse.hpp"

#include <utility>

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
std::pair<NTL::ZZ, NTL::ZZ> basic_elgamal_encrypt(NTL::ZZ public_key, NTL::ZZ message, NTL::ZZ ephemeral_key, NTL::ZZ generator, NTL::ZZ prime_field)
{
    NTL::ZZ c1 = fast_exponent(generator, ephemeral_key, prime_field);
    NTL::ZZ public_key_ephemeral = fast_exponent(public_key, ephemeral_key, prime_field);
    NTL::ZZ cipher_component = message * public_key_ephemeral;
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
NTL::ZZ basic_elgamal_decrypt(NTL::ZZ c1, NTL::ZZ c2, NTL::ZZ private_key, NTL::ZZ prime_field)
{
    NTL::ZZ shared_secret = fast_exponent(c1, private_key, prime_field);
    NTL::ZZ shared_secret_inverse = get_multiplicative_inverse(shared_secret, prime_field);
    NTL::ZZ message = shared_secret_inverse * c2;
    return mod(message, prime_field);
}

#endif
