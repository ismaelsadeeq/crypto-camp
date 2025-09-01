### Secp256k1 Sessions by sipa

#### Background
- [Abstract Algebra (Cyclic groups, Isomorphism, Homomorphism, Endomorphism, Finite Fields, and Scalars)](https://cryptocamp.website/t/secp256k1-part-1-abstract-algebra-overview/32)
- [Secp256k1 curve, Weierstrass Equation, Point Operations in the group, Scalar Multiplication, Isomorphism between scalars and points, The GLV endomorphism, Jacobian coordinates](https://cryptocamp.website/t/secp256k1-part-2-the-secp256k1-group/66)
- [Prep Exercises](https://github.com/ismaelsadeeq/crypto-camp-prep)

#### Session

- Point addition: $a=(x_1, y_1) + b=(x_2, y_2)$
  The Weierstrass equation is simplified to create a formula for calculating the sum of two points. When the points are identical, the point doubling formula is used. These points are known as affine coordinates, with affine point addition and point doubling operations.

- Dividing two coordinates $a/b$
  This is equivalent to $a \cdot b^{-1}$, requiring calculation of the multiplicative inverse. The naive approach using Fermat's Little Theorem (FLT) is slow.

- Converting affine coordinates to Jacobian coordinates to avoid division
  Affine coordinates convert to Jacobian coordinates represented as $(X, Y, Z)$, which corresponds to affine $(X/Z^2, Y/Z^3)$. This avoids expensive inversions, providing 50x multiplication improvement for variable-time operations and 80x for constant-time operations.

- GLV endomorphism
  $\lambda \in GF(n)$, $\beta \in GF(p)$ are cube roots of 1 in respective fields.
  For any point on the curve, $\lambda(x,y) = (\beta x, y)$.
  This means you don't need to find the y-coordinate.
  This is a property of the secp256k1 elliptic curve because $n \equiv p \equiv 1 \pmod{3}$, but not applicable to all curves.

- Batch Inversion avoids division when you want the inverse of multiple points
  If you want to get 2 inverses, batch inverse $(a,b)^{-1} = (a^{-1}, b^{-1})$.
  To get $a^{-1}$ and $b^{-1}$, compute $(a \times b)^{-1}$, then $a^{-1} = (a \times b)^{-1} \times b$ and $b^{-1} = (a \times b)^{-1} \times a$.
  For a batch of $n$ inverses, only need 1 inverse and $3(n-1)$ multiplications.
  This means you can convert many Jacobian coordinates to affine ones together.
  Toy example use case (don't do it): looking for a vanity address by trying many coordinates in batches.

- X-point only operations (ECDH): You can avoid sending the y-coordinate because you can derive it using the Weierstrass equation $(\pm y) = \sqrt{x^3 + 7}$. On our curve, the square root of $x$ is $x^{(p+1)/4}$ because $p \equiv 3 \pmod{4}$. However, you want to avoid square root calculations. This uses only 64 bytes and doesn't reveal as much about the connection.
  See Elligator Swift.

- Detecting whether a given x-coordinate is on the curve
- Avoiding square root computation to determine if an x-coordinate belongs to the secp256k1 curve
- Quadratic twist of secp256k1 curve

- Speeding operations by isomorphic mapping: Map secp256k1 curve points to another curve, perform operations, and map back

- BIP 340 Tie-breaker rule
  If we only want the x-coordinate of $a \cdot P$ and only know the x-coordinate of $P$, we never need to calculate the y-coordinate (avoiding square root).
  
  Key: Define isomorphism $f_R(X, Y) = (R^2X, R^3Y)$ mapping to curve $Y^2 = X^3 + 7R^6$
  
  Process: Do isomorphism with $R =$ y-coordinate → multiply by $a$ → convert back without y-coordinate: $X/(Z^2 \cdot X^3 + 7)$
  
  CAVEAT: Must verify $X$ is on the curve, otherwise you compute on $Y^2 = -X^3 + 7$ (quadratic twist). Relevant in BIP 324.

- Verification methods:
  - Calculate y-coordinate (requires square root)
  - Use Jacobi symbol: Legendre symbol $(x/p) = x^{(p-1)/2} \bmod p$, Jacobi symbol extends to composite numbers, computed with GCD algorithms

- BIP 340: X-only signatures need tie-breaker (evenness) to determine which y-coordinate to use

- History: Initially wanted Jacobi for performance, later realized benchmark was flawed

- [Binary GCD](https://en.wikipedia.org/wiki/Binary_GCD_algorithm)
    ```python3
     def bingcd(f, g):
        assert (f & 1)
        while g != 0:
            if (f & 1) == 0:
                f >>= 1
                continue
            if (g & 1) == 0:
                g >>= 1
                continue
            if g > f:
                f, g = g, f
            f -= g
        return f
    ```
- [Safe GCD](https://github.com/bitcoin-core/secp256k1/blob/master/doc/safegcd_implementation.md)

    ```python3
        def safegcd(f, g):
            assert (f & 1)
            delta = 1
            while g != 0:
                if delta > 0 and g & 1:
                    delta, f, g = 1 - delta, g, (g - f) >> 1
                elif g & 1:
                    delta, f, g = 1 + delta, f, (g + f) >> 1
                else:
                    delta, f, g = 1 + delta, f, (g    ) >> 1
        return abs(f)
    ```
- [Safe GCD correctness validation by sipa (AKA safeGCD bounds)](https://github.com/sipa/safegcd-bounds)

- Jacobi was mistakenly thought to be way faster than inverses. That was wrong, even in GMP it is quite similar. With safegcd, inverses became faster than Jacobi.

- Point multiplication: When a scalar is multiplied by a point, e.g., $a \cdot B$ where $a = 5$ means $5 \cdot B = B + B + B + B + B$
  Example: Schnorr verification checks $s \cdot G = R + H(m,R) \cdot P$, same as $R = s \cdot G - H(m,R) \cdot P$
  More general form: $a_1 \cdot P_1 + a_2 \cdot P_2 + \ldots + a_m \cdot P_m$
  The larger the scalar, the more operations required.

- Bos-Coster point multiplication algorithm: Reduce scalar factors by rewriting the equation to subtract them from each other
  $a_1P_1 + a_2P_2 + \ldots + a_mP_m$
  $= a_1P_1 - a_2P_1 + a_2P_1 + a_2P_2 + \ldots + a_mP_m$
  $= (a_1 - a_2)P_1 + a_2(P_1 + P_2) + \ldots + a_mP_m$
  After this step, calculate new scalars and points, delete if 0, sort, and repeat.
  Weakness: When $a_2$ is much smaller than $a_1$, subtracting doesn't reduce significantly. Detect this and reduce by doing some multiplication first.

- Single point multiplication using double-and-add algorithm
  For $P = G$, precompute all values at compile time.
  Number of doubling operations equals number of bits in the scalar.
  Number of addition operations equals number of 1s in the scalar (128 on average).

- Straus's method reduces additions: Group 2 bits (or more) and precompute $T_1 = [0, P_1, 2P_1, 3P_1]$, $T_2 = [0, P_2, 2P_2, 3P_2]$
  At each iteration, double the result 2 times, add $T_1[a_1[i]]$ and $T_2[a_2[i]]$.

- Non-adjacent form (NAF) for writing numbers: Each digit is $\{1, 0, -1\}$ with no adjacent nonzero numbers
  $7$ in NAF $= 1\ 0\ 0\ -1$ (since $8 - 1 = 7$)
  Converting from binary right to left: $0$ if $0$, $1$ if $1$. If the $1$ makes adjacent nonzeros, change the last bit to $-1$ and add $1\ 0$.
  Example: $1\ 1$ becomes $1\ 0\ -1$, so $3 = 4 - 1$
  3-bit window NAF requires at least 3 zeros between any nonzero bits.

- Effective affine: When adding Jacobian coordinates with the same $Z$, map to an isomorphic curve where affine/Jacobian conversion is trivial
  $(X, Y, Z_1) + (X_2, Y_2, Z_2)$ such that $Z_1 = Z_2$
  Define $f$ such that $R = Z$: $f_R(X, Y, Z) = (X, Y, Z/R)$
  Now we have $(X, Y, 1)$ and $(X_2, Y_2, 1)$ in isomorphic curve $Y^2 = X^3 + 7Z^6$
  Affine coordinates can be calculated without division.
  Use before addition to add these coordinates, get resulting coordinate $(X_3, Y_3, Z_3)$
  Convert back by applying $f_R^{-1}$ which gives $(X_3, Y_3, Z_3 \cdot Z)$
  Implementation: `secp256k1_gej_add_zinv_var`

- Pippenger's algorithm: When you have many points, create buckets and distribute additions so the number of additions doesn't scale with the number of points
  Only need $2 * numBuckets * numGroups$ additions
  Ecmult uses this when multiplying with more than 100 points
  Relevant for batch verification and MuSig key aggregation.

- Benchmarking the point multiplication techniques in `ec_mult`
  - `./bench` for high-level benchmarking
  - `./bench_internal` for low-level benchmarking

- C programming language is used in secp256k1 for better compatibility

- Secp256k1 vs other curves: ed25519 is another curve sometimes used in place of secp256k1 because it has one uniform point addition formula. Secp256k1 point addition has 3 forms that need to be made constant-time to avoid side-channel attacks. This means evaluating all branches, which is inefficient. However, ed25519 has a drawback where the curve order has a cofactor, so there can be subgroups with small order. To maintain security, points must be verified to be in the correct subgroup.

- The drawback of secp256k1 curve: Big number representation in C is tedious, and libsecp256k1 has 8 scalar representations. The magnitude of the carry is what other components of the codebase need to be aware of.

### Provable Cryptography for Bitcoin by Jonas

- [Prework Materials](https://github.com/cryptography-camp/prework)
- [Workbook](https://github.com/cryptography-camp/workbook)

### Bitcoin for Privacy: Bulletproofs by Waxwing

- [Ideas behind Bulletproofs](https://github.com/AdamISZ/ideas-about-bulletproofs/blob/master/bulletproofs.pdf)
- [ZKP MOOC Lecture Series: Introduction and History](https://www.youtube.com/watch?v=uchjTIlPzFo)
- [From Zero (Knowledge) to Bulletproofs](https://github.com/AdamISZ/from0k2bp/blob/master/from0k2bp.pdf)
- [Bulletproofs Paper](https://eprint.iacr.org/2017/1066.pdf)
- [Efficient Zero-Knowledge Arguments for Arithmetic Circuits in the Discrete Log Setting](https://eprint.iacr.org/2016/263.pdf)
- [Confidential Transactions](https://gnusha.org/pi/bitcoindev/?q=confidential+transactions)
- [Recent Work OP_ZKP](https://gnusha.org/pi/bitcoindev/?q=bulletproofs)

#### Resources

- [Guide to ECC (Hankerson, Menezes, Vanstone)](http://tomlr.free.fr/Math%E9matiques/Math%20Complete/Cryptography/Guide%20to%20Elliptic%20Curve%20Cryptography%20-%20D.%20Hankerson,%20A.%20Menezes,%20S.%20Vanstone.pdf)
- [Guide to the secp256k1 library](https://cryptocamp.website/t/secp256k1-part-3-guide-to-the-library/174)
- [Advances in Cryptology](https://link.springer.com/book/10.1007/3-540-69053-0)
- [Sipa class Plan](https://gist.github.com/sipa/eee05b5eec6fe0a58ddd2d3b92d49c4e)
- [glozow notes](https://docs.google.com/document/d/1i5SwViHNfJtZAY4u0ekDpVo9mirkXW3JIih8vdA0DVM/edit?usp=sharing)
