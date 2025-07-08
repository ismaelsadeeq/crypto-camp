#include "common/big_int.hpp"
#include "common/fast_exp.hpp"
#include "common/multiplicative_inverse.hpp"
#include "elgamal/elgamal.hpp"
#include "secp256k1/secp256k1.hpp"

#include "cassert"
#include "iostream"

void fast_exp_tests() {
    // ================================
    // Test without modulo (small range)
    assert(fast_exponent(big_int(3), big_int(218), big_int(1000)) == big_int(489));
    assert(fast_exponent(big_int(2), big_int(10)) == big_int(1024));
    assert(fast_exponent(big_int(3), big_int(3)) == big_int(27));
    assert(fast_exponent(big_int(5), big_int(0)) == big_int(1));
    assert(fast_exponent(big_int(0), big_int(5)) == big_int(0));
    assert(fast_exponent(big_int(0), big_int(0)) == big_int(1));  // defined as 1 in your implementation

    // ================================
    // With modulo
    assert(fast_exponent(big_int(3), big_int(2), big_int(4)) == big_int(1));
    assert(fast_exponent(big_int(2), big_int(10), big_int(1000)) == big_int(24)); // 1024 % 1000 = 24
    assert(fast_exponent(big_int(7), big_int(4), big_int(5)) == big_int(1));      // 2401 % 5 = 1
    assert(fast_exponent(big_int(10), big_int(0), big_int(6)) == big_int(1));     // 1 % 6 = 1
    assert(fast_exponent(big_int(0), big_int(5), big_int(7)) == big_int(0));      // Base zero, result 0

    // Negative base, modulo
    assert(fast_exponent(big_int(-3), big_int(7), big_int(7)) == big_int(4));     // (-3)^7 = -2187 % 7 = 4

    // ================================
    // Test with huge numbers
    big_int huge_base = conv<big_int>("84265675725482892459719348378630146162719620409152809167814480007059199482163");
    big_int huge_exp  = conv<big_int>("123456789123456789");
    big_int huge_mod  = conv<big_int>("115792089237316195423570985008687907853269984665640564039457584007913129639937"); // 2^256 - 1

    big_int huge_result = conv<big_int>("27255166557789855811328668216188941042231292868742423663488474786644472725517");
    assert(fast_exponent(huge_base, huge_exp, huge_mod) == huge_result);

    std::cout << "All fast exponent test vectors passed!" << std::endl;
}

void test_multiplicative_inverse() {
    // ============================================
    // Well-known small primes
    assert(get_multiplicative_inverse(big_int(2), big_int(7)) == big_int(4)); // 2 * 4 ≡ 1 mod 7
    assert(get_multiplicative_inverse(big_int(3), big_int(7)) == big_int(5)); // 3 * 5 ≡ 1 mod 7
    assert(get_multiplicative_inverse(big_int(4), big_int(7)) == big_int(2));
    assert(get_multiplicative_inverse(big_int(5), big_int(7)) == big_int(3));
    assert(get_multiplicative_inverse(big_int(6), big_int(7)) == big_int(6)); // self-inverse

    // ============================================
    // Edge cases
    assert(get_multiplicative_inverse(big_int(1), big_int(7)) == big_int(1));
    assert(get_multiplicative_inverse(big_int(6), big_int(7)) == big_int(6));

    // ============================================
    // Medium primes
    assert(get_multiplicative_inverse(big_int(10), big_int(17)) == big_int(12));
    assert(get_multiplicative_inverse(big_int(4), big_int(17)) == big_int(13));

    // ============================================
    // Large primes
    assert(get_multiplicative_inverse(big_int(123), big_int(1009)) == big_int(484));
    assert(get_multiplicative_inverse(big_int(456), big_int(10007)) == big_int(1602));

    // ============================================
    // Negative base
    assert(get_multiplicative_inverse(big_int(-3), big_int(7)) == big_int(2)); // -3 ≡ 4 mod 7; inverse of 4 is 2

    std::cout << "All multiplicative inverse test vectors passed!" << std::endl;
}

void test_elgamal() {
    struct Vector {
        big_int p, g, x, k, m, a, b;
    };

    std::vector<Vector> test_vectors = {
        {big_int(71), big_int(33), big_int(62), big_int(31), big_int(15), big_int(62), big_int(18)},
        {big_int(23), big_int(11), big_int(6), big_int(3), big_int(10), big_int(20), big_int(22)},
        {big_int(809), big_int(3), big_int(68), big_int(89), big_int(100), big_int(345), big_int(517)},
        {big_int(17), big_int(6), big_int(5), big_int(10), big_int(13), big_int(15), big_int(9)},
        {
            conv<big_int>("84265675725482892459719348378630146162719620409152809167814480007059199482163"),
            big_int(5),
            conv<big_int>("2799014790424892046701478888900891009403869701173893426"),
            conv<big_int>("23517683968368899022119256606644551548285683288848885921"),
            conv<big_int>("87521618088882658227876453"),
            conv<big_int>("22954586883013884818653063688294540134886732496160582262267014428782771199687"),
            conv<big_int>("56046128113101346099694619669629128017849277484825379502821514323706183544424")
        },
        {
            conv<big_int>("12658517083168187407924345155971956101250996576825115113297001855799796437288935576230034157578333666497170430505565580165565829633685607504706642034926119"),
            conv<big_int>(7),
            conv<big_int>("2001688878140630728014209681954697141876038523595247208"),
            conv<big_int>("5446024688717452254835115775456957961297236108858862823"),
            conv<big_int>("87521618088882658227876453"),
            conv<big_int>("2150519483988769855483983776372336742288374425191291528256965705108393490638750082340115568718132372731853110762124400441550538499580316268601341087676203"),
            conv<big_int>("1540471266850557563382406324432354117072109094950140952195099581066490559252112349492583688225692526496193879919152401794896907007394565292272866724291488")
        }
    };

    for (size_t i = 0; i < test_vectors.size(); ++i) {
        const auto& v = test_vectors[i];

        big_int public_key = fast_exponent(v.g, v.x, v.p);
        auto [c1, c2] = basic_elgamal_encrypt(public_key, v.m, v.k, v.g, v.p);

        // Check against known ciphertexts
        assert(c1 == v.a);
        assert(c2 == v.b);

        // Now decrypt and check the original message
        big_int decrypted = basic_elgamal_decrypt(c1, c2, v.x, v.p);
        assert(decrypted == v.m);

    }

    std::cout << "All ElGamal test vectors passed!" << std::endl;
}

void test_affine_addition() {
    // Test 1
    point P1 = {conv<big_int>("67021774492365321256634043516869791044054964063002935266026048760627130221114"),
                conv<big_int>("22817883221438079958217963063610327523693969913024717835557258242342029550595")};
    point Q1 = {conv<big_int>("61124938217888369397608518626468079588341162087856379517664485009963441753645"),
                conv<big_int>("5723382937169086635766392599511664586625983027860520036338464885987365575658")};
    point R1 = affine_point_addition(P1, Q1);
    assert(R1.first == conv<big_int>("78518484088348927894279633941273782106215956054783044881924083038087974375069"));
    assert(R1.second == conv<big_int>("18400956471605157290158330638123206056219981947313880254846397293938760781200"));

    // Test 2 (reflection — result is 0 point)
    point P2 = {conv<big_int>("44797955726860071483167773525787460171685721903803276437396496681708013097206"),
                conv<big_int>("112878323467240798018200025047246733779416351939079609883282945822975931592141")};
    point Q2 = {conv<big_int>("44797955726860071483167773525787460171685721903803276437396496681708013097206"),
                conv<big_int>("2913765770075397405370959961441174073853632726560954156174638184932903079522")};
    point R2 = affine_point_addition(P2, Q2);
    assert(R2.first == 0);
    assert(R2.second == 0);

    // Test 3 (point doubling)
    point P3 = {conv<big_int>("95200151225387174391707134980196577229773167465894787919263504089948495725202"),
                conv<big_int>("94213123740092242124032541289267941722641721980066680728855126898974205181980")};
    point R3 = affine_point_addition(P3, P3);
    assert(R3.first == conv<big_int>("5909177817561749019375996132097716007690336893057112295739767175467136927121"));
    assert(R3.second == conv<big_int>("32162989297956602751967132742255814558956860587998309119003795115938320862381"));

    // Test 4
    point P4 = {conv<big_int>("24050370140998638157368766089090079788245793492514664296883668741529047882113"),
                conv<big_int>("14478882322437672032054487172211630444001167135141445302555096737662467817571")};
    point Q4 = {conv<big_int>("15045863282447234231848775263852322721143017336655001075698483887751182719636"),
                conv<big_int>("14478882322437672032054487172211630444001167135141445302555096737662467817571")};
    point R4 = affine_point_addition(P4, Q4);
    assert(R4.first == conv<big_int>("76695855813870323034353443655745505343881173836470898666875431378628604069914"));
    assert(R4.second == conv<big_int>("101313206914878523391516497836476277409268817530499118736902487270246366854092"));

    // Test 5
    point P5 = {conv<big_int>("14256779447437936128616290794341059890063336098474125854681710102809814868320"),
                conv<big_int>("90566103014364716248988921534849031279541603477816641946022463390335657035131")};
    point Q5 = {conv<big_int>("2303067510121489830312323422056091166740725427601969990117485452141659178613"),
                conv<big_int>("25225986222951479174582063473838876573728381187823922093435120617573177636532")};
    point R5 = affine_point_addition(P5, Q5);
    assert(R5.first == conv<big_int>("95430772898311369787541983276504378677140303663720683940530878996024106515165"));
    assert(R5.second == conv<big_int>("48068184564993462938397020947826677061288691733511084479824032705110581338856"));

    std::cout << "All affine addition test vectors passed!\n";
}

int main() {
    fast_exp_tests();
    test_multiplicative_inverse();
    test_elgamal();
    test_affine_addition();
    return 0;
}