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
    assert(point_are_equal(R2, POINT_AT_INFINITY));

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

void test_jacobian_addition() {

    // Test 1
    jacobian_point J1 = make_tuple(
        conv<big_int>("61168739479711927142764658335960185139044138470269152817362835609619277248733"),
        conv<big_int>("21365265259791813296359020025112135293342760115353080382870338561918313862807"),
        conv<big_int>("37064183328797598544560694959943799168750358913858865780091974718018553562419")
    );

    jacobian_point J2 = make_tuple(
        conv<big_int>("75776791705958340557958402430698975706422201066979121642449913138944604425660"),
        conv<big_int>("66383280047496136929271400526347103822935621943780462161181840552194350141564"),
        conv<big_int>("75975606300704613123930174557625573844043347281105167940536468038500802717509")
    );

    point R1 = convert_jacobian_to_affine(jacobian_point_addition(J1, J2));
    assert(R1.first == conv<big_int>("72863032945283280953636129059545959529634042357753453132026174732744194676931"));
    assert(R1.second == conv<big_int>("111529132148508388427246132585785101600429639308058372390604751264868469767543"));

    // Test 2 — result is point at infinity
    jacobian_point J3 = make_tuple(
        conv<big_int>("89959325059742944430358451400705002920926825355225869621717936807494095714290"),
        conv<big_int>("96093053924735119484524007701924861311484651710593769022900107977673928960245"),
        conv<big_int>("66142611799578950251083409575885695298839488135797694779041885661190727675299")
    );

    jacobian_point J4 = make_tuple(
        conv<big_int>("61152793683249667605361745755257610395039301799655537107480658643593848781730"),
        conv<big_int>("108824838086741573141078213715633247883899533027170274847878148878014138167046"),
        conv<big_int>("20026567909062914103680712539641599080083135680565932483453732406779235372092")
    );

    point R2 = convert_jacobian_to_affine(jacobian_point_addition(J3, J4));
    assert(point_are_equal(R2, POINT_AT_INFINITY));

    // Test 3 — normal result
    jacobian_point J5 = make_tuple(
        conv<big_int>("1547568827951595983041825486208171785819741431893371520256763714464258127790"),
        conv<big_int>("87153109579099129796596751254693228766379983077346253255841414029284516911078"),
        conv<big_int>("105104885998309941273615701006706417602105584887217436384728254947105995740715")
    );

    jacobian_point J6 = make_tuple(
        conv<big_int>("102754269592907928248165438489539780821724369832426272173645274109108284691770"),
        conv<big_int>("38298190034438650883752719589335983487411860447931052099125319988280170002045"),
        conv<big_int>("56745928453254477537417735654158445415425453625586007664329168279192608303666")
    );

    point R3 = convert_jacobian_to_affine(jacobian_point_addition(J5, J6));
    assert(R3.first == conv<big_int>("21324256287414615615026299379536579336529998865990184416926039607504524853626"));
    assert(R3.second == conv<big_int>("96719670966356830360698314514227297774284915420887284954650836535688914930874"));

    std::cout << "All Jacobian addition test vectors passed!\n";
}



void test_scalar_muliplication()
{ 
    big_int z = conv<big_int>(1);
    // Test 1
    big_int s1 = conv<big_int>("23529072936145521956642440150769408702836782170707519110832596096096916532363");
    point P1 = {conv<big_int>("94777218176490725267733209794395406270863807953747235979017564313980479098344"), conv<big_int>("53121120406880321033414824968851949358991212541220678285657788880408683486672")};
    point R1 = affine_scalar_multiplication(s1, P1);    

    assert(R1.first == conv<big_int>("81492582484984365721511233996054540050314813088236204730182464710703690737195"));
    assert(R1.second == conv<big_int>("84165397430175583340352582740254662715932722835371860159802475562062898918484"));

    jacobian_point R1_j = jacobian_scalar_multiplication(s1, convert_affine_to_jacobian(P1));
    point R1_j_p = convert_jacobian_to_affine(R1_j);
    assert(point_are_equal(R1, R1_j_p));

    // Test 2
    big_int s2 = conv<big_int>("77770687059601253501098075906318324640585620643934538062621691587089455400301");
    point Q2 = {conv<big_int>("5187380010089560191829928600869675928625207216422014112981972591844926771008"), conv<big_int>("75026050083095897004323393777174635055491620440662638678606562665317466685019")};
    point R2 = affine_scalar_multiplication(s2, Q2);
    assert(R2.first == conv<big_int>("76999255841974189685876230118581110410155956505185745130247574937430232984638"));
    assert(R2.second == conv<big_int>("87571171775685157828750403037960903210473289232782306139148947195874900187006"));

    jacobian_point R2_j = jacobian_scalar_multiplication(s2, convert_affine_to_jacobian(Q2));
    point R2_j_p = convert_jacobian_to_affine(R2_j);
    assert(point_are_equal(R2, R2_j_p));

    // Test 3
    big_int s3 = conv<big_int>("3747619523960563074315083315669137577217731866086110333821423552891044218266");
    point Q3 = {conv<big_int>("66371586610273545144505648512343824229224003523952192165787799288317344396675"), conv<big_int>("6489011411151914877089190610663845093649879070897583530615192453262848111419")};
    point R3 = affine_scalar_multiplication(s3, Q3);
    assert(R3.first == conv<big_int>("109441138145498884726545575659592733193661671281368885246963601136369148387669"));
    assert(R3.second == conv<big_int>("83708880322787879701338478937074052809697986569225329829504559758598509123336"));

    jacobian_point R3_j = jacobian_scalar_multiplication(s3, convert_affine_to_jacobian(Q3));
    point R3_j_p = convert_jacobian_to_affine(R3_j);
    assert(point_are_equal(R3, R3_j_p));
    std::cout << "All scalar multiplication test vectors passed!\n";
}

int main() {
    fast_exp_tests();
    test_multiplicative_inverse();
    test_elgamal();
    test_affine_addition();
    test_jacobian_addition();
    test_scalar_muliplication();
    return 0;
}