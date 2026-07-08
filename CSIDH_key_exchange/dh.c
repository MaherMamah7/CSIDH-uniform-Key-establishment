/* Commutative group-action Diffie-Hellman (CSIDH / qt-Pegasis style).
 *
 *   E0 = starting curve (j = 1728).
 *   Alice: secret ideal [a] (uniform Cayley sampler) -> public E_A = [a] E0.
 *   Bob:   secret ideal [b] (uniform Cayley sampler) -> public E_B = [b] E0.
 *   They exchange ONLY the Montgomery coefficient A of their public curve.
 *   Alice computes [a] E_B, Bob computes [b] E_A.
 *   By commutativity of the class-group action, [a][b]E0 = [b][a]E0 = [ab]E0,
 *   so both land on the SAME curve.
 *
 * qt-Pegasis' 2^e/4D-split action determines that curve only up to its
 * Montgomery model and quadratic twist (the internal twist selection is keyed on
 * the ideal, not the target curve), so the raw coefficients agree only about half
 * the time. Canonicalisation is therefore applied here, in the protocol layer, to
 * the FINAL shared curve only (never to the intermediate public keys, whose
 * orientation the second action still needs): canonical_A returns the
 * deterministic minimum coefficient over all F_p Montgomery models of the shared
 * curve and its twist, so both parties derive ONE identical shared secret in a
 * single field element.
 *
 * The order is the qt-Pegasis Frobenius order O = Z[(1+sqrt(-p))/2], D = -p,
 * shared by the sampler and the action. The secret ideal's norm-equation
 * solution (qt_output_t) is curve-independent, so it is applied both to E0 and
 * to the peer's curve.
 */
#include <ec.h>
#include <qt_pegasis.h>
#include <qt_normeq.h>
#include "qform.h"
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

static double now_ms(void){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return 1e3*t.tv_sec+1e-6*t.tv_nsec; }

static int sample_split_prime(mpz_t ell_out, const mpz_t L, const mpz_t D, gmp_randstate_t rng)
{
    mpz_t span, ell; mpz_inits(span, ell, NULL); int ok = 0;
    mpz_sub_ui(span, L, 2); if (mpz_sgn(span) <= 0) mpz_set_ui(span, 1);
    for (int i = 0; i < 1000000; i++) {
        mpz_urandomm(ell, rng, span); mpz_add_ui(ell, ell, 3); mpz_nextprime(ell, ell);
        if (mpz_cmp(ell, L) > 0) continue;
        if (mpz_kronecker(D, ell) != 1) continue;
        mpz_set(ell_out, ell); ok = 1; break;
    }
    mpz_clears(span, ell, NULL); return ok;
}

static void sample_uniform_class(qform_t *x, const mpz_t D, double A, gmp_randstate_t rng)
{
    double lnp = mpz_sizeinbase(D, 2) * log(2.0), lnlnp = log(lnp);
    mpz_t L; mpz_init(L); mpz_set_d(L, pow(lnp, A)); if (mpz_cmp_ui(L, 5) < 0) mpz_set_ui(L, 5);
    unsigned long t = (unsigned long)llround(0.75 * lnp / ((A/2.0 - 1.0) * lnlnp)); if (t < 1) t = 1;
    qform_t pf, step; qform_init(&pf); qform_init(&step); mpz_t ell; mpz_init(ell);
    qform_principal(x, D);
    for (unsigned long i = 0; i < t; i++) {
        sample_split_prime(ell, L, D, rng);
        unsigned long bit; { mpz_t rb; mpz_init(rb); mpz_urandomb(rb, rng, 1); bit = mpz_get_ui(rb); mpz_clear(rb); }
        qform_prime(&pf, ell, D);
        if (bit) qform_set(&step, &pf); else qform_inverse(&step, &pf);
        qform_compose(x, x, &step, D);
    }
    qform_reduce(x);
    mpz_clear(ell); qform_clear(&pf); qform_clear(&step); mpz_clear(L);
}

/* secret key = norm-equation solution of a uniform-sampled ideal [a] */
typedef struct { qt_output_t sol; qform_t cls; } secret_key;

/* sample a uniform ideal whose (N, alpha) is qt_normeq-solvable; fill sk */
static int keygen(secret_key *sk, const mpz_t D, double A, gmp_randstate_t rng)
{
    ibz_t N, aa, ab; ibz_init(&N); ibz_init(&aa); ibz_init(&ab);
    qform_init(&sk->cls);
    int ok = 0;
    for (int tries = 0; tries < 200 && !ok; tries++) {
        sample_uniform_class(&sk->cls, D, A, rng);
        if (mpz_fdiv_ui(sk->cls.a, 4) == 2) continue;
        mpz_t g; mpz_init(g); mpz_gcd(g, sk->cls.b, sk->cls.a);
        int coprime = (mpz_cmp_ui(g, 1) == 0); mpz_clear(g);
        if (!coprime) continue;
        mpz_set(N, sk->cls.a);
        mpz_sub_ui(aa, sk->cls.b, 1); mpz_fdiv_q_2exp(aa, aa, 1); mpz_set_ui(ab, 1);
        ok = qt_normeq(&sk->sol, &N, &aa, &ab);
    }
    ibz_finalize(&N); ibz_finalize(&aa); ibz_finalize(&ab);
    return ok;
}

/* act by the secret ideal on the curve given by Montgomery coefficient inA;
 * write the resulting curve's Montgomery coefficient to outA. */
static void act_on_coeff(fp_t *outA, const qt_output_t *sol, const fp_t *inA)
{
    ec_curve_t E, Ea, Etmp;
    ec_curve_init_from_A(&E, inA);
    compute_action(&Ea, &Etmp, &E, sol);
    fp_copy(outA, &Ea.A);
}

/* model transmitting one Montgomery coefficient over the wire (encode/decode) */
static void wire_send_recv(fp_t *dst, const fp_t *src)
{
    uint8_t buf[FP_ENCODED_BYTES];
    fp_encode(&buf, src);
    fp_decode(dst, &buf);
}

/* --- Canonical, twist-pair Montgomery coefficient -------------------------
 *
 * The action returns the shared curve in an arbitrary Montgomery model and up to
 * its quadratic twist. E_A: y^2 = x^3 + A x^2 + x has up to three rational
 * 2-torsion points (x = 0 and the roots of x^2 + Ax + 1); moving any one to the
 * origin gives another Montgomery coefficient of the SAME curve, and the twist is
 * E_{-A}. canonical_A returns the deterministic minimum over every F_p Montgomery
 * model of the curve and its twist, so both parties -- who land on the same curve
 * up to model/twist -- agree on one coefficient. */

/* MSB-first compare of canonical byte encodings; returns 1 iff a < b.
 * Buffers are zero-initialised: modexp writes only the field byte length
 * (<= FP_ENCODED_BYTES), so the unused high bytes must be cleared. */
static int fp_lt(const fp_t *a, const fp_t *b)
{
    uint8_t ba[FP_ENCODED_BYTES] = {0}, bb[FP_ENCODED_BYTES] = {0};
    fp_encode(&ba, a); fp_encode(&bb, b);
    for (int i = FP_ENCODED_BYTES - 1; i >= 0; i--)
        if (ba[i] != bb[i]) return ba[i] < bb[i];
    return 0;
}

/* F_p Montgomery models A' of E_A: for each rational 2-torsion x-coord t,
 * A' = (3t + A)/u^2 with u^2 the square root of den = 3t^2 + 2At + 1 that is
 * itself a square (den must be a square, else that translation is on the twist). */
static int same_orient_models(const fp_t *Ain, fp_t out[3])
{
    fp_t A; fp_copy(&A, Ain);
    fp_t roots[3]; int nr = 0;
    fp_set_zero(&roots[nr++]);
    fp_t A2, four, disc;
    fp_sqr(&A2, &A); fp_set_small(&four, 4); fp_sub(&disc, &A2, &four);
    if (fp_is_square(&disc)) {
        fp_t sq, negA, r, s, tmp;
        fp_copy(&sq, &disc); fp_sqrt(&sq);
        fp_neg(&negA, &A);
        fp_add(&tmp, &negA, &sq); fp_half(&r, &tmp);
        fp_sub(&tmp, &negA, &sq); fp_half(&s, &tmp);
        fp_copy(&roots[nr++], &r); fp_copy(&roots[nr++], &s);
    }
    fp_t one, three, two; fp_set_one(&one); fp_set_small(&three, 3); fp_set_small(&two, 2);
    int nc = 0;
    for (int i = 0; i < nr; i++) {
        fp_t t, t2, den, num, tmp, u;
        fp_copy(&t, &roots[i]); fp_sqr(&t2, &t);
        fp_mul(&den, &three, &t2);
        fp_mul(&tmp, &two, &A); fp_mul(&tmp, &tmp, &t); fp_add(&den, &den, &tmp);
        fp_add(&den, &den, &one);
        if (!fp_is_square(&den)) continue;
        fp_mul(&num, &three, &t); fp_add(&num, &num, &A);
        fp_copy(&u, &den); fp_sqrt(&u);
        if (!fp_is_square(&u)) fp_neg(&u, &u);   /* pick the square root that is a square */
        fp_inv(&u);
        fp_mul(&out[nc], &num, &u);
        nc++;
    }
    return nc;
}

/* Deterministic minimum coefficient over the whole twist pair (curve + twist). */
static void canonical_A(fp_t *out, const fp_t *Ain)
{
    fp_t all[6]; int n = 0;
    fp_t negA; fp_neg(&negA, Ain);
    n += same_orient_models(Ain,  all + n);
    n += same_orient_models(&negA, all + n);
    int mi = 0;
    for (int i = 1; i < n; i++) if (fp_lt(&all[i], &all[mi])) mi = i;
    fp_copy(out, &all[mi]);
}

int main(int argc, char **argv)
{
    unsigned long seed = (argc > 1) ? strtoul(argv[1], NULL, 10) : 1UL;
    double A = (argc > 2) ? atof(argv[2]) : 0.0;

    gmp_randstate_t rng; gmp_randinit_default(rng); gmp_randseed_ui(rng, seed);
    mpz_t p, D; mpz_init(p); mpz_init(D);
    mpz_set(p, CHARACTERISTIC); mpz_neg(D, p);
    if (A == 0.0) A = (mpz_sizeinbase(p, 2) > 2000) ? 6.0 : 5.0;

    printf("Group-action DH over qt-Pegasis Frobenius order (disc -p), p has %zu bits, A=%.3g, seed=%lu\n\n",
           mpz_sizeinbase(p, 2), A, seed);

    /* public starting curve E0 (j=1728) */
    ec_curve_t E0; ec_curve_init(&E0);
    fp_t A0c; fp_copy(&A0c, &E0.A);

    /* --- key generation --- */
    double tk = now_ms();
    secret_key ska, skb;
    if (!keygen(&ska, D, A, rng)) { printf("Alice keygen failed\n"); return 1; }
    if (!keygen(&skb, D, A, rng)) { printf("Bob keygen failed\n"); return 1; }
    double t_keygen = now_ms() - tk;
    gmp_printf("Alice secret [a]: reduced form (a=%Zd, ...)\n", ska.cls.a);
    gmp_printf("Bob   secret [b]: reduced form (a=%Zd, ...)\n\n", skb.cls.a);

    /* --- public keys: E_A = [a]E0, E_B = [b]E0 --- */
    fp_t A_pub_alice, A_pub_bob;
    double tp = now_ms();
    act_on_coeff(&A_pub_alice, &ska.sol, &A0c);
    act_on_coeff(&A_pub_bob,   &skb.sol, &A0c);
    double t_pub = now_ms() - tp;

    /* --- exchange: send ONLY the Montgomery A coefficient --- */
    fp_t A_from_alice, A_from_bob;
    wire_send_recv(&A_from_alice, &A_pub_alice);  /* Alice -> Bob */
    wire_send_recv(&A_from_bob,   &A_pub_bob);    /* Bob   -> Alice */

    /* --- shared secret --- */
    fp_t A_shared_alice, A_shared_bob;
    double ts = now_ms();
    act_on_coeff(&A_shared_alice, &ska.sol, &A_from_bob);   /* Alice: [a]E_B */
    act_on_coeff(&A_shared_bob,   &skb.sol, &A_from_alice); /* Bob:   [b]E_A */
    double t_shared = now_ms() - ts;

    /* Canonicalise each shared curve to the deterministic representative of its
     * twist pair, so both parties agree on ONE coefficient. */
    fp_t key_alice, key_bob;
    canonical_A(&key_alice, &A_shared_alice);
    canonical_A(&key_bob,   &A_shared_bob);

    int raw_match   = fp_is_equal(&A_shared_alice, &A_shared_bob);
    int canon_match = fp_is_equal(&key_alice, &key_bob);

    fp_print("Alice raw   A([a]E_B)     = ", &A_shared_alice);
    fp_print("Bob   raw   A([b]E_A)     = ", &A_shared_bob);
    fp_print("shared secret (canonical) = ", &key_alice);
    printf("\nraw coefficient matches            : %s\n", raw_match ? "YES" : "no (action returned a twin/other model)");
    printf("canonical shared secret matches    : %s\n", canon_match ? "YES" : "NO");
    printf("KEY EXCHANGE: %s\n", canon_match ? "SUCCESS" : "FAIL");

    printf("\ntiming:\n");
    printf("  keygen (2 parties, sample+normeq): %.3f ms\n", t_keygen);
    printf("  public keys (2 actions)          : %.3f ms\n", t_pub);
    printf("  shared secret (2 actions)        : %.3f ms\n", t_shared);
    printf("  per group action                 : %.3f ms\n", t_shared / 2.0);

    qform_clear(&ska.cls); qform_clear(&skb.cls);
    mpz_clears(p, D, NULL); gmp_randclear(rng);
    return canon_match ? 0 : 1;
}
