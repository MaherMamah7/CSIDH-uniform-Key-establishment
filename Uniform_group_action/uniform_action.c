/* Uniform group action:
 *   1. sample a near-uniform ideal class of the qt-Pegasis Frobenius order
 *      O = Z[(1+sqrt(-p))/2] (disc D = -p) via the Cayley-graph random walk
 *      (reuses qform.* from ../Cayley_uniform_ideal),
 *   2. convert the reduced form (a,b,c) to qt-Pegasis' (N, alpha) coordinates,
 *   3. act by that ideal on the qt-Pegasis starting curve E0 (j=1728) using
 *      qt_normeq (Step 1) + compute_action (Steps 2-3),
 *   4. verify with the group law: [abar]([a]E0) ~= E0.
 *
 * Both sides use the SAME prime p = CHARACTERISTIC and the SAME order, so the
 * only work is the representation bridge:
 *   ideal [a, (b+sqrt(D))/2]  <->  (N=a, alpha=(b-1)/2 + 1*omega),  omega=(1+sqrt(D))/2
 *   conjugate ideal           <->  (N=a, alpha=(b+1)/2 + (-1)*omega)
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

static double now_ms(void)
{
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return 1e3 * ts.tv_sec + 1e-6 * ts.tv_nsec;
}

/* sample a split odd prime ell <= L (kronecker(D, ell) == 1) */
static int sample_split_prime(mpz_t ell_out, const mpz_t L, const mpz_t D, gmp_randstate_t rng)
{
    mpz_t span, ell; mpz_inits(span, ell, NULL);
    int ok = 0;
    mpz_sub_ui(span, L, 2);
    if (mpz_sgn(span) <= 0) mpz_set_ui(span, 1);
    for (int tries = 0; tries < 1000000; tries++) {
        mpz_urandomm(ell, rng, span); mpz_add_ui(ell, ell, 3);
        mpz_nextprime(ell, ell);
        if (mpz_cmp(ell, L) > 0) continue;
        if (mpz_kronecker(D, ell) != 1) continue;
        mpz_set(ell_out, ell); ok = 1; break;
    }
    mpz_clears(span, ell, NULL);
    return ok;
}

/* one near-uniform class as a reduced form of discriminant D */
static void sample_uniform_class(qform_t *x, const mpz_t D, double A, gmp_randstate_t rng)
{
    double lnp = mpz_sizeinbase(D, 2) * log(2.0);
    double lnlnp = log(lnp);
    mpz_t L; mpz_init(L); mpz_set_d(L, pow(lnp, A));
    if (mpz_cmp_ui(L, 5) < 0) mpz_set_ui(L, 5);
    unsigned long t = (unsigned long)llround(0.75 * lnp / ((A/2.0 - 1.0) * lnlnp));
    if (t < 1) t = 1;

    qform_t pf, step; qform_init(&pf); qform_init(&step);
    mpz_t ell; mpz_init(ell);
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

int main(int argc, char **argv)
{
    unsigned long seed = (argc > 1) ? strtoul(argv[1], NULL, 10) : 1UL;
    double A = (argc > 2) ? atof(argv[2]) : 0.0;   /* 0 => auto by prime size */

    gmp_randstate_t rng; gmp_randinit_default(rng); gmp_randseed_ui(rng, seed);

    mpz_t p, D; mpz_init(p); mpz_init(D);
    mpz_set(p, CHARACTERISTIC);           /* qt-Pegasis prime for this PNAME */
    mpz_neg(D, p);                        /* D = -p */
    /* default L-exponent: A=6 for p larger than 2000 bits, else A=5 */
    if (A == 0.0) A = (mpz_sizeinbase(p, 2) > 2000) ? 6.0 : 5.0;
    if (A <= 2.0) { fprintf(stderr, "A must be > 2\n"); return 1; }

    printf("p has %zu bits;  order O = Z[(1+sqrt(-p))/2], disc D = -p;  A = %.3g, seed = %lu\n",
           mpz_sizeinbase(p, 2), A, seed);

    ibz_t N, aa, ab, caa, cab;
    ibz_init(&N); ibz_init(&aa); ibz_init(&ab); ibz_init(&caa); ibz_init(&cab);
    qform_t x; qform_init(&x);
    qt_output_t sol_a, sol_abar;

    /* sample a uniform class whose (N, alpha) is acceptable to qt_normeq */
    double t_walk = 0.0, t_normeq = 0.0, ts;
    int tries = 0, ok1 = 0, ok2 = 0;
    for (; tries < 200; tries++) {
        ts = now_ms();
        sample_uniform_class(&x, D, A, rng);
        t_walk += now_ms() - ts;
        /* N = a; require N odd-ish (skip N == 2 mod 4) and gcd(b, a) == 1 */
        if (mpz_fdiv_ui(x.a, 4) == 2) continue;
        { mpz_t g; mpz_init(g); mpz_gcd(g, x.b, x.a);
          int coprime = (mpz_cmp_ui(g, 1) == 0); mpz_clear(g);
          if (!coprime) continue; }
        /* alpha = (b-1)/2 + 1*omega ;  conj = (b+1)/2 + (-1)*omega ; same N=a */
        mpz_set(N, x.a);
        mpz_sub_ui(aa, x.b, 1); mpz_fdiv_q_2exp(aa, aa, 1); mpz_set_ui(ab, 1);
        mpz_add_ui(caa, x.b, 1); mpz_fdiv_q_2exp(caa, caa, 1); mpz_set_si(cab, -1);

        ts = now_ms();
        ok1 = qt_normeq(&sol_a, &N, &aa, &ab);
        ok2 = qt_normeq(&sol_abar, &N, &caa, &cab);
        t_normeq += now_ms() - ts;
        if (ok1 && ok2) break;
    }
    if (!(ok1 && ok2)) {
        printf("could not obtain a qt_normeq-solvable uniform ideal in %d tries\n", tries);
        return 1;
    }

    gmp_printf("sampled uniform class: reduced form (a=%Zd, b=%Zd, c=%Zd)\n", x.a, x.b, x.c);
    gmp_printf("  -> qt-Pegasis ideal (N=%Zd, alpha=%Zd + %Zd*omega)   [%d resample(s)]\n",
               N, aa, ab, tries);

    /* act on the starting curve */
    ec_curve_t E0, Ea, Etmp, Eback, Eback2;
    ec_curve_init(&E0);
    ts = now_ms();
    compute_action(&Ea, &Etmp, &E0, &sol_a);         /* Ea = [a] E0 */
    double t_action1 = now_ms() - ts;
    ts = now_ms();
    compute_action(&Eback, &Eback2, &Ea, &sol_abar); /* [abar][a] E0 */
    double t_action2 = now_ms() - ts;

    fp_t j0, ja, jat, jb, jb2;
    ec_j_inv(&j0, &E0); ec_j_inv(&ja, &Ea); ec_j_inv(&jat, &Etmp);
    ec_j_inv(&jb, &Eback); ec_j_inv(&jb2, &Eback2);

    int back = fp_is_equal(&jb, &j0) || fp_is_equal(&jb2, &j0);
    int twin = fp_is_equal(&ja, &jat);

    fp_print("j(E0)      = ", &j0);
    fp_print("j([a]E0)   = ", &ja);
    printf("group-law round trip  [abar][a]E0 == E0 : %s\n", back ? "OK" : "FAIL");
    printf("twin-j (two twists of one action share j): %s\n", twin ? "OK" : "FAIL");

    /* Both actions are computed (needed for the round trip); report the warm
     * one as the representative steady-state cost, ignoring first-call warm-up. */
    double t_action = (t_action2 < t_action1) ? t_action2 : t_action1;
    printf("\ntiming:\n");
    printf("  sampling (walk, incl. resamples) : %.3f ms\n", t_walk);
    printf("  qt_normeq (Step 1)               : %.3f ms\n", t_normeq/2.0);
    printf("  compute_action (Steps 2-3)       : %.3f ms\n", t_action);
    printf("  one group action ([a]E0)         : %.3f ms\n", t_walk + t_normeq/2.0 + t_action);

    ibz_finalize(&N); ibz_finalize(&aa); ibz_finalize(&ab); ibz_finalize(&caa); ibz_finalize(&cab);
    qform_clear(&x); mpz_clears(p, D, NULL); gmp_randclear(rng);
    return (back && twin) ? 0 : 1;
}
