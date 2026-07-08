/* Near-uniform ideal-class sampler via a random walk on the Cayley graph of the
 * class group Cl(O), where O = Z[(1+sqrt(-p))/2] is the qt-Pegasis Frobenius
 * order (discriminant D = -p, p = 3 mod 4).
 *
 * Algorithm:
 *   x_0 = principal class
 *   for i = 1..t:  sample a split prime ideal p_i of norm <= L,  eps_i in {-1,+1}
 *                  x_i = x_{i-1} * p_i^{eps_i}   (kept reduced)
 *   output Gaussian-reduced canonical representative of x_t.
 * with  L = (ln p)^2  and  t = round( (3/4) * ln p / ln ln p ).
 *
 * Ideal classes are binary quadratic forms of discriminant D (see qform.*).
 */
#include "qform.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

/* qt-Pegasis primes: p = odd * 2^e - 1 */
static int qtp_prime(mpz_t p, int pname)
{
    unsigned long odd, e;
    switch (pname) {
        case 500:  odd = 27; e = 500;  break;
        case 1000: odd = 15; e = 1004; break;
        case 1500: odd = 5;  e = 1522; break;
        case 2000: odd = 5;  e = 2014; break;
        case 4000: odd = 45; e = 4024; break;
        default: return 0;
    }
    mpz_set_ui(p, odd);
    mpz_mul_2exp(p, p, e);
    mpz_sub_ui(p, p, 1);
    return 1;
}

/* sample a split odd prime ell <= L (kronecker(D, ell) == 1) into ell_out.
 * returns 1 on success, 0 on failure. */
static int sample_split_prime(mpz_t ell_out, const mpz_t L, const mpz_t D, gmp_randstate_t rng)
{
    mpz_t span, ell; mpz_inits(span, ell, NULL);
    int ok = 0;
    mpz_sub_ui(span, L, 2);                 /* [0, L-2) -> +3 gives [3, L+1) */
    if (mpz_sgn(span) <= 0) mpz_set_ui(span, 1);
    for (int tries = 0; tries < 1000000; tries++) {
        mpz_urandomm(ell, rng, span);
        mpz_add_ui(ell, ell, 3);
        mpz_nextprime(ell, ell);
        if (mpz_cmp(ell, L) > 0) continue;
        if (mpz_kronecker(D, ell) != 1) continue;   /* must split */
        mpz_set(ell_out, ell);
        ok = 1;
        break;
    }
    mpz_clears(span, ell, NULL);
    return ok;
}

int main(int argc, char **argv)
{
    int pname = 500;
    unsigned long seed = 1;
    int verbose = 0;
    double A = 5.0;                 /* L = (ln p)^A ; requires A > 2.
                                     * A~5-6 is fastest: since x is reduced each
                                     * step, per-step cost is ~constant, so larger
                                     * A (fewer steps) wins until L=(ln p)^A exceeds
                                     * ~2^64 and prime sampling gets costly. */
    long reps = 1;                  /* timing: repeat the walk this many times */
    mpz_t p; mpz_init(p);
    int have_custom = 0;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-p") && i+1 < argc) { mpz_set_str(p, argv[++i], 10); have_custom = 1; }
        else if (!strcmp(argv[i], "-s") && i+1 < argc) seed = strtoul(argv[++i], NULL, 10);
        else if (!strcmp(argv[i], "-A") && i+1 < argc) A = atof(argv[++i]);
        else if (!strcmp(argv[i], "-r") && i+1 < argc) reps = atol(argv[++i]);
        else if (!strcmp(argv[i], "-v")) verbose = 1;
        else pname = atoi(argv[i]);
    }
    if (A <= 2.0) {
        fprintf(stderr, "A must be > 2 (t = (3/4) log p / ((A/2-1) log log p) diverges at A=2)\n");
        return 1;
    }
    if (!have_custom && !qtp_prime(p, pname)) {
        fprintf(stderr, "unknown PNAME %d (use 500/1000/1500/2000/4000 or -p <decimal>)\n", pname);
        return 1;
    }
    if (mpz_fdiv_ui(p, 4) != 3) { fprintf(stderr, "p must be 3 mod 4\n"); return 1; }

    mpz_t D; mpz_init(D); mpz_neg(D, p);            /* D = -p */

    /* parameters: L = (ln p)^A,  t = round( (3/4) ln p / ((A/2 - 1) ln ln p) ) */
    double lnp = mpz_sizeinbase(p, 2) * log(2.0);   /* ~ ln p */
    double lnlnp = log(lnp);
    double Ld = pow(lnp, A);
    mpz_t L; mpz_init(L);
    mpz_set_d(L, Ld);                               /* L = round((ln p)^A) */
    if (mpz_cmp_ui(L, 5) < 0) mpz_set_ui(L, 5);
    unsigned long t = (unsigned long)llround(0.75 * lnp / ((A/2.0 - 1.0) * lnlnp));
    if (t < 1) t = 1;

    gmp_randstate_t rng; gmp_randinit_default(rng); gmp_randseed_ui(rng, seed);

    printf("p has %zu bits (p = %s)\n", mpz_sizeinbase(p, 2), have_custom ? "custom" : "qt-Pegasis prime");
    gmp_printf("A = %.4g,  L = (ln p)^A = %Zd,  t = round(3/4 * ln p / ((A/2-1) ln ln p)) = %lu,  seed = %lu\n",
               A, L, t, seed);

    /* random walk (repeated `reps` times for timing) */
    qform_t x, pf, step; qform_init(&x); qform_init(&pf); qform_init(&step);
    mpz_t ell; mpz_init(ell);
    struct timespec ts0, ts1;
    clock_gettime(CLOCK_MONOTONIC, &ts0);
    for (long rep = 0; rep < reps; rep++) {
        qform_principal(&x, D);
        for (unsigned long i = 0; i < t; i++) {
            if (!sample_split_prime(ell, L, D, rng)) { fprintf(stderr, "failed to sample split prime\n"); return 1; }
            unsigned long bit; { mpz_t rb; mpz_init(rb); mpz_urandomb(rb, rng, 1); bit = mpz_get_ui(rb); mpz_clear(rb); }
            int eps = bit ? 1 : -1;         /* epsilon_i in {-1,+1} */

            qform_prime(&pf, ell, D);
            if (eps < 0) qform_inverse(&step, &pf); else qform_set(&step, &pf);
            qform_compose(&x, &x, &step, D);
            if (verbose && reps == 1) gmp_printf("  step %2lu: ell=%Zd eps=%+d -> a=%Zd\n", i+1, ell, eps, x.a);
        }
        qform_reduce(&x);
    }
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    double walk_ms = 1e3*((ts1.tv_sec-ts0.tv_sec) + (ts1.tv_nsec-ts0.tv_nsec)*1e-9) / (double)reps;
    mpz_clear(ell);
    printf("avg walk time: %.4f ms/class over %ld reps  (%.5f ms/step)\n", walk_ms, reps, walk_ms/(double)t);

    /* sanity: discriminant preserved */
    mpz_t disc; mpz_init(disc); qform_disc(disc, &x);
    int disc_ok = (mpz_cmp(disc, D) == 0);

    printf("\nnear-uniform ideal class (reduced binary quadratic form of disc -p):\n");
    qform_print("  ", &x);
    gmp_printf("  norm N(class) = a = %Zd  (%zu bits)\n", x.a, mpz_sizeinbase(x.a, 2));
    printf("  discriminant check: %s\n", disc_ok ? "OK" : "MISMATCH");

    mpz_clears(disc, D, p, NULL);
    qform_clear(&x); qform_clear(&pf); qform_clear(&step);
    gmp_randclear(rng);
    return disc_ok ? 0 : 1;
}
