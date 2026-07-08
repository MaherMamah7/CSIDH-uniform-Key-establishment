/* End-to-end self-consistency test of the full C pipeline:
 *   Step 1 (our qt_normeq)  ->  Steps 2-3 (compute_action).
 *
 * Group-action law: for an ideal a with conjugate abar, a*abar = (N) is
 * principal, so [abar]([a]E) is isomorphic to E. We sample a, solve the norm
 * equation for both a and abar, act in sequence, and check j([abar][a]E)==j(E).
 * This needs no SageMath and no start-curve matching.
 */
#include <ec.h>
#include <qt_pegasis.h>
#include <qt_normeq.h>
#include <gmp.h>
#include <stdio.h>

int main(int argc, char **argv)
{
    unsigned long seed = (argc > 1) ? strtoul(argv[1], NULL, 10) : 20260705UL;
    int ncases = (argc > 2) ? atoi(argv[2]) : 8;

    gmp_randstate_t rng;
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, seed);

    int passed = 0, total = 0;

    for (int i = 0; i < ncases; i++) {
        ibz_t N, a, b, ca, cb;
        ibz_init(&N); ibz_init(&a); ibz_init(&b); ibz_init(&ca); ibz_init(&cb);

        int sampled = 0;
        while (!sampled) {
            if (!qt_sample_ideal(&N, &a, &b, rng)) continue;
            if (mpz_fdiv_ui(N, 4) == 2) continue;
            sampled = 1;
        }

        /* conjugate ideal abar: conj(a + b*w) = (a+b) - b*w */
        mpz_add(ca, a, b);
        mpz_neg(cb, b);

        qt_output_t sol_a, sol_abar;
        int ok1 = qt_normeq(&sol_a, &N, &a, &b);
        int ok2 = qt_normeq(&sol_abar, &N, &ca, &cb);

        total++;
        if (!ok1 || !ok2) {
            printf("case %d: SOLVE_FAIL (a=%d abar=%d)\n", i, ok1, ok2);
            goto next;
        }

        ec_curve_t E0, Ea, Etmp, Eback, Eback2;
        ec_curve_init(&E0);

        /* Ea = [a] E0 */
        compute_action(&Ea, &Etmp, &E0, &sol_a);
        /* Eback = [abar] Ea */
        compute_action(&Eback, &Eback2, &Ea, &sol_abar);

        fp_t j0, jb, jb2, ja, jat;
        ec_j_inv(&j0, &E0);
        ec_j_inv(&ja, &Ea);
        ec_j_inv(&jat, &Etmp);
        ec_j_inv(&jb, &Eback);
        ec_j_inv(&jb2, &Eback2);

        /* round trip should land back on j(E0) (either output branch) */
        int back = fp_is_equal(&jb, &j0) || fp_is_equal(&jb2, &j0);
        /* the two twists produced by a single action share j */
        int twin = fp_is_equal(&ja, &jat);

        if (back && twin) {
            passed++;
            printf("case %d: OK  (roundtrip=yes, twin-j=yes)\n", i);
        } else {
            printf("case %d: FAIL (roundtrip=%d, twin-j=%d)\n", i, back, twin);
        }

    next:
        ibz_finalize(&N); ibz_finalize(&a); ibz_finalize(&b);
        ibz_finalize(&ca); ibz_finalize(&cb);
    }

    gmp_randclear(rng);
    printf("\nqt_normeq end-to-end: %d/%d passed\n", passed, total);
    return passed == total ? 0 : 1;
}
