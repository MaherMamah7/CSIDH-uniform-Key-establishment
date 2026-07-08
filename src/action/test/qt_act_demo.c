/* Minimal demo: sample a random ideal [a], act on the starting curve E0
 * (j-invariant 1728), and print j(E0) and j([a]E0).
 *
 * Build (after `make PNAME=500 run` has produced the libs in build/):
 *   gcc -I/opt/homebrew/include -L/opt/homebrew/lib -DP_500 -O3 -DNDEB UG \
 *       -Isrc/action/include -Isrc/qt_interface/include -Iinclude \
 *       -Isrc/dim4/include -Isrc/ec/include -Isrc/gf/include -Isrc/mp/include \
 *       -Isrc/params -Isrc/common/include -Isrc/quadratic/include \
 *       src/action/test/qt_act_demo.c -o build/qt_act_demo_500 \
 *       -L./build -laction -lqt -ldim4 -lec -lgf_test_utils_500 -lgf_500 \
 *       -lparams_500 -lquadratic -lcommon -lmp -lgmp -lm
 *   ./build/qt_act_demo_500 [seed]
 */
#include <ec.h>
#include <qt_pegasis.h>
#include <qt_normeq.h>
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
    unsigned long seed = (argc > 1) ? strtoul(argv[1], NULL, 10) : 20260705UL;

    gmp_randstate_t rng;
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, seed);

    /* Step 1: sample a random ideal a = (N, alpha) and solve its norm eq. */
    ibz_t N, a, b;
    ibz_init(&N); ibz_init(&a); ibz_init(&b);

    qt_output_t sol;
    int ok = 0;
    while (!ok) {
        if (!qt_sample_ideal(&N, &a, &b, rng)) continue;
        if (mpz_fdiv_ui(N, 4) == 2) continue;         /* skip N = 2 (mod 4) */
        ok = qt_normeq(&sol, &N, &a, &b);
    }

    gmp_printf("sampled ideal a: N = %Zd,  alpha = %Zd + %Zd*omega\n", N, a, b);

    /* Starting curve E0 with j-invariant 1728. */
    ec_curve_t E0, Ea, Eabar;
    ec_curve_init(&E0);

    /* Steps 2-3: act by [a] on E0. */
    compute_action(&Ea, &Eabar, &E0, &sol);

    fp_t j0, ja, jabar;
    ec_j_inv(&j0, &E0);
    ec_j_inv(&ja, &Ea);
    ec_j_inv(&jabar, &Eabar);

    fp_print("j(E0)      ", &j0);
    fp_print("j([a]E0)   ", &ja);
    fp_print("j([abar]E0)", &jabar);

    ibz_finalize(&N); ibz_finalize(&a); ibz_finalize(&b);
    gmp_randclear(rng);
    return 0;
}
