/* Fine-grained timing of the full-C action from A0:
 *   Step 1 = qt_normeq,  Step 2 = kernel points,  Step 3 = isogeny chain.
 * Build with -DQTP_TIMING and compile qt_pegasis.c directly (not libaction). */
#include <ec.h>
#include <qt_pegasis.h>
#include <qt_normeq.h>
#include <gmp.h>
#include <stdio.h>
#include <time.h>

extern double qtp_time_step2;   /* accumulated inside compute_action */
extern double qtp_time_step3;

static double now_s(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

int main(int argc, char **argv)
{
    unsigned long seed = (argc > 1) ? strtoul(argv[1], NULL, 10) : 1UL;
    int ncases = (argc > 2) ? atoi(argv[2]) : 50;

    gmp_randstate_t rng;
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, seed);

    qtp_time_step2 = 0.0;
    qtp_time_step3 = 0.0;
    double t_step1 = 0.0;
    int done = 0;

    for (int i = 0; i < ncases; i++) {
        ibz_t N, a, b;
        ibz_init(&N); ibz_init(&a); ibz_init(&b);
        int ok = 0;
        while (!ok) {
            if (!qt_sample_ideal(&N, &a, &b, rng)) continue;
            if (mpz_fdiv_ui(N, 4) == 2) continue;
            ok = 1;
        }
        qt_output_t sol;
        double t0 = now_s();
        int solved = qt_normeq(&sol, &N, &a, &b);
        double t1 = now_s();
        if (!solved) { ibz_finalize(&N); ibz_finalize(&a); ibz_finalize(&b); continue; }
        t_step1 += (t1 - t0);

        ec_curve_t E, Ea, Eabar;
        ec_curve_init(&E);
        compute_action(&Ea, &Eabar, &E, &sol);   /* accumulates step2/step3 */
        done++;

        ibz_finalize(&N); ibz_finalize(&a); ibz_finalize(&b);
    }
    gmp_randclear(rng);

    double s1 = 1e3 * t_step1 / done;
    double s2 = 1e3 * qtp_time_step2 / done;
    double s3 = 1e3 * qtp_time_step3 / done;
    printf("Full-C timing over %d actions (from A0, 500-bit):\n", done);
    printf("  Step 1  norm equation   : %.3f ms\n", s1);
    printf("  Step 2  kernel points   : %.3f ms\n", s2);
    printf("  Step 3  isogeny chain   : %.3f ms\n", s3);
    printf("  Steps 2+3 (action)      : %.3f ms\n", s2 + s3);
    printf("  total (1+2+3)           : %.3f ms\n", s1 + s2 + s3);
    return 0;
}
