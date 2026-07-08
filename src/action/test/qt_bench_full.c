/* Times the two C components of one qt-Pegasis action, from the base curve A0:
 *   - qt_normeq   (Step 1, our norm-equation solver)
 *   - compute_action (Steps 2-3, the 4D isogeny chain)
 * Reports per-call averages so the three implementations can be compared. */
#include <ec.h>
#include <qt_pegasis.h>
#include <qt_normeq.h>
#include <gmp.h>
#include <stdio.h>
#include <time.h>

static double now_s(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

int main(int argc, char **argv)
{
    unsigned long seed = (argc > 1) ? strtoul(argv[1], NULL, 10) : 1UL;
    int ncases = (argc > 2) ? atoi(argv[2]) : 20;

    gmp_randstate_t rng;
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, seed);

    double t_normeq = 0, t_action = 0;
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

        ec_curve_t E, Ea, Eabar;
        ec_curve_init(&E);
        double t2 = now_s();
        compute_action(&Ea, &Eabar, &E, &sol);
        double t3 = now_s();

        t_normeq += (t1 - t0);
        t_action += (t3 - t2);
        done++;

        ibz_finalize(&N); ibz_finalize(&a); ibz_finalize(&b);
    }

    gmp_randclear(rng);
    printf("C timing over %d actions (from A0):\n", done);
    printf("  qt_normeq (Step 1)      : %.3f ms/action\n", 1e3 * t_normeq / done);
    printf("  compute_action (2-3, 4D): %.3f ms/action\n", 1e3 * t_action / done);
    printf("  full C total            : %.3f ms/action\n", 1e3 * (t_normeq + t_action) / done);
    return 0;
}
