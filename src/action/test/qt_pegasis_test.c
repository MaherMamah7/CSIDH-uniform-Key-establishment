#include <ec.h>
#include <qt_pegasis.h>
#include <solutions_small.h>
#include <rng.h>
#include <bench_test_arguments.h>

int
test_qt_pegasis_from_sol()
{
    ec_curve_t E, Ea, Eabar;

    ec_curve_init(&E);

    for (int i = 0; i < 10; i++) {

        #ifndef NDEBUG
            printf("testing case %d\n", i);
        #endif
        compute_action(&Ea, &Eabar, &E, &SOL[i]);
        copy_curve(&E, &Ea);
    }

    printf("qt_pegasis tests from precomputed norm equation solutions.............. PASSED\n");
    return 0;
}


int
main(int argc, char *argv[])
{
    uint32_t seed[12] = { 0 };
    int iterations = 10;
    int help = 0;
    int seed_set = 0;
    int res = 0;

    for (int i = 1; i < argc; i++) {
        if (!help && strcmp(argv[i], "--help") == 0) {
            help = 1;
            continue;
        }

        if (!seed_set && !parse_seed(argv[i], seed)) {
            seed_set = 1;
            continue;
        }

        if (sscanf(argv[i], "--iterations=%d", &iterations) == 1) {
            continue;
        }
    }

    if (help || iterations <= 0) {
        printf("Usage: %s [--iterations=<iterations>] [--seed=<seed>]\n", argv[0]);
        printf("Where <iterations> is the number of iterations used for testing; if not "
               "present, uses the default: %d)\n",
               iterations);
        printf("Where <seed> is the random seed to be used; if not present, a random seed is "
               "generated\n");
        return 1;
    }

    if (!seed_set) {
        randombytes_select((unsigned char *)seed, sizeof(seed));
    }

    printf("--------------------------------------------------------------------------------\n\n");
    printf("Testing qt_pegasis executions:\n\n");
    print_seed(seed);

    res = test_qt_pegasis_from_sol();

    return res;
}
