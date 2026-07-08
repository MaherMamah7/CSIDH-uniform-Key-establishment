#include <bench_test_arguments.h>
#include <bench.h>
#include <rng.h>
#include <constants.h>
#include <stdlib.h>
#include <qt_pegasis.h>
#include <qt_normeq.h>
#include <solutions.h>

// for median computation
int
uint64_comparison(const void *a, const void *b)
{
    return (*(uint64_t *)a < *(uint64_t *)b);
}

int
bench_qt_pegasis_from_sol(int iterations)
{
    uint64_t start, end, median;
    uint64_t sum = 0;
    uint64_t *timings;
    float dt_ns = 0, t, ns;

    ec_curve_t E, Ea, Eabar;
    timings = malloc(iterations * sizeof(uint64_t));
    ec_curve_init(&E);
    printf("\n\nRunning qt_pegasis benchmarks with %d iterations and the %d-bit prime\n",
           iterations,
           ibz_bitsize(&(CHARACTERISTIC)));
    // Bug on entries

    for (int i = 0; i < iterations; i++) {
        
        #ifndef NDEBUG
        printf("Generated input %d\n", i);
        #endif
        t = clock();
        start = cpucycles();
        compute_action(&Ea, &Eabar, &E, &SOL[i]);
        end = cpucycles();
        ns = (1000000000. * (float)(clock()-t) / CLOCKS_PER_SEC);
        dt_ns = dt_ns+ns; 

        timings[i] = end-start;
        sum = sum + timings[i];
        copy_curve(&E, &Ea);

        //#ifndef NDEBUG
        printf("Benchmark %d passed\n", i);
        //#endif

    }
    qsort(timings, iterations, sizeof(uint64_t), uint64_comparison);
    
    //for (int i = 0; i < iterations; i++) {
    //    ibz_printf("%d %" PRIu64 "\n", i, timings[i]);
    //}
    median = timings[(int)(iterations / 2)];
    printf("qt_pegasis (without normeq) took %" PRIu64 " cycles [%.3f ms] on average per iteration\n", sum / iterations, (dt_ns / (1000000 * iterations)));
    printf("                             and %" PRIu64 " cycles was the median over %d iterations\n", median, iterations);
    free(timings);
    return 0;
}

int
main(int argc, char *argv[])
{
    uint32_t seed[12];
    int iterations = 100;
    int help = 0;
    int seed_set = 100;
    int tests = 1;
    int invalid = 0;



#ifndef NDEBUG
    fprintf(stderr,
            "\x1b[31mIt looks like the code was compiled with assertions enabled.\n"
            "This will severely impact performance measurements.\x1b[0m\n");
#endif

    for (int i = 1; i < argc; i++) {
        if (!tests && strcmp(argv[i], "--sol") == 0) {
            tests = 1;
            continue;
        }

        if (!seed_set && !parse_seed(argv[i], seed)) {
            seed_set = 1;
            continue;
        }

        if (!help && strcmp(argv[i], "--help") == 0) {
            help = 1;
            continue;
        }

        if (sscanf(argv[i], "--iterations=%d", &iterations) == 1) {
            continue;
        }
    }

    invalid = invalid || (argc > 3);
    invalid = invalid || (iterations < 0) || ((tests) && (iterations > 1000));

    if (help || invalid) {
        if (invalid) {
            printf("Invalid input\n");
        }
        printf("Usage: %s [--bitsize=<bitsize>] [--iterations=<iterations>] [--sol] [--seed=<seed>]\n", argv[0]);
        printf("Where <iterations> is the number of iterations used for benchmarking (below 10 if --sol is set); if "
               "not present, uses the "
               "default: %d)\n",
               iterations);
        printf("Where <seed> is the random seed to be used; if not present, a random seed is generated\n");
        printf("If --sol is passed, precomputed solutions from solutions.h are used to test the HD part\n");
        printf("Output has last bit set if tests failed, second-to-last if randomness failed\n");
        return 1;
    }

    if (!seed_set) {
        randombytes_select((unsigned char *)seed, sizeof(seed));
    }

    print_seed(seed);

#if defined(TARGET_BIG_ENDIAN)
    for (int i = 0; i < 12; i++) {
        seed[i] = BSWAP32(seed[i]);
    }
#endif

    randombytes_init((unsigned char *)seed, NULL, 256);
    cpucycles_init();
    
    
    return bench_qt_pegasis_from_sol(iterations);
}