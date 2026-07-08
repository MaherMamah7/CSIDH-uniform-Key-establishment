# Near-uniform ideal-class sampler (Cayley-graph random walk)

Generates a near-uniform random ideal class in the class group `Cl(O)` of the
qt-Pegasis Frobenius order

    O = Z[(1 + sqrt(-p))/2],   discriminant  D = -p,   p = 3 (mod 4)

by a random walk on the Cayley graph of `Cl(O)`, then returns the Gaussian-
reduced (canonical, smallest-norm) representative of the resulting class.

This is standalone: it depends only on **GMP**, not on the qt-Pegasis build. It
uses the same primes `p` as qt-Pegasis (`p = odd * 2^e - 1` for each PNAME).

## Algorithm

    x_0 = principal class (identity)
    for i = 1 .. t:
        sample a split prime ideal p_i  with norm N(p_i) = ell_i <= L
        pick eps_i in {-1, +1} uniformly
        x_i = x_{i-1} * p_i^{eps_i}          (kept Gaussian-reduced)
    output reduce(x_t)

with (parameter A > 2, chosen via `-A`, default 5)

    L = (ln p)^A                                 # norm bound on the sampled prime ideals
    t = round( (3/4) * ln p / ((A/2 - 1) ln ln p) )   # number of steps

Larger A means a bigger generating set (larger L) and fewer steps t; as A -> 2+
the step count diverges, so A must be strictly greater than 2. (`ln` is the
natural logarithm.)

### Choosing A for speed

Because the class element is Gaussian-reduced after every step, the per-step
cost is essentially constant in A (one composition on ~log p-bit numbers). So
total time is roughly proportional to t and *decreases* as A grows -- until the
sampled primes ell <= L = (ln p)^A stop fitting in a 64-bit machine word
(L ~ 2^64), where per-step cost jumps. The empirical optimum is therefore

    A_opt  ~  64*ln2 / ln(ln p)  ~  44 / ln(ln p)

i.e. about 7.6 at 500-bit, 6.1 at 2000-bit, 5.6 at 4000-bit. Measured walk time
(ms/class): 500-bit 0.63/0.48/0.38 for A=4/5/6; 4000-bit 6.4/5.3/7.8 for
A=4/5/6. The default A=5 is a robust choice: fast and below the word-size cliff
at every supported size. Use A=6 for extra speed at <=2000-bit.

Ideal classes are represented as primitive positive-definite binary quadratic
forms `(a, b, c)` with `b^2 - 4ac = D`. Class multiplication is Dirichlet
composition; the class inverse is the conjugate form `(a, -b, c)`; the canonical
representative is the Gaussian-reduced form (`-a < b <= a <= c`). A split prime
`ell` (`kronecker(D, ell) = 1`) gives the prime form `(ell, b_ell, c_ell)`.

The reduced form's `a` is the least norm of an integral ideal in the class.

## Build

    make            # builds `cayley` and `test_qform`
    make test       # runs correctness tests (class groups C3, C5, C7)

Requires a C11 compiler and GMP. On macOS a Homebrew GMP is auto-detected;
otherwise pass `GMP_PREFIX=/path/to/gmp`.

## Run

    ./cayley <PNAME> [-A a] [-s SEED] [-v]   # PNAME in {500,1000,1500,2000,4000}
    ./cayley -p <decimal-prime> [-A a] [-s SEED]

Flags:
    -A a   set the exponent in L = (ln p)^A (real number > 2; default 5)
    -s n   PRNG seed (default 1)
    -r n   repeat the walk n times and report average time (timing mode)
    -v     print each walk step (ell_i, eps_i, running norm)

Examples:

    ./cayley 500              # 505-bit qt-Pegasis prime, A=5, seed 1
    ./cayley 500 -A 6         # larger generating set, fewer steps
    ./cayley 500 -A 2.5       # smaller L, more steps
    ./cayley 2000 -s 7        # 2017-bit prime, seed 7
    ./cayley 500 -v           # print each walk step
    ./cayley -p 23 -s 3       # custom prime (must be 3 mod 4)

Or via make:

    make run PNAME=1000 SEED=5

Output is the reduced form `(a, b, c)` of discriminant `-p` and its norm `a`
(~ half the bit-length of `p`, as expected for a reduced class). The result is
deterministic in the seed; different seeds give (near-)independent uniform
classes.

## Timing

Use `-r n` to repeat the walk `n` times and report the average; it prints
`ms/class` and `ms/step`:

    ./cayley 500 -r 1000            # average one 505-bit walk over 1000 reps
    ./cayley 2000 -A 6 -r 200       # 2017-bit, A=6, 200 reps

To time a single walk from the shell instead:

    make
    time ./cayley 2000 -s 1

## Files

- `qform.h` / `qform.c` — binary quadratic form arithmetic: reduction,
  Dirichlet composition, inverse, prime-ideal forms (GMP).
- `cayley.c`            — parameter setup (`p`, `L`, `t`), prime-ideal sampling,
  the random walk, and final reduction.
- `test_qform.c`       — validates composition/reduction against C3/C5/C7.

## Notes / caveats

- Only odd split primes are used as generators (the standard generating set;
  under GRH the small split primes generate `Cl(O)`).
- The walk is non-constant-time (it is a rejection sampler over primes) and is
  intended for research/experimentation, matching the qt-Pegasis Step-1 style.
- "Near-uniform" is the standard mixing result for such walks with a
  sub-exponential generating set; this code implements the walk, it does not
  prove the mixing bound.
