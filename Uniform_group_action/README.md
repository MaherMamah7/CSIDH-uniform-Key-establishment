# Uniform group action

Takes a near-uniform ideal class produced by the Cayley-graph sampler
(`../Cayley_uniform_ideal`) and **acts with it on the qt-Pegasis starting curve
E0 (j = 1728)** using the qt-Pegasis pipeline (`qt_normeq` + `compute_action`),
then checks the group law.

Both sides use the *same* prime `p = CHARACTERISTIC` and the *same* order
`O = Z[(1+sqrt(-p))/2]` (discriminant `D = -p`), so the only glue needed is a
representation bridge.

## The representation bridge

The sampler returns an ideal class as a reduced binary quadratic form `(a,b,c)`
of discriminant `D`, i.e. the ideal `[a, (b+sqrt(D))/2]`. qt-Pegasis wants the
ideal as `(N, alpha = alpha_a + alpha_b * omega)` with `omega = (1+sqrt(D))/2`.
Since `(b+sqrt(D))/2 = (b-1)/2 + omega`:

    N       = a
    alpha_a = (b-1)/2 ,  alpha_b = 1            # the ideal [a]
    conj:  alpha_a = (b+1)/2 ,  alpha_b = -1    # its conjugate class, same N

`qt_normeq` then solves the norm equation and `compute_action` evaluates the
dimension-4 isogeny, exactly as in `../src/action/test/qt_normeq_e2e.c`.

### Acceptance constraints (handled by resampling)

`qt_normeq` needs (as in the qt-Pegasis reference):
- `N != 2 (mod 4)`, and
- `gcd(trace(alpha)=b, N=a) = 1`  (so the trace is invertible mod N).

The driver samples a fresh uniform class if these fail or if the norm equation
has no solution. In practice this needs 0-1 resamples.

## Build

    make PNAME=500        # or 1000 / 1500 / 2000 / 4000

Self-contained: it compiles the required qt-Pegasis sources plus
`../Cayley_uniform_ideal/qform.c` directly, with `-DP_<PNAME>`. Needs GMP.
No prior qt-Pegasis build is required. **Rebuild for each PNAME.**

## Run

    ./uniform_action_<PNAME> [seed] [A]

- `seed` : PRNG seed (default 1)
- `A`    : L-exponent for the sampler. If omitted it is chosen automatically:
           **A = 6 when p has more than 2000 bits, else A = 5.**

Examples:

    ./uniform_action_500 1        # auto A=5
    ./uniform_action_2000 7       # auto A=6 (2017-bit)
    ./uniform_action_500 1 6      # force A=6

Output: the sampled reduced class, its `(N, alpha)` form, `j(E0)`, `j([a]E0)`,
and two checks:
- **round trip** `[abar]([a]E0) == E0` (the group law, since `a * abar = (N)` is principal),
- **twin-j** (the two twists from one action share their j-invariant).

## Timing

The program prints a timing breakdown for the run:

    timing:
      sampling (walk, incl. resamples) : ...
      qt_normeq (Step 1)               : ...
      compute_action (Steps 2-3)       : ...
      one group action ([a]E0)         : ...   # walk + norm eq + action

`compute_action` dominates and is the steady-state per-action cost (it runs a
fixed-length 4D chain independent of the ideal): ~47 ms (500-bit), ~2.3 s
(2000-bit), ~18 s (4000-bit); the sampler is sub-ms to a few ms.

For a benchmark averaged over many actions, use the repo-root
`make PNAME=500 bench_full`.

To time the whole process from the shell instead:

    time ./uniform_action_2000 1

## Status

Validated: round trip and twin-j pass at 500-bit (seeds 1-8, 0-1 resamples each)
and at 2000-bit. Because the orders match exactly, the uniform-sampled ideal
acts correctly through the qt-Pegasis 4D machinery.

## Notes / caveats

- The action is by the reduced ideal of *composite* norm `a ~ sqrt(p)`; this is
  within the qt-Pegasis norm-equation regime (their own sampled primes are also
  ~sqrt(p)), and `qt_normeq` accepts composite `N` as long as the two
  constraints above hold.
- This is a research/experimentation tool (non-constant-time), matching the
  qt-Pegasis Step-1 style.
