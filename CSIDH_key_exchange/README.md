# Group-action Diffie–Hellman (CSIDH / qt-Pegasis style)

A non-interactive key exchange built on the qt-Pegasis effective class-group
action, in the CSIDH style.

    E0 = starting curve (j = 1728)
    Alice: secret ideal [a] (uniform Cayley sampler) -> public E_A = [a] E0
    Bob:   secret ideal [b] (uniform Cayley sampler) -> public E_B = [b] E0
    exchange only the Montgomery coefficient A of E_A, E_B
    Alice: [a] E_B = [ab] E0        Bob: [b] E_A = [ab] E0

By commutativity of the class-group action `[a][b]E0 = [b][a]E0 = [ab]E0`, both
parties land on the same curve and derive a shared secret from a single field
element (one Montgomery coefficient).

Everything shares the qt-Pegasis Frobenius order `O = Z[(1+sqrt(-p))/2]`
(discriminant `D = -p`, `p = 3 mod 4`): the sampler walks `Cl(O)`, and the
action evaluates the ideal on curves over the same `p`.

## Canonical shared coefficient

The qt-Pegasis `2^e / 4D`-split action fixes the target curve only up to its
Montgomery *model* and its quadratic *twist* (`A <-> -A`), and its internal twist
selection is keyed on the ideal rather than on the target curve. So the two
parties' raw coefficients agree only about half the time (their `j`-invariants
always agree).

To recover a single shared value, canonicalisation is applied **in this protocol
layer, to the final shared curve only** (never to the intermediate public keys,
whose orientation the second action still needs). `canonical_A` returns the
deterministic minimum coefficient over every `F_p` Montgomery model of the shared
curve *and its twist*, so both parties agree on one coefficient. This keeps the
qt-Pegasis and Cayley-sampler code untouched.

## Build

    make PNAME=500        # or 1000 / 1500 / 2000 / 4000

Self-contained: it compiles the required qt-Pegasis sources plus
`../Cayley_uniform_ideal/qform.c` directly, with `-DP_<PNAME>`. Needs a C11
compiler and GMP (auto-detected on macOS Homebrew; else pass
`GMP_PREFIX=/path/to/gmp`). **Rebuild for each PNAME** (`make clean` first).

## Run

    ./dh_<PNAME> [seed] [A]

- `seed` : PRNG seed for both parties' ideal sampling (default 1)
- `A`    : L-exponent for the sampler. If omitted it is chosen automatically:
           **A = 6 when p has more than 2000 bits, else A = 5.**

Examples:

    ./dh_500 1            # 505-bit prime, auto A=5, seed 1
    ./dh_2000 7           # 2017-bit prime, auto A=6
    ./dh_500 3 6          # force A=6

Output: Alice's and Bob's raw shared coefficients, the canonical shared secret,
and `KEY EXCHANGE: SUCCESS` when the canonical coefficients match.

## Timing

The program prints a per-phase breakdown for the run:

    timing:
      keygen (2 parties, sample+normeq): ...   # walk + norm equation, x2
      public keys (2 actions)          : ...   # [a]E0 and [b]E0
      shared secret (2 actions)        : ...   # [a]E_B and [b]E_A
      per group action                 : ...   # one compute_action

`compute_action` dominates and is the steady-state per-action cost (a
fixed-length 4D chain, independent of the ideal): ~47 ms (500-bit), ~2.3 s
(2000-bit), ~18 s (4000-bit); the sampler is sub-ms to a few ms.

### Time the whole exchange from the shell

Build once, then wrap the run in `time`:

    make PNAME=500 dh
    time ./dh_500 1

Or in a single line (build + timed run):

    make PNAME=500 dh && time ./dh_500 1

`time` reports the wall-clock `real` for the complete exchange (keygen + public
keys + shared secret + canonicalisation); the program's own breakdown attributes
that time to each phase.

## Files

- `dh.c`      — the protocol: ideal sampling (via the Cayley sampler), key
  generation (`qt_normeq`), public keys and shared secret (`compute_action`),
  the canonical-coefficient reduction, and timing.
- `Makefile`  — per-PNAME self-contained build.

## Notes / caveats

- Secret ideals are the reduced classes from the Cayley sampler; a fresh class
  is drawn if `qt_normeq`'s constraints (`N != 2 mod 4`, `gcd(b, a) = 1`,
  solvable norm equation) fail. In practice 0–1 resamples.
- Research/experimentation tool (non-constant-time). It demonstrates that the
  qt-Pegasis action supports a commutative-group-action key exchange when the
  shared curve is reduced to a canonical Montgomery coefficient.
