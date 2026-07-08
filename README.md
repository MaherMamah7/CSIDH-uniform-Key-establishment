# CSIDH-Uniform Key Establishment using qt-Pegasis — effective class-group action in C

A standalone **C** implementation of CSIDH-Key exchange using the qt-Pegasis effective class-group action,
together a uniform ideal sampler on Cayley graph.
The core solves the **norm equation in C** and thus complements the earlier work
of Dartois and Duparc, so the entire algorithm runs end-to-end in C with **no
SageMath in the loop**.

Maintaind by *Maher Mamah*
---

## Repository layout

Each component below is self-contained and **has its own `README.md`** with build,
run, and timing instructions.

| Path                    | What it is                                                                                   | README |
|-------------------------|----------------------------------------------------------------------------------------------|--------|
| `src/` (+ root Makefile)| **qt-Pegasis core** — the effective class-group action: norm equation (Step 1, in C) plus the dimension-4 isogeny chain (Steps 2–3). | this file |
| `Cayley_uniform_ideal/` | **Near-uniform ideal-class sampler** — a random walk on the Cayley graph of `Cl(O)`, returning a Gaussian-reduced ideal class. Standalone (GMP only). | [`Cayley_uniform_ideal/README.md`](Cayley_uniform_ideal/README.md) |
| `Uniform_group_action/` | **Uniform group action** — takes a sampled ideal class and acts with it on the starting curve `E0` (`j = 1728`) via the qt-Pegasis pipeline; checks the group law. | [`Uniform_group_action/README.md`](Uniform_group_action/README.md) |
| `CSIDH_key_exchange/`   | **Group-action Diffie–Hellman** — a CSIDH-style non-interactive key exchange: two parties sample secret ideals, act on `E0`, and agree on one canonical Montgomery coefficient. | [`CSIDH_key_exchange/README.md`](CSIDH_key_exchange/README.md) |

All four share the same Frobenius order `O = Z[(1+sqrt(-p))/2]` (discriminant
`D = -p`, `p = 3 mod 4`) and the same per-prime constants, so ideals sampled by
one component act correctly through another.

---

## Quick start (qt-Pegasis core)

Requirements: a C11 compiler and **GMP** (6.0.0+). On macOS the Makefile
auto-detects a Homebrew GMP (`/opt/homebrew` or `/usr/local`); on Linux the
system GMP is used. If GMP lives elsewhere, pass `GMP_PREFIX=/path/to/gmp`.

```sh
# Pick a prime size: 500, 1000, 1500, 2000, or 4000.
make PNAME=500 run          # end-to-end demo: sample ideal -> norm eq (C) ->
                            # 4D action (C); verifies the group-law round trip
make PNAME=500 test         # correctness tests (fp, ec, dim4, action)
make PNAME=500 bench_full   # time the full-C action: Step 1 + Steps 2-3
make PNAME=500 bench        # time fp arithmetic and the 4D chain

make clean                  # remove build artifacts
```

Always `make clean` before switching `PNAME` (the build tree is prime-specific).

Example `make PNAME=500 run` output:

```
case 0: OK  (roundtrip=yes, twin-j=yes)
...
qt_normeq end-to-end: 8/8 passed
```

To act on your own ideal from C, call `qt_normeq(&sol, &N, &alpha_a, &alpha_b)`
(declared in `src/qt_interface/include/qt_normeq.h`) to fill a solution, then
`compute_action(&Ea, &Eabar, &E, &sol)`. See `src/action/test/qt_normeq_e2e.c`.

---

## Timing

Each component prints its own timing breakdown and documents how to time it in
its README. Quick reference:

```sh
# qt-Pegasis core: averaged full-C action (Step 1 + Steps 2-3)
make PNAME=500 bench_full

# CSIDH key exchange: time the whole exchange from the shell
cd CSIDH_key_exchange && make PNAME=500 dh && time ./dh_500 1

# uniform group action: time one sample + norm eq + action
cd Uniform_group_action && make PNAME=500 && time ./uniform_action_500 1

# Cayley sampler: time a batch of walks
cd Cayley_uniform_ideal && make && ./cayley 500 -r 1000
```

`compute_action` is the steady-state per-action cost (a fixed-length 4D chain,
independent of the ideal): ~47 ms (500-bit), ~2.3 s (2000-bit), ~18 s (4000-bit).

---

## What this project adds

The three algorithmic steps of a qt-Pegasis action are:

1. **Norm equation** — solve `N1 + N2 = 2^e` for a sampled ideal.
2. **Kernel points** — build the dimension-4 kernel torsion basis.
3. **Isogeny chain** — the dimension-4 `(2,2)`-isogeny chain and splitting.

The upstream C library implements Steps 2–3 only. This project adds **Step 1 in
C** — ideal sampling, sum-of-two-squares, and the `qlapoti` norm-equation
solver (`src/qt_interface/qt_normeq_solve.c`) — wired into the action via
`qt_normeq()`. The results were validated against the SageMath reference:
identical output j-invariants from the same starting curve, and the group-law
round trip `[ā][a]E ≅ E`.

On top of the core, `Cayley_uniform_ideal`, `Uniform_group_action`, and
`CSIDH_key_exchange` add near-uniform sampling of `Cl(O)`, a bridge that applies a
sampled class to `E0`, and a commutative-group-action key exchange, respectively.

---

## Credits / provenance

This work builds directly on the following, with thanks:

- **qt-Pegasis** — P. Dartois, J. Komada Eriksen, R. Invernizzi, F. Vercauteren,
  *"qt-Pegasis: Simpler and Faster Effective Class Group Actions"*
  (Eurocrypt 2026, [eprint 2025/1859](https://eprint.iacr.org/2025/1859)),
  and its SageMath reference implementation.
- **C 4D-isogeny implementation of qt-Pegasis** — the constant-time
  dimension-4 `(2,2)`-isogeny chain (Steps 2–3), on which the `ec`, `dim4`,
  `action`, `params`, `gf`, and `mp` modules here are based.
- **SQIsign** ([the-sqisign](https://github.com/SQIsign/the-sqisign)) — the
  `quadratic`/`intbig` big-integer module was copied and adapted from it.

---

## Code structure (core `src/`)

- `common`     — utilities (PRNG, benchmarking, mod-4 matrices).
- `params`     — per-prime constants and precomputed values.
- `mp`         — multiprecision arithmetic.
- `gf`         — `GF(p)` arithmetic (generated with `modarith`).
- `ec`         — elliptic-curve arithmetic.
- `dim4`       — dimension-4 `(2,2)`-isogenies.
- `quadratic`  — arbitrary-precision integers (adapted from SQIsign).
- `qt_interface` — **norm-equation solver (Step 1)** and its interface.
- `action`     — the qt-Pegasis group action tying the steps together.

Generating parameters for a new prime uses SageMath + `modarith`; see
`scripts/generate_fp.py`.

---

## License

See [`LICENSE`](LICENSE). Portions derive from SQIsign (LGPL) and the upstream
qt-Pegasis C implementation; their licenses and attribution are retained.
