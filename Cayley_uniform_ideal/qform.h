/* Binary quadratic forms of discriminant D < 0, D = 1 (mod 4).
 *
 * A primitive positive-definite form (a, b, c) with b^2 - 4ac = D represents an
 * ideal class of the imaginary quadratic order O of discriminant D. For the
 * qt-Pegasis Frobenius order O = Z[(1+sqrt(-p))/2] we have D = -p (p = 3 mod 4).
 *
 * Provides: reduction (canonical representative), Dirichlet composition (class
 * group multiplication), inverse (conjugate class), and prime-ideal forms.
 */
#ifndef QFORM_H
#define QFORM_H

#include <gmp.h>

typedef struct {
    mpz_t a, b, c;
} qform_t;

void qform_init(qform_t *f);
void qform_clear(qform_t *f);
void qform_set(qform_t *dst, const qform_t *src);
void qform_print(const char *name, const qform_t *f);

/* discriminant b^2 - 4ac into out */
void qform_disc(mpz_t out, const qform_t *f);

/* principal (identity) form of discriminant D: (1, 1, (1-D)/4) */
void qform_principal(qform_t *f, const mpz_t D);

/* Gaussian reduction: replace f by the unique reduced form in its class
 * (-a < b <= a <= c, and b >= 0 if a == c or a == b). */
void qform_reduce(qform_t *f);

/* class-group inverse: conjugate form (a, -b, c), then reduced */
void qform_inverse(qform_t *dst, const qform_t *src);

/* Dirichlet composition dst = f1 * f2 (then reduced). dst may alias inputs. */
void qform_compose(qform_t *dst, const qform_t *f1, const qform_t *f2, const mpz_t D);

/* Build the prime-ideal form of norm ell (an odd prime) for discriminant D.
 * Requires ell to split, i.e. kronecker(D, ell) == 1. Returns 1 on success,
 * 0 if ell does not split (form left untouched). */
int qform_prime(qform_t *f, const mpz_t ell, const mpz_t D);

#endif
