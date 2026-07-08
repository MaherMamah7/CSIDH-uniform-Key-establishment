#include "qform.h"
#include <stdlib.h>

void qform_init(qform_t *f)  { mpz_inits(f->a, f->b, f->c, NULL); }
void qform_clear(qform_t *f) { mpz_clears(f->a, f->b, f->c, NULL); }
void qform_set(qform_t *d, const qform_t *s) { mpz_set(d->a, s->a); mpz_set(d->b, s->b); mpz_set(d->c, s->c); }

void qform_print(const char *name, const qform_t *f)
{
    gmp_printf("%s(a=%Zd, b=%Zd, c=%Zd)\n", name, f->a, f->b, f->c);
}

void qform_disc(mpz_t out, const qform_t *f)
{
    mpz_t t; mpz_init(t);
    mpz_mul(out, f->b, f->b);       /* b^2 */
    mpz_mul(t, f->a, f->c); mpz_mul_ui(t, t, 4);
    mpz_sub(out, out, t);           /* b^2 - 4ac */
    mpz_clear(t);
}

void qform_principal(qform_t *f, const mpz_t D)
{
    mpz_set_ui(f->a, 1);
    mpz_set_ui(f->b, 1);
    mpz_sub_ui(f->c, D, 1);         /* D - 1 */
    mpz_neg(f->c, f->c);            /* 1 - D */
    mpz_fdiv_q_2exp(f->c, f->c, 2); /* (1 - D)/4 */
}

/* Tonelli-Shanks: square root of n mod odd prime p. 1 on success. */
static int sqrtmod(mpz_t rop, const mpz_t n_in, const mpz_t p)
{
    mpz_t n, tmp, q, z, c, t, r, b, e, pm1;
    mpz_inits(n, tmp, q, z, c, t, r, b, e, pm1, NULL);
    mpz_mod(n, n_in, p);
    if (mpz_sgn(n) == 0) { mpz_set_ui(rop, 0); goto done_ok; }
    if (mpz_legendre(n, p) != 1) { mpz_clears(n,tmp,q,z,c,t,r,b,e,pm1,NULL); return 0; }
    if (mpz_fdiv_ui(p, 4) == 3) {
        mpz_add_ui(e, p, 1); mpz_fdiv_q_2exp(e, e, 2);
        mpz_powm(rop, n, e, p);
        goto done_ok;
    }
    mpz_sub_ui(pm1, p, 1); mpz_set(q, pm1);
    unsigned long s = 0;
    while (mpz_even_p(q)) { mpz_fdiv_q_2exp(q, q, 1); s++; }
    mpz_set_ui(z, 2);
    while (mpz_legendre(z, p) != -1) mpz_add_ui(z, z, 1);
    mpz_powm(c, z, q, p);
    mpz_powm(t, n, q, p);
    mpz_add_ui(tmp, q, 1); mpz_fdiv_q_2exp(tmp, tmp, 1);
    mpz_powm(r, n, tmp, p);
    unsigned long m = s;
    while (mpz_cmp_ui(t, 1) != 0) {
        unsigned long i = 0; mpz_set(tmp, t);
        while (mpz_cmp_ui(tmp, 1) != 0) { mpz_powm_ui(tmp, tmp, 2, p); i++; if (i == m) break; }
        mpz_set(b, c);
        for (unsigned long k = 0; k + i + 1 < m; k++) mpz_powm_ui(b, b, 2, p);
        mpz_powm_ui(c, b, 2, p);
        mpz_mul(t, t, c); mpz_mod(t, t, p);
        mpz_mul(r, r, b); mpz_mod(r, r, p);
        m = i;
    }
    mpz_set(rop, r);
done_ok:
    mpz_clears(n, tmp, q, z, c, t, r, b, e, pm1, NULL);
    return 1;
}

/* normalize: shift b into (-a, a] and recompute c */
static void normalize(qform_t *f, const mpz_t D)
{
    mpz_t q, twoa, t; mpz_inits(q, twoa, t, NULL);
    mpz_mul_2exp(twoa, f->a, 1);                 /* 2a */
    /* q = round((a - b) / (2a)) via floor((a - b) / (2a))? we want b' in (-a,a].
       b' = b + 2a*k with b' in (-a, a]. k = floor((a - b)/(2a)) gives b' in (-a, a]. */
    mpz_sub(t, f->a, f->b);                      /* a - b */
    mpz_fdiv_q(q, t, twoa);                       /* k = floor((a-b)/2a) */
    mpz_mul(t, q, twoa); mpz_add(f->b, f->b, t);  /* b += 2a*k */
    /* c = (b^2 - D)/(4a) */
    mpz_mul(t, f->b, f->b); mpz_sub(t, t, D);
    mpz_mul_2exp(f->a, f->a, 0); /* no-op keep */
    { mpz_t fa4; mpz_init(fa4); mpz_mul_2exp(fa4, f->a, 2); mpz_divexact(f->c, t, fa4); mpz_clear(fa4); }
    mpz_clears(q, twoa, t, NULL);
}

void qform_reduce(qform_t *f)
{
    mpz_t D, t; mpz_inits(D, t, NULL);
    qform_disc(D, f);
    normalize(f, D);
    while (mpz_cmp(f->a, f->c) > 0) {
        /* rho step: (a,b,c) -> (c, -b, a) then normalize */
        mpz_swap(f->a, f->c);
        mpz_neg(f->b, f->b);
        normalize(f, D);
    }
    /* canonical sign: if a == c or a == b then b >= 0 */
    if (mpz_sgn(f->b) < 0) {
        if (mpz_cmp(f->a, f->c) == 0) {          /* a == c -> b = |b| */
            mpz_neg(f->b, f->b);
            mpz_mul(t, f->b, f->b); mpz_sub(t, t, D);
            { mpz_t fa4; mpz_init(fa4); mpz_mul_2exp(fa4, f->a, 2); mpz_divexact(f->c, t, fa4); mpz_clear(fa4); }
        } else {
            mpz_t nb; mpz_init(nb); mpz_neg(nb, f->b);
            if (mpz_cmp(nb, f->a) == 0) {         /* b == -a -> b = a */
                mpz_set(f->b, f->a);
                mpz_mul(t, f->b, f->b); mpz_sub(t, t, D);
                { mpz_t fa4; mpz_init(fa4); mpz_mul_2exp(fa4, f->a, 2); mpz_divexact(f->c, t, fa4); mpz_clear(fa4); }
            }
            mpz_clear(nb);
        }
    }
    mpz_clears(D, t, NULL);
}

void qform_inverse(qform_t *dst, const qform_t *src)
{
    qform_set(dst, src);
    mpz_neg(dst->b, dst->b);
    qform_reduce(dst);
}

/* Dirichlet composition. dst = reduce(f1 * f2). Uses d = gcd(a1,a2,s) via two
 * extended gcds. Validated against known class groups and disc preservation
 * on random compositions (see test_qform.c). */
void qform_compose(qform_t *dst, const qform_t *f1, const qform_t *f2, const mpz_t D)
{
    mpz_t a1,b1,c1,a2,b2,c2;
    mpz_inits(a1,b1,c1,a2,b2,c2,NULL);
    if (mpz_cmp(f1->a, f2->a) <= 0) {                   /* ensure a1 <= a2 */
        mpz_set(a1,f1->a); mpz_set(b1,f1->b); mpz_set(c1,f1->c);
        mpz_set(a2,f2->a); mpz_set(b2,f2->b); mpz_set(c2,f2->c);
    } else {
        mpz_set(a1,f2->a); mpz_set(b1,f2->b); mpz_set(c1,f2->c);
        mpz_set(a2,f1->a); mpz_set(b2,f1->b); mpz_set(c2,f1->c);
    }

    mpz_t s, m, d0, x, y, d, u, w, A, a1d, a2d, t, b3, c3, tmp;
    mpz_inits(s, m, d0, x, y, d, u, w, A, a1d, a2d, t, b3, c3, tmp, NULL);

    mpz_add(s, b1, b2); mpz_fdiv_q_2exp(s, s, 1);       /* s = (b1+b2)/2 */
    mpz_sub(m, b1, b2); mpz_fdiv_q_2exp(m, m, 1);       /* m = (b1-b2)/2 */

    mpz_gcdext(d0, x, y, a2, a1);                       /* x*a2 + y*a1 = d0 = gcd(a1,a2) */
    mpz_gcdext(d, u, w, d0, s);                         /* u*d0 + w*s = d = gcd(a1,a2,s) */

    mpz_mul(A, a1, a2);
    mpz_divexact(A, A, d); mpz_divexact(A, A, d);       /* A = a1*a2/d^2 */
    mpz_divexact(a1d, a1, d);                           /* a1/d */
    mpz_divexact(a2d, a2, d);                           /* a2/d */

    /* t = ( u*x*m - w*c2 ) mod (a1/d) */
    mpz_mul(t, u, x); mpz_mul(t, t, m);
    mpz_mul(tmp, w, c2);
    mpz_sub(t, t, tmp);
    mpz_mod(t, t, a1d);

    /* b3 = b2 + 2*(a2/d)*t */
    mpz_mul(b3, a2d, t); mpz_mul_2exp(b3, b3, 1); mpz_add(b3, b3, b2);

    /* center b3 into (-A, A] mod 2A */
    mpz_mul_2exp(tmp, A, 1);
    mpz_mod(b3, b3, tmp);
    if (mpz_cmp(b3, A) > 0) mpz_sub(b3, b3, tmp);

    /* c3 = (b3^2 - D)/(4A) */
    mpz_mul(c3, b3, b3); mpz_sub(c3, c3, D);
    mpz_mul_2exp(tmp, A, 2);
    mpz_divexact(c3, c3, tmp);

    mpz_set(dst->a, A); mpz_set(dst->b, b3); mpz_set(dst->c, c3);
    qform_reduce(dst);

    mpz_clears(s, m, d0, x, y, d, u, w, A, a1d, a2d, t, b3, c3, tmp, NULL);
    mpz_clears(a1,b1,c1,a2,b2,c2,NULL);
}

int qform_prime(qform_t *f, const mpz_t ell, const mpz_t D)
{
    mpz_t Dm, r, b, c, tmp; mpz_inits(Dm, r, b, c, tmp, NULL);
    if (mpz_kronecker(D, ell) != 1) { mpz_clears(Dm,r,b,c,tmp,NULL); return 0; }
    mpz_mod(Dm, D, ell);
    if (!sqrtmod(r, Dm, ell)) { mpz_clears(Dm,r,b,c,tmp,NULL); return 0; }
    /* choose b odd with b == r (mod ell) */
    mpz_set(b, r);
    if (mpz_even_p(b)) mpz_add(b, b, ell);              /* r + ell is odd (ell odd) */
    /* c = (b^2 - D) / (4*ell) */
    mpz_mul(c, b, b); mpz_sub(c, c, D);
    mpz_mul_ui(tmp, ell, 4);
    mpz_divexact(c, c, tmp);
    mpz_set(f->a, ell); mpz_set(f->b, b); mpz_set(f->c, c);
    qform_reduce(f);
    mpz_clears(Dm, r, b, c, tmp, NULL);
    return 1;
}
