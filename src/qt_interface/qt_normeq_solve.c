/** @file
 *
 * @brief qt-Pegasis Step 1 (norm equation) in C.
 *
 * Supplies the piece the 4D library is missing: sampling a reduced ideal and
 * solving the norm equation N1 + N2 = 2^{e_sol}, producing a qt_output_t that
 * compute_action() consumes. Ported from the SageMath reference
 * (sos.py / norm_eq.py / qt_pegasis.py) and validated against it.
 *
 * ibz_t == mpz_t, so this uses GMP directly. Non-constant-time (matches the
 * current milestone: the norm equation is inherently a rejection sampler).
 */

#include <qt_normeq.h>
#include <string.h>
#include <stdlib.h>

extern const char *QTP_GOOD_PRIME_PROD;
extern const char *QTP_BAD_PRIME_PROD;

#define SIEVE_BOUND 100000UL

/* =========================== small helpers =========================== */

static unsigned long modu(const mpz_t x, unsigned long m) { return mpz_fdiv_ui(x, m); }

/* Tonelli-Shanks sqrt of n mod odd prime p. Returns 1 on success. */
static int sqrtmod(mpz_t rop, const mpz_t n_in, const mpz_t p)
{
    mpz_t n, tmp;
    mpz_inits(n, tmp, NULL);
    mpz_mod(n, n_in, p);
    if (mpz_sgn(n) == 0) { mpz_set_ui(rop, 0); mpz_clears(n, tmp, NULL); return 1; }
    if (mpz_legendre(n, p) != 1) { mpz_clears(n, tmp, NULL); return 0; }
    if (mpz_fdiv_ui(p, 4) == 3) {
        mpz_t e; mpz_init(e);
        mpz_add_ui(e, p, 1); mpz_fdiv_q_2exp(e, e, 2);
        mpz_powm(rop, n, e, p);
        mpz_clears(n, tmp, e, NULL);
        return 1;
    }
    mpz_t q, z, c, t, r, b, pm1;
    mpz_inits(q, z, c, t, r, b, pm1, NULL);
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
    mpz_clears(n, tmp, q, z, c, t, r, b, pm1, NULL);
    return 1;
}

/* =========================== quadratic order =========================== */
/* Element a + b*w, w = (1+sqrt(-p))/2. Context holds p and wnorm=(p+1)/4. */

typedef struct { mpz_t p, wnorm; } qo_ctx;
typedef struct { mpz_t a, b; } qelem;

static void qo_ctx_init(qo_ctx *c, const mpz_t p)
{
    mpz_init_set(c->p, p);
    mpz_init(c->wnorm);
    mpz_add_ui(c->wnorm, p, 1);
    mpz_fdiv_q_2exp(c->wnorm, c->wnorm, 2);
}
static void qo_ctx_clear(qo_ctx *c) { mpz_clears(c->p, c->wnorm, NULL); }
static void qe_init(qelem *x) { mpz_inits(x->a, x->b, NULL); }
static void qe_clear(qelem *x) { mpz_clears(x->a, x->b, NULL); }
static void qe_set(qelem *d, const qelem *s) { mpz_set(d->a, s->a); mpz_set(d->b, s->b); }

static void qo_norm(const qo_ctx *c, mpz_t out, const qelem *x)
{
    mpz_t t; mpz_init(t);
    mpz_mul(out, x->a, x->a);
    mpz_mul(t, x->a, x->b); mpz_add(out, out, t);
    mpz_mul(t, x->b, x->b); mpz_mul(t, t, c->wnorm); mpz_add(out, out, t);
    mpz_clear(t);
}
static void qo_trace(mpz_t out, const qelem *x) { mpz_mul_2exp(out, x->a, 1); mpz_add(out, out, x->b); }
static void qo_conj(qelem *out, const qelem *x)
{
    mpz_t na; mpz_init(na);
    mpz_add(na, x->a, x->b);
    mpz_neg(out->b, x->b);
    mpz_set(out->a, na);
    mpz_clear(na);
}
static void qo_mul(const qo_ctx *c, qelem *out, const qelem *x, const qelem *y)
{
    mpz_t ra, rb, t; mpz_inits(ra, rb, t, NULL);
    mpz_mul(ra, x->a, y->a); mpz_mul(t, x->b, y->b); mpz_mul(t, t, c->wnorm); mpz_sub(ra, ra, t);
    mpz_mul(rb, x->a, y->b); mpz_mul(t, x->b, y->a); mpz_add(rb, rb, t);
    mpz_mul(t, x->b, y->b); mpz_add(rb, rb, t);
    mpz_set(out->a, ra); mpz_set(out->b, rb);
    mpz_clears(ra, rb, t, NULL);
}

/* =========================== sum of two squares =========================== */

typedef struct { mpz_t p; unsigned long e; } qfac;
typedef struct { qfac *v; size_t n, cap; } qfaclist;
static void fl_init(qfaclist *f) { f->v = NULL; f->n = 0; f->cap = 0; }
static void fl_clear(qfaclist *f){ for (size_t i=0;i<f->n;i++) mpz_clear(f->v[i].p); free(f->v); f->v=NULL; f->n=f->cap=0; }
static void fl_push(qfaclist *f, const mpz_t p, unsigned long e)
{
    if (f->n == f->cap) { f->cap = f->cap ? f->cap*2 : 8; f->v = realloc(f->v, f->cap*sizeof(qfac)); }
    mpz_init_set(f->v[f->n].p, p); f->v[f->n].e = e; f->n++;
}
static void fl_push_ui(qfaclist *f, unsigned long p, unsigned long e)
{ mpz_t t; mpz_init_set_ui(t, p); fl_push(f, t, e); mpz_clear(t); }

static void rep_gcd(mpz_t a, const mpz_t b, mpz_t g_out)
{
    mpz_t g; mpz_init(g); mpz_set_ui(g_out, 1);
    mpz_gcd(g, a, b);
    while (mpz_cmp_ui(g, 1) != 0) { mpz_mul(g_out, g_out, g); mpz_divexact(a, a, g); mpz_gcd(g, a, g); }
    mpz_clear(g);
}

static unsigned long *g_primes = NULL;
static size_t g_nprimes = 0;
static void ensure_sieve(void)
{
    if (g_primes) return;
    char *comp = calloc(SIEVE_BOUND + 1, 1);
    for (unsigned long i = 2; i*i <= SIEVE_BOUND; i++)
        if (!comp[i]) for (unsigned long j = i*i; j <= SIEVE_BOUND; j += i) comp[j] = 1;
    size_t cnt = 0;
    for (unsigned long i = 2; i <= SIEVE_BOUND; i++) if (!comp[i]) cnt++;
    g_primes = malloc(cnt * sizeof(unsigned long)); g_nprimes = 0;
    for (unsigned long i = 2; i <= SIEVE_BOUND; i++) if (!comp[i]) g_primes[g_nprimes++] = i;
    free(comp);
}
static void factor_smooth(const mpz_t m_in, qfaclist *f)
{
    ensure_sieve();
    mpz_t m; mpz_init_set(m, m_in);
    for (size_t i = 0; i < g_nprimes && mpz_cmp_ui(m, 1) != 0; i++) {
        unsigned long p = g_primes[i];
        if (!mpz_divisible_ui_p(m, p)) continue;
        unsigned long e = 0;
        while (mpz_divisible_ui_p(m, p)) { mpz_divexact_ui(m, m, p); e++; }
        fl_push_ui(f, p, e);
    }
    if (mpz_cmp_ui(m, 1) != 0) fl_push(f, m, 1);
    mpz_clear(m);
}
static int cornacchia_prime(mpz_t r, mpz_t s, const mpz_t p)
{
    mpz_t neg1, x, t; mpz_inits(neg1, x, t, NULL);
    mpz_set_si(neg1, -1); mpz_mod(neg1, neg1, p);
    if (!sqrtmod(x, neg1, p)) { mpz_clears(neg1, x, t, NULL); return 0; }
    mpz_set(r, p); mpz_set(s, x);
    mpz_mul(t, s, s);
    while (mpz_cmp(t, p) > 0) { mpz_mod(x, r, s); mpz_set(r, s); mpz_set(s, x); mpz_mul(t, s, s); }
    mpz_mod(r, r, s);
    mpz_clears(neg1, x, t, NULL);
    return 1;
}
static int two_squares_from_factors(mpz_t a, mpz_t b, const qfaclist *f)
{
    mpz_t r, s, na, nb, m, tmp; mpz_inits(r, s, na, nb, m, tmp, NULL);
    mpz_set_ui(a, 1); mpz_set_ui(b, 0);
    int ok = 1;
    for (size_t i = 0; i < f->n && ok; i++) {
        const mpz_t *p = (const mpz_t *)&f->v[i].p;
        unsigned long e = f->v[i].e;
        if (mpz_cmp_ui(*p, 1) == 0) continue;
        if ((e & 1) && mpz_fdiv_ui(*p, 4) == 3) { ok = 0; break; }
        if (e >= 2) { mpz_pow_ui(m, *p, e/2); mpz_mul(a, a, m); mpz_mul(b, b, m); }
        if (e & 1) {
            if (mpz_cmp_ui(*p, 2) == 0) { mpz_sub(na, a, b); mpz_add(nb, a, b); mpz_set(a, na); mpz_set(b, nb); }
            else {
                if (!cornacchia_prime(r, s, *p)) { ok = 0; break; }
                mpz_mul(na, a, r); mpz_mul(tmp, b, s); mpz_sub(na, na, tmp);
                mpz_mul(nb, b, r); mpz_mul(tmp, a, s); mpz_add(nb, nb, tmp);
                mpz_set(a, na); mpz_set(b, nb);
            }
        }
    }
    if (ok) { mpz_abs(a, a); mpz_abs(b, b); if (mpz_cmp(a, b) > 0) mpz_swap(a, b); }
    mpz_clears(r, s, na, nb, m, tmp, NULL);
    return ok;
}
/* The good/bad prime products are fixed ~2150-digit constants. Parse them from
 * their decimal strings once and reuse, rather than on every sum_of_squares
 * call (base-10 parsing of thousands of digits is O(d^2) and was a dominant
 * cost of Step 1). */
static mpz_t g_good_prod, g_bad_prod;
static int g_prime_prods_ready = 0;
static void ensure_prime_prods(void)
{
    if (g_prime_prods_ready) return;
    mpz_init_set_str(g_good_prod, QTP_GOOD_PRIME_PROD, 10);
    mpz_init_set_str(g_bad_prod, QTP_BAD_PRIME_PROD, 10);
    g_prime_prods_ready = 1;
}

static int sum_of_squares(mpz_t b1, mpz_t b2, const mpz_t n_in)
{
    if (mpz_sgn(n_in) < 0) return 0;
    if (mpz_sgn(n_in) == 0) { mpz_set_ui(b1, 0); mpz_set_ui(b2, 0); return 1; }
    ensure_prime_prods();
    const mpz_t *good = (const mpz_t *)&g_good_prod;
    const mpz_t *bad = (const mpz_t *)&g_bad_prod;
    mpz_t n_odd, bad_cof, good_cof, sbf, tmp; mpz_inits(n_odd, bad_cof, good_cof, sbf, tmp, NULL);
    int result = 0;
    qfaclist f; fl_init(&f);
    unsigned long v2 = mpz_scan1(n_in, 0);
    mpz_fdiv_q_2exp(n_odd, n_in, v2);
    if (mpz_fdiv_ui(n_odd, 4) == 3) goto done;
    rep_gcd(n_odd, *bad, bad_cof);
    if (mpz_cmp_ui(bad_cof, 1) != 0) { if (!mpz_root(sbf, bad_cof, 2)) goto done; }
    rep_gcd(n_odd, *good, good_cof);
    fl_push_ui(&f, 2, v2);
    if (mpz_cmp_ui(bad_cof, 1) != 0) fl_push(&f, sbf, 2);
    if (mpz_cmp_ui(good_cof, 1) != 0) factor_smooth(good_cof, &f);
    if (mpz_cmp_ui(n_odd, 1) != 0) {
        if (!mpz_probab_prime_p(n_odd, 25)) goto done;
        if (mpz_fdiv_ui(n_odd, 4) == 3) goto done;
        fl_push(&f, n_odd, 1);
    }
    if (two_squares_from_factors(b1, b2, &f)) {
        mpz_mul(tmp, b1, b1); mpz_addmul(tmp, b2, b2);
        result = (mpz_cmp(tmp, n_in) == 0);
    }
done:
    fl_clear(&f);
    mpz_clears(n_odd, bad_cof, good_cof, sbf, tmp, NULL);
    return result;
}

/* =========================== norm equation (qlapoti) =========================== */

#define MAX_ITERS 10000000UL

static void round_div(mpz_t out, const mpz_t a, const mpz_t b)
{
    mpz_t q, r, twor; mpz_inits(q, r, twor, NULL);
    mpz_fdiv_qr(q, r, a, b);
    mpz_mul_2exp(twor, r, 1);
    if (mpz_cmp(twor, b) >= 0) mpz_add_ui(q, q, 1);
    mpz_set(out, q);
    mpz_clears(q, r, twor, NULL);
}

typedef struct {
    mpz_t A1, A2, B1, B2, C1, C2, D1, D2, E1, E2, N1, N2, Nb1;
} neq_res;
static void neq_init(neq_res *r){ mpz_inits(r->A1,r->A2,r->B1,r->B2,r->C1,r->C2,r->D1,r->D2,r->E1,r->E2,r->N1,r->N2,r->Nb1,NULL); }
static void neq_clear(neq_res *r){ mpz_clears(r->A1,r->A2,r->B1,r->B2,r->C1,r->C2,r->D1,r->D2,r->E1,r->E2,r->N1,r->N2,r->Nb1,NULL); }

static int fix_coeff(const qo_ctx *ctx, const qelem *d1, const qelem *d2,
                     const mpz_t N, const mpz_t z, int *state,
                     mpz_t xC1, mpz_t xC2, mpz_t D1, mpz_t D2, mpz_t xE1, mpz_t xE2)
{
    qelem cj2, cj1, prod; qe_init(&cj2); qe_init(&cj1); qe_init(&prod);
    mpz_t t; mpz_init(t);
    qo_conj(&cj2, d2); qo_conj(&cj1, d1);
    qo_mul(ctx, &prod, d1, &cj2);
    mpz_divexact(D1, prod.a, N); mpz_divexact(D2, prod.b, N);
    mpz_set(xC1, cj2.a); mpz_set(xC2, cj2.b);
    mpz_set(xE1, cj1.a); mpz_set(xE2, cj1.b);
    int ret = 1; *state = 0;
    if (mpz_even_p(N)) { *state = 1; goto out; }
    unsigned long z8 = modu(z, 8);
    if (z8 == 1) {
        if (modu(D1,4)==2){ret=0;goto out;} mpz_add(t,D1,D2); if(modu(t,4)==2){ret=0;goto out;} *state=0;
    } else if (z8 == 5) {
        if (modu(D1,4)==0){ret=0;goto out;} mpz_add(t,D1,D2); if(modu(t,4)==0){ret=0;goto out;} *state=0;
    } else {
        int has2 = 0;
        mpz_add(t,N,D1); if(modu(t,4)==2)has2=1;
        mpz_add(t,N,D1); mpz_add(t,t,D2); if(modu(t,4)==2)has2=1;
        mpz_sub(t,xC1,xE1); if(modu(t,4)==2)has2=1;
        mpz_add(t,xC1,xC2); mpz_sub(t,t,xE1); mpz_sub(t,t,xE2); if(modu(t,4)==2)has2=1;
        if (!has2) { *state=1; ret=1; goto out; }
        has2 = 0;
        mpz_sub(t,D1,N); if(modu(t,4)==2)has2=1;
        mpz_sub(t,D1,N); mpz_add(t,t,D2); if(modu(t,4)==2)has2=1;
        mpz_add(t,xC1,xE1); if(modu(t,4)==2)has2=1;
        mpz_add(t,xC1,xC2); mpz_add(t,t,xE1); mpz_add(t,t,xE2); if(modu(t,4)==2)has2=1;
        if (!has2) { *state=2; ret=1; goto out; }
        ret = 0;
    }
out:
    qe_clear(&cj2); qe_clear(&cj1); qe_clear(&prod); mpz_clear(t);
    return ret;
}

static int qlapoti(const qo_ctx *ctx, const mpz_t N, const qelem *alpha_in,
                   unsigned long e, neq_res *result)
{
    if (mpz_fdiv_ui(N, 4) == 2) return 0;
    int success = 0;
    qelem alpha; qe_init(&alpha); qe_set(&alpha, alpha_in);
    mpz_t tr_alpha, n_alpha, tr_alpha_inv, r, M, max_k;
    mpz_inits(tr_alpha, n_alpha, tr_alpha_inv, r, M, max_k, NULL);
    qo_trace(tr_alpha, &alpha);
    if (mpz_sgn(tr_alpha) < 0) { mpz_neg(alpha.a, alpha.a); mpz_neg(alpha.b, alpha.b); mpz_neg(tr_alpha, tr_alpha); }
    qo_norm(ctx, n_alpha, &alpha);
    if (!mpz_invert(tr_alpha_inv, tr_alpha, N)) goto cleanup0;
    mpz_divexact(r, n_alpha, N);
    mpz_set_ui(M, 1); mpz_mul_2exp(M, M, e);
    mpz_mul_2exp(max_k, N, e - 2); mpz_fdiv_q(max_k, max_k, ctx->p); mpz_sqrt(max_k, max_k);

    mpz_t k, ksq, bound, v, t1, r_tli, two_r_tli, tmp, tmp2;
    mpz_inits(k, ksq, bound, v, t1, r_tli, two_r_tli, tmp, tmp2, NULL);
    mpz_set(k, max_k); mpz_mul(ksq, k, k);
    mpz_mul_2exp(tmp, M, 2); mpz_mul(tmp2, ksq, N); mpz_sub(tmp, tmp, tmp2);
    if (mpz_sgn(tmp) <= 0) mpz_set_ui(bound, 0);
    else { mpz_mul(tmp, tmp, max_k); mpz_mul(tmp, tmp, max_k); mpz_mul_2exp(tmp2, N, 2); mpz_fdiv_q(tmp, tmp, tmp2); mpz_sqrt(bound, tmp); }
    mpz_add_ui(tmp, ksq, 1); mpz_mul(tmp, tmp, r); mpz_sub(tmp, M, tmp);
    mpz_mul(tmp, tmp, tr_alpha_inv); mpz_neg(tmp, tmp); mpz_mod(v, tmp, N);
    mpz_mul(r_tli, r, tr_alpha_inv);
    mpz_mul_2exp(tmp, k, 1); mpz_sub_ui(tmp, tmp, 1); mpz_mul(tmp, tmp, r_tli); mpz_neg(tmp, tmp); mpz_mod(t1, tmp, N);
    mpz_mul_2exp(two_r_tli, r_tli, 1);
    mpz_mod(two_r_tli, two_r_tli, N);   /* reduce increment once; v,t1 stay in [0,N) via cond. subtract */
    int N_even = mpz_even_p(N);

    qelem delta_1, delta_2, kalpha; qe_init(&delta_1); qe_init(&delta_2); qe_init(&kalpha);
    mpz_t s, t, z, _z, d1, d2, b1, b2, xC1, xC2, D1, D2, xE1, xE2, xV, acc;
    mpz_inits(s, t, z, _z, d1, d2, b1, b2, xC1, xC2, D1, D2, xE1, xE2, xV, acc, NULL);

    for (unsigned long iter = 0; iter < MAX_ITERS; iter++) {
        mpz_sub_ui(k, k, 1);
        /* v, t1, two_r_tli are all in [0,N); each update lands in [0,2N), so a
         * single conditional subtraction replaces a full big-integer division. */
        mpz_add(v, v, t1);         if (mpz_cmp(v, N) >= 0)  mpz_sub(v, v, N);
        mpz_add(t1, t1, two_r_tli); if (mpz_cmp(t1, N) >= 0) mpz_sub(t1, t1, N);
        if (N_even && mpz_even_p(k)) continue;
        if (mpz_cmp(v, bound) >= 0) continue;
        mpz_mul(ksq, k, k);
        mpz_mul(tmp, k, v); mpz_neg(tmp, tmp); mpz_add_ui(tmp2, ksq, 1); round_div(xV, tmp, tmp2);
        mpz_mul(s, xV, k); mpz_neg(s, s); mpz_sub(s, s, v);
        mpz_set(t, xV);
        if (N_even) {
            if (mpz_even_p(k)) continue;
            if (modu(s,2)==modu(r,2)) continue;
            if (modu(t,2)==modu(r,2)) continue;
        }
        mpz_add_ui(tmp, ksq, 1); mpz_mul(tmp, tmp, r); mpz_sub(_z, M, tmp);
        mpz_mul(tmp, tr_alpha, v); mpz_add(_z, _z, tmp);
        mpz_fdiv_q(z, _z, N);
        mpz_mul(tmp, s, s); mpz_sub(z, z, tmp); mpz_mul(tmp, t, t); mpz_sub(z, z, tmp);
        if (mpz_sgn(z) < 0) continue;
        { unsigned long z4 = modu(z,4), z8 = modu(z,8); if (z4==3 || z8==6) continue; }
        mpz_mul(tmp, s, N); mpz_add(delta_1.a, alpha.a, tmp); mpz_set(delta_1.b, alpha.b);
        qo_norm(ctx, tmp, &delta_1);
        if (!mpz_divisible_p(tmp, N)) continue;
        mpz_divexact(d1, tmp, N);
        if (modu(z,4)==0) continue;
        mpz_mul(kalpha.a, alpha.a, k); mpz_mul(kalpha.b, alpha.b, k);
        mpz_mul(tmp, t, N); mpz_add(delta_2.a, kalpha.a, tmp); mpz_set(delta_2.b, kalpha.b);
        if (modu(z,4)==2) {
            mpz_add(tmp, delta_1.a, delta_2.a); mpz_add(tmp2, delta_1.b, delta_2.b); mpz_add(tmp, tmp, tmp2);
            if (mpz_even_p(tmp)) continue;
        }
        int state;
        if (!fix_coeff(ctx, &delta_1, &delta_2, N, z, &state, xC1, xC2, D1, D2, xE1, xE2)) continue;
        if (!sum_of_squares(b1, b2, z)) continue;

        unsigned long z4 = modu(z,4), z8 = modu(z,8);
        if (z4==1 && mpz_even_p(b1)) mpz_swap(b1, b2);
        if (z8==2) {
            if (modu(b1,4)==3) mpz_neg(b1, b1);
            if ((state==1 && modu(b2,4)==3) || (state==2 && modu(b2,4)==1)) mpz_neg(b2, b2);
        }
        qo_norm(ctx, tmp, &delta_2); mpz_divexact(d2, tmp, N);
        mpz_mul(acc, b1, b1); mpz_mul(acc, acc, N);
        mpz_set(result->Nb1, acc); mpz_add(result->N1, d1, acc);
        mpz_mul(acc, b2, b2); mpz_mul(acc, acc, N); mpz_add(result->N2, d2, acc);
        mpz_mul(result->A1, b1, delta_1.a); mpz_mul(result->A2, b1, delta_1.b);
        mpz_mul(result->B1, b1, b2); mpz_mul(result->B1, result->B1, N); mpz_set_ui(result->B2, 0);
        mpz_mul(result->C1, b1, xC1); mpz_mul(result->C2, b1, xC2);
        mpz_mul(result->E1, b2, xE1); mpz_mul(result->E2, b2, xE2);
        mpz_set(result->D1, D1); mpz_set(result->D2, D2);
        if (N_even) {
            mpz_set(result->Nb1, d1);
            mpz_swap(result->E1, result->C1); mpz_swap(result->E2, result->C2);
            mpz_swap(result->B1, result->D1); mpz_swap(result->B2, result->D2);
        }
        success = 1;
        break;
    }
    qe_clear(&delta_1); qe_clear(&delta_2); qe_clear(&kalpha);
    mpz_clears(s, t, z, _z, d1, d2, b1, b2, xC1, xC2, D1, D2, xE1, xE2, xV, acc, NULL);
    mpz_clears(k, ksq, bound, v, t1, r_tli, two_r_tli, tmp, tmp2, NULL);
cleanup0:
    qe_clear(&alpha);
    mpz_clears(tr_alpha, n_alpha, tr_alpha_inv, r, M, max_k, NULL);
    return success;
}

/* =========================== public interface =========================== */

int qt_sample_ideal(ibz_t *N_out, ibz_t *alpha_a, ibz_t *alpha_b, gmp_randstate_t rng)
{
    qo_ctx ctx; qo_ctx_init(&ctx, CHARACTERISTIC);
    mpz_t ub, ell, negp, negp_mod, r, inv2, lam, two, tmp, n0, n1;
    mpz_inits(ub, ell, negp, negp_mod, r, inv2, lam, two, tmp, n0, n1, NULL);
    int ok = 0;

    mpz_sqrt(ub, ctx.p); mpz_fdiv_q_2exp(ub, ub, 1);
    mpz_neg(negp, ctx.p);
    mpz_add_ui(tmp, ub, 1);
    mpz_urandomm(ell, rng, tmp);
    mpz_nextprime(ell, ell);
    while (mpz_kronecker(negp, ell) != 1) mpz_nextprime(ell, ell);

    mpz_mod(negp_mod, negp, ell);
    if (!sqrtmod(r, negp_mod, ell)) goto out;
    mpz_set_ui(two, 2);
    if (!mpz_invert(inv2, two, ell)) goto out;
    mpz_sub_ui(tmp, r, 1); mpz_mul(tmp, tmp, inv2); mpz_mod(lam, tmp, ell);

    mpz_add_ui(*alpha_a, lam, 1);
    mpz_set_si(*alpha_b, -1);

    { qelem al; qe_init(&al); mpz_set(al.a, *alpha_a); mpz_set(al.b, *alpha_b);
      qo_norm(&ctx, n0, &al);
      mpz_sub(al.a, *alpha_a, ell); qo_norm(&ctx, n1, &al);
      if (mpz_cmp(n1, n0) < 0) mpz_sub(*alpha_a, *alpha_a, ell);
      qe_clear(&al); }

    mpz_set(*N_out, ell);
    ok = 1;
out:
    mpz_clears(ub, ell, negp, negp_mod, r, inv2, lam, two, tmp, n0, n1, NULL);
    qo_ctx_clear(&ctx);
    return ok;
}

int qt_normeq(qt_output_t *sol, const ibz_t *N, const ibz_t *alpha_a, const ibz_t *alpha_b)
{
    unsigned long e_sol = (unsigned long)TORSION_EVEN_POWER - 3;
    qo_ctx ctx; qo_ctx_init(&ctx, CHARACTERISTIC);
    qelem alpha; qe_init(&alpha); mpz_set(alpha.a, *alpha_a); mpz_set(alpha.b, *alpha_b);

    neq_res res; neq_init(&res);
    int ok = qlapoti(&ctx, *N, &alpha, e_sol, &res);
    if (!ok) { neq_clear(&res); qe_clear(&alpha); qo_ctx_clear(&ctx); return 0; }

    /* sanity: N1 + N2 == 2^{e_sol} */
    mpz_t MOD, chk; mpz_inits(MOD, chk, NULL);
    mpz_set_ui(MOD, 1); mpz_mul_2exp(MOD, MOD, e_sol);
    mpz_add(chk, res.N1, res.N2);
    int id_ok = (mpz_cmp(chk, MOD) == 0);

    /* Reduce every field to [0, 2^{e_sol+2}) and export to digit arrays. */
    mpz_set_ui(MOD, 1); mpz_mul_2exp(MOD, MOD, e_sol + 2);
    mpz_t *fields[12] = { &res.N1, &res.Nb1, &res.A1, &res.A2, &res.B1, &res.B2,
                          &res.C1, &res.C2, &res.D1, &res.D2, &res.E1, &res.E2 };
    digit_t *outs[12] = { sol->N1, sol->Nb1, sol->A1, sol->A2, sol->B1, sol->B2,
                          sol->C1, sol->C2, sol->D1, sol->D2, sol->E1, sol->E2 };
    for (int i = 0; i < 12; i++) {
        mpz_mod(chk, *fields[i], MOD);
        memset(outs[i], 0, NWORDS_ORDER * sizeof(digit_t));
        ibz_to_digits(outs[i], (const ibz_t *)&chk);
    }

    sol->Nmod4 = (uint8_t)mpz_fdiv_ui(*N, 4);
    /* two_isog_choice = (a + b) mod 2 of the ideal generator alpha = a + b*w */
    mpz_add(chk, *alpha_a, *alpha_b);
    sol->two_isog_choice = (mpz_odd_p(chk) != 0);

    mpz_clears(MOD, chk, NULL);
    neq_clear(&res); qe_clear(&alpha); qo_ctx_clear(&ctx);
    return id_ok;
}

int qt_sample_and_normeq(qt_output_t *sol, gmp_randstate_t rng)
{
    ibz_t N, a, b; ibz_init(&N); ibz_init(&a); ibz_init(&b);
    int ok = 0;
    for (int tries = 0; tries < 64 && !ok; tries++) {
        if (!qt_sample_ideal(&N, &a, &b, rng)) continue;
        if (mpz_fdiv_ui(N, 4) == 2) continue;
        ok = qt_normeq(sol, &N, &a, &b);
    }
    ibz_finalize(&N); ibz_finalize(&a); ibz_finalize(&b);
    return ok;
}
