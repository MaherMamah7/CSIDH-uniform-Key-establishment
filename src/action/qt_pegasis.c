#include <qt_pegasis.h>
#include <inttypes.h>

#ifdef QTP_TIMING
#include <time.h>
double qtp_time_step2 = 0.0;   /* kernel points */
double qtp_time_step3 = 0.0;   /* isogeny chain + splitting */
static double qtp_now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}
#endif

static inline void
correct_A_coeff(ec_curve_t *E)
{
    // A-2 should not be a square
    // If it is the case, take the other Fp-isomorphic A'
    // If A = -2(a^4+b^4)/(a^4-b^4),
    // then A' = -2((a+b)^4+(a-b)^4)/((a+b)^4-(a-b)^4)
    // where a/b = ((A-2)/(A+2))^{1/4}

    fp_t t0, t1, t2;
    uint32_t sel;

    fp_set_small(&t2, 2);    // t2 = 2
    fp_sub(&t0, &E->A, &t2); // t0 = A-2
    sel = fp_is_square(&t0);
    fp_add(&t1, &E->A, &t2); // t1 = A+2
    fp_inv(&t1);             // t1 = 1/(A+2)
    fp_mul(&t0, &t0, &t1);   // t0 = (A-2)/(A+2)
    fp_ftrt(&t0, &t0);       // t0 = ((A-2)/(A+2))^{1/4} = a/b

    fp_set_one(&t2);       // t2 = 1
    fp_add(&t1, &t0, &t2); // t1 = a/b+1
    fp_sub(&t2, &t0, &t2); // t2 = a/b-1
    fp_sqr(&t1, &t1);      // t1 = (a/b+1)^2
    fp_sqr(&t1, &t1);      // t1 = (a/b+1)^4
    fp_sqr(&t2, &t2);      // t2 = (a/b-1)^2
    fp_sqr(&t2, &t2);      // t2 = (a/b-1)^4
    fp_add(&t0, &t1, &t2); // t0 = (a/b+1)^4 + (a/b-1)^4
    fp_sub(&t1, &t2, &t1); // t1 = (a/b-1)^4 - (a/b+1)^4
    fp_inv(&t1);           // t1 = 1/((a/b-1)^4 - (a/b+1)^4)
    fp_mul(&t0, &t0, &t1); // t0 = -((a/b+1)^4 + (a/b-1)^4)/((a/b+1)^4 - (a/b-1)^4)
    fp_add(&t0, &t0, &t0); // t0 = -2((a/b+1)^4 + (a/b-1)^4)/((a/b+1)^4 - (a/b-1)^4) = A'

    // If A-2 is a square, A = A'
    // else, A stays A
    fp_select(&E->A, &E->A, &t0, sel);

    E->is_ws_model_computed = false;
    E->is_A24_computed = false;

    // Check that either A-2 or A'-2 is not a square
    if (true) {
        fp_set_small(&t2, 2);
        fp_sub(&t1, &t0, &t2); // t1 = A'-2
        assert(~sel || ~fp_is_square(&t1));
    }
}

static inline void
two_isogenous_curve(ec_curve_t *E2, ec_curve_t *E, const bool two_isog_choice)
{
    // Implementation based on https://joostrenes.nl/publications/pqc-eprint.pdf,
    // Proposition 2 and Remark 6, pp. 10-11.
    ec_point_t T0, T1, Tm1, P0, P1;
    uint32_t sel1, sel2;

    two_torsion(&T0, &T1, &Tm1, E);
    sel1 = fp_is_zero(&T1.x);

    // If T1 == (0,0), P0 = T1 and P1 = Tm1
    // else, P0 = Tm1 and P1 = T1
    select_ec_point(&P0, &Tm1, &T1, sel1);
    select_ec_point(&P1, &T1, &Tm1, sel1);

    fp_t A0, A1, t0, t1;

    // Compute A0, whose kernel is <(0,0)> (P0 = either T1 or Tm1)
    fp_add(&A0, &E->A, &E->A); // A0 = 2*A
    fp_add(&t0, &P1.x, &P1.x); // t0 = 2*x1
    fp_add(&t0, &t0, &E->A);   // t0 = 2*x1+A = sqrt(A^2-4)
    fp_inv(&t0);               // t0 = 1/sqrt(A^2-4)
    fp_mul(&A0, &A0, &t0);     // A0 = 2*A/sqrt(A^2-4)
    fp_neg(&t1, &A0);          // t1 = -2*A/sqrt(A^2-4)

    // If sqrt(A^2-4) is a square, A0 = -2*A/sqrt(A^2-4)
    // else, A0 = 2*A/sqrt(A^2-4)
    sel2 = fp_is_square(&t0);
    fp_select(&A0, &A0, &t1, sel2);

    // Compute A1, whose kernel is not <(0,0)> (P1 = either Tm1 or T1)
    fp_set_small(&t0, 2);      // t0 = 2
    fp_add(&A1, &P1.x, &P1.x); // A1 = 2*x1
    fp_sqr(&A1, &A1);          // A1 = 4*x1^2
    fp_sub(&A1, &t0, &A1);     // A1 = 2-4*x1^2 = 2(1-2*x1^2)
    fp_neg(&t1, &A1);          // t1 = -2(1-2*x1^2)

    // If x1 is a square, A1 = 2(1-2*x1^2)
    // else, A1 = -2(1-2*x1^2)
    sel2 = fp_is_square(&P1.x);
    fp_select(&A1, &t1, &A1, sel2);

    // Choose which coefficient correpsonds to T1 and Tm1
    // A0 should correspond to T1 and A1 to Tm1
    fp_cswap(&A0, &A1, ~sel1);

    // If a = 1 mod 2, then two_isog_choice == True
    // and we act by the ideal is (2,omega-1) of kernel <T1>
    // else, two_isog_choice == False and
    // we act by the ideal is (2,omega) of kernel <Tm1>
    sel2 = -(uint32_t)two_isog_choice;

    fp_select(&E2->A, &A1, &A0, sel2);

    // Already done in correct_A_coeff
    // E2->is_ws_model_computed=false;
    // E2->is_A24_computed=false;

    correct_A_coeff(E2);
}

void
compute_action(ec_curve_t *Ea, ec_curve_t *Eabar, ec_curve_t *E, const norm_eq_sol_t *sol)
{
    ec_curve_t E2, domain_EC[2];
    jac_point_t P, Q, TP, TQ, OEstart, s5TP, s6TP, as5TQ, as6TQ, P16, Q16, maQ4, S, Ptmp;
    quad_point_t T[4], T16[4];
    first_isogenies_t first_isogenies;
    theta_struct_dim4_t codomain; // Ea x Eabar x A
    digit_t POW2[NWORDS_ORDER];
    digit_t alpha[NWORDS_ORDER], s1[NWORDS_ORDER], s2[NWORDS_ORDER], s3[NWORDS_ORDER], s4[NWORDS_ORDER],
        s5[NWORDS_ORDER], s6[NWORDS_ORDER], as1[NWORDS_ORDER], as2[NWORDS_ORDER], as3[NWORDS_ORDER], as4[NWORDS_ORDER],
        tmp[NWORDS_ORDER];
    digit_t N1mod16, s1mod16, s2mod16, s3mod16, s4mod16, as1mod16, as2mod16, as3mod16, as4mod16;
    uint32_t sel, ctls5, ctls6;
    int kbits;
    uint8_t N1mod4, Nb1mod4, A1mod4, A2mod4, position_codomain,
        position_theta_structure, position_theta_change_splitting;
    basis_components_16_t compo;

#ifdef QTP_TIMING
    double _qtp_t0 = qtp_now();
#endif

    two_isogenous_curve(&E2, E, sol->two_isog_choice);

    // If N = 2 mod 4, choose the 2-isogenous curve E2
    sel = -(uint32_t)(sol->Nmod4 == 2);
    select_curve(&domain_EC[0], E, &E2, sel);
    ec_twist(&domain_EC[1], &domain_EC[0]);

    ec_compute_A24(&domain_EC[0]);
    ec_compute_A24(&domain_EC[1]);

    // Generating Fp-rational basis of domain_EC[0][2^f]
    jac_rational_basis(&P, &Q, &TP, &TQ, &domain_EC[1], &domain_EC[0]);

    // Computation of 16-torsion and 4-torsion basis induced by (P,Q)
    jac_dbl_iter(&P16, &P, TORSION_EVEN_POWER - 5, &domain_EC[0]);
    jac_dbl_iter(&Q16, &Q, TORSION_EVEN_POWER - 5, &domain_EC[1]);

    // init 2^{e}
    POW2[0] = 1;
    for (uint8_t i = 1; i < NWORDS_ORDER; i++) {
        POW2[i] = 0;
    }
    multiple_mp_shiftl(POW2, TORSION_EVEN_POWER - 1, NWORDS_ORDER);

    /* Computing the 2^{e+2}-torsion above the kernel */

    // alpha = 1/N1 mod 2^{e+2}
    mp_inv_2e(alpha, sol->N1, TORSION_EVEN_POWER - 1, NWORDS_ORDER);

    // s1 = C1-E1 mod 2^{e+2}
    mp_sub(tmp, POW2, sol->E1, NWORDS_ORDER);
    mp_add(s1, sol->C1, tmp, NWORDS_ORDER);
    mp_mod_2exp(s1, TORSION_EVEN_POWER - 1, NWORDS_ORDER);
    mp_add(s1, s1, POW2, NWORDS_ORDER); // constant time

    // s5 = C2 - E2
    mp_sub(tmp, POW2, sol->E2, NWORDS_ORDER);
    mp_add(s5, sol->C2, tmp, NWORDS_ORDER); // constant time

    // s2 = C1 + C2 - E1 - E2 mod 2^{e+2} = s1 + s5 mod 2^{e+2}
    mp_add(s2, s1, s5, NWORDS_ORDER);
    mp_mod_2exp(s2, TORSION_EVEN_POWER - 1, NWORDS_ORDER);
    mp_add(s2, s2, POW2, NWORDS_ORDER); // constant time

    // s4 = B1 + D1 mod 2^{e+2}
    mp_add(s4, sol->B1, sol->D1, NWORDS_ORDER);
    mp_mod_2exp(s4, TORSION_EVEN_POWER - 1, NWORDS_ORDER);
    mp_add(s4, s4, POW2, NWORDS_ORDER); // constant time

    // s6 = B2 + D2
    mp_add(s6, sol->B2, sol->D2, NWORDS_ORDER);

    // s3 = B1 + B2 + D1 + D2 mod 2^{e+2} = s4 + s6 mod 2^{e+2}
    mp_add(s3, s4, s6, NWORDS_ORDER);
    mp_mod_2exp(s3, TORSION_EVEN_POWER - 1, NWORDS_ORDER);
    mp_add(s3, s3, POW2, NWORDS_ORDER); // constant time

    // as1 = alpha*s1 mod 2^{e+2}
    mp_mul(as1, alpha, s1, NWORDS_ORDER);
    mp_mod_2exp(as1, TORSION_EVEN_POWER - 1, NWORDS_ORDER);
    mp_add(as1, as1, POW2, NWORDS_ORDER); // constant time

    // as2 = alpha*s2 mod 2^{e+2}
    mp_mul(as2, alpha, s2, NWORDS_ORDER);
    mp_mod_2exp(as2, TORSION_EVEN_POWER - 1, NWORDS_ORDER);
    mp_add(as2, as2, POW2, NWORDS_ORDER); // constant time

    // as3 = alpha*s3 mod 2^{e+2}
    mp_mul(as3, alpha, s3, NWORDS_ORDER);
    mp_mod_2exp(as3, TORSION_EVEN_POWER - 1, NWORDS_ORDER);
    mp_add(as3, as3, POW2, NWORDS_ORDER); // constant time

    // as4 = alpha*s4 mod 2^{e+2}
    mp_mul(as4, alpha, s4, NWORDS_ORDER);
    mp_mod_2exp(as4, TORSION_EVEN_POWER - 1, NWORDS_ORDER);
    mp_add(as4, as4, POW2, NWORDS_ORDER); // constant time

    ctls5 = -(uint32_t)(s5[0] & 1);
    ctls6 = -(uint32_t)(s6[0] & 1);

    jac_point_init(&OEstart);

    select_jac_point(&s5TP, &OEstart, &TP, ctls5);
    select_jac_point(&s6TP, &OEstart, &TP, ctls6);
    select_jac_point(&as5TQ, &OEstart, &TQ, ctls5);
    select_jac_point(&as6TQ, &OEstart, &TQ, ctls6);

    // if alpha == 1 mod 4, then maQ4 = -Q4
    // else, alpha == 3 mod 4 and maQ4 = Q4
    sel = -(uint32_t)((alpha[0] & 3) == 1);
    jac_dbl_iter(&maQ4, &Q16, 2, &domain_EC[1]);
    jac_neg(&S, &maQ4);
    select_jac_point(&maQ4, &maQ4, &S, sel);

    // T1 = (N1*P, 0, s3*P+s6TP,-s1*P+s5TP)
    kbits = n_bits(sol->N1, NWORDS_ORDER);
    jac_MUL(&T[0][0], &P, sol->N1, kbits, &domain_EC[0]);

    copy_jac_point(&T[0][1], &OEstart);

    kbits = n_bits(s3, NWORDS_ORDER);
    jac_MUL(&Ptmp, &P, s3, kbits, &domain_EC[0]);
    ADD(&T[0][2], &Ptmp, &s6TP, &domain_EC[0]);

    kbits = n_bits(s1, NWORDS_ORDER);
    jac_MUL(&T[0][3], &P, s1, kbits, &domain_EC[0]);
    jac_neg(&Ptmp, &T[0][3]);
    ADD(&T[0][3], &Ptmp, &s5TP, &domain_EC[0]);

    // T2 = (0, N1*P, s2*P+s5TP,s4*P-s6TP)
    copy_jac_point(&T[1][0], &OEstart);
    copy_jac_point(&T[1][1], &T[0][0]);

    kbits = n_bits(s2, NWORDS_ORDER);
    jac_MUL(&Ptmp, &P, s2, kbits, &domain_EC[0]);
    ADD(&T[1][2], &Ptmp, &s5TP, &domain_EC[0]);

    kbits = n_bits(s4, NWORDS_ORDER);
    jac_MUL(&Ptmp, &P, s4, kbits, &domain_EC[0]);
    ADD(&T[1][3], &Ptmp, &s6TP, &domain_EC[0]); // TP = -TP

    // T3 = ([1-alpha*2^e]Q,0,as4*Q+as6TQ,-as2*Q+as5TQ)
    // with [-alpha*2^e]Q = maQ4
    ADD(&T[2][0], &Q, &maQ4, &domain_EC[1]);

    copy_jac_point(&T[2][1], &OEstart);

    kbits = n_bits(as4, NWORDS_ORDER);
    jac_MUL(&Ptmp, &Q, as4, kbits, &domain_EC[1]);
    ADD(&T[2][2], &Ptmp, &as6TQ, &domain_EC[1]);

    kbits = n_bits(as2, NWORDS_ORDER);
    jac_MUL(&T[2][3], &Q, as2, kbits, &domain_EC[1]);
    jac_neg(&Ptmp, &T[2][3]);
    ADD(&T[2][3], &Ptmp, &as5TQ, &domain_EC[1]);

    // T4 = (0,[1-alpha*2^e]Q,as1*Q+as5TQ,as3*Q-as6TQ)
    copy_jac_point(&T[3][0], &OEstart);
    copy_jac_point(&T[3][1], &T[2][0]);

    kbits = n_bits(as1, NWORDS_ORDER);
    jac_MUL(&Ptmp, &Q, as1, kbits, &domain_EC[1]);
    ADD(&T[3][2], &Ptmp, &as5TQ, &domain_EC[1]);

    kbits = n_bits(as3, NWORDS_ORDER);
    jac_MUL(&Ptmp, &Q, as3, kbits, &domain_EC[1]);
    ADD(&T[3][3], &Ptmp, &as6TQ, &domain_EC[1]); // TQ = -TQ

    /* Computing the 16-torsion above the kernel of f2*f1 */

    N1mod16 = ((sol->N1)[0] & 15) ^ 16;
    s1mod16 = (s1[0] & 15) ^ 16;
    s2mod16 = (s2[0] & 15) ^ 16;
    s3mod16 = (s3[0] & 15) ^ 16;
    s4mod16 = (s4[0] & 15) ^ 16;
    as1mod16 = (as1[0] & 15) ^ 16;
    as2mod16 = (as2[0] & 15) ^ 16;
    as3mod16 = (as3[0] & 15) ^ 16;
    as4mod16 = (as4[0] & 15) ^ 16;

    // [2^{e-2}]T1 = (N1*P16, 0, s3*P16,-s1*P16)
    kbits = n_bits(&N1mod16, 1);

    // printf("k_bits is %d\n",kbits);
    jac_MUL(&T16[0][0], &P16, &N1mod16, kbits, &domain_EC[0]);

    copy_jac_point(&T16[0][1], &OEstart);

    kbits = n_bits(&s3mod16, 1);
    jac_MUL(&T16[0][2], &P16, &s3mod16, kbits, &domain_EC[0]);

    kbits = n_bits(&s1mod16, 1);
    jac_MUL(&T16[0][3], &P16, &s1mod16, kbits, &domain_EC[0]);
    jac_neg(&T16[0][3], &T16[0][3]);

    // [2^{e-2}]T2 = (0, N1*P16, s2*P16, s4*P16)
    copy_jac_point(&T16[1][0], &OEstart);
    copy_jac_point(&T16[1][1], &T16[0][0]);

    kbits = n_bits(&s2mod16, 1);
    jac_MUL(&T16[1][2], &P16, &s2mod16, kbits, &domain_EC[0]);

    kbits = n_bits(&s4mod16, 1);
    jac_MUL(&T16[1][3], &P16, &s4mod16, kbits, &domain_EC[0]);

    // [2^{e-2}]T3 = (Q16,0,as4*Q16,-as2*Q16)
    copy_jac_point(&T16[2][0], &Q16);
    copy_jac_point(&T16[2][1], &OEstart);

    kbits = n_bits(&as4mod16, 1);
    jac_MUL(&T16[2][2], &Q16, &as4mod16, kbits, &domain_EC[1]);

    kbits = n_bits(&as2mod16, 1);
    jac_MUL(&T16[2][3], &Q16, &as2mod16, kbits, &domain_EC[1]);
    jac_neg(&T16[2][3], &T16[2][3]);

    // [2^{e-2}]T4 = (0,Q16,as1*Q16,as3*Q16)
    copy_jac_point(&T16[3][0], &OEstart);
    copy_jac_point(&T16[3][1], &T16[2][0]);

    kbits = n_bits(&as1mod16, 1);
    jac_MUL(&T16[3][2], &Q16, &as1mod16, kbits, &domain_EC[1]);

    kbits = n_bits(&as3mod16, 1);
    jac_MUL(&T16[3][3], &Q16, &as3mod16, kbits, &domain_EC[1]);
#ifndef NDEBUG
    printf("End kernel\n");
#endif
    /* Computation of f2*f1 */
    N1mod4 = (sol->N1)[0] & 3;
    Nb1mod4 = (sol->Nb1)[0] & 3;
    A1mod4 = (sol->A1)[0] & 3;
    A2mod4 = (sol->A2)[0] & 3;
    compo.N = N1mod16;
    compo.s1 = s1mod16;
    compo.s2 = s2mod16;
    compo.s3 = s3mod16;
    compo.s4 = s4mod16;

#ifndef NDEBUG

    quad_point_t T16_2[4];

    quad_point_dbl_iter(&T16_2[0], &T[0], TORSION_EVEN_POWER - 5, &domain_EC[0]);
    quad_point_dbl_iter(&T16_2[1], &T[1], TORSION_EVEN_POWER - 5, &domain_EC[0]);
    quad_point_dbl_iter(&T16_2[2], &T[2], TORSION_EVEN_POWER - 5, &domain_EC[1]);
    quad_point_dbl_iter(&T16_2[3], &T[3], TORSION_EVEN_POWER - 5, &domain_EC[1]);

    uint32_t equ = -1;
    for (uint8_t q = 0; q < 4; q++) {
        for (uint8_t j = 0; j < 4; j++) {
            equ &= jac_is_equal(&T16_2[q][j], &T16[q][j]);
        }
    }
    printf("Both basis are coherent %d\n", equ);


    printf("N = %d, s1 = %d, s2 = %d, s3 = %d, s4 = %d \n",
           compo.N & 15,
           compo.s1 & 15,
           compo.s2 & 15,
           compo.s3 & 15,
           compo.s4 & 15);

#endif
    // Bug on entries
    // SOL[4] leads to position_codomain and position_theta_structure = 255
    position_codomain = compute_codomain_position(compo);
    position_theta_structure = compute_theta_position(compo);

#ifndef NDEBUG
    printf("position_codomain = %i\n", position_codomain);
    printf("position_theta_structure = %i\n", position_theta_structure);
    printf("After position\n");
#endif

#ifdef QTP_TIMING
    double _qtp_t1 = qtp_now();
    qtp_time_step2 += _qtp_t1 - _qtp_t0;
#endif

    compute_first_isogenies(&first_isogenies, domain_EC, T16, &P16, &Q16, position_codomain, position_theta_structure);

#ifndef NDEBUG
    printf("End first isogenies\n");
#endif
    /* Isogeny chain computation */

    compute_chain_codomain(&codomain, T, &first_isogenies);

#ifndef NDEBUG
    printf("End chain\n");
#endif
    /* Splitting */
    position_theta_change_splitting = compute_theta_position_splitting(N1mod4, Nb1mod4, A1mod4, A2mod4);

#ifndef NDEBUG
    printf("position_splitting = %i\n", position_theta_change_splitting);
#endif

    codomain_to_curves(Ea, Eabar, &codomain, position_theta_change_splitting);

#ifndef NDEBUG
    printf("End splitting\n");
#endif

    // If N is even, we have acted with the conjugate
    sel = -(uint32_t)(((sol->Nmod4) & 1) == 0);
    cswap_curve(Ea, Eabar, sel);

#ifdef QTP_TIMING
    qtp_time_step3 += qtp_now() - _qtp_t1;
#endif
}
