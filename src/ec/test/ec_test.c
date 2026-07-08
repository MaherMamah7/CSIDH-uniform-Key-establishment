#include <assert.h>
#include <stdio.h>
#include <inttypes.h>

#include "test_extras.h"
#include <ec.h>
//#include <isog.h>
#include <rng.h>
#include <bench_test_arguments.h>

/******************************
Test functions
******************************/

void
jac_normalize(jac_point_t *P){
    fp_t t0, t1;
    
    fp_copy(&t0,&P->z); // t0=z
    fp_inv(&t0);// t0=1/z
    fp_sqr(&t1,&t0);// t1=1/z^2
    fp_mul(&P->x,&P->x,&t1);// x=x/z^2
    fp_mul(&t1,&t1,&t0);// t1=1/z^3
    fp_mul(&P->y,&P->y,&t1);// y=y/z^3
    fp_set_one(&P->z);// z=1
}

void
ADD_test(jac_point_t *R, const jac_point_t *P, const jac_point_t *Q, const ec_curve_t *E){
    fp_t t0, t1, t2, t3, t4, t5, t6, t7, t8, dx, dy;

    fp_sqr(&t0,&P->z); // t0 = z1^2
    fp_sqr(&t1,&Q->z); // t1 = z2^2
    fp_mul(&t2,&P->x,&t1); // t2 = x1*z2^2
    fp_mul(&t3,&Q->x,&t0); // t3 = x2*z1^2
    fp_sub(&dx,&t3,&t2); // dx = x2*z1^2 - x1*z2^2

    fp_mul(&t0,&t0,&P->z); // t0 = z1^3
    fp_mul(&t1,&t1,&Q->z); // t1 = z2^3
    fp_mul(&t4,&P->y,&t1); // t4 = y1*z2^3
    fp_mul(&t5,&Q->y,&t0); // t5 = y2*z1^3
    fp_sub(&dy,&t5,&t4); // dy = y2*z1^3 - y1*z2^3

    fp_mul(&t6, &P->z, &Q->z); // t6 = z1*z2
    fp_mul(&R->z, &t6, &dx); // z3 = z1*z2*dx

    fp_sqr(&t6,&t6); // t6 = (z1*z2)^2
    fp_mul(&t6,&E->A,&t6); // t6 = A*(z1*z2)^2
    fp_add(&t6,&t6,&t2); // t6 = A*(z1*z2)^2 + x1*z2^2
    fp_add(&t6,&t6,&t3); // t6 = A*(z1*z2)^2 + x1*z2^2 + x2*z1^2
    fp_sqr(&t7,&dx); // t7 = (dx)^2
    fp_sqr(&t8,&dy); // t8 = (dy)^2
    fp_mul(&t6,&t7,&t6); // t6 = (dx)^2*(A*(z1*z2)^2 + x1*z2^2 + x2*z1^2)
    fp_sub(&R->x,&t8,&t6); // x3 = (dy)^2 - (dx)^2*(A*(z1*z2)^2 + x1*z2^2 + x2*z1^2)

    fp_mul(&t6,&t2,&t7); // t6 = x1*z2^2*(dx)^2
    fp_sub(&t6,&t6,&R->x); // t6 = x1*z2^2*(dx)^2 - x3
    fp_mul(&t6,&t6,&dy); // t6 = dy*(x1*z2^2*(dx)^2 - x3)
    fp_mul(&t7,&t7,&dx); // t7 = (dx)^3
    fp_mul(&t8,&t4,&t7); // t8 = y1*z2^3*(dx)^3
    fp_sub(&R->y, &t6, &t8); // y3 = dy*(x1*z2^2*(dx)^2 - x3) - y1*z2^3*(dx)^3
}

uint32_t jac_is_on_curve(const jac_point_t *P, const ec_curve_t *E){
    fp_t t0, t1, t2, x, y2;

    fp_copy(&t0,&P->z);
    fp_inv(&t0); // t0 =1/Z
    fp_sqr(&t1,&t0); // t1 = 1/Z^2
    fp_mul(&x,&P->x,&t1); // x = X/Z^2
    fp_mul(&t1, &t1, &t0); // t1 = 1/Z^3
    fp_mul(&y2,&P->y,&t1); // y = Y/Z^3
    fp_sqr(&y2,&y2); // y2 = y^2

    fp_sqr(&t0,&x); // t0 = x^2
    fp_mul(&t1, &E->A, &x); // t1 = A*x
    fp_set_one(&t2); // t2 = 1
    fp_add(&t2,&t2,&t1); //t2 = A*x + 1
    fp_add(&t2,&t2,&t0); // t2 = x^2 + A*x +1
    fp_mul(&t2,&t2,&x); // t2 = x*(x^2 + A*x +1)

    return fp_is_equal(&t2,&y2);
}

uint32_t ws_is_on_curve(const jac_ws_point_t *P, const ec_curve_t *E){
    fp_t t0, t1, t2, x, y2, b;

    fp_copy(&t0,&P->z);
    fp_inv(&t0); // t0 =1/Z
    fp_sqr(&t1,&t0); // t1 = 1/Z^2
    fp_mul(&x,&P->x,&t1); // x = X/Z^2
    fp_mul(&t1, &t1, &t0); // t1 = 1/Z^3
    fp_mul(&y2,&P->y,&t1); // y = Y/Z^3
    fp_sqr(&y2,&y2); // y2 = y^2

    fp_set_one(&t0); // t0 = 1
    fp_sqr(&t1,&E->ao3); // t1 = A^2/9
    fp_add(&t1,&t1,&t1); // t1 = 2*A^2/9
    fp_sub(&t1,&t1,&t0); // t1 = 2*A^2/9-1
    fp_mul(&b,&E->ao3,&t1); // b = A/3*(2*A^2/9-1)

    //fp_t A, a, one;
    //fp_add(&A,&E->ao3,&E->ao3);
    //fp_add(&A,&A,&E->ao3); 
    //printf("A == 3*A/3: %u\n",fp_is_equal(&A,&E->A));

    //fp_set_one(&one);
    //fp_mul(&a,&E->ao3,&E->A);
    //fp_sub(&a,&one,&a);
    //printf("a == 1-A^2/3: %u\n",fp_is_equal(&a,&E->a));

    fp_sqr(&t0,&x); // t0 = x^2
    fp_mul(&t0,&t0,&x); // t0 = x^3
    fp_mul(&t1, &E->a, &x); // t1 = a*x
    fp_add(&t2,&t1,&b); //t2 = a*x + b
    fp_add(&t2,&t0,&t2); // t2 = x^3 + a*x + b

    return fp_is_equal(&t2,&y2)|(fp_is_zero(&P->x)&fp_is_zero(&P->z));
}

int
test_jacobian(ec_curve_t *curve, unsigned int Ntest)
{
    unsigned int i;

    ec_point_t P, Q, X, P1;
    jac_point_t R, S, T, U, jac_zero;
    jac_ws_point_t Rw;

    jac_point_init(&jac_zero);

    for (i = 0; i < Ntest; i++) {
        ec_random_test(&P, curve);
        ec_random_test(&Q, curve);
        ec_random_test(&X, curve);

        // Convert to Jacobian coordinates.
        ec_to_jac(&R, &P);
        ec_to_jac(&S, &Q);
        ec_to_jac(&T, &X);

        // Convert back to (X:Y:Z)
        jac_to_ec(&P1,&R);
        if(!ec_is_equal(&P1,&P)){
            printf("Failed Montgomery (X:Y:Z) <--> Jacobian conversion\n");
            return 1;
        }

        ADD(&R, &jac_zero, &jac_zero, curve);
        if (!jac_is_equal(&R, &jac_zero)) {
            printf("Failed 0 + 0 = 0 in jac\n");
            return 1;
        }

        DBL(&R, &jac_zero, curve);
        if (!jac_is_equal(&R, &jac_zero)) {
            printf("Failed 2*0 = 0 in jac\n");
            return 1;
        }

        jac_neg(&R, &S);
        ADD(&R, &S, &R, curve);
        if (!jac_is_equal(&R, &jac_zero)) {
            printf("Failed P - P = 0 in jac\n");
            return 1;
        }

        ADD(&R, &S, &jac_zero, curve);
        if (!jac_is_equal(&R, &S)) {
            printf("Failed P + 0 = P in jac\n");
            return 1;
        }

        ADD(&R, &jac_zero, &S, curve);
        if (!jac_is_equal(&R, &S)) {
            printf("Failed P + 0 = P in jac\n");
            return 1;
        }

        ADD(&R, &S, &jac_zero, curve);
        if (!jac_is_equal(&R, &S)) {
            printf("Failed 0 + P = P in jac\n");
            return 1;
        }

        DBL(&R, &S, curve);
        ADD(&U, &S, &S, curve);
        if (!jac_is_equal(&R, &U)) {
            printf("Failed P + P = 2*P in jac\n");
            return 1;
        }

        ADD(&R, &T, &S, curve);
        ADD(&T, &S, &T, curve);
        if (!jac_is_equal(&R, &T)) {
            printf("Failed P + Q = Q + P in jac\n");
            return 1;
        }

        ADD(&R, &T, &S, curve);
        ADD(&U, &S, &T, curve);
        if (!jac_is_equal(&R, &U)) {
            printf("Failed P + Q = Q + P in jac\n");
            return 1;
        }

        // Double R to make it different than (T + S).
        DBL(&R, &R, curve);
        ADD(&U, &S, &T, curve);
        ADD(&U, &U, &R, curve);
        ADD(&R, &R, &T, curve);
        ADD(&R, &R, &S, curve);
        if (!jac_is_equal(&R, &U)) {
            printf("Failed (P + Q) + R = P + (Q + R) in jac\n");
            return 1;
        }

        ADD(&R, &S, &T, curve);
        DBL(&R, &R, curve);
        DBL(&S, &S, curve);
        DBL(&T, &T, curve);
        ADD(&U, &S, &T, curve);
        if (!jac_is_equal(&R, &U)) {
            printf("Failed 2*(P + Q) = 2*P + 2*Q in jac\n");
            return 1;
        }

        
        jac_to_ws(&Rw, &jac_zero, curve);
        ws_to_jac(&R, &Rw, curve);
        if (!jac_is_equal(&R, &jac_zero)) {
            printf("Failed converting to Weierstrass 0\n");
            return 1;
        }


        jac_to_ws(&Rw, &S, curve);
        ws_to_jac(&R, &Rw, curve);
        if (!jac_is_equal(&S, &R)) {
            printf("Failed converting to Weierstrass S\n");
            return 1;
        }

        
        DBL(&S, &S, curve);
        jac_to_ws(&Rw, &S, curve);
        ws_to_jac(&R, &Rw, curve);
        if (!jac_is_equal(&S, &R)) {
            printf("Failed converting to Weierstrass 2*S\n");
            return 1;
        }

        
        jac_to_ws(&Rw, &jac_zero, curve);
        DBL_ws(&Rw, &Rw);
        ws_to_jac(&R, &Rw, curve);
        if (!jac_is_equal(&R, &jac_zero)) {
            printf("Failed 2*0 = 0 in Weierstrass\n");
            return 1;
        }

        jac_to_ws(&Rw, &S, curve);
        DBL_ws(&Rw, &Rw);
        ws_to_jac(&R, &Rw, curve);
        DBL(&S, &S, curve);
        if (!jac_is_equal(&S, &R)) {
            printf("Failed doubling in Weierstrass\n");
            return 1;
        }
    }
    printf("Jacobian tests....................................................... PASSED\n");
    return 0;
}

int
test_xDBL_xADD(ec_curve_t *curve, unsigned int Ntest)
{
    ec_point_t P, Q;
    jac_point_t jP, jQ, jPQ;
    xz_point_t xP, xQ, xPQ, xR1, xR2;
    for (unsigned int i = 0; i < Ntest; i++) {
        ec_random_test(&P, curve);
        ec_random_test(&Q, curve);

        ec_to_jac(&jP, &P);
        ec_to_jac(&jQ, &Q);
        jac_neg(&jQ,&jQ);

        ADD(&jPQ, &jP, &jQ, curve);
        jac_neg(&jQ,&jQ);

        ec_to_xz(&xP, &P);
        ec_to_xz(&xQ, &Q);
        jac_to_xz(&xPQ, &jPQ);

        // 2(P + Q) = 2P + 2Q
        xADD(&xR1, &xP, &xQ, &xPQ);
        xDBL(&xR1, &xR1, curve);
        xDBL(&xP, &xP, curve);
        xDBL(&xQ, &xQ, curve);
        xDBL(&xPQ, &xPQ, curve);
        xADD(&xR2, &xP, &xQ, &xPQ);
        if (!xz_is_equal(&xR1, &xR2)) {
            printf("Failed 2(P + Q) = 2P + 2Q\n");
            return 1;
        }

        // (P+Q) + (P-Q) = 2P
        xADD(&xR1, &xP, &xQ, &xPQ);
        xDBL(&xQ, &xQ, curve);
        xADD(&xR1, &xR1, &xPQ, &xQ);
        xDBL(&xP, &xP, curve);
        if (!xz_is_equal(&xR1, &xP)) {
            printf("Failed (P+Q) + (P-Q) = 2P\n");
            return 1;
        }
    }
    printf("xDBL and xADD tests.................................................. PASSED\n");
    return 0;
}

int
test_xDBLADD(ec_curve_t *curve, unsigned int Ntest)
{
    ec_point_t P, Q;
    jac_point_t jP, jQ, jPQ;
    xz_point_t xP, xQ, xPQ, xR1, xR2;

    for (unsigned int i = 0; i < Ntest; i++) {
        ec_random_test(&P, curve);
        ec_random_test(&Q, curve);

        ec_to_jac(&jP, &P);
        ec_to_jac(&jQ, &Q);
        jac_neg(&jQ,&jQ);

        ADD(&jPQ, &jP, &jQ, curve);
        jac_neg(&jQ,&jQ);

        ec_to_xz(&xP, &P);
        ec_to_xz(&xQ, &Q);
        jac_to_xz(&xPQ, &jPQ);

        xDBLADD(&xR1, &xR2, &xP, &xQ, &xPQ, curve);
        xADD(&xPQ, &xP, &xQ, &xPQ);
        if (!xz_is_equal(&xR2, &xPQ)) {
            printf("Failed addition in xDBLADD\n");
            return 1;
        }
        xDBL(&xP, &xP, curve);
        if (!xz_is_equal(&xR1, &xP)) {
            printf("Failed doubling in xDBLADD\n");
            return 1;
        }
    }
    printf("xDBLADD tests........................................................ PASSED\n");
    return 0;
}

int
test_zero_identities(ec_curve_t *curve, unsigned int Ntest)
{
    unsigned int i;

    xz_point_t P, Q, R, ec_zero;

    xz_point_init(&ec_zero);

    assert(xz_is_zero(&ec_zero));

    for (i = 0; i < Ntest; i++) {
        xz_random_test(&P, curve);

        xADD(&R, &ec_zero, &ec_zero, &ec_zero);
        if (!xz_is_zero(&R)) {
            printf("Failed 0 + 0 = 0\n");
            return 1;
        }

        xDBL(&R, &P, curve);
        xADD(&R, &P, &P, &R);
        if (!xz_is_zero(&R)) {
            printf("Failed P - P = 0\n");
            return 1;
        }

        xDBL(&R, &ec_zero, curve);
        if (!xz_is_zero(&R)) {
            printf("Failed 2*0 = 0\n");
            return 1;
        }

        xADD(&R, &P, &ec_zero, &P);
        if (!xz_is_equal(&R, &P)) {
            printf("Failed P + 0 = P\n");
            return 1;
        }
        xADD(&R, &ec_zero, &P, &P);
        if (!xz_is_equal(&R, &P)) {
            printf("Failed P + 0 = P\n");
            return 1;
        }

        xDBLADD(&R, &Q, &P, &ec_zero, &P, curve);
        if (!xz_is_equal(&Q, &P)) {
            printf("Failed P + 0 = P in xDBLADD\n");
            return 1;
        }
        xDBLADD(&R, &Q, &ec_zero, &P, &P, curve);
        if (!xz_is_equal(&Q, &P)) {
            printf("Failed P + 0 = P in xDBLADD\n");
            return 1;
        }
        if (!xz_is_zero(&R)) {
            printf("Failed 2*0 = 0 in xDBLADD\n");
            return 1;
        }
    }
    printf("xz zero identities................................................... PASSED\n");
    return 0;
}

int
test_mul(ec_curve_t *curve, unsigned int Ntest){
    ibz_t one_ibz, pp1_ibz;
    digit_t pp1[NWORDS_ORDER], c;

    ibz_init(&one_ibz);
    ibz_init(&pp1_ibz);
    ibz_set(&one_ibz, 1);
    ibz_add(&pp1_ibz, &CHARACTERISTIC, &one_ibz);
    ibz_to_digits(pp1, &pp1_ibz);
    ibz_to_digits(&c, &TORSION_ODD);
    ibz_finalize(&one_ibz);
    ibz_finalize(&pp1_ibz);

    int kbits = n_bits(pp1, NWORDS_ORDER);
    int kbits_c = n_bits(&c, 1);
    //printf("kbits = %i\n",kbits);

    xz_point_t xP, xQ, xR;
    ec_point_t P, Q;
    jac_point_t jP, jQ;

    for (unsigned int i = 0; i < Ntest; i++) {
        xz_random_test(&xP, curve);
        ec_random_test(&P, curve);

        xMUL(&xQ, &xR, &xP, pp1, kbits, curve);

        if(!xz_is_zero(&xQ)){
            printf("Failed (p+1)*P = 0 in xMUL\n");
            return 1;
        }
        if(!xz_is_equal(&xR,&xP)){
            printf("Failed (p+2)*P = P in xMUL\n");
            return 1;
        }

        ec_MUL(&Q, &P, pp1, kbits, curve);
        if(!ec_is_zero(&Q)){
            printf("Failed (p+1)*P = 0 in ec_MUL\n");
            return 1;
        }

        ec_MUL(&Q, &P, &c, kbits_c, curve);
        ec_to_jac(&jQ,&Q);
        jac_dbl_iter(&jQ, &jQ, TORSION_EVEN_POWER, curve);
        if(!jac_is_zero(&jQ)){
            printf("Failed (p+1)*P = 0 in ec_MUL and jac_dbl_iter\n");
            return 1;
        }

        ec_to_jac(&jP,&P);

        jac_MUL(&jQ,&jP,pp1,kbits,curve);
        if(!jac_is_zero(&jQ)){
            printf("Failed (p+1)*P = 0 in jac_MUL\n");
            return 1;
        }
    }
    printf("Multiplication tests................................................. PASSED\n");
    return 0;
}

int 
test_basis(ec_curve_t *curve){

    jac_point_t P, Q, TP, TQ, R;
    ec_curve_t twist;

    jac_rational_basis(&P, &Q, &TP, &TQ, &twist, curve);

    if(!jac_is_on_curve(&P,curve)){
        printf("P not on curve\n");
        return 1;
    }
    if(!jac_is_on_curve(&Q,&twist)){
        printf("Q not on twist\n");
        return 1;
    }
    if(!jac_is_on_curve(&TP,curve)){
        printf("TP not on curve\n");
        return 1;
    }
    if(!jac_is_on_curve(&TQ,&twist)){
        printf("TQ not on twist\n");
        return 1;
    }
        
    jac_dbl_iter(&P, &P, TORSION_EVEN_POWER-2, curve);
    if(jac_is_zero(&P)){
        printf("Wrong order: 2^(e-2)*P = 0\n");
        return 1;
    }
    DBL(&P,&P,curve);
    if(!jac_is_zero(&P)){
        printf("Wrong order: 2^(e-1)*P != 0\n");
        return 1;
    }

    jac_dbl_iter(&Q, &Q, TORSION_EVEN_POWER-2, &twist);
    if(jac_is_zero(&Q)){
        printf("Wrong order: 2^(e-2)*Q = 0\n");
        return 1;
    }
    DBL(&Q,&Q,&twist);
    if(!jac_is_zero(&Q)){
        printf("Wrong order: 2^(e-1)*Q != 0\n");
        return 1;
    }

    DBL(&R,&TP,curve);
    if(!jac_is_zero(&R)){
        printf("TP is not a point of 2-torsion\n");
        return 1;
    }

    DBL(&R,&TQ,&twist);
    if(!jac_is_zero(&R)){
        printf("TQ is not a point of 2-torsion\n");
        return 1;
    }

    printf("Basis generation tests............................................... PASSED\n");
    return 0;

}

int
test_add_components(ec_curve_t *curve, unsigned int Ntest){
    ec_curve_t twist;
    ec_point_t P, Q;
    jac_point_t jP, jQ, jPpQ, jPmQ;
    xz_point_t xPpQ, xPmQ, xR, xS;
    add_components_t uvw;

    for(unsigned int i=0; i<Ntest; i++){
        ec_random_test(&P,curve);
        ec_random_test(&Q,curve);

        ec_to_jac(&jP,&P);
        ec_to_jac(&jQ,&Q);

        ADD(&jPpQ,&jP,&jQ,curve); // P+Q
        jac_to_xz(&xPpQ,&jPpQ);
        jac_neg(&jQ,&jQ);
        ADD(&jPmQ,&jP,&jQ,curve); // P-Q
        jac_neg(&jQ,&jQ);
        jac_to_xz(&xPmQ,&jPmQ);

        jac_to_xz_add_components(&uvw, &jP, &jQ, curve, false, false);

        fp_sub(&xR.x,&uvw.u,&uvw.v);
        fp_copy(&xR.z,&uvw.w);
        fp_add(&xS.x,&uvw.u,&uvw.v);
        fp_copy(&xS.z,&uvw.w);

        if(!xz_is_equal(&xR,&xPpQ)){
            printf("From curve: P+Q != (u-v:*:w)\n");
            return 1;
        }
        if(!xz_is_equal(&xS,&xPmQ)){
            printf("From curve: P-Q != (u+v:*:w)\n");
            return 1;
        }

        // Tests on the twist
        ec_twist(&twist, curve);
        ec_compute_A24(&twist);

        ec_random_test(&P,&twist);
        ec_random_test(&Q,&twist);

        ec_to_jac(&jP,&P);
        ec_to_jac(&jQ,&Q);

        ADD(&jPpQ,&jP,&jQ,&twist); // P+Q
        jac_to_xz(&xPpQ,&jPpQ);
        fp_neg(&xPpQ.x,&xPpQ.x);
        jac_neg(&jQ,&jQ);
        ADD(&jPmQ,&jP,&jQ,&twist); // P-Q
        jac_neg(&jQ,&jQ);
        jac_to_xz(&xPmQ,&jPmQ);
        fp_neg(&xPmQ.x,&xPmQ.x);

        jac_to_xz_add_components(&uvw, &jP, &jQ, curve, true, true);

        fp_add(&xR.x,&uvw.u,&uvw.v);
        fp_copy(&xR.z,&uvw.w);
        fp_sub(&xS.x,&uvw.u,&uvw.v);
        fp_copy(&xS.z,&uvw.w);

        if(!xz_is_equal(&xR,&xPpQ)){
            printf("From twist: P+Q != (u+v:*:w)\n");
            return 1;
        }
        if(!xz_is_equal(&xS,&xPmQ)){
            printf("From twist: P-Q != (u-v:*:w)\n");
            return 1;
        }

        // Test sum zeros
        ec_random_test(&P,curve);
        ec_to_jac(&jP,&P);

        jac_point_t zero;
        xz_point_t xP,xP2;
        jac_point_init(&zero);  
        jac_to_xz(&xP,&jP);
        jac_to_xz_add_components(&uvw, &jP, &zero, curve, false, false);
        fp_copy(&xP2.x, &uvw.u);
        fp_copy(&xP2.z, &uvw.w);
        if(!xz_is_equal(&xP,&xP2)){
            printf("From curve: P+0 != (xP:*:zP)\n");
            return 1;
        }
        jac_to_xz_add_components(&uvw, &zero, &jP ,curve, false, false);
        fp_copy(&xP2.x, &uvw.u);
        fp_copy(&xP2.z, &uvw.w);
        if(!xz_is_equal(&xP,&xP2)){
            printf("From curve: 0+P != (xP:*:zP)\n");
            return 1;
        }
        
        // Test double on curve
        xz_point_t x2P, x2P2;
        jac_point_t mP;
        jac_to_xz(&x2P,&jP);
        xDBL(&x2P,&x2P,curve);
        fp_add(&x2P.x,&x2P.x,&x2P.x);
        fp_add(&x2P.z,&x2P.z,&x2P.z);
        jac_to_xz_add_components(&uvw, &jP, &jP, curve, false, false);
        fp_copy(&x2P2.x, &uvw.u);
        fp_sub(&x2P2.z,&uvw.w, &uvw.v);
        if(!xz_is_equal(&x2P,&x2P2)){
            printf("From curve: P+P != (x_{2P}:*:z_{2P})\n");
            return 1;
        }
        fp_add(&x2P2.z,&uvw.w, &uvw.v);
        if(!xz_is_zero(&x2P2)){
            printf("From curve: P-P != (1:0)\n");
            return 1;
        }
        jac_neg(&mP, &jP);
        jac_to_xz_add_components(&uvw, &jP, &mP, curve, false, false);
        fp_copy(&x2P2.x, &uvw.u);
        fp_sub(&x2P2.z,&uvw.w, &uvw.v);
        if(!xz_is_zero(&x2P2)){
            printf("From curve: P+-P != (1:0)\n");
            return 1;
        }
        fp_add(&x2P2.z,&uvw.w, &uvw.v);
        if(!xz_is_equal(&x2P,&x2P2)){
            printf("From curve: P--P != (x_{2P}:*:z_{2P})\n");
            return 1;
        }

        // Test double on twist
        ec_random_test(&Q,&twist);
        ec_to_jac(&jQ,&Q);

        //xz_point_t x2P, x2P2;
        jac_to_xz(&x2P,&jQ);
        fp_neg(&x2P.x,&x2P.x);
        xDBL(&x2P,&x2P,curve);
        jac_to_xz_add_components(&uvw, &jQ, &jQ, curve, true, true);
        fp_copy(&x2P2.x, &uvw.u);
        fp_add(&x2P2.z,&uvw.w, &uvw.v);
        if(!xz_is_equal(&x2P,&x2P2)){
            printf("From twist: P+P != (x_{2P}:*:z_{2P})\n");
            return 1;
        }
        fp_sub(&x2P2.z,&uvw.w, &uvw.v);
        if(!xz_is_zero(&x2P2)){
            printf("From twist: P-P != (1:0)\n");
            return 1;
        }
        jac_neg(&mP, &jQ);
        jac_to_xz_add_components(&uvw, &jQ, &mP, curve, true, true);
        fp_copy(&x2P2.x, &uvw.u);
        fp_add(&x2P2.z,&uvw.w, &uvw.v);
        if(!xz_is_zero(&x2P2)){
            printf("From twist: P+-P != (1:0)\n");
            return 1;
        }
        fp_sub(&x2P2.z,&uvw.w, &uvw.v);
        if(!xz_is_equal(&x2P,&x2P2)){
            printf("From twist: P--P != (x_{2P}:*:z_{2P})\n");
            return 1;
        }

        
    }
    printf("Barycentric coordinates tests........................................ PASSED\n");
    return 0;
}

int
main(int argc, char *argv[])
{
    uint32_t seed[12] = { 0 };
    int iterations = 100;
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
    printf("Testing elliptic curve arithmetic over GF(p):\n\n");
    print_seed(seed);

//#if defined(TARGET_BIG_ENDIAN)
    //for (int i = 0; i < 12; i++) {
        //seed[i] = BSWAP32(seed[i]);
    //}
//#endif

    randombytes_init((unsigned char *)seed, NULL, 256);

    ec_curve_t curve;
    ec_curve_init(&curve);

    res |= test_jacobian(&curve, iterations);
    res |= test_xDBL_xADD(&curve, iterations);
    res |= test_xDBLADD(&curve, iterations);
    res |= test_zero_identities(&curve, iterations);
    res |= test_mul(&curve, iterations);
    res |= test_basis(&curve);
    res |= test_add_components(&curve, iterations);
    //res |= test_xDBL_variants(&curve, iterations);
    //res |= test_zero_identities(&curve, iterations);
    //res |= test_jacobian(&curve, iterations);

    //fp2_random_test(&(curve.C));
    //fp2_mul(&(curve.A), &(curve.A), &(curve.C));
    //ec_curve_normalize_A24(&curve);

    //res |= test_xDBL_xADD(&curve, iterations);
    //res |= test_xDBLADD(&curve, iterations);
    //res |= test_xDBL_variants(&curve, iterations);
    //res |= test_zero_identities(&curve, iterations);
    //res |= test_jacobian(&curve, iterations);

    if (res) {
        printf("Tests failed!\n");
    } else {
        printf("All ec arithmetic tests passed.\n");
    }

    return res;
}
