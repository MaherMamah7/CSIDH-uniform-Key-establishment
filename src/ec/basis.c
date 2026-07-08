#include <ec.h>
#include <mp.h>
#include <gmp.h>
#include <assert.h>

uint32_t
x_is_on_curve(fp_t *y2, const fp_t *x, const ec_curve_t *E){
	fp_t t0, t1;

	fp_sqr(&t0,x); // t0 = x^2
	fp_mul(&t1,&E->A,x); // t1 = A*x
	fp_add(&t0,&t0,&t1); // t0 = x^2+A*x
	fp_set_one(&t1); // t1 =1
	fp_add(&t0,&t0,&t1); // t0 = x^2+A*x+1
	fp_mul(y2,&t0,x); // y2 = x*(x^2+A*x+1)
	return fp_is_square(y2);
}

uint32_t ec_is_on_curve(const ec_point_t *P, const ec_curve_t *E){
    fp_t t0, t1, t2, x, y2;

    fp_copy(&t0,&P->z);
    fp_inv(&t0); // t0 =1/Z
    fp_mul(&x,&P->x,&t0); // x = X/Z
    fp_mul(&y2,&P->y,&t0); // y = Y/Z
    fp_sqr(&y2,&y2); // y2 = y^2

    fp_sqr(&t0,&x); // t0 = x^2
    fp_mul(&t1, &E->A, &x); // t1 = A*x
    fp_set_one(&t2); // t2 = 1
    fp_add(&t2,&t2,&t1); //t2 = A*x + 1
    fp_add(&t2,&t2,&t0); // t2 = x^2 + A*x +1
    fp_mul(&t2,&t2,&x); // t2 = x*(x^2 + A*x +1)

    return fp_is_equal(&t2,&y2)|(fp_is_zero(&P->x)&fp_is_zero(&P->z));
}

uint32_t
tate_2(const ec_point_t *P, const ec_point_t *Q, const bool normalized){
	// Takes P a 2-torsion point and Q any point Q!=P
	// Returns 0 if tate_2(P,Q)=1 and 0xF..F if tate_2(P,Q)=-1
	// Cost: 1I + 1Exp + 4M + 1a
	if(normalized){// Assumes z(P)==z(Q)==1
		fp_t t0;
		fp_sub(&t0,&Q->x,&P->x);
		return ~fp_is_square(&t0);
	}
	else{
		fp_t t0, t1;
		fp_mul(&t0,&Q->x,&P->z); // t0 = x*z_P
		fp_mul(&t1,&Q->z,&P->x); // t1 = z*x_P
		fp_sub(&t0,&t0,&t1); // t0 = x*z_P - z*x_P
		fp_mul(&t1,&Q->z,&P->z); // t1 = z*z_P
		fp_inv(&t1); // t1 = 1/(z*z_P)
		fp_mul(&t0,&t0,&t1); // t0 = (x*z_P - z*x_P)/(z*z_P)
		return ~fp_is_square(&t0);
	}
	
}

uint32_t
self_tate_2(ec_point_t *P, const ec_curve_t *E){
	// Takes P a 2-torsion point and returns tate_2(P,P)
	// Returns 0 if tate_2(P,P)=1 and 0xF..F if tate_2(P,P)=-1
	// Uses the formula t_2(P,P)=f_{2,P}(D_P)^{(p-1)/2}
	// with f_{2,P}(D_P)=lt_P(f_{2,P})=lt_P(x-x_P)
	// =1/(3*x_P^2+2*A*x_P+1)=-1/(A*x_P+2) when x_P\neq 0
	// (with x_P normalised) 
	// Cost: 1Exp + 1M + 2a (+1 I)
	// Not constant time but used deterministically independently
	// of the secret.
	fp_t t0, t1;
	uint32_t ret;

	// Normalize
	if(!fp_is_one(&P->z)){
		fp_inv(&P->z);
		fp_mul(&P->x,&P->x,&P->z);
		fp_set_one(&P->z);
	}

	fp_mul(&t0,&P->x,&E->A); // t0 = A*x_P
	fp_set_small(&t1,2); // t1 = 2
	fp_add(&t0,&t0,&t1); // t0 = A*x_P+2 = -(3*x_P^2+2*A*x_P+1)

	ret = fp_is_square(&t0); // -1 is not a square in Fp
	
	return ret^(ret & fp_is_zero(&P->x)); // if x_P==0, tate_2(P,P) must be 1 (return 0)
}

void
two_torsion(ec_point_t *T0, ec_point_t *T1, ec_point_t *Tm1,const ec_curve_t *E){
	// T0, T1 and Tm1 are normalized (z=1)
	// Uses [PEGASIS https://eprint.iacr.org/2025/401.pdf, Lemma D.1]
	fp_t delta, four;
	ec_point_t P0, P1, Pm1;
	uint32_t sel0, sel1;

	// Point (0:0:1)
	fp_set_zero(&P0.x);
	fp_set_zero(&P0.y);
	fp_set_one(&P0.z);

	// Point (alpha:0:1) (alpha = (-A + sqrt(A^2-4))/2)
	fp_sqr(&delta,&E->A); // A^2
	fp_set_small(&four,4);
	fp_sub(&delta,&delta,&four); // A^2-4
	fp_sqrt(&delta); // sqrt(A^2-4)

	fp_sub(&P1.x,&delta,&E->A); // -A+sqrt(A^2-4)
	fp_half(&P1.x,&P1.x); // (-A+sqrt(A^2-4))/2

	fp_set_zero(&P1.y);
	fp_set_one(&P1.z);

	// Point (1/alpha:0:1)
	fp_sub(&Pm1.x,&P1.x,&delta);// 1/alpha=(-A-sqrt(A^2-4))/2
	fp_set_zero(&Pm1.y);
	fp_set_one(&Pm1.z);

	// Find T0 (among P1 and Pm1 because tate_2(P0,P0)=1)
	sel0 = self_tate_2(&P1,E);
	select_ec_point(T0,&Pm1,&P1,sel0); // tate_2(T0,T0) = -1
	select_ec_point(T1,&Pm1,&P1,~sel0); // Store T1 or Tm1 (to be determined)

	// Distinguish T1 and Tm1
	copy_point(Tm1,&P0);
	sel1 = tate_2(T1,T0,true);
	cswap_points(Tm1,T1,~sel1); // tate_2(T1,T0) = -1 and tate_2(Tm1,T0) = 1

	//assert(ec_is_on_curve(T0, E));
	//assert(ec_is_on_curve(T1, E));
	//assert(ec_is_on_curve(Tm1, E));
	//printf("T0 == P0: %u\n",ec_is_equal(T0,&P0));
	//printf("T0 == P1: %u\n",ec_is_equal(T0,&P1));
	//printf("T0 == Pm1: %u\n",ec_is_equal(T0,&Pm1));
	//printf("T1 == P0: %u\n",ec_is_equal(T1,&P0));
	//printf("T1 == P1: %u\n",ec_is_equal(T1,&P1));
	//printf("T1 == Pm1: %u\n",ec_is_equal(T1,&Pm1));
	//printf("Tm1 == P0: %u\n",ec_is_equal(Tm1,&P0));
	//printf("Tm1 == P1: %u\n",ec_is_equal(Tm1,&P1));
	//printf("Tm1 == Pm1: %u\n",ec_is_equal(Tm1,&Pm1));
}

void
ec_rational_basis(ec_point_t *P, ec_point_t *Q, ec_point_t *TP, ec_point_t *TQ, 
	ec_curve_t *Et, ec_curve_t *E){
	// Write p=c*2^f-1
	// Returns P\in E(Fp) and Q\in E^t(Fp) each of order 2^{f-1}
	// where E: y^2 = x^3+Ax^2+x and E^t: y^2 = x^3-Ax^2+x
	
	ec_point_t T0, T1, Tm1, OE, P1, P2, T0t, T1t, Tm1t;
	fp_t xP, yP, xQ, yQ, x, y;
	uint32_t sel, mu1, mu2;
	digit_t c;
	ibz_to_digits(&c, &TORSION_ODD);
	int kbits = n_bits(&c,1);

	two_torsion(&T0, &T1, &Tm1, E);

	// Choice of P
	fp_set_zero(&xP);
	fp_set_zero(&yP);
	for(int i=0; i<20; i++){
		fp_add(&x,&NQR_TABLE[i],&Tm1.x);
		sel = x_is_on_curve(&y,&x,E);
		fp_select(&xP,&xP,&x,sel);
		fp_select(&yP,&yP,&y,sel);
	}

	fp_copy(&P->x,&xP);
	//assert(x_is_on_curve(&y,&P->x,E));
	//assert(fp_is_equal(&y,&yP));
	//assert(fp_is_square(&yP));
	fp_sqrt(&yP);
	fp_copy(&P->y,&yP);
	fp_set_one(&P->z);

	//assert(ec_is_on_curve(P, E));

	ec_MUL(P, P, &c, kbits, E);

	//assert(ec_is_on_curve(P, E));

	// Choice of Q
	ec_twist(Et,E);

	fp_set_zero(&xQ);
	fp_set_zero(&yQ);
	for(int i=0; i<20; i++){
		// Recall that xQ in Et is -xQ in E
		// x_Q-x(T1) is a square iff -x_Q+x(T1) is not a square
		fp_sub(&x,&NQR_TABLE[i],&T1.x); // So QR_TABLE is useless (and always was...)
		sel = x_is_on_curve(&y,&x,Et);
		fp_select(&xQ,&xQ,&x,sel);
		fp_select(&yQ,&yQ,&y,sel);
	}

	fp_copy(&Q->x,&xQ);
	fp_sqrt(&yQ);
	fp_copy(&Q->y,&yQ);
	fp_set_one(&Q->z);

	//assert(ec_is_on_curve(Q, Et));

	ec_MUL(Q, Q, &c, kbits, Et);

	//assert(ec_is_on_curve(Q, Et));

	// Computing TP = omega(P)-P
	mu1 = tate_2(&Tm1,P,false);
	mu2 = tate_2(&T1,P,false);

	ec_point_init(&OE);

	// If mu1 == 0 and mu2 == 0 then TP = 0_E
	// If mu1 == 0 and mu2 == 0xF..F then TP = Tm1
	// If mu1 == 0xF..F and mu2 == 0 then TP = T1
	// If mu1 == 0xF..F and mu2 == 0xF..F then TP = T0  
	select_ec_point(&P1,&OE,&T1,mu1);
	select_ec_point(&P2,&Tm1,&T0,mu1);
	select_ec_point(TP,&P1,&P2,mu2);

	//assert(ec_is_on_curve(TP, E));

	// Computing TQ = omega(Q)
	fp_neg(&T0t.x,&T0.x);
	fp_set_zero(&T0t.y);
	fp_copy(&T0t.z,&T0.z);

	fp_neg(&T1t.x,&T1.x);
	fp_set_zero(&T1t.y);
	fp_copy(&T1t.z,&T1.z);

	fp_neg(&Tm1t.x,&Tm1.x);
	fp_set_zero(&Tm1t.y);
	fp_copy(&Tm1t.z,&Tm1.z);

	mu1 = tate_2(&Tm1t,Q,false);
	mu2 = tate_2(&T1t,Q,false);

	// If mu1 == 0 and mu2 == 0 then TQ = 0_E
	// If mu1 == 0 and mu2 == 0xF..F then TQ = T1t
	// If mu1 == 0xF..F and mu2 == 0 then TQ = Tm1t
	// If mu1 == 0xF..F and mu2 == 0xF..F then TQ = T0t  
	select_ec_point(&P1,&OE,&Tm1t,mu1);
	select_ec_point(&P2,&T1t,&T0t,mu1);
	select_ec_point(TQ,&P1,&P2,mu2);// TQ already on the twist

	//assert(ec_is_on_curve(TQ, Et));
}

void jac_rational_basis(jac_point_t *P, jac_point_t *Q, jac_point_t *TP, jac_point_t *TQ, 
	ec_curve_t *Et, ec_curve_t *E){

	ec_point_t Pec, Qec, TPec, TQec;

	ec_rational_basis(&Pec, &Qec,&TPec, &TQec, Et, E);

	ec_to_jac(P, &Pec);
	ec_to_jac(Q, &Qec);
	ec_to_jac(TP, &TPec);
	ec_to_jac(TQ, &TQec);
}
