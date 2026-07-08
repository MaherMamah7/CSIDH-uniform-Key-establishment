#include <assert.h>
#include <stdio.h>
#include <mp.h>
#include <ec.h>

void 
ec_curve_init(ec_curve_t *E){
    fp_copy(&E->A,&A0);
    E->is_ws_model_computed=false;
    E->is_A24_computed=false;
}

void
jac_point_init(jac_point_t *P)
{ // Initialize point as identity element (0:1:0)
	fp_set_zero(&(P->x));
    fp_set_one(&(P->y));
    fp_set_zero(&(P->z));
}

void
xz_point_init(xz_point_t *P)
{ // Initialize point as identity element (1:0)
    fp_set_one(&(P->x));
    fp_set_zero(&(P->z));
}

void
ec_point_init(ec_point_t *P)
{ // Initialize point as identity element (0:1:0)
    fp_set_zero(&(P->x));
    fp_set_one(&(P->y));
    fp_set_zero(&(P->z));
}

int 
ec_curve_verify_A(const fp_t *A){
	// Verify the Montgomery coefficient A is valid (A^2-4 \ne 0)
	// Return 1 if curve is valid, 0 otherwise
    fp_t t;
    fp_set_small(&t,2);// t=2
    if (fp_is_equal(A, &t))
        return 0;
    fp_neg(&t, &t); // t=-2
    if (fp_is_equal(A, &t))
        return 0;
    return 1;
}

int ec_curve_init_from_A(ec_curve_t *E, const fp_t *A){
	fp_copy(&E->A,A);
	E->is_ws_model_computed=false;
    E->is_A24_computed=false;
	return ec_curve_verify_A(A);
}

void ec_compute_ws(ec_curve_t *E){
	// Cost 2M + 1a
	// ao3 = A/3, a = 1-A^2/3
	if(!(E->is_ws_model_computed)){
		fp_t one;
		fp_set_one(&one);
		fp_div3(&E->ao3,&E->A);
		fp_mul(&E->a,&E->ao3,&E->A);
		fp_sub(&E->a,&one,&E->a);
        E->is_ws_model_computed = true;
	}
}

void ec_compute_A24(ec_curve_t *E){
    // Cost: 2M + a
    if(!(E->is_A24_computed)){
        fp_t two;
        fp_set_small(&two,2);
        fp_add(&E->A24,&E->A,&two);
        fp_half(&E->A24,&E->A24);
        fp_half(&E->A24,&E->A24);
        E->is_A24_computed = true;
    }
}

void ec_twist(ec_curve_t *Et, const ec_curve_t *E){
    fp_neg(&Et->A,&E->A);
    Et->is_ws_model_computed = false;
    Et->is_A24_computed = false;
}

void ec_j_inv(fp_t *j_inv, const ec_curve_t *curve){
	// j-invariant computation for Montgommery coefficient A
	fp_t t0, t1;

	fp_sqr(&t1,&curve->A); // t1 = A^2
	fp_set_small(&t0,3); // t0 = 3
	fp_sub(&t1,&t1,&t0); // t1 = t1-t0 = A^2-3
	fp_sqr(j_inv,&t1); // (A^2-3)^2
	fp_mul(j_inv,j_inv,&t1); // (A^2-3)^3
	fp_mul_small(j_inv,j_inv,256); // 256(A^2-3)^3
	fp_set_one(&t0); // t0 = 1
	fp_sub(&t1,&t1,&t0); // t1 = t1-t0 = t1 - 1 = A^2-4
	fp_inv(&t1); // t1 = 1/(A^2-4)
	fp_mul(j_inv,j_inv,&t1); // j_inv = 256(A^2-3)^3/(A^2-4)
}

void ec_to_jac(jac_point_t *Q, const ec_point_t *P){
	// x = X_mont/Z_mont = X_mont*Z_mont/Z_mont^2
	// y = Y_mont/Z_mont = Y_mont*Z_mont^2/Z_mont^3
	// (X_jac: Y_jac: Z_jac) = (X_mont*Z_mont: Y_mont*Z_mont^2: Z_mont)

	fp_t z2, one;

	fp_sqr(&z2,&P->z);
	fp_mul(&Q->x,&P->x,&P->z);
	fp_mul(&Q->y,&P->y,&z2);
	fp_copy(&Q->z,&P->z);

	// If Z == 0, return (0:1:0)
	uint32_t sel = fp_is_zero(&P->z);
	fp_set_one(&one);
	fp_select(&Q->y,&Q->y,&one,sel);
}

void jac_to_ec(ec_point_t *Q, const jac_point_t *P){
	// x = X_jac/Z_jac^2 = X_jac*Z_jac/Z_jac^3
	// y = Y_jac/Z_jac^3
	// (X_mont: Y_mont: Z_mont) = (X_jac*Z_jac: Y_jac: Z_jac^3)

	fp_t z2;

	fp_sqr(&z2,&P->z);
	fp_mul(&Q->z,&z2,&P->z);
	fp_mul(&Q->x,&P->x,&P->z);
	fp_copy(&Q->y,&P->y);
}



void jac_to_ws(jac_ws_point_t *Q, const jac_point_t *P, ec_curve_t *E){
	// X_ws = X_jac+(A*Z_jac^2)/3, Y_ws = Y_jac, Z_ws = Z_jac, T_ws = a*Z_ws^4
	// Cost of 2M + 2S + a + (precomputation ec_compute_ws(E)).
	ec_compute_ws(E);

	fp_sqr(&Q->t, &P->z);
	fp_mul(&Q->x, &E->ao3, &Q->t);
	fp_add(&Q->x, &Q->x, &P->x);
	fp_sqr(&Q->t, &Q->t);
	fp_mul(&Q->t, &Q->t, &E->a);
	fp_copy(&Q->y,&P->y);
	fp_copy(&Q->z,&P->z);
}



void ws_to_jac(jac_point_t *Q, const jac_ws_point_t *P, ec_curve_t *E){
	// X_jac = X_ws-(A*Z_ws^2)/3, Y_jac = Y_ws, Z_jac = Z_ws
	// Cost of 1M + 1S + a + (precomputation ec_compute_ws(E)).
	ec_compute_ws(E);

	fp_sqr(&Q->x, &P->z);
	fp_mul(&Q->x, &E->ao3, &Q->x);
	fp_sub(&Q->x, &P->x, &Q->x);
	fp_copy(&Q->y,&P->y);
	fp_copy(&Q->z,&P->z);
}

void ec_to_xz(xz_point_t *Q, const ec_point_t *P){
    // Maps infty = (0:1:0) to (1:0) and not (0:0)

    fp_copy(&Q->x,&P->x);
    fp_copy(&Q->z,&P->z);

    fp_t one;
    uint32_t ctl = ec_is_zero(P);
    fp_set_one(&one);

    fp_select(&Q->x,&Q->x,&one,ctl);
}

void jac_to_xz(xz_point_t *Q, const jac_point_t *P){
    // Maps infty = (0:1:0) to (1:0) and not (0:0)
    // Cost 1S

    fp_copy(&Q->x,&P->x);
    fp_sqr(&Q->z,&P->z);

    fp_t one;
    uint32_t ctl = jac_is_zero(P);
    fp_set_one(&one);

    fp_select(&Q->x,&Q->x,&one,ctl);
}

uint32_t
ec_is_zero(const ec_point_t *P)
{
    return fp_is_zero(&P->z);
}

uint32_t
ec_is_equal(const ec_point_t *P, const ec_point_t *Q)
{ // Evaluate if two points in standard Montgomery coordinates (X:Y:Z) are equal
  // Returns 0xFFFFFFFF (true) if P=Q, 0 (false) otherwise
    fp_t t0, t1, t2, t3;

    // Check if P, Q are the points at infinity
    uint32_t l_zero = ec_is_zero(P);
    uint32_t r_zero = ec_is_zero(Q);

    // Check if PX * QZ = QX * PZ AND PY * QZ = QY * PZ
    fp_mul(&t0, &P->x, &Q->z);
    fp_mul(&t1, &P->z, &Q->x);
    fp_mul(&t2, &P->y, &Q->z);
    fp_mul(&t3, &P->z, &Q->y);
    uint32_t lr_equal = fp_is_equal(&t0, &t1) & fp_is_equal(&t2, &t3);

    // Points are equal if
    // - Both are zero, or
    // - neither are zero AND PX * QZ = QX * PZ AND PY * QZ = QY * PZ
    return (l_zero & r_zero) | (~l_zero & ~r_zero & lr_equal);
}

uint32_t
jac_is_zero(const jac_point_t *P)
{
    return fp_is_zero(&P->z);
}

uint32_t
jac_is_equal(const jac_point_t *P, const jac_point_t *Q)
{ // Evaluate if two points in Jacobian Montgomery coordinates (X:Y:Z) are equal
  // Returns 0xFFFFFFFF (true) if P=Q, 0 (false) otherwise
    fp_t u0, u1, t0, t1, t2, t3;

    // Check if P, Q are the points at infinity
    uint32_t l_zero = jac_is_zero(P);
    uint32_t r_zero = jac_is_zero(Q);

    // Check if PX * QZ^2 = QX * PZ^2 AND PY * QZ^3 = QY * PZ^3
    fp_sqr(&u0, &Q->z);// u0 = QZ^2
    fp_sqr(&u1, &P->z);// u1 = PZ^2
    fp_mul(&t0, &P->x, &u0);// t0 = PX*QZ^2
    fp_mul(&t1, &u1, &Q->x);// t1 = QX*PZ^2
    fp_mul(&u0, &u0, &Q->z);// u0 = QZ^3
    fp_mul(&u1, &u1, &P->z);// u1 = PZ^3
    fp_mul(&t2, &P->y, &u0);// t2 = PY*QZ^3
    fp_mul(&t3, &u1, &Q->y);// t3 = QY*PZ^3
    uint32_t lr_equal = fp_is_equal(&t0, &t1) & fp_is_equal(&t2, &t3);
    //printf("PX * QZ^2 == QX * PZ^2 : %u\n",fp_is_equal(&t0, &t1));
    //printf("PY * QZ^3 = QY * PZ^3 : %u\n",fp_is_equal(&t2, &t3));

    // Points are equal if
    // - Both are zero, or
    // - neither are zero AND PX * QZ = QX * PZ AND PY * QZ = QY * PZ
    return (l_zero & r_zero) | (~l_zero & ~r_zero & lr_equal);
}

uint32_t
xz_is_equal(const xz_point_t *P, const xz_point_t *Q)
{ // Evaluate if two points in Jacobian Montgomery coordinates (X:Z) are equal
  // Returns 0xFFFFFFFF (true) if P=Q, 0 (false) otherwise
    fp_t t0, t1;

    // Check if P, Q are the points at infinity
    uint32_t l_zero = xz_is_zero(P);
    uint32_t r_zero = xz_is_zero(Q);

    // Check if PX * QZ = QX * PZ
    fp_mul(&t0, &P->x, &Q->z);
    fp_mul(&t1, &P->z, &Q->x);
    uint32_t lr_equal = fp_is_equal(&t0, &t1);

    // Points are equal if
    // - Both are zero, or
    // - neither are zero AND PX * QZ = QX * PZ
    return (l_zero & r_zero) | (~l_zero & ~r_zero & lr_equal);
}

uint32_t
xz_is_zero(const xz_point_t *P){
    return fp_is_zero(&P->z);
}

void jac_neg(jac_point_t *Q, const jac_point_t *P){
	fp_copy(&Q->x,&P->x);
	fp_neg(&Q->y,&P->y);
	fp_copy(&Q->z,&P->z);
}

void ADD(jac_point_t *R, const jac_point_t *P, const jac_point_t *Q, const ec_curve_t *E){
	// Addition on a Montgomery curve, representation in Jacobian coordinates (X:Y:Z) corresponding
    // to (x,y) = (X/Z^2,Y/Z^3) This version receives the coefficient value A
    //
    // Complete routine, to handle all edge cases:
    //   if ZP == 0:            # P == inf
    //       return Q
    //   if ZQ == 0:            # Q == inf
    //       return P
    //   dy <- YQ*ZP**3 - YP*ZQ**3
    //   dx <- XQ*ZP**2 - XP*ZQ**2
    //   if dx == 0:             # x1 == x2
    //       if dy == 0:         # ... and y1 == y2: doubling case
    //           dy <- ZP*ZQ * (3*XP^2 + ZP^2 * (2*A*XP + ZP^2))
    //           dx <- 2*YP*ZP
    //       else:              # ... but y1 != y2, thus P = -Q
    //           return inf
    //   XR <- dy**2 - dx**2 * (A*ZP^2*ZQ^2 + XP*ZQ^2 + XQ*ZP^2)
    //   YR <- dy * (XP*ZQ^2 * dx^2 - XR) - YP*ZQ^3 * dx^3
    //   ZR <- dx * ZP * ZQ

    // Constant time processing:
    // - The case for P == 0 or Q == 0 is handled at the end with conditional select
    // - dy and dx are computed for both the normal and doubling cases, we switch when
    //   dx == dy == 0 for the normal case.
    // - If we have that P = -Q then dx = 0 and so ZR will be zero, giving us the point
    //   at infinity for "free".
    //
    // These current formula are expensive and I'm probably missing some tricks...
    // Thought I'd get the ball rolling.
    // Cost 17M + 6S + 13a
    fp_t t0, t1, t2, t3, u1, u2, v1, dx, dy;

    // If P is zero or Q is zero we will conditionally swap before returning.
    uint32_t ctl1 = fp_is_zero(&P->z);
    uint32_t ctl2 = fp_is_zero(&Q->z);

    // Precompute some values
    fp_sqr(&t0, &P->z); // t0 = z1^2
    fp_sqr(&t1, &Q->z); // t1 = z2^2

    // Compute dy and dx for ordinary case
    fp_mul(&v1, &t1, &Q->z); // v1 = z2^3
    fp_mul(&t2, &t0, &P->z); // t2 = z1^3
    fp_mul(&v1, &v1, &P->y); // v1 = y1z2^3
    fp_mul(&t2, &t2, &Q->y); // t2 = y2z1^3
    fp_sub(&dy, &t2, &v1);   // dy = y2z1^3 - y1z2^3
    fp_mul(&u2, &t0, &Q->x); // u2 = x2z1^2
    fp_mul(&u1, &t1, &P->x); // u1 = x1z2^2
    fp_sub(&dx, &u2, &u1);   // dx = x2z1^2 - x1z2^2

    // Compute dy and dx for doubling case
    fp_add(&t1, &P->y, &P->y);   // dx_dbl = t1 = 2y1
    fp_add(&t2, &E->A, &E->A); // t2 = 2A
    fp_mul(&t2, &t2, &P->x);     // t2 = 2Ax1
    fp_add(&t2, &t2, &t0);       // t2 = 2Ax1 + z1^2
    fp_mul(&t2, &t2, &t0);       // t2 = z1^2 * (2Ax1 + z1^2)
    fp_sqr(&t0, &P->x);          // t0 = x1^2
    fp_add(&t2, &t2, &t0);       // t2 = x1^2 + z1^2 * (2Ax1 + z1^2)
    fp_add(&t2, &t2, &t0);       // t2 = 2*x1^2 + z1^2 * (2Ax1 + z1^2)
    fp_add(&t2, &t2, &t0);       // t2 = 3*x1^2 + z1^2 * (2Ax1 + z1^2)
    fp_mul(&t2, &t2, &Q->z);     // dy_dbl = t2 = z2 * (3*x1^2 + z1^2 * (2Ax1 + z1^2))

    // If dx is zero and dy is zero swap with double variables
    uint32_t ctl = fp_is_zero(&dx) & fp_is_zero(&dy);
    fp_select(&dx, &dx, &t1, ctl);
    fp_select(&dy, &dy, &t2, ctl);

    // Some more precomputations 
    fp_mul(&t0, &P->z, &Q->z); // t0 = z1z2
    fp_sqr(&t1, &t0);          // t1 = (z1z2)^2
    fp_sqr(&t2, &dx);          // t2 = dx^2
    fp_sqr(&t3, &dy);          // t3 = dy^2

    // Compute x3 = dy**2 - dx**2 * (A*ZP^2*ZQ^2 + XP*ZQ^2 + XQ*ZP^2)
    fp_mul(&R->x, &E->A, &t1); // x3 = A*(z1z2)^2
    fp_add(&R->x, &R->x, &u1);  // x3 = A*(z1z2)^2 + u1
    fp_add(&R->x, &R->x, &u2);  // x3 = A*(z1z2)^2 + u1 + u2
    fp_mul(&R->x, &R->x, &t2);  // x3 = dx^2 * (A*(z1z2)^2 + u1 + u2)
    fp_sub(&R->x, &t3, &R->x);  // x3 = dy^2 - dx^2 * (A*(z1z2)^2 + u1 + u2)

    // Compute y3 = dy * (XP*ZQ^2 * dx^2 - XR) - YP*ZQ^3 * dx^3
    fp_mul(&R->y, &u1, &t2);     // y3 = u1 * dx^2
    fp_sub(&R->y, &R->y, &R->x); // y3 = u1 * dx^2 - x3
    fp_mul(&R->y, &R->y, &dy);   // y3 = dy * (u1 * dx^2 - x3)
    fp_mul(&t3, &t2, &dx);       // t3 = dx^3
    fp_mul(&t3, &t3, &v1);       // t3 = v1 * dx^3
    fp_sub(&R->y, &R->y, &t3);   // y3 = dy * (u1 * dx^2 - x3) - v1 * dx^3

    // Compute z3 = dx * z1 * z2
    fp_mul(&R->z, &dx, &t0);

    // Finally, we need to set R = P is Q.Z = 0 and R = Q if P.Z = 0
    select_jac_point(R, R, Q, ctl1);
    select_jac_point(R, R, P, ctl2);
}

void DBL(jac_point_t *Q, const jac_point_t *P, const ec_curve_t *E){
	// Cost of 6M + 6S.
	// Doubling on a Montgomery curve, representation in Jacobian coordinates (X:Y:Z) corresponding to
	// (X/Z^2,Y/Z^3) This version receives the coefficient value A
    fp_t t0, t1, t2, t3;

    uint32_t flag = fp_is_zero(&P->x) & fp_is_zero(&P->z);

    fp_sqr(&t0, &P->x); // t0 = x1^2
    fp_add(&t1, &t0, &t0);
    fp_add(&t0, &t0, &t1); // t0 = 3x1^2
    fp_sqr(&t1, &P->z);    // t1 = z1^2
    fp_mul(&t2, &P->x, &E->A);
    fp_add(&t2, &t2, &t2); // t2 = 2Ax1
    fp_add(&t2, &t1, &t2); // t2 = 2Ax1+z1^2
    fp_mul(&t2, &t1, &t2); // t2 = z1^2(2Ax1+z1^2)
    fp_add(&t2, &t0, &t2); // t2 = alpha = 3x1^2 + z1^2(2Ax1+z1^2)
    fp_mul(&Q->z, &P->y, &P->z);
    fp_add(&Q->z, &Q->z, &Q->z); // z2 = 2y1z1
    fp_sqr(&t0, &Q->z);
    fp_mul(&t0, &t0, &E->A); // t0 = 4Ay1^2z1^2
    fp_sqr(&t1, &P->y);
    fp_add(&t1, &t1, &t1);     // t1 = 2y1^2
    fp_add(&t3, &P->x, &P->x); // t3 = 2x1
    fp_mul(&t3, &t1, &t3);     // t3 = 4x1y1^2
    fp_sqr(&Q->x, &t2);        // x2 = alpha^2
    fp_sub(&Q->x, &Q->x, &t0); // x2 = alpha^2 - 4Ay1^2z1^2
    fp_sub(&Q->x, &Q->x, &t3);
    fp_sub(&Q->x, &Q->x, &t3); // x2 = alpha^2 - 4Ay1^2z1^2 - 8x1y1^2
    fp_sub(&Q->y, &t3, &Q->x); // y2 = 4x1y1^2 - x2
    fp_mul(&Q->y, &Q->y, &t2); // y2 = alpha(4x1y1^2 - x2)
    fp_sqr(&t1, &t1);          // t1 = 4y1^4
    fp_sub(&Q->y, &Q->y, &t1);
    fp_sub(&Q->y, &Q->y, &t1); // y2 = alpha(4x1y1^2 - x2) - 8y1^4

    fp_select(&Q->x, &Q->x, &P->x, -flag); // What kind of magic is that?
    fp_select(&Q->z, &Q->z, &P->z, -flag);
}

void DBL_ws(jac_ws_point_t *Q, const jac_ws_point_t *P){
	// Cost of 3M + 5S + 14a.
	// Doubling on a Weierstrass curve, representation in modified Jacobian coordinates
	// (X:Y:Z:T=a*Z^4) corresponding to (X/Z^2,Y/Z^3), where a is the curve coefficient.
	// Formula from https://hyperelliptic.org/EFD/g1p/auto-shortw-modified.html

    //uint32_t flag = fp_is_zero(&P->x) & fp_is_zero(&P->z);
    //printf("flag = %u\n",flag);

    fp_t xx, c, cc, r, s, m;
    // XX = X^2
    fp_sqr(&xx, &P->x);
    // A = 2*Y^2
    fp_sqr(&c, &P->y);
    fp_add(&c, &c, &c);
    // AA = A^2
    fp_sqr(&cc, &c);
    // R = 2*AA
    fp_add(&r, &cc, &cc);
    // S = (X+A)^2-XX-AA
    fp_add(&s, &P->x, &c);
    fp_sqr(&s, &s);
    fp_sub(&s, &s, &xx);
    fp_sub(&s, &s, &cc);
    // M = 3*XX+T1
    fp_add(&m, &xx, &xx);
    fp_add(&m, &m, &xx);
    fp_add(&m, &m, &P->t);
    // X3 = M^2-2*S
    fp_sqr(&Q->x, &m);
    fp_sub(&Q->x, &Q->x, &s);
    fp_sub(&Q->x, &Q->x, &s);
    // Z3 = 2*Y*Z
    fp_mul(&Q->z, &P->y, &P->z);
    fp_add(&Q->z, &Q->z, &Q->z);
    // Y3 = M*(S-X3)-R
    fp_sub(&Q->y, &s, &Q->x);
    fp_mul(&Q->y, &Q->y, &m);
    fp_sub(&Q->y, &Q->y, &r);
    // T3 = 2*R*T1
    fp_mul(&Q->t, &P->t, &r);
    fp_add(&Q->t, &Q->t, &Q->t);

    //fp_t QX, QZ;

    //fp_copy(&QX,&Q->x);
    //fp_copy(&QZ,&Q->z);

    //fp_select(&Q->x, &Q->x, &P->x, -flag);
    //fp_select(&Q->z, &Q->z, &P->z, -flag);

    //printf("QX == QX: %u\n",fp_is_equal(&QX,&Q->x));
    //printf("QZ == QZ: %u\n",fp_is_equal(&QZ,&Q->z));
}

void jac_dbl_iter(jac_point_t *Q, const jac_point_t *P, const unsigned int n, ec_curve_t *E){
	if (n == 0) {
        copy_jac_point(Q, P);
        return;
    }

    if (n<=2){
    	DBL(Q, P, E);
    	for(unsigned int i=1; i<n; i++){
    		DBL(Q, Q, E);
    	}
    }
    else{
    	jac_ws_point_t R;
    	
    	jac_to_ws(&R, P, E);

    	for(unsigned int i=0; i<n; i++){
    		DBL_ws(&R, &R);
    	}

    	ws_to_jac(Q,&R,E);
    }
}

void jac_MUL(jac_point_t *Q, const jac_point_t *P, const digit_t *k, const int kbits, 
    ec_curve_t *curve){

    xz_point_t P0, Q0, R0;
    ec_point_t P_ec, Q_ec;

    jac_to_ec(&P_ec,P);
    ec_to_xz(&P0,&P_ec);

    xMUL(&Q0, &R0, &P0, k, kbits, curve);

    recover_y(&Q_ec, &P_ec, &Q0, &R0, curve);
    ec_to_jac(Q,&Q_ec);
}

uint8_t
jac_to_xz_add_components(add_components_t *add_comp, const jac_point_t *P, const jac_point_t *Q, 
    const ec_curve_t *E,  const bool twistP, const bool twistQ)
{
    // Take P and Q in E distinct, two jac_point_t, return three components u,v and w in the base fp such
    // that the xz coordinates of P+Q are (u-v:w) and of P-Q are (u+v:w)
    // Now also works when P == 0 or Q == 0 (but not both)
    // Cost 11M + 5S + 7a

    fp_t t0, t1, t2, t3, t4, t5, t6, xP, xQ, mt4, mt5, zP2, zQ2;
    uint32_t ctl1 = jac_is_zero(P), ctl2 = jac_is_zero(Q), ctlt1 = -(int32_t)twistP, ctlt2 = -(int32_t)twistQ, ctl_dbl, ctl_same_point;

    //ec_compute_A24(E);

    // Change x into -x if points are on the twist
    fp_neg(&zP2,&P->x);
    fp_select(&xP,&P->x,&zP2,ctlt1);
    fp_neg(&zQ2,&Q->x);
    fp_select(&xQ,&Q->x,&zQ2,ctlt2);

    fp_sqr(&zP2, &P->z);             // zP2 = z1^2
    fp_sqr(&zQ2, &Q->z);             // zQ2 = z2^2

    fp_mul(&t2, &xP, &zQ2);          // t2 = x1z2^2
    fp_mul(&t3, &zP2, &xQ);          // t3 = z1^2x2
    fp_mul(&t4, &P->y, &Q->z);       // t4 = y1z2
    fp_mul(&t4, &t4, &zQ2);          // t4 = y1z2^3
    fp_mul(&t5, &P->z, &Q->y);       // t5 = z1y2
    fp_mul(&t5, &t5, &zP2);          // t5 = z1^3y2

    ctl_same_point = fp_is_equal(&t4, &t5);

    fp_mul(&t0, &zP2, &zQ2);         // t0 = (z1z2)^2
    fp_mul(&t6, &t4, &t5);           // t6 = (z1z_2)^3y1y2
    fp_add(&add_comp->v, &t6, &t6);  // v  = 2(z1z_2)^3y1y2
    fp_sqr(&t4, &t4);                // t4 = y1^2z2^6
    fp_sqr(&t5, &t5);                // t5 = z1^6y_2^2

    // If P is on the twist, y1 <- i*y1 so t4 <- -t4 
    fp_neg(&mt4,&t4);
    fp_select(&t4,&t4,&mt4,ctlt1);

    // If Q is on the twist, y2 <- i*y2 so t5 <- -t5 
    fp_neg(&mt5,&t5);
    fp_select(&t5,&t5,&mt5,ctlt2);

    fp_add(&t4, &t4, &t5);          // t4 = z1^6y_2^2 + y1^2z2^6
    fp_add(&t5, &t2, &t3);          // t5 = x1z2^2 + z_1^2x2
    fp_add(&t6, &t3, &t3);          // t6 = 2z_1^2x2
    fp_sub(&t6, &t5, &t6);          // t6 = lambda = x1z2^2 - z_1^2x2
    fp_sqr(&t6, &t6);               // t6 = lambda^2 = (x1z2^2 - z_1^2x2)^2
    fp_mul(&t1, &E->A, &t0);        // t1 = A*(z1z2)^2
    fp_add(&t1, &t5, &t1);          // t1 = gamma =A*(z1z2)^2 + x1z2^2 +z_1^2x2
    fp_mul(&t1, &t1, &t6);          // t1 = gamma*lambda^2
    fp_sub(&add_comp->u, &t4, &t1); // u  = z1^6y_2^2 + y1^2z2^6 - gamma*lambda^2
    fp_mul(&add_comp->w, &t6, &t0); // w  = (z1z2)^2(lambda)^2
    
    // Deal with P == ±Q
    ctl_dbl = fp_is_zero(&t6);      // true if P == ±Q
    fp_add(&t0, &xP, &zP2);         // t0 = x1 + z1^2
    fp_sqr(&t0, &t0);               // t0 = (x1 + z1^2)^2
    fp_sub(&t1, &xP, &zP2);         // t1 = (x1 - z1^2)
    fp_sqr(&t1, &t1);               // t1 = (x1 - z1^2)^2
    fp_sub(&t2, &t0, &t1);          // t2 = 4x1z1^2
    fp_mul(&t5, &t0, &t1);          // t5 = (x1 + z1^2)^2(x1 - z1^2)^2 = X_{2P}
    fp_add(&t5, &t5, &t5);          // t5 = 2(x1 + z1^2)^2(x1 - z1^2)^2 = 2X_{2P}
    fp_mul(&t0, &t2, &E->A24);      // t0 = ((A+2)/4)*4x1z1^2
    fp_add(&t0, &t0, &t1);          // t0 = (x1 - z1^2)^2 + ((A+2)/4)*4x1z1^2 
    fp_mul(&t6, &t0, &t2);          // t6 = 4x1z1^2[(x1 - z1^2)^2 + ((A+2)/4)*4x1z1^2] = Z_{2P}
    fp_neg(&t4, &t6);
 
    fp_select(&t4,&t6,&t4, ctl_same_point ^ ctlt1 );  // t4 = Z_{2P} or -Z_{2P}

    // P ==  Q, then (u:v:w) = (2X_(2P): -Z_(2P): Z_{2P}) with v position-swapped.
    // P == -Q, then (u:v:w) = (2X_(2P):  Z_(2P): Z_{2P}) with v position-swapped.
    fp_select(&add_comp->u,&add_comp->u,&t5,ctl_dbl);
    fp_select(&add_comp->v,&add_comp->v,&t4,ctl_dbl);
    fp_select(&add_comp->w,&add_comp->w,&t6,ctl_dbl);

    fp_set_zero(&t5);
    // If P == 0, (u:v:w)=(x_Q:0:z_Q^2)
    fp_select(&add_comp->u,&add_comp->u,&xQ,ctl1);
    fp_select(&add_comp->v,&add_comp->v,&t5,ctl1|ctl2);
    fp_select(&add_comp->w,&add_comp->w,&zQ2,ctl1);

    // If Q == 0, (u:v:w)=(x_P:0:z_P^2)
    fp_select(&add_comp->u,&add_comp->u,&xP, ctl2);
    fp_select(&add_comp->w,&add_comp->w,&zP2,ctl2);


    return (ctl_dbl & 1) & ~(ctl1|ctl2); // 1 if double but not when zero happend.  
}

void
xDBL(xz_point_t *Q, const xz_point_t *P, ec_curve_t *E)
{
    // Doubling of a Montgomery point in projective coordinates (X:Z).
    // Input: projective Montgomery x-coordinates P = (XP:ZP), where xP=XP/ZP, and
    //        the Montgomery curve constants A24 = (A+2)/4 is (pre)computed.
    // Output: projective Montgomery x-coordinates Q <- 2*P = (XQ:ZQ) such that x(2P)=XQ/ZQ.
    // Cost 3M + 2S + 4a
    fp_t t0, t1, t2;

    ec_compute_A24(E);

    fp_add(&t0, &P->x, &P->z);
    fp_sqr(&t0, &t0);
    fp_sub(&t1, &P->x, &P->z);
    fp_sqr(&t1, &t1);
    fp_sub(&t2, &t0, &t1);
    fp_mul(&Q->x, &t0, &t1);
    fp_mul(&t0, &t2, &E->A24);
    fp_add(&t0, &t0, &t1);
    fp_mul(&Q->z, &t0, &t2);
}

void
xADD(xz_point_t *R, const xz_point_t *P, const xz_point_t *Q, const xz_point_t *PQ)
{
    // Differential addition of Montgomery points in projective coordinates (X:Z).
    // Input: projective Montgomery points P=(XP:ZP) and Q=(XQ:ZQ) such that xP=XP/ZP and xQ=XQ/ZQ, and difference
    //        PQ=P-Q=(XPQ:ZPQ).
    // Output: projective Montgomery point R <- P+Q = (XR:ZR) such that x(P+Q)=XR/ZR.
    // Cost: 4M + 2S + 6a
    fp_t t0, t1, t2, t3;

    fp_add(&t0, &P->x, &P->z);
    fp_sub(&t1, &P->x, &P->z);
    fp_add(&t2, &Q->x, &Q->z);
    fp_sub(&t3, &Q->x, &Q->z);
    fp_mul(&t0, &t0, &t3);
    fp_mul(&t1, &t1, &t2);
    fp_add(&t2, &t0, &t1);
    fp_sub(&t3, &t0, &t1);
    fp_sqr(&t2, &t2);
    fp_sqr(&t3, &t3);
    fp_mul(&t2, &PQ->z, &t2);
    fp_mul(&R->z, &PQ->x, &t3);
    fp_copy(&R->x, &t2);
}

void
xDBLADD(xz_point_t *R,
        xz_point_t *S,
        const xz_point_t *P,
        const xz_point_t *Q,
        const xz_point_t *PQ,
        ec_curve_t *E)
{ 
    // Simultaneous doubling and differential addition.
    // Input:  projective Montgomery points P=(XP:ZP) and Q=(XQ:ZQ) such that xP=XP/ZP and xQ=XQ/ZQ, the difference
    //         PQ=P-Q=(XPQ:ZPQ), and the Montgomery curve constants A24 = (A+2)/4 is (pre)computed.
    // Output: projective Montgomery points R <- 2*P = (XR:ZR) such that x(2P)=XR/ZR, and S <- P+Q = (XS:ZS) such that =
    //         x(Q+P)=XS/ZS.
    // Cost: 7M + 4S + 8a (gain: it ain't much, but that's still something...)
    fp_t t0, t1, t2;

    ec_compute_A24(E);

    fp_add(&t0, &P->x, &P->z);
    fp_sub(&t1, &P->x, &P->z);
    fp_sqr(&R->x, &t0);
    fp_sub(&t2, &Q->x, &Q->z);
    fp_add(&S->x, &Q->x, &Q->z);
    fp_mul(&t0, &t0, &t2);
    fp_sqr(&R->z, &t1);
    fp_mul(&t1, &t1, &S->x);
    fp_sub(&t2, &R->x, &R->z);
    fp_mul(&R->x, &R->x, &R->z);
    fp_mul(&S->x, &E->A24, &t2);
    fp_sub(&S->z, &t0, &t1);
    fp_add(&R->z, &R->z, &S->x);
    fp_add(&S->x, &t0, &t1);
    fp_mul(&R->z, &R->z, &t2);
    fp_sqr(&S->z, &S->z);
    fp_sqr(&S->x, &S->x);
    fp_mul(&S->z, &S->z, &PQ->x);
    fp_mul(&S->x, &S->x, &PQ->z);
}


void recover_y(ec_point_t *R, const ec_point_t *P, const xz_point_t *Q, 
    const xz_point_t *PQ, const ec_curve_t *E){
    // Cost: 13M + 1S + 9a
    fp_t v1, v2, v3, v4, v5, v6, one, t;

    fp_mul(&v1,&P->x,&Q->z); // v1 = X_P*Z_Q
    fp_mul(&v5,&Q->x,&P->z); // v5 = X_Q*z_P
    fp_add(&v2,&v5,&v1); // v2 = X_Q*z_P + X_P*Z_Q 
    fp_sub(&v3,&v5,&v1); // v3 = X_Q*z_P - X_P*Z_Q 
    fp_sqr(&v3,&v3); // v3 = (X_Q*z_P - X_P*Z_Q)^2
    fp_mul(&v3,&v3,&PQ->x); // v3 = X_PQ(X_Q*z_P - X_P*Z_Q)^2
    fp_mul(&v6,&P->z,&Q->z); // v6 = z_P*z_Q
    fp_mul(&v1,&E->A,&v6); // v1 = A*z_P*_Q
    fp_add(&v1,&v1,&v1); // v1 = 2*A*z_P*z_Q
    fp_add(&v2,&v1,&v2); // v2 = X_Q*z_P + X_P*Z_Q + 2*A*z_P*z_Q
    fp_mul(&v4,&P->x,&Q->x); // v4 = x_P*x_Q
    fp_add(&v4,&v4,&v6); // v4 = x_P*x_Q + z_P*z_Q
    fp_mul(&v2,&v2,&v4); // v2 = (x_P*x_Q + z_P*z_Q)*(X_Q*z_P + X_P*Z_Q + 2*A*z_P*z_Q)
    fp_mul(&v1,&v1,&v6); // v1 = 2*A*z_P^2*z_Q^2
    fp_sub(&v2,&v2,&v1); // v2 = (x_P*x_Q + z_P*z_Q)*(X_Q*z_P + X_P*Z_Q + 2*A*z_P*z_Q) 
    // -2*A*z_P^2*z_Q^2
    fp_mul(&v2,&v2,&PQ->z); // v2 = z_PQ*((x_P*x_Q + z_P*z_Q)*(X_Q*z_P + X_P*Z_Q + 2*A*z_P*z_Q)
    // -2*A*z_P^2*z_Q^2)
    fp_sub(&R->y,&v2,&v3); // Y_Q = z_PQ*((x_P*x_Q + z_P*z_Q)*(X_Q*z_P + X_P*Z_Q + 2*A*z_P*z_Q)
    // -2*A*z_P^2*z_Q^2) - X_PQ(X_Q*z_P - X_P*Z_Q)^2
    fp_add(&v1,&P->y,&P->y); // v1 = 2*y_P
    fp_mul(&v1,&v1,&v6); // v1 = 2*y_P*z_P*z_Q
    fp_mul(&v1,&v1,&PQ->z); // v1 = 2*y_P*z_P*z_Q*z_PQ
    fp_mul(&R->x,&v1,&Q->x); // v1 = 2*y_P*z_P*z_Q*z_PQ*x_Q
    fp_mul(&R->z,&v1,&Q->z); // v1 = 2*y_P*z_P*z_Q^2*z_PQ

    // If P == 0 or Q == 0, return R = 0
    uint32_t ctl = ec_is_zero(P)|xz_is_zero(Q);
    fp_set_one(&one);
    fp_select(&R->y,&R->y,&one,ctl);

    // If P+Q == 0, return R = -P
    ctl = xz_is_zero(PQ);
    fp_neg(&t,&P->y);
    fp_select(&R->x,&R->x,&P->x,ctl);
    fp_select(&R->y,&R->y,&t,ctl);
    fp_select(&R->z,&R->z,&P->z,ctl);
}

void
xMUL(xz_point_t *Q, xz_point_t *R, const xz_point_t *P, const digit_t *k, const int kbits, 
    ec_curve_t *curve)
{ // The Montgomery ladder
  // Input: projective Montgomery point P=(XP:ZP) such that xP=XP/ZP, a scalar k of bitlength kbits, and
  //        the Montgomery curve constants (A:C) (or A24 = (A+2C/4C:1) if normalized).
  // Output: projective Montgomery points Q <- k*P = (XQ:ZQ) such that x(k*P)=XQ/ZQ
    // and R <- [k+1]*P = (XR:ZR) such that x([k+1]*P)=XR/ZR
    // Uses https://eprint.iacr.org/2017/212, Algorithm 6.
    xz_point_t R0, R1;
    uint32_t mask;
    unsigned int bit, prevbit = 0, swap;

    // R0 <- P, R1 <- 2*P
    copy_xz_point(&R0,P);
    xDBL(&R1,&R0,curve);

    // Main loop
    for (int i = kbits - 1; i >= 0; i--) {
        bit = (k[i >> LOG2RADIX] >> (i & (RADIX - 1))) & 1;
        swap = bit ^ prevbit;
        prevbit = bit;
        mask = 0 - (uint32_t)swap;

        cswap_xz_points(&R0, &R1, mask);
        xDBLADD(&R0, &R1, &R0, &R1, P, curve);
    }

    swap = 0 ^ prevbit;
    mask = 0 - (uint32_t)swap;

    cswap_xz_points(&R0, &R1, mask);

    copy_xz_point(Q,&R0);
    copy_xz_point(R,&R1);
}

void ec_MUL(ec_point_t *Q, const ec_point_t *P, const digit_t *k, const int kbits, 
    ec_curve_t *curve){
    
    xz_point_t P0, Q0, R0;
    ec_point_t Q_copy;

    ec_to_xz(&P0,P);

    xMUL(&Q0, &R0, &P0, k, kbits, curve);

    recover_y(&Q_copy, P, &Q0, &R0, curve);
    copy_point(Q,&Q_copy);
}

