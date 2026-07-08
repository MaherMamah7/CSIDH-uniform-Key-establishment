/** @file
 *
 * @authors Pierrick Dartois, based on 
 * work by Luca De Feo and Francisco RH
 * for the SQIsign v2.0 NIST submission
 *
 * @brief Elliptic curve stuff
 */

#ifndef EC_H
#define EC_H

#include <fp.h>
//#include <tools.h>
#include <stdio.h>
#include <constants.h>
#include <fp_constants.h>
#include <mp.h>
#include <gmp.h>

/** @defgroup ec Elliptic curves
 * @{
 */

/** @defgroup ec_t Data structures
 * @{
 */

/** @brief Projective point in (X:Z) Montgomery 
 * coordinates.
 *
 * @typedef xz_point_t
 *
 * @struct xz_point_t
 *
 * A projective point (X:Z) representing a point
 * on the Kummer line.
 */
typedef struct xz_point_t
{
    fp_t x;
    fp_t z;
} xz_point_t;

/** @brief Projective point in usual coordinates
 * in the Montgomery model y^2 = x^3 + Ax^2 + x.
 * May not be necessary (see pairings and basis 
 * generation).
 *
 * @typedef ec_point_t
 *
 * @struct ec_point_t
 *
 * A projective point in (X:Y:Z) coordinates such that
 * x = X/Z and y = Y/Z.
 */
typedef struct ec_point_t
{
    fp_t x;
    fp_t y;
    fp_t z;
} ec_point_t;

/** @brief Projective point in Jacobian coordinates
 * in the Montgomery model y^2 = x^3 + Ax^2 + x
 *
 * @typedef jac_point_t
 *
 * @struct jac_point_t
 *
 * A projective point in (X:Y:Z) coordinates such that
 * x = X/Z^2 and y = X/Z^3.
 */
typedef struct jac_point_t
{
    fp_t x;
    fp_t y;
    fp_t z;
} jac_point_t;

/** @brief Projective point in Jacobian coordinates
 * in the Weierstrass model y^2 = x^3 + ax + b
 * with a = 1-A^2/3 and b = A/3(2A^2/9-1)
 *
 * @typedef jac_ws_point_t
 *
 * @struct jac_ws_point_t
 *
 * A projective point in (X:Y:Z:T) coordinates such that
 * x = X/Z^2, y = X/Z^2 and T=a*Z^4
 */
typedef struct jac_ws_point_t
{
    fp_t x;
    fp_t y;
    fp_t z;
    fp_t t;
} jac_ws_point_t;

/** @brief Addition components
 *
 * @typedef add_components_t
 *
 * @struct add_components_t
 *
 * 3 components u,v,w that define the (X:Z) coordinates of both
 * addition and substraction of two distinct points with
 * P+Q =(u-v:w) and P-Q = (u+v:w)
 */
typedef struct add_components_t
{
    fp_t u;
    fp_t v;
    fp_t w;
} add_components_t;

/** @brief An elliptic curve
 *
 * @typedef ec_curve_t
 *
 * @struct ec_curve_t
 *
 * An elliptic curve in projective Montgomery form
 */
typedef struct ec_curve_t
{
    fp_t A;
    // Coefficients for the Weierstrass Jacobian model
    fp_t ao3;// ao3 = A/3
    fp_t a; // a = 1-A^2/3
    bool is_ws_model_computed; // says if ao3 and a have been computed
    // Coefficient for (X:Z) coordinates
    fp_t A24; // A24 = (A+2)/4
    bool is_A24_computed;
} ec_curve_t;

// end ec_t
/** @}
 */

/** @defgroup ec_curve_t Curves and isomorphisms
 * @{
 */

// Initalisation for curves and points
void ec_curve_init(ec_curve_t *E);
void jac_point_init(jac_point_t *P);
void xz_point_init(xz_point_t *P);
void ec_point_init(ec_point_t *P);

/**
 * @brief Verify that a Montgomery coefficient is valid
 *
 * @param A an fp_t
 *
 * @return 0  if curve is invalid, 1 otherwise
 */
int ec_curve_verify_A(const fp_t *A);

/**
 * @brief Initialize an elliptic curve from a coefficient
 *
 * @param A an fp_t
 * @param E the elliptic curve to initialize
 *
 * @return 0  if curve is invalid, 1 otherwise
 */
int ec_curve_init_from_A(ec_curve_t *E, const fp_t *A);

// Copying points and curves
static inline void
copy_point(ec_point_t *P, const ec_point_t *Q)
{
    fp_copy(&P->x, &Q->x);
    fp_copy(&P->y, &Q->y);
    fp_copy(&P->z, &Q->z);
}

static inline void
copy_xz_point(xz_point_t *P, const xz_point_t *Q)
{
    fp_copy(&P->x, &Q->x);
    fp_copy(&P->z, &Q->z);
}

static inline void
copy_jac_point(jac_point_t *P, const jac_point_t *Q)
{
    fp_copy(&P->x, &Q->x);
    fp_copy(&P->y, &Q->y);
    fp_copy(&P->z, &Q->z);
}

static inline void
copy_jac_ws_point(jac_ws_point_t *P, const jac_ws_point_t *Q)
{
    fp_copy(&P->x, &Q->x);
    fp_copy(&P->y, &Q->y);
    fp_copy(&P->z, &Q->z);
    fp_copy(&P->t, &Q->t);
}

static inline void
copy_curve(ec_curve_t *E1, const ec_curve_t *E2)
{
    fp_copy(&(E1->A), &(E2->A));
    E1->is_ws_model_computed = E2->is_ws_model_computed;
    fp_copy(&(E1->ao3), &(E2->ao3));
    fp_copy(&(E1->a), &(E2->a));
    E1->is_A24_computed = E2->is_A24_computed;
    fp_copy(&(E1->A24),&(E2->A24));
}

static inline void
select_curve(ec_curve_t *E,const ec_curve_t *E1, const ec_curve_t *E2, const uint32_t option)
{
    // If option = 0 then E <- E1, else if option = 0xFF...FF then E <- E2
    fp_select(&(E->A), &(E1->A), &(E2->A), option);
    E->is_ws_model_computed = (E1->is_ws_model_computed&(~option))||(E2->is_ws_model_computed&option);
    fp_select(&(E->ao3), &(E1->ao3), &(E2->ao3), option);
    fp_select(&(E->a), &(E1->a), &(E2->a), option);
    fp_select(&(E->A24), &(E1->A24), &(E2->A24), option);
    E->is_A24_computed = (E1->is_A24_computed&(~option))||(E2->is_A24_computed&option);
}

static inline void
cswap_curve(ec_curve_t *E1, ec_curve_t *E2, const uint32_t option)
{
    // If option = 0 we do not swap, else if option = 0xFF...FF then E1, E2 <- E2, E1
    bool b1, b2;

    fp_cswap(&(E1->A), &(E2->A), option);
    b1 = (E1->is_ws_model_computed&(~option))||(E2->is_ws_model_computed&option);
    b2 = (E2->is_ws_model_computed&(~option))||(E1->is_ws_model_computed&option);
    E1->is_ws_model_computed = b1;
    E2->is_ws_model_computed = b2;
    fp_cswap(&(E1->ao3), &(E2->ao3), option);
    fp_cswap(&(E1->a), &(E2->a), option);
    fp_cswap(&(E1->A24), &(E2->A24), option);
    b1 = (E1->is_A24_computed&(~option))||(E2->is_A24_computed&option);
    b2 = (E2->is_A24_computed&(~option))||(E1->is_A24_computed&option);
    E1->is_A24_computed = b1;
    E2->is_A24_computed = b2;
}

// Selecting points and curves
static inline void
select_jac_point(jac_point_t *Q, const jac_point_t *P1, const jac_point_t *P2, const uint32_t option)
{ // Select points in constant time
  // If option = 0 then Q <- P1, else if option = 0xFF...FF then Q <- P2
    fp_select(&(Q->x), &(P1->x), &(P2->x), option);
    fp_select(&(Q->y), &(P1->y), &(P2->y), option);
    fp_select(&(Q->z), &(P1->z), &(P2->z), option);
}

static inline void
select_ec_point(ec_point_t *Q, const ec_point_t *P1, const ec_point_t *P2, const uint32_t option)
{ // Select points in constant time
  // If option = 0 then Q <- P1, else if option = 0xFF...FF then Q <- P2
    fp_select(&(Q->x), &(P1->x), &(P2->x), option);
    fp_select(&(Q->y), &(P1->y), &(P2->y), option);
    fp_select(&(Q->z), &(P1->z), &(P2->z), option);
}

static inline void
cswap_points(ec_point_t *P, ec_point_t *Q, const uint32_t option)
{ // Swap points in constant time
  // If option = 0 then P <- P and Q <- Q, else if option = 0xFF...FF then P <- Q and Q <- P
    fp_cswap(&(P->x), &(Q->x), option);
    fp_cswap(&(P->y), &(Q->y), option);
    fp_cswap(&(P->z), &(Q->z), option);
}

static inline void
cswap_jac_points(jac_point_t *P, jac_point_t *Q, const uint32_t option)
{ // Swap points in constant time
  // If option = 0 then P <- P and Q <- Q, else if option = 0xFF...FF then P <- Q and Q <- P
    fp_cswap(&(P->x), &(Q->x), option);
    fp_cswap(&(P->y), &(Q->y), option);
    fp_cswap(&(P->z), &(Q->z), option);
}

static inline void
cswap_xz_points(xz_point_t *P, xz_point_t *Q, const uint32_t option)
{ // Swap points in constant time
  // If option = 0 then P <- P and Q <- Q, else if option = 0xFF...FF then P <- Q and Q <- P
    fp_cswap(&(P->x), &(Q->x), option);
    fp_cswap(&(P->z), &(Q->z), option);
}

/**
 * @brief Compute Jacobian Weierstrass coefficients ao3 and a
 *
 * @param E a curve
 */
void ec_compute_ws(ec_curve_t *E);

/**
 * @brief Compute Montgomery coefficients A24 = (A+2)/4
 *
 * @param E a curve
 */
void ec_compute_A24(ec_curve_t *E);

/**
 * @brief Compute the twist of a Montgomery curve
 *
 * @param Et output twisted curve
 * @param E input curve
 * 
 * By convention E: y^2 = x^3 + A*x^2 + x and
 * Et: y^2 = x^3 - A*x^2 + x via the ismorphism
 * tau: (x,y) --> (-x,iy) 
 */
void ec_twist(ec_curve_t *Et, const ec_curve_t *E);

/**
 * @brief j-invariant.
 *
 * @param j_inv computed j_invariant
 * @param curve input curve
 */
void ec_j_inv(fp_t *j_inv, const ec_curve_t *curve);

// Conversion functions

/**
 * @brief Converts a point from standard projective coordinates
 * to Jacobian coordinates in the Montgomery model
 *
 * @param Q output converted point
 * @param P input point to convert
 */
void ec_to_jac(jac_point_t *Q, const ec_point_t *P);

/**
 * @brief Converts a point from Jacobian coordinates to standard
 * projective coordinates in the Montgomery model
 *
 * @param Q output converted point
 * @param P input point to convert
 */
void jac_to_ec(ec_point_t *Q, const jac_point_t *P);

/**
 * @brief Converts points from Montgomery to Weierstrass Jacobian
 * coordinates
 *
 * @param Q output converted point
 * @param P input point to convert
 * @param E parent curve
 * 
 * Computes the additional coordinate T = a*Z^4
 */
void jac_to_ws(jac_ws_point_t *Q, const jac_point_t *P, ec_curve_t *E);

/**
 * @brief Converts points from Weierstrass to Montgomery Jacobian
 * coordinates
 *
 * @param Q output converted point
 * @param P input point to convert
 * @param E parent curve
 * 
 * Forgets the T coordinate.
 */
void ws_to_jac(jac_point_t *Q, const jac_ws_point_t *P, ec_curve_t *E);

/**
 * @brief Converts points from (X:Y:Z) to (X:Z) Montgomery coordinates
 *
 * @param Q output converted point
 * @param P input point to convert
 * 
 * Forgets the Y coordinate.
 */
void ec_to_xz(xz_point_t *Q, const ec_point_t *P);

/**
 * @brief Converts points from  Jacobian (X:Y:Z) to (X:Z) Montgomery 
 * coordinates
 *
 * @param Q output converted point
 * @param P input point to convert
 * 
 * Forgets the Y coordinate.
 */
void jac_to_xz(xz_point_t *Q, const jac_point_t *P);


/** @}
 */

/** @defgroup ec_point_t Point operations
 * @{
 */

/**
 * @brief Point equality
 *
 * @param P a point in standard projective coordinates
 * @param Q a point in standard projective coordinates
 * @return 0xFFFFFFFF if equal, zero otherwise
 */
uint32_t ec_is_equal(const ec_point_t *P, const ec_point_t *Q);

/**
 * @brief Point at infinity equality
 *
 * @param P a point in standard projective coordinates
 * @return 0xFFFFFFFF if point at infinity, zero otherwise
 */
uint32_t ec_is_zero(const ec_point_t *P);

/**
 * @brief Point equality
 *
 * @param P a point in Jacobian coordinates
 * @param Q a point in Jacobian coordinates
 * @return 0xFFFFFFFF if equal, zero otherwise
 */
uint32_t jac_is_equal(const jac_point_t *P, const jac_point_t *Q);

/**
 * @brief Point at infinity equality
 *
 * @param P a point in Jacobian coordinates
 * @return 0xFFFFFFFF if point at infinity, zero otherwise
 */
uint32_t jac_is_zero(const jac_point_t *P);

/**
 * @brief Point equality
 *
 * @param P a point in (X:Z) Montgomery coordinates
 * @param Q a point in (X:Z) Montgomery coordinates
 * @return 0xFFFFFFFF if equal, zero otherwise
 */
uint32_t xz_is_equal(const xz_point_t *P, const xz_point_t *Q);

/**
 * @brief Point at infinity equality
 *
 * @param P a point in (X:Z) Montgomery coordinates
 * @return 0xFFFFFFFF if point at infinity, zero otherwise
 */
uint32_t xz_is_zero(const xz_point_t *P);

// Jacobian arithmetic

/**
 * @brief Negation of a point in Jacobian coordinates in the 
 * Montgomery model.
 *
 * @param Q output negated point
 * @param P input point to negate
 */
void jac_neg(jac_point_t *Q, const jac_point_t *P);

/**
 * @brief Sums two points in Jacobian coordinates in the 
 * Montgomery model.
 *
 * @param R output sum of points P+Q
 * @param P input point to sum
 * @param Q input point to sum
 * @param E parent elliptic curve
 * 
 * Complete algorithm that can handle all edge cases (P==0, Q==0, P=+/- Q...)
 */
void ADD(jac_point_t *R, const jac_point_t *P, const jac_point_t *Q, const ec_curve_t *E);

/**
 * @brief Doubles a point in Jacobian coordinates in the Montgomery model.
 *
 * @param Q output point duplication 2P
 * @param P input point to double
 * @param E parent elliptic curve
 * 
 * Complete algorithm that can handle the case P == 0.
 */
void DBL(jac_point_t *Q, const jac_point_t *P, const ec_curve_t *E);

/**
 * @brief Doubles a point in Jacobian coordinates in the Weierestrass model.
 *
 * @param Q output point duplication 2P
 * @param P input point to double
 * 
 * Complete algorithm that can handle the case P == 0.
 */
void DBL_ws(jac_ws_point_t *Q, const jac_ws_point_t *P);

/**
 * @brief Doubles a point multiple time in Jacobian coordinates in 
 * the Jacobian model.
 *
 * @param Q output point duplication [2^n]P
 * @param P input point to double n times
 * @param n number of times to double
 * @param E parent elliptic curve
 * 
 * Complete algorithm that can handle the case P == 0. Optimised 
 * using the Weierstrass model when multiple iterations are needed.
 */
void jac_dbl_iter(jac_point_t *Q, const jac_point_t *P, const unsigned int n, ec_curve_t *E);

/** 
 * @brief The Montgomery ladder with Jacobian projective 
 * coordinates in the Montgomery model
 * 
 * @param Q output point [k]P
 * @param P input point
 * @param k scalar multiple
 * @param kbits number of bits of k
 * @param curve parent curve
 * 
 * Uses xMUL and y reconstruction (revover_y) and converts back
 * to Jacobian coordinates.
 */
void jac_MUL(jac_point_t *Q, const jac_point_t *P, const digit_t *k, const int kbits, 
    ec_curve_t *curve);

/**
 * @brief Take P and Q in E distinct, two jac_point_t, and returns 
 * the barycentric coordinates i.e. three components u,v and w in 
 * the base fp such that the (X:Z) coordinates of P+Q are (u-d_P*d_Q*v:w) 
 * and of P-Q are (u+d_P*d_Q*v:w) on E, where d_P = i if P is on the twist of E
 * and 1 otherwise and d_Q = i if Q is on the twist of E and 1 otherwise 
 * (see Lemma 1 and Proposition 2 of https://eprint.iacr.org/2025/736.pdf).
 *
 * @param uvw output barycentric coordinates
 * @param P input point
 * @param Q input point
 * @param E parent elliptic curve (or its twist)
 * @param twistP true if P is on the twist of E, false otherwise
 * @param twistQ true if Q is on the twist of E, false otherwise
 * 
 * Now also works when P == 0 or Q == 0 (but not both).
 */
uint8_t jac_to_xz_add_components(add_components_t *uvw, const jac_point_t *P, const jac_point_t *Q, 
	const ec_curve_t *E, const bool twistP, const bool twistQ);

// xz-aritmetic

/**
 * @brief Doubles a point in (X:Z) coordinates in the Montgomery model.
 *
 * @param Q output point duplication 2P
 * @param P input point to double
 * @param E parent elliptic curve
 * 
 * Complete algorithm that can handle the case P == 0.
 */
void xDBL(xz_point_t *Q, const xz_point_t *P, ec_curve_t *E);

/**
 * @brief Differential addition of points P, Q, P-Q in (X:Z) coordinates
 * in the Montgomery model.
 * 
 * @param R output sum of points P+Q
 * @param P input point to sum
 * @param Q input point to sum
 * @param PQ difference P-Q
 * 
 * Correct in the case P == 0 or Q == 0 but not P == Q (PQ == 0).
 */
void xADD(xz_point_t *R, const xz_point_t *P, const xz_point_t *Q, const xz_point_t *PQ);

/**
 * @brief Simultaneous doubling and differential addition of
 * points in (X:Z) coordinates in the Montgomery model.
 * 
 * @param R output point duplication 2P
 * @param S sum of points P+Q
 * @param P input point to sum
 * @param Q input point to sum
 * @param PQ difference P-Q
 * @param E parent curve
 * 
 * Correct in the case P == 0 or Q == 0 but not P == Q (PQ == 0).
 */
void xDBLADD(xz_point_t *R, xz_point_t *S, const xz_point_t *P, const xz_point_t *Q, 
    const xz_point_t *PQ, ec_curve_t *E);

/** 
 * @brief The Montgomery ladder with (X:Z) Montgomery coordinates
 * 
 * @param Q output point [k]P
 * @param R output point [k+1]P
 * @param P input point
 * @param k scalar multiple
 * @param kbits number of bits of k
 * @param curve parent curve
 * 
 * Uses https://eprint.iacr.org/2017/212, Algorithm 6.
 */
void
xMUL(xz_point_t *Q, xz_point_t *R, const xz_point_t *P, const digit_t *k, const int kbits, 
    ec_curve_t *curve);

/** 
 * @brief The Montgomery ladder with (X:Y:Z) standard projective 
 * coordinates in the Montgomery model
 * 
 * @param Q output point [k]P
 * @param P input point
 * @param k scalar multiple
 * @param kbits number of bits of k
 * @param curve parent curve
 * 
 * Uses xMUL and y reconstruction (revover_y).
 */
void ec_MUL(ec_point_t *Q, const ec_point_t *P, const digit_t *k, const int kbits, 
    ec_curve_t *curve);

/**
 * @brief Computes a full projective point Q=(X_Q:Y_Q:Z_Q) in the Montgomery model
 * from the Montgomery (X:Z) coordinates of Q and P+Q and the full projective 
 * coordinates (X_P:Y_P:Z_P) of P.
 * 
 * @param R output point to recover
 * @param P input point in Montgomery (X:Y:Z) coordinates
 * @param Q input point in Montgomery (X:Z) coordinates
 * @param PQ sum P+Q in Montgomery (X:Z) coorfinates
 * @param E parent curve
 * 
 * Uses Okeya–Sakurai's algorithm (DOI 10.1007/3-540-44709-1_12,
 * see also https://eprint.iacr.org/2017/212, Algorithm 5). Complete
 * algorithm that handles edge cases (P, Q or PQ ==0).
 */
void recover_y(ec_point_t *R, const ec_point_t *P, const xz_point_t *Q, 
    const xz_point_t *PQ, const ec_curve_t *E);

/** @}
 */

/** @defgroup Basis basis generation
 * @{
 */

/**
 * @brief Determines if the x-coordinate defines an Fp-rational point.
 * 
 * @param y2 output value of x^3 + A*x^2 + x (= y^2)
 * @param x input affine x-coordinate
 * @param E Montgomery curve y^2 = x^3 + A*x^2 + x
 * 
 * @return 0xF...F if x defines an Fp-rational point and 0 otherwise.
 * 
 * Also computes y2 to compute square roots afterwards.
 */
uint32_t x_is_on_curve(fp_t *y2, const fp_t *x, const ec_curve_t *E);

uint32_t ec_is_on_curve(const ec_point_t *P, const ec_curve_t *E);

/**
 * @brief Computes the 2-Tate pairing of DISTINCT points P!=Q.
 * 
 * @param P input point
 * @param Q input point
 * @param normalized true if zP==zQ==1 and false otherwise
 * 
 * @return 0 if tate_2(P,Q)==1 and 0xF...F if tate_2(P,Q)==-1
 */
uint32_t tate_2(const ec_point_t *P, const ec_point_t *Q, const bool normalized);

/**
 * @brief Computes the 2-Tate self pairing of a point.
 * 
 * @param P input point
 * @param E parent elliptic curve
 * 
 * @return 0 if tate_2(P,P)==1 and 0xF...F if tate_2(P,P)==-1
 */
uint32_t self_tate_2(ec_point_t *P, const ec_curve_t *E);

/**
 * @brief Computes the the three non-zero 2-torsion points of a Montgomery 
 * curve defined over Fp on the surface of the Volcano as follows:
 * - T0 is the point corresponding to a descending 2-isogeny;
 * - T1 corresponds to a horizontal 2-isogeny that lifts to 
 * an Fp-rational 4-torsion point;
 * - Tm1 corresponds to a horizontal 2-isogeny that does not lift to 
 * an Fp-rational 4-torsion point.
 * 
 * @param T0 output point
 * @param T1 output point
 * @param Tm1 output point
 * @param E parent curve
 */
void two_torsion(ec_point_t *T0, ec_point_t *T1, ec_point_t *Tm1,const ec_curve_t *E);

/**
 * @brief Computes a basis (P,Q) of the 2^{f-1}-torsion of E (where p=c*2^f-1) 
 * such that P is Fp-rational on E, Q is Fp-rational on the twist of E, along
 * with TP = omega(P)-P and TQ = omega(Q) (where omega= (pi+1)/2 defines the
 * orientation of E). P and TP are defined over E and Q and TQ are defined on
 * the twist of E (which is also computed).
 * 
 * @param P output point
 * @param Q output point on the twist
 * @param TP output point
 * @param TQ output point on the twist
 * @param Et output twist of E
 * @param E input curve
 * 
 * The points are returned in standard projective coordinates.
 */
void ec_rational_basis(ec_point_t *P, ec_point_t *Q, ec_point_t *TP, ec_point_t *TQ, 
    ec_curve_t *Et, ec_curve_t *E);

/**
 * @brief See ec_rational_basis function. The returned points are in Jacobian
 * coordinates.
 */
void jac_rational_basis(jac_point_t *P, jac_point_t *Q, jac_point_t *TP, jac_point_t *TQ, 
    ec_curve_t *Et, ec_curve_t *E);

/** @}
 */

#endif


