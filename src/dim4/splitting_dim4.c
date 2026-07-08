#include <isog_dim4.h>
#include <assert.h>

// const for splitting. 
const uint8_t splitting_automorphism[12] = { 0, 1, 2, 3, 0, 3, 2, 1, 0, 1, 3, 2 };


// Given  N, Nb1, A1, A2 are mod 4. 
//Define which change of theta structure matrix we have to use. 
uint8_t compute_theta_position_splitting(const uint8_t N1,const uint8_t Nb1,const uint8_t A1,const uint8_t A2)
{
uint8_t unique_id = ((N1 & 3)<< 6)| ((Nb1 & 3) << 4)| ((A1 & 3) << 2)| (A2 & 3); 
uint8_t ans; 
    switch (unique_id) {
        case 80: case 82: case 88: case 90: case 240: case 242: case 248: case 250: ans = 0; break;
        case 81: case 83: case 121: case 123: case 217: case 219: case 241: case 243: ans = 1; break;
        case 87: case 93: case 117: case 127: case 213: case 223: case 247: case 253: ans = 2; break; 
        default: ans = 255; break; 
    }

    return ans; 
}

#ifndef NDEBUG

void
print_theta_point(char *message, const theta_point_dim4_t *P)
{
    fp_t tmp, tmp2;
    fp_copy(&tmp, &(*P)[0]);
    fp_inv(&tmp);

    for (uint8_t i = 0; i < 16; i++) {
        fp_mul(&tmp2, &(*P)[i], &tmp);
        fp_print(message, &tmp2);
    }
}

static uint8_t
scalprod(uint8_t i, uint8_t j)
{
    uint8_t ret = 0;
    ret ^= (i & 1) & (j & 1);
    ret ^= ((i >> 1) & 1) & ((j >> 1) & 1);
    ret ^= ((i >> 2) & 1) & ((j >> 2) & 1);
    ret ^= ((i >> 3) & 1) & ((j >> 3) & 1);
    return ret;
}

uint16_t
count_even_coeff_dim4(theta_point_dim4_t *theta_null_point)
{
    uint16_t nbr_even = 0;
    fp_t counter, tmp;

    for (uint8_t i = 0; i < 16; i++) {
        for (uint8_t j = 0; j < 16; j++) {
            if (scalprod(i, j) == 0) {
                fp_set_zero(&counter);
                for (uint8_t k = 0; k < 16; k++) {
                    fp_mul(&tmp, &(*theta_null_point)[k], &(*theta_null_point)[k ^ j]);
                    if (scalprod(i, k) == 1) {
                        fp_sub(&counter, &counter, &tmp);
                    } else {
                        fp_add(&counter, &counter, &tmp);
                    }
                }
                if (fp_is_zero(&counter) != (uint32_t)-1) {
                    nbr_even += 1;
                }
            }
        }
    }

    return nbr_even;
}


#endif


void
apply_change_coords_splitting(theta_point_dim4_t *out, const uint8_t splitting_position_theta_change_variable)
{

    uint8_t i;
    theta_point_dim4_t tmp;
    fp_t tmp1[2], tmp2[2];

    for (i = 0; i < 4; i++) {
        fp_add(&tmp1[0], &(*out)[(i << 2)], &(*out)[(i << 2) ^ 1]);
        fp_sub(&tmp1[1], &(*out)[(i << 2)], &(*out)[(i << 2) ^ 1]);
        fp_add(&tmp2[0], &(*out)[(i << 2) ^ 2], &(*out)[(i << 2) ^ 3]);
        fp_sub(&tmp2[1], &(*out)[(i << 2) ^ 2], &(*out)[(i << 2) ^ 3]);

        fp_add(&tmp[(i << 2) ^ (0 ^ i)], &tmp1[0], &tmp2[0]);
        fp_add(&tmp[(i << 2) ^ (1 ^ i)], &tmp1[1], &tmp2[1]);
        fp_sub(&tmp[(i << 2) ^ (2 ^ i)], &tmp1[0], &tmp2[0]);
        fp_sub(&tmp[(i << 2) ^ (3 ^ i)], &tmp1[1], &tmp2[1]);
    }

    for (i = 0; i < 4; i++) {
        fp_copy(&(*out)[(i << 2) ^ splitting_automorphism[(splitting_position_theta_change_variable << 2) ^ 0]],
                &tmp[(i << 2) ^ 0]);
        fp_copy(&(*out)[(i << 2) ^ splitting_automorphism[(splitting_position_theta_change_variable << 2) ^ 1]],
                &tmp[(i << 2) ^ 1]);
        fp_copy(&(*out)[(i << 2) ^ splitting_automorphism[(splitting_position_theta_change_variable << 2) ^ 2]],
                &tmp[(i << 2) ^ 2]);
        fp_copy(&(*out)[(i << 2) ^ splitting_automorphism[(splitting_position_theta_change_variable << 2) ^ 3]],
                &tmp[(i << 2) ^ 3]);
    }
}

uint32_t
is_product(const theta_point_dim4_t *null_point, const theta_point_dim1_t *theta1, const theta_point_dim1_t *theta2)
{
    fp_t null_point_dim2[4], t0, t1;
    theta_point_dim4_t test;
    uint32_t res = -(uint32_t)1;

    for (uint8_t i = 0; i < 4; i++) {
        fp_copy(&null_point_dim2[i], &(*null_point)[i << 2]);
    }

    for (uint8_t i = 0; i < 2; i++) {
        for (uint8_t j = 0; j < 2; j++) {
            for (uint8_t k = 0; k < 4; k++) {
                fp_mul(&test[i + (j << 1) + (k << 2)], &(*theta1)[i], &(*theta2)[j]);
                fp_mul(&test[i + (j << 1) + (k << 2)], &test[i + (j << 1) + (k << 2)], &null_point_dim2[k]);
            }
        }
    }

    for (uint8_t i = 0; i < 16; i++) {
        for (uint8_t j = 0; j < 16; j++) {
            if (i < j) {
                fp_mul(&t0, &(*null_point)[i], &test[j]);
                fp_mul(&t1, &(*null_point)[j], &test[i]);
                res = res & fp_is_equal(&t0, &t1);
            }
        }
    }

    return res;
}

static inline void
theta_point_to_curve(ec_curve_t *E, const theta_point_dim1_t *theta)
{
    // We use a formula that stays over Fp to avoid the twist.
    fp_t t0, t1, t2, t3, inv[2], A1, A2;
    uint32_t sel;

    fp_sqr(&t0, &(*theta)[0]); // t0 = a^2
    fp_sqr(&t0, &t0);          // t0 = a^4
    fp_sqr(&t1, &(*theta)[1]); // t1 = b^2
    fp_sqr(&t1, &t1);          // t1 = b^4

    fp_add(&A1, &t0, &t1);     // A1 = a^4+b^4
    fp_sub(&inv[0], &t1, &t0); // inv[0] = b^4-a^4

    fp_add(&t2, &(*theta)[0], &(*theta)[1]); // t2 = a+b
    fp_sqr(&t2, &t2);                        // t2 = (a+b)^2
    fp_sqr(&t2, &t2);                        // t2 = (a+b)^4
    fp_sub(&t3, &(*theta)[0], &(*theta)[1]); // t3 = a-b
    fp_sqr(&t3, &t3);                        // t3 = (a-b)^2
    fp_sqr(&t3, &t3);                        // t3 = (a-b)^4

    fp_add(&A2, &t2, &t3);     // A2 = (a+b)^4+(a-b)^4
    fp_sub(&inv[1], &t3, &t2); // inv[1] = (a-b)^4-(a+b)^4

    fp_batched_inv(inv, 2);
    fp_mul(&A1, &A1, &inv[0]); // A1 = -(a^4+b^4)/(a^4-b^4)
    fp_add(&A1, &A1, &A1);     // A1 = -2(a^4+b^4)/(a^4-b^4)
    fp_mul(&A2, &A2, &inv[1]); // A2 = -((a+b)^4+(a-b)^4)/((a+b)^4-(a-b)^4)
    fp_add(&A2, &A2, &A2);     // A2 = -2((a+b)^4+(a-b)^4)/((a+b)^4-(a-b)^4)

    fp_set_small(&t0, 2);  // t0 = 2
    fp_sub(&t1, &A1, &t0); // t1 = A1-2
    sel = fp_is_square(&t1);

    // Test of correctness
    // It can be proved that A1-2 and A2-2 cannot be simultaneously squares
    // TODO: set false if it works
    if (true) {
        fp_sub(&t2, &A2, &t0); // t2 = A2-2
        assert((~sel) || (~fp_is_square(&t2)));
    }

    // If (A1-2) is a not a square, A <- A1
    // else ((A1-2) is a square and (A2-2) is not) A <- A2
    fp_select(&E->A, &A1, &A2, sel);

    E->is_ws_model_computed = false;
    E->is_A24_computed = false;
}

void
codomain_to_curves(ec_curve_t *Ea,
                   ec_curve_t *Eabar,
                   const theta_struct_dim4_t *codomain,
                   const uint8_t position_theta_change_variable)
{
    theta_point_dim4_t null_point;
    theta_point_dim1_t theta_Ea, theta_Eabar;

    // print_theta_point("inv_dual_null_point = ", &codomain->inv_dual_null_point);

    invert_theta_point(&null_point, &codomain->inv_dual_null_point);
    inline_hadamard(&null_point);

#ifndef NDEBUG
    assert(count_even_coeff_dim4(&null_point) == 90);
    printf("codomain type before change of basis is good\n");
#endif

    apply_change_coords_splitting(&null_point, position_theta_change_variable);

#ifndef NDEBUG
    assert(count_even_coeff_dim4(&null_point) == 90);

#endif

    // Select i0 such that null_point[4*i0]!=0
    uint8_t i0 = 0;
    uint32_t sel;
    for (uint8_t j = 0; j < 4; j++) {
        sel = ~fp_is_zero(&null_point[j << 2]);
        i0 = i0 ^ ((i0 ^ j) & sel);
    }

    // theta_Ea = (theta_{00i_0}, theta_{10i_0})
    fp_copy(&theta_Ea[0], &null_point[i0 << 2]);
    fp_copy(&theta_Ea[1], &null_point[1 + (i0 << 2)]);
    // theta_Eabar = (theta_{00i_0}, theta_{01i_0})
    fp_copy(&theta_Eabar[0], &null_point[i0 << 2]);
    fp_copy(&theta_Eabar[1], &null_point[2 + (i0 << 2)]);

// To check the theta null point is a product
// TODO: set false if it works
#ifndef NDEBUG
    if (true) {
        assert(is_product(&null_point, &theta_Ea, &theta_Eabar));
    }
#endif
    // Computing Ea
    theta_point_to_curve(Ea, &theta_Ea);
    // Computing Eabar
    theta_point_to_curve(Eabar, &theta_Eabar);
}
