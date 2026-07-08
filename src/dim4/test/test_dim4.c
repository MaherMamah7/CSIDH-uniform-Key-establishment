#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include "test_utils.h"
#include <rng.h>
#include <bench_test_arguments.h>
#include <time.h>

#include <isog_dim4.h>
#include <mp.h>
#include <ec.h>
#include <ec.h>

#include "gluing_test_cases.c"

const uint8_t test_int_to_ker[16] = { 255, 0, 1, 255, 2, 255, 255, 255, 3, 255, 255, 255, 255, 255, 255, 255 };
const uint8_t Sigmas_test[24][4] = {
    { 0, 1, 2, 3 },
    { 0, 1, 3, 2 },
    { 0, 2, 1, 3 },
    { 0, 2, 3, 1 },
    { 0, 3, 1, 2 },
    { 0, 3, 2, 1 },
    { 1, 0, 2, 3 },
    { 1, 0, 3, 2 },
    { 1, 2, 0, 3 },
    { 1, 2, 3, 0 },
    { 1, 3, 0, 2 },
    { 1, 3, 2, 0 },
    { 2, 0, 1, 3 },
    { 2, 0, 3, 1 },
    { 2, 1, 0, 3 },
    { 2, 1, 3, 0 },
    { 2, 3, 0, 1 },
    { 2, 3, 1, 0 },
    { 3, 0, 1, 2 },
    { 3, 0, 2, 1 },
    { 3, 1, 0, 2 },
    { 3, 1, 2, 0 },
    { 3, 2, 0, 1 },
    { 3, 2, 1, 0 }
};

// function used to test the product theta structure
uint8_t
bitreverse4_test(uint8_t x)
{
    x = ((x & 0x1) << 3) | ((x & 0x2) << 1) | ((x & 0x4) >> 1) | ((x & 0x8) >> 3);
    return x & 0xF;
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

static inline bool
theta_is_equal(theta_point_dim4_t *P, theta_point_dim4_t *Q)
{
    fp_t tij, tji;
    for (uint8_t i = 0; i < 16; i++) {
        for (uint8_t j = 0; j < i; j++) {
            fp_mul(&tij, &(*P)[i], &(*Q)[j]);
            fp_mul(&tji, &(*P)[j], &(*Q)[i]);
            if (!(fp_is_equal(&tij, &tji))) {
                return false;
            }
        }
    }
    return true;
}

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

bool
is_correct_gen_isog_test(theta_struct_dim4_t *domain, const theta_point_dim4_t *ker4)
{
    theta_point_dim4_t ker2[4], ker2_bis[4];
    uint8_t s = 0;

    for (uint8_t i = 0; i < 4; i++) {
        theta_double_iter(&ker2[i], &ker4[i], domain, 1);
    }

    for (uint8_t i = 0; i < 4; i++) {
        for (uint8_t j = 0; j < 16; j++) {
            s = scalprod((1 << i), j);
            if (s) {
                fp_neg(&ker2_bis[i][j], &domain->null_point[j]);
            } else {
                fp_copy(&ker2_bis[i][j], &domain->null_point[j]);
            }
        }
        if (!(theta_is_equal(&ker2_bis[i], &ker2[i]))) {
            return false;
        }
    }

    for (uint8_t j = 0; j < 16; j++) {
        // This is where it fails. Zeros everywhere...
        if (fp_is_zero(&domain->null_point[j])) {
            return false;
        }
    }
    return true;
}

int
test_precomputation()
{
    bool test_theta = false;
    bool test_codomain = false;
    uint8_t position_theta;
    uint8_t position_codomain;

    for (int i = 0; i < 8192; i++) {
        basis_components_16_t compo_16;
        compo_16.N = all_valid_gluing_16_cases[i][0];
        compo_16.s1 = all_valid_gluing_16_cases[i][1];
        compo_16.s2 = all_valid_gluing_16_cases[i][2];
        compo_16.s3 = all_valid_gluing_16_cases[i][3];
        compo_16.s4 = all_valid_gluing_16_cases[i][4];

        position_theta = compute_theta_position(compo_16);
        test_theta |= (position_theta == 255);
        position_codomain = compute_codomain_position(compo_16);
        test_codomain |= (position_codomain == 255);
    }

    if (test_theta) {
        printf("Problem with precomputed theta structures\n");
        return 1;
    }
    if (test_codomain) {
        printf("Problem with precomputed codomain\n");
        return 1;
    }
    printf("Gluing precomputation tests............................................ PASSED\n");

    return 0;
}

int
test_gluing_domain_shape(unsigned int Ntest)
{
    for (unsigned int repetition = 0; repetition < Ntest; repetition++) {

        jac_point_t P, Q, TP, TQ, P16, Q16, P8j, P4j, Q8j, Q4j;
        ec_curve_t E, Et;
        ec_curve_init_from_A(&E, &A0);
        ec_twist(&Et, &E);

        uint64_t iter_to_16 = TORSION_EVEN_POWER - 5;

        jac_rational_basis(&P, &Q, &TP, &TQ, &Et, &E);

        jac_dbl_iter(&P16, &P, iter_to_16, &E);
        jac_dbl_iter(&Q16, &Q, iter_to_16, &Et);

        DBL(&P8j, &P16, &E);
        DBL(&P4j, &P8j, &E);
        DBL(&Q8j, &Q16, &Et);
        DBL(&Q4j, &Q8j, &Et);

        theta_point_dim4_t product_theta_point, theta_domain, square_theta_codomain;
        ec_point_t P4, Q4;
        jac_to_ec(&P4, &P4j);
        jac_to_ec(&Q4, &Q4j);

        uint8_t vv[2];
        randombytes(vv, 2);
        uint16_t case_action = ((uint16_t)(vv[0] << 8) | vv[1]) & 8191;
        basis_components_16_t compo_16;
        compo_16.N = all_valid_gluing_16_cases[case_action][0];
        compo_16.s1 = all_valid_gluing_16_cases[case_action][1];
        compo_16.s2 = all_valid_gluing_16_cases[case_action][2];
        compo_16.s3 = all_valid_gluing_16_cases[case_action][3];
        compo_16.s4 = all_valid_gluing_16_cases[case_action][4];

        uint8_t theta_position = compute_theta_position(compo_16);
        uint8_t basis_swap = compute_product_theta_coord(&product_theta_point, &P4, &Q4);
        theta_copy(&theta_domain, &product_theta_point);
        apply_change_coords_theta(&theta_domain, theta_position);

        theta_copy(&square_theta_codomain, &theta_domain);

        // Count number of zeros of the codomain:
        squared(&square_theta_codomain, &square_theta_codomain);
        hadamard(&square_theta_codomain, &square_theta_codomain);

        uint8_t number_zeros = 0;

        for (int i = 0; i < 16; i++) {
            if (fp_is_zero(&square_theta_codomain[i])) {
                number_zeros += 1;
            }
        }

        if (number_zeros != 8) {
            printf("codomain has not enough zeros\n");
            return 1;
        }

        theta_point_dim4_t zero_point;
        quad_point_t zero_vector;
        jac_point_t U;
        jac_point_init(&U);

        for (int i = 0; i < 4; i++) {
            copy_jac_point(&zero_vector[i], &U);
        }

        quad_point_to_dim4_vec(&zero_point, &zero_vector, false);
        dim4_vector_to_theta_point(&zero_point, &product_theta_point, basis_swap, theta_position);

        uint32_t same_point = fp_is_zero(&ZERO);
        for (int i = 0; i < 16; i++) {
            same_point &= fp_is_equal(&zero_point[i], &theta_domain[i]);
        }

        if (!same_point) {
            printf("Problem with codomain and image of zero as a theta point\n");
            return 1;
        }
    }

    printf("Gluing domain shape tests.............................................. PASSED\n");
    return 0;
}

int
test_gluing_product_theta_domain(unsigned int Ntest)
{
    for (unsigned int repetition = 0; repetition < Ntest; repetition++) {

        jac_point_t P, Q, TP, TQ, P16, Q16, P8j, P4j, Q8j, Q4j;
        ec_curve_t E, Et;
        ec_curve_init_from_A(&E, &A0);
        ec_twist(&Et, &E);

        uint64_t iter_to_16 = TORSION_EVEN_POWER - 5;

        jac_rational_basis(&P, &Q, &TP, &TQ, &Et, &E);

        jac_dbl_iter(&P16, &P, iter_to_16, &E);
        jac_dbl_iter(&Q16, &Q, iter_to_16, &Et);

        DBL(&P8j, &P16, &E);
        DBL(&P4j, &P8j, &E);
        DBL(&Q8j, &Q16, &Et);
        DBL(&Q4j, &Q8j, &Et);

        if (rand() % 2 == 1) {
            fp_neg(&Q4j.x, &Q4j.x);
        }

        theta_point_dim4_t product_theta_point;
        ec_point_t P4, Q4;
        jac_to_ec(&P4, &P4j);
        jac_to_ec(&Q4, &Q4j);

        // fp_print("Q4.x = ", &Q4.x);
        // fp_print("Q4.z = ", &Q4.z);

        uint8_t vv[2];
        randombytes(vv, 2);
        //uint16_t case_action = ((uint16_t)vv[0] << 8 | vv[1]) & 8191;
        //basis_components_16_t compo_16;
        //compo_16.N = all_valid_gluing_16_cases[case_action][0];
        //compo_16.s1 = all_valid_gluing_16_cases[case_action][1];
        //compo_16.s2 = all_valid_gluing_16_cases[case_action][2];
        //compo_16.s3 = all_valid_gluing_16_cases[case_action][3];
        //compo_16.s4 = all_valid_gluing_16_cases[case_action][4];

        uint8_t basis_swap = compute_product_theta_coord(&product_theta_point, &P4, &Q4);

        quad_point_t product_T4_points[4], product_S4_points[4];

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                if (i == j) {
                    copy_jac_point(&product_S4_points[i][j], &P4j);
                    copy_jac_point(&product_T4_points[i][j], &Q4j);

                } else {
                    jac_point_init(&product_S4_points[i][j]);
                    jac_point_init(&product_T4_points[i][j]);
                }
            }
        }

        theta_point_dim4_t T4[4], S4[4];

        for (int i = 0; i < 4; i++) {
            quad_point_to_dim4_vec(&T4[i], &product_T4_points[i], true);
            quad_point_to_dim4_vec(&S4[i], &product_S4_points[i], false);
        }

        theta_point_dim4_t tmp;
        for (int i = 0; i < 4; i++) {
            hadamard(&tmp, &T4[i]);
            for (uint8_t j = 0; j < 16; j++) {
                fp_copy(&T4[i][j], &tmp[bitreverse4_test(basis_swap ^ j)]);
            }
            dot_prod(&T4[i], &T4[i], &product_theta_point);

            hadamard(&tmp, &S4[i]);
            for (uint8_t j = 0; j < 16; j++) {
                fp_copy(&S4[i][j], &tmp[bitreverse4_test(basis_swap ^ j)]);
            }
            dot_prod(&S4[i], &S4[i], &product_theta_point);
        }

        uint32_t test_T, test_S;
        test_T = fp_is_zero(&ZERO);
        test_S = fp_is_zero(&ZERO);

        for (int g = 0; g < 4; g++) {
            for (int i = 0; i < 16; i++) {
                test_S &= fp_is_equal(&S4[g][i], &S4[g][i ^ (1 << g)]);

                if ((i >> g) & 1) {
                    test_T &= fp_is_zero(&T4[g][i]);
                }
            }
        }

        if (!test_S) {
            printf("Problem with the S points of the product theta structure\n");
            return 1;
        }
        if (!test_T) {
            printf("Problem with the T points of the product theta structure\n");
            return 1;
        }
    }

    printf("Gluing product theta basis tests ...................................... PASSED\n");
    return 0;
}

int
test_gluing_isotropic_basis(unsigned int Ntest)
{

    for (unsigned int repetition = 0; repetition < Ntest; repetition++) {

        jac_point_t P, Q, TP, TQ, P16, Q16, P8j, P4j, Q8j, Q4j;
        ec_curve_t E, Et;
        ec_curve_init_from_A(&E, &A0);
        ec_twist(&Et, &E);

        uint64_t iter_to_16 = TORSION_EVEN_POWER - 5;

        jac_rational_basis(&P, &Q, &TP, &TQ, &Et, &E);

        jac_dbl_iter(&P16, &P, iter_to_16, &E);
        jac_dbl_iter(&Q16, &Q, iter_to_16, &Et);

        DBL(&P8j, &P16, &E);
        DBL(&P4j, &P8j, &E);
        DBL(&Q8j, &Q16, &Et);
        DBL(&Q4j, &Q8j, &Et);

        theta_point_dim4_t product_theta_point;
        ec_point_t P4, Q4;
        jac_to_ec(&P4, &P4j);
        jac_to_ec(&Q4, &Q4j);

        uint8_t vv[2];
        randombytes(vv, 2);
        uint16_t case_action = ((uint16_t)vv[0] << 8 | vv[1]) & 8191;
        basis_components_16_t compo_16;
        compo_16.N = all_valid_gluing_16_cases[case_action][0];
        compo_16.s1 = all_valid_gluing_16_cases[case_action][1];
        compo_16.s2 = all_valid_gluing_16_cases[case_action][2];
        compo_16.s3 = all_valid_gluing_16_cases[case_action][3];
        compo_16.s4 = all_valid_gluing_16_cases[case_action][4];

        uint8_t theta_position = compute_theta_position(compo_16);
        uint8_t basis_swap = compute_product_theta_coord(&product_theta_point, &P4, &Q4);

        quad_point_t symplectic_basis_4[8];

        construct_symplectic_basis_4(symplectic_basis_4, &P4j, &Q4j, &E, &Et, compo_16, true);

        theta_point_dim4_t T4[4], S4[4];

        bool twist_list[4] = { false, false, true, true };
        for (int i = 0; i < 4; i++) {
            quad_point_to_dim4_vec(&T4[i], &symplectic_basis_4[i], twist_list[i]);
            quad_point_to_dim4_vec(&S4[i], &symplectic_basis_4[4 + i], !twist_list[i]);

            dim4_vector_to_theta_point(&T4[i], &product_theta_point, basis_swap, theta_position);
            dim4_vector_to_theta_point(&S4[i], &product_theta_point, basis_swap, theta_position);
        }

        uint32_t test_T, test_S;
        test_T = fp_is_zero(&ZERO);
        test_S = fp_is_zero(&ZERO);

        for (int g = 0; g < 4; g++) {
            for (int i = 0; i < 16; i++) {
                test_S &= fp_is_equal(&S4[g][i], &S4[g][i ^ (1 << g)]);

                if ((i >> g) & 1) {
                    test_T &= fp_is_zero(&T4[g][i]);
                }
            }
        }

        if (!test_S) {
            printf("Problem with the S points of the product theta structure\n");
            return 1;
        }

        if (!test_T) {
            printf("Problem with the T points of the product theta structure\n");
            return 1;
        }
    }

    printf("Gluing symplectic 4 torsion basis tests ............................... PASSED\n");
    return 0;
}

int
test_gluing_barycentric_lifting(unsigned int Ntest)
{

    for (unsigned int repetition = 0; repetition < Ntest; repetition++) {

        jac_point_t P, Q, TP, TQ, P16, Q16, P8j, P4j, Q8j, Q4j;
        ec_curve_t E, Et;
        ec_curve_init_from_A(&E, &A0);
        ec_twist(&Et, &E);
        ec_compute_A24(&E);
        ec_compute_A24(&Et);

        uint64_t iter_to_16 = TORSION_EVEN_POWER - 5;

        jac_rational_basis(&P, &Q, &TP, &TQ, &Et, &E);

        jac_dbl_iter(&P16, &P, iter_to_16, &E);
        jac_dbl_iter(&Q16, &Q, iter_to_16, &Et);

        DBL(&P8j, &P16, &E);
        DBL(&P4j, &P8j, &E);
        DBL(&Q8j, &Q16, &Et);
        DBL(&Q4j, &Q8j, &Et);

        theta_point_dim4_t product_theta_point, domain_theta;
        ec_point_t P4, Q4;
        jac_to_ec(&P4, &P4j);
        jac_to_ec(&Q4, &Q4j);

        uint8_t vv[2];
        randombytes(vv, 2);
        uint16_t case_action = ((uint16_t)vv[0] << 8 | vv[1]) & 8191;
        basis_components_16_t compo_16;
        compo_16.N = all_valid_gluing_16_cases[case_action][0];
        compo_16.s1 = all_valid_gluing_16_cases[case_action][1];
        compo_16.s2 = all_valid_gluing_16_cases[case_action][2];
        compo_16.s3 = all_valid_gluing_16_cases[case_action][3];
        compo_16.s4 = all_valid_gluing_16_cases[case_action][4];

        uint8_t theta_position = compute_theta_position(compo_16);

        uint8_t basis_swap = compute_product_theta_coord(&product_theta_point, &P4, &Q4);
        theta_copy(&domain_theta, &product_theta_point);
        apply_change_coords_theta(&domain_theta, theta_position);

        quad_point_t symplectic_basis_16[4], symplectic_basis_8[4], symplectic_basis_4[4], symplectic_basis_4_ver[4];

        construct_symplectic_basis_4(symplectic_basis_4_ver, &P4j, &Q4j, &E, &Et, compo_16, false);
        construct_symplectic_basis_16(symplectic_basis_16, &P16, &Q16, &E, &Et, compo_16);

        for (int i = 0; i < 2; i++) {
            quad_point_DBL(&symplectic_basis_8[i], &symplectic_basis_16[i], &E);
            quad_point_DBL(&symplectic_basis_8[i + 2], &symplectic_basis_16[i + 2], &Et);
            quad_point_DBL(&symplectic_basis_4[i], &symplectic_basis_8[i], &E);
            quad_point_DBL(&symplectic_basis_4[i + 2], &symplectic_basis_8[i + 2], &Et);
        }

        uint32_t test_quad_DBL;
        test_quad_DBL = fp_is_zero(&ZERO);

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                test_quad_DBL &= jac_is_equal(&(symplectic_basis_4)[i][j], &symplectic_basis_4_ver[i][j]);
            }
        }

        if (!test_quad_DBL) {
            printf("Problem with the construction of 16 torsion basis\n");
            return 1;
        }

        // Test on curve.
        quad_point_t P1, Q1, P1_p_Q1, P1_m_Q1;
        for (int k = 0; k < 4; k++) {
            // TEST T0 + 2T0
            for (int i = 0; i < 4; i++) {
                copy_jac_point(&P1[i], &symplectic_basis_16[(k & 2) >> 1][i]);
                copy_jac_point(&Q1[i], &symplectic_basis_16[k & 1][i]);
            }
            quad_point_DBL(&Q1, &Q1, &E); // To avoid P +- P.
            quad_point_ADD(&P1_p_Q1, &P1, &Q1, &E);
            for (int i = 0; i < 4; i++) {
                jac_neg(&Q1[i], &Q1[i]);
            }
            quad_point_ADD(&P1_m_Q1, &P1, &Q1, &E);
            for (int i = 0; i < 4; i++) {
                jac_neg(&Q1[i], &Q1[i]);
            }

            theta_point_dim4_t real_vec, im_vec;
            ec_couplepoints_to_dim4_vectors(&real_vec, &im_vec, &P1, &Q1, &E, false, false);

            theta_point_dim4_t P_p_Q_bary, P_p_Q_stand, P_m_Q_bary, P_m_Q_stand;
            for (int j = 0; j < 16; j++) {
                fp_sub(&P_p_Q_bary[j], &real_vec[j], &im_vec[j]);
                fp_add(&P_m_Q_bary[j], &real_vec[j], &im_vec[j]);
            }

            quad_point_to_dim4_vec(&P_p_Q_stand, &P1_p_Q1, false);
            quad_point_to_dim4_vec(&P_m_Q_stand, &P1_m_Q1, false);

            uint32_t plus = fp_is_zero(&ZERO);
            uint32_t minus = fp_is_zero(&ZERO);

            fp_t inv_el_bary, inv_el_stand;
            fp_copy(&inv_el_bary, &P_p_Q_bary[0]);
            fp_copy(&inv_el_stand, &P_p_Q_stand[0]);
            fp_inv(&inv_el_stand);
            fp_inv(&inv_el_bary);
            for (int j = 0; j < 16; j++) {
                fp_mul(&P_p_Q_bary[j], &P_p_Q_bary[j], &inv_el_bary);
                fp_mul(&P_p_Q_stand[j], &P_p_Q_stand[j], &inv_el_stand);
                //printf("j = %d\n", j);
                //fp_print("P + Q bary = ", &P_p_Q_bary[j]);
                //fp_print("P + Q stand = ", &P_p_Q_stand[j]);
                //printf("plus equal %d\n", fp_is_equal(&P_p_Q_bary[j],&P_p_Q_stand[j]));
                plus &= fp_is_equal(&P_p_Q_bary[j], &P_p_Q_stand[j]);
            }
            fp_copy(&inv_el_bary, &P_m_Q_bary[0]);
            fp_copy(&inv_el_stand, &P_m_Q_stand[0]);
            fp_inv(&inv_el_stand);
            fp_inv(&inv_el_bary);
            for (int j = 0; j < 16; j++) {
                fp_mul(&P_m_Q_bary[j], &P_m_Q_bary[j], &inv_el_bary);
                fp_mul(&P_m_Q_stand[j], &P_m_Q_stand[j], &inv_el_stand);
                //printf("j = %d\n", j);
                //fp_print("P - Q bary = ", &P_m_Q_bary[j]);
                //fp_print("P - Q stand = ", &P_m_Q_stand[j]);
                //printf("minus equal %d\n", fp_is_equal(&P_m_Q_stand[j],&P_m_Q_bary[j]));
                minus &= fp_is_equal(&P_m_Q_bary[j], &P_m_Q_stand[j]);
            }

            if (!plus) {
                printf("Problem with the T16_%d + T8_%d points\n", (k & 2) >> 1, k & 1);
                return 1;
            }
            if (!minus) {
                printf("Problem with the T16_%d - T8_%d points\n", (k & 2) >> 1, k & 1);
                return 1;
            }
        }

        // Test on twist.
        for (int k = 0; k < 4; k++) {
            for (int i = 0; i < 4; i++) {
                copy_jac_point(&P1[i], &symplectic_basis_16[2 + ((k & 2) >> 1)][i]);
                copy_jac_point(&Q1[i], &symplectic_basis_16[2 + (k & 1)][i]);
            }
            quad_point_DBL(&Q1, &Q1, &E); // To avoid P +- P.
            quad_point_ADD(&P1_p_Q1, &P1, &Q1, &Et);
            for (int i = 0; i < 4; i++) {
                jac_neg(&Q1[i], &Q1[i]);
            }
            quad_point_ADD(&P1_m_Q1, &P1, &Q1, &Et);
            for (int i = 0; i < 4; i++) {
                jac_neg(&Q1[i], &Q1[i]);
            }

            theta_point_dim4_t real_vec, im_vec;
            ec_couplepoints_to_dim4_vectors(&real_vec, &im_vec, &P1, &Q1, &E, true, true);

            theta_point_dim4_t P_p_Q_bary, P_p_Q_stand, P_m_Q_bary, P_m_Q_stand;
            for (int j = 0; j < 16; j++) {
                fp_add(&P_p_Q_bary[j], &real_vec[j], &im_vec[j]);
                fp_sub(&P_m_Q_bary[j], &real_vec[j], &im_vec[j]);
            }

            quad_point_to_dim4_vec(&P_p_Q_stand, &P1_p_Q1, true);
            quad_point_to_dim4_vec(&P_m_Q_stand, &P1_m_Q1, true);

            uint32_t plus = fp_is_zero(&ZERO);
            uint32_t minus = fp_is_zero(&ZERO);

            fp_t inv_el_bary, inv_el_stand;
            fp_copy(&inv_el_bary, &P_p_Q_bary[0]);
            fp_copy(&inv_el_stand, &P_p_Q_stand[0]);
            fp_inv(&inv_el_stand);
            fp_inv(&inv_el_bary);
            for (int j = 0; j < 16; j++) {
                fp_mul(&P_p_Q_bary[j], &P_p_Q_bary[j], &inv_el_bary);
                fp_mul(&P_p_Q_stand[j], &P_p_Q_stand[j], &inv_el_stand);
                //printf("j = %d\n", j);
                //fp_print("P + Q bary = ", &P_p_Q_bary[j]);
                //fp_print("P + Q stand = ", &P_p_Q_stand[j]);
                //printf("plus equal %d\n", fp_is_equal(&P_p_Q_bary[j],&P_p_Q_stand[j]));
                plus &= fp_is_equal(&P_p_Q_bary[j], &P_p_Q_stand[j]);
            }
            fp_copy(&inv_el_bary, &P_m_Q_bary[0]);
            fp_copy(&inv_el_stand, &P_m_Q_stand[0]);
            fp_inv(&inv_el_stand);
            fp_inv(&inv_el_bary);
            for (int j = 0; j < 16; j++) {
                fp_mul(&P_m_Q_bary[j], &P_m_Q_bary[j], &inv_el_bary);
                fp_mul(&P_m_Q_stand[j], &P_m_Q_stand[j], &inv_el_stand);
                //fp_print("P - Q bary = ", &P_m_Q_bary[j]);
                //fp_print("P - Q stand = ", &P_m_Q_stand[j]);
                //printf("minus equal %d\n", fp_is_equal(&P_m_Q_stand[j],&P_m_Q_bary[j]));
                minus &= fp_is_equal(&P_m_Q_bary[j], &P_m_Q_stand[j]);
            }

            if (!plus) {
                printf("Problem with the T16_%d + T8_%d points\n", 2 + ((k & 2) >> 1), 2 + (k & 1));
                return 1;
            }
            if (!minus) {
                printf("Problem with the T16_%d - T8_%d points\n", 2 + ((k & 2) >> 1), 2 + (k & 1));
                return 1;
            }
        }



        // Test double on curve
        for (int k = 0; k < 4; k++) {
            // TEST T0 + 2T0
            for (int i = 0; i < 4; i++) {
                copy_jac_point(&P1[i], &symplectic_basis_16[(k & 2) >> 1][i]);
                copy_jac_point(&Q1[i], &symplectic_basis_16[k & 1][i]);
            }
            //quad_point_DBL(&Q1, &Q1, &E); // To avoid P +- P.
            quad_point_ADD(&P1_p_Q1, &P1, &Q1, &E);
            for (int i = 0; i < 4; i++) {
                jac_neg(&Q1[i], &Q1[i]);
            }
            quad_point_ADD(&P1_m_Q1, &P1, &Q1, &E);
            for (int i = 0; i < 4; i++) {
                jac_neg(&Q1[i], &Q1[i]);
            }

            theta_point_dim4_t real_vec, im_vec;
            ec_couplepoints_to_dim4_vectors(&real_vec, &im_vec, &P1, &Q1, &E, false, false);

            theta_point_dim4_t P_p_Q_bary, P_p_Q_stand, P_m_Q_bary, P_m_Q_stand;
            for (int j = 0; j < 16; j++) {
                fp_sub(&P_p_Q_bary[j], &real_vec[j], &im_vec[j]);
                fp_add(&P_m_Q_bary[j], &real_vec[j], &im_vec[j]);
            }

            quad_point_to_dim4_vec(&P_p_Q_stand, &P1_p_Q1, false);
            quad_point_to_dim4_vec(&P_m_Q_stand, &P1_m_Q1, false);

            uint32_t plus = fp_is_zero(&ZERO);
            uint32_t minus = fp_is_zero(&ZERO);

            fp_t inv_el_bary, inv_el_stand;
            fp_copy(&inv_el_bary, &P_p_Q_bary[0]);
            fp_copy(&inv_el_stand, &P_p_Q_stand[0]);
            fp_inv(&inv_el_stand);
            fp_inv(&inv_el_bary);
            //plus &= -(uint32_t)theta_is_equal(&P_p_Q_bary, &P_p_Q_bary);

            for (int j = 0; j < 16; j++) {
                fp_mul(&P_p_Q_bary[j], &P_p_Q_bary[j], &inv_el_bary);
                fp_mul(&P_p_Q_stand[j], &P_p_Q_stand[j], &inv_el_stand);
                //printf("j = %d\n", j);
                //fp_print("P + Q bary = ", &P_p_Q_bary[j]);
                //fp_print("P + Q stand = ", &P_p_Q_stand[j]);
                //printf("plus equal %d\n", fp_is_equal(&P_p_Q_bary[j],&P_p_Q_stand[j]));
                plus &= -(uint32_t)fp_is_equal(&P_p_Q_bary[j], &P_p_Q_stand[j]);
            }
            fp_copy(&inv_el_bary, &P_m_Q_bary[0]);
            fp_copy(&inv_el_stand, &P_m_Q_stand[0]);
            fp_inv(&inv_el_stand);
            fp_inv(&inv_el_bary);
            //minus &= -(uint32_t)theta_is_equal(&P_m_Q_bary, &P_m_Q_bary);

            for (int j = 0; j < 16; j++) {
                fp_mul(&P_m_Q_bary[j], &P_m_Q_bary[j], &inv_el_bary);
                fp_mul(&P_m_Q_stand[j], &P_m_Q_stand[j], &inv_el_stand);
                //printf("j = %d\n", j);
                //fp_print("P - Q bary = ", &P_m_Q_bary[j]);
                //fp_print("P - Q stand = ", &P_m_Q_stand[j]);
                //printf("minus equal %d\n", fp_is_equal(&P_m_Q_stand[j],&P_m_Q_bary[j]));
                minus &= -(uint32_t)fp_is_equal(&P_m_Q_bary[j], &P_m_Q_stand[j]);
            }

            if (!plus) {
                printf("Problem with the T16_%d + T16_%d points\n", (k & 2) >> 1, k & 1);
                return 1;
            }
            if (!minus) {
                printf("Problem with the T16_%d - T16_%d points\n", (k & 2) >> 1, k & 1);
                return 1;
            }
        }

        
        // Test double on twist
        //printf("Test double on twist\n");
        for (int k = 0; k < 2; k++) {
            // TEST T0 + T0
            for (int i = 0; i < 4; i++) {
                copy_jac_point(&P1[i], &symplectic_basis_16[2 + k][i]);
                copy_jac_point(&Q1[i], &symplectic_basis_16[2 + k][i]);
            }
            //quad_point_DBL(&Q1, &Q1, &E); // To avoid P +- P.
            quad_point_ADD(&P1_p_Q1, &P1, &Q1, &Et);
            for (int i = 0; i < 4; i++) {
                jac_neg(&Q1[i], &Q1[i]);
            }
            quad_point_ADD(&P1_m_Q1, &P1, &Q1, &Et);
            for (int i = 0; i < 4; i++) {
                jac_neg(&Q1[i], &Q1[i]);
            }

            theta_point_dim4_t real_vec, im_vec;
            ec_couplepoints_to_dim4_vectors(&real_vec, &im_vec, &P1, &Q1, &E, true, true);

            theta_point_dim4_t P_p_Q_bary, P_p_Q_stand, P_m_Q_bary, P_m_Q_stand;
            for (int j = 0; j < 16; j++) {
                fp_sub(&P_p_Q_bary[j], &real_vec[j], &im_vec[j]);
                fp_add(&P_m_Q_bary[j], &real_vec[j], &im_vec[j]);
            }

            quad_point_to_dim4_vec(&P_p_Q_stand, &P1_p_Q1, true);
            quad_point_to_dim4_vec(&P_m_Q_stand, &P1_m_Q1, true);

            uint32_t plus = fp_is_zero(&ZERO);
            uint32_t minus = fp_is_zero(&ZERO);

            fp_t inv_el_bary, inv_el_stand;
            fp_copy(&inv_el_bary, &P_p_Q_bary[0]);
            fp_copy(&inv_el_stand, &P_p_Q_stand[0]);
            fp_inv(&inv_el_stand);
            fp_inv(&inv_el_bary);
            plus &= -(uint32_t)theta_is_equal(&P_p_Q_bary, &P_p_Q_bary);

            fp_copy(&inv_el_bary, &P_m_Q_bary[0]);
            fp_copy(&inv_el_stand, &P_m_Q_stand[0]);
            fp_inv(&inv_el_stand);
            fp_inv(&inv_el_bary);
            minus &= -(uint32_t)theta_is_equal(&P_m_Q_bary, &P_m_Q_bary);

            if (!plus) {
                printf("Problem with the T16_%d + T16_%d points\n", 2 + k , 2 + k);
                return 1;
            }
            if (!minus) {
                printf("Problem with the T16_%d - T16_%d points\n", 2 + k , 2 + k);
                return 1;
            }
        }
        

        // Test P in E; Q in Et
        for (int k = 0; k < 4; k++) {
            for (int i = 0; i < 4; i++) {
                copy_jac_point(&P1[i], &symplectic_basis_16[((k & 2) >> 1)][i]);
                copy_jac_point(&Q1[i], &symplectic_basis_16[2 + (k & 1)][i]);
            }
            quad_point_DBL(&Q1, &Q1, &Et);
            quad_point_DBL(&Q1, &Q1, &Et);
            quad_point_DBL(&Q1, &Q1, &Et);
            theta_point_dim4_t P_p_Q_bary, PP;
            couplepoint_to_product_of_alt_sum_theta_point(
                &P_p_Q_bary, &P1, &Q1, &E, &product_theta_point, false, true, basis_swap, theta_position);

            quad_point_to_dim4_vec(&PP, &P1, false);
            dim4_vector_to_theta_point(&PP, &product_theta_point, basis_swap, theta_position);
            squared(&PP, &PP);
            // quad_point_to_dim4_vec(&P_m_Q_stand,&P1_m_Q1,true);

            uint32_t plus = fp_is_zero(&ZERO);

            fp_t inv_el_bary, inv_el_PP;
            fp_copy(&inv_el_bary, &P_p_Q_bary[0]);
            fp_copy(&inv_el_PP, &PP[0]);
            fp_inv(&inv_el_PP);
            fp_inv(&inv_el_bary);
            for (int j = 0; j < 16; j++) {
                fp_mul(&P_p_Q_bary[j], &P_p_Q_bary[j], &inv_el_bary);
                fp_mul(&PP[j], &PP[j], &inv_el_PP);
                plus &= fp_is_equal(&P_p_Q_bary[j], &PP[j]);
            }

            if (!plus) {
                printf("Problem with the T16_%d + T2_%d points\n", ((k & 2) >> 1), 2 + (k & 1));
                return 1;
            }
        }
        // Test P in Et; Q in E
        for (int k = 0; k < 4; k++) {
            for (int i = 0; i < 4; i++) {
                copy_jac_point(&P1[i], &symplectic_basis_16[2 + ((k & 2) >> 1)][i]);
                copy_jac_point(&Q1[i], &symplectic_basis_16[(k & 1)][i]);
            }
            quad_point_DBL(&Q1, &Q1, &E);
            quad_point_DBL(&Q1, &Q1, &E);
            quad_point_DBL(&Q1, &Q1, &E);
            theta_point_dim4_t P_p_Q_bary, PP;
            couplepoint_to_product_of_alt_sum_theta_point(
                &P_p_Q_bary, &P1, &Q1, &E, &product_theta_point, true, false, basis_swap, theta_position);

            quad_point_to_dim4_vec(&PP, &P1, true);
            dim4_vector_to_theta_point(&PP, &product_theta_point, basis_swap, theta_position);
            squared(&PP, &PP);
            // quad_point_to_dim4_vec(&P_m_Q_stand,&P1_m_Q1,true);

            uint32_t plus = fp_is_zero(&ZERO);

            fp_t inv_el_bary, inv_el_PP;
            fp_copy(&inv_el_bary, &P_p_Q_bary[0]);
            fp_copy(&inv_el_PP, &PP[0]);
            fp_inv(&inv_el_PP);
            fp_inv(&inv_el_bary);
            for (int j = 0; j < 16; j++) {
                fp_mul(&P_p_Q_bary[j], &P_p_Q_bary[j], &inv_el_bary);
                fp_mul(&PP[j], &PP[j], &inv_el_PP);
                plus &= fp_is_equal(&P_p_Q_bary[j], &PP[j]);
            }

            if (!plus) {
                printf("Problem with the T16_%d + T2_%d points\n", 2 + ((k & 2) >> 1), (k & 1));
                return 1;
            }
        }
    }

    printf("Gluing barycentric coeff tests ........................................ PASSED\n");
    return 0;
}

int
test_gluing_16_torsion_points(unsigned int Ntest)
{

    // for(int repetition =0; repetition < Ntest; repetition++)
    //{

    jac_point_t P, Q, TP, TQ, P16, Q16, P8j, P4j, Q8j, Q4j;
    ec_curve_t E, Et;
    ec_curve_init_from_A(&E, &A0);
    ec_twist(&Et, &E);

    uint64_t iter_to_16 = TORSION_EVEN_POWER - 5;

    jac_rational_basis(&P, &Q, &TP, &TQ, &Et, &E);

    jac_dbl_iter(&P16, &P, iter_to_16, &E);
    jac_dbl_iter(&Q16, &Q, iter_to_16, &Et);

    DBL(&P8j, &P16, &E);
    DBL(&P4j, &P8j, &E);
    DBL(&Q8j, &Q16, &Et);
    DBL(&Q4j, &Q8j, &Et);

    theta_point_dim4_t product_theta_point, domain_theta;
    ec_point_t P4, Q4;
    jac_to_ec(&P4, &P4j);
    jac_to_ec(&Q4, &Q4j);

    uint16_t case_action = Ntest & 8191;
    basis_components_16_t compo_16;
    compo_16.N = all_valid_gluing_16_cases[case_action][0];
    compo_16.s1 = all_valid_gluing_16_cases[case_action][1];
    compo_16.s2 = all_valid_gluing_16_cases[case_action][2];
    compo_16.s3 = all_valid_gluing_16_cases[case_action][3];
    compo_16.s4 = all_valid_gluing_16_cases[case_action][4];

    uint8_t theta_position = compute_theta_position(compo_16);

    uint8_t basis_swap = compute_product_theta_coord(&product_theta_point, &P4, &Q4);
    theta_copy(&domain_theta, &product_theta_point);
    apply_change_coords_theta(&domain_theta, theta_position);

    quad_point_t symplectic_basis_16[4], symplectic_basis_8[4];
    construct_symplectic_basis_16(symplectic_basis_16, &P16, &Q16, &E, &Et, compo_16);

    for (int i = 0; i < 2; i++) {
        quad_point_DBL(&symplectic_basis_8[i], &symplectic_basis_16[i], &E);
        quad_point_DBL(&symplectic_basis_8[i + 2], &symplectic_basis_16[i + 2], &Et);
    }

    theta_point_dim4_t domain;
    theta_copy(&domain, &product_theta_point);
    apply_change_coords_theta(&domain, theta_position);

    // Test on curve.
    quad_point_t PP;
    theta_point_dim4_t HSK;

    for (int k = 0; k < 4; k++) {
        for (int i = 0; i < 4; i++) {
            copy_jac_point(&PP[i], &symplectic_basis_8[k][i]);
        }
        quad_point_to_dim4_vec(&HSK, &PP, k >> 1);
        dim4_vector_to_theta_point(&HSK, &product_theta_point, basis_swap, theta_position);

        dot_prod(&HSK, &HSK, &domain);
        hadamard(&HSK, &HSK);

        printf("This is  T16_%d \n", k);

        for (int i = 0; i < 16; i++) {
            fp_print("16 torsion points", &HSK[i]);
        }
        printf("\n\n");
    }

    quad_point_ADD(&PP, &symplectic_basis_8[0], &symplectic_basis_8[1], &E);
    quad_point_to_dim4_vec(&HSK, &PP, false);
    dim4_vector_to_theta_point(&HSK, &product_theta_point, basis_swap, theta_position);

    dot_prod(&HSK, &HSK, &domain);
    hadamard(&HSK, &HSK);

    printf("This is T16_3 \n");

    for (int i = 0; i < 16; i++) {
        fp_print("16 torsion points ", &HSK[i]);
    }
    printf("\n\n");

    quad_point_ADD(&PP, &symplectic_basis_8[2], &symplectic_basis_8[3], &Et);
    quad_point_to_dim4_vec(&HSK, &PP, true);
    dim4_vector_to_theta_point(&HSK, &product_theta_point, basis_swap, theta_position);

    dot_prod(&HSK, &HSK, &domain);
    hadamard(&HSK, &HSK);

    printf("This is T16_12 \n");

    for (int i = 0; i < 16; i++) {
        fp_print("16 torsion points ", &HSK[i]);
    }
    printf("\n\n");

    printf("Gluing 16 torsion points.................................... PASSED\n");
    return 0;
}

int
test_gluing_compute(unsigned int Ntest)
{

    for (unsigned int repetition = 0; repetition < Ntest; repetition++) {

        jac_point_t P, Q, TP, TQ, P16, Q16;
        ec_curve_t EC_domain[2];
        ec_curve_init_from_A(&EC_domain[0], &A0);
        ec_twist(&EC_domain[1], &EC_domain[0]);

        uint64_t iter_to_16 = TORSION_EVEN_POWER - 5;

        jac_rational_basis(&P, &Q, &TP, &TQ, &EC_domain[1], &EC_domain[0]);

        jac_dbl_iter(&P16, &P, iter_to_16, &EC_domain[0]);
        jac_dbl_iter(&Q16, &Q, iter_to_16, &EC_domain[1]);

        uint8_t vv[2];
        randombytes(vv, 2);
        uint16_t case_action = ((uint16_t)vv[0] << 8 | vv[1]) & 8191;
        basis_components_16_t compo_16;
        compo_16.N = all_valid_gluing_16_cases[case_action][0];
        compo_16.s1 = all_valid_gluing_16_cases[case_action][1];
        compo_16.s2 = all_valid_gluing_16_cases[case_action][2];
        compo_16.s3 = all_valid_gluing_16_cases[case_action][3];
        compo_16.s4 = all_valid_gluing_16_cases[case_action][4];

        quad_point_t T16[4];

        construct_symplectic_basis_16(T16, &P16, &Q16, &EC_domain[0], &EC_domain[1], compo_16);

        gluing_isog_dim4_t gluing_iso;

        uint8_t position_theta = compute_theta_position(compo_16);
        uint8_t position_codomain = compute_codomain_position(compo_16);

        gluing_compute(&gluing_iso, T16, &P16, &Q16, EC_domain, position_codomain, position_theta);
        theta_point_dim4_t auxiliary_point_theta, domain;
        theta_copy(&domain, &gluing_iso.theta_product_domain);
        apply_change_coords_theta(&domain, gluing_iso.position_theta_structure);

        quad_point_t auxiliary_point;

        for (int i = 0; i < 4; i++) {
            copy_jac_point(&auxiliary_point[i], &gluing_iso.auxiliary_point[i]);
        }

        quad_point_DBL(&auxiliary_point, &auxiliary_point, &EC_domain[(position_codomain > 1)]);
        quad_point_to_dim4_vec(&auxiliary_point_theta, &auxiliary_point, (position_codomain > 1));
        dim4_vector_to_theta_point(
            &auxiliary_point_theta, &gluing_iso.theta_product_domain, gluing_iso.basis_swap, position_theta);
        dot_prod(&auxiliary_point_theta, &auxiliary_point_theta, &domain);

        hadamard(&auxiliary_point_theta, &auxiliary_point_theta);

        uint32_t are_the_same = fp_is_zero(&ZERO);
        fp_t valu, tempo;
        fp_sqr(&valu, &gluing_iso.inv_codomain_auxiliary_point[0]);
        fp_mul(&valu, &valu, &auxiliary_point_theta[0]);

        for (int i = 0; i < 16; i++) {
            fp_sqr(&tempo, &gluing_iso.inv_codomain_auxiliary_point[i]);
            fp_mul(&tempo, &tempo, &auxiliary_point_theta[i]);
            are_the_same &= fp_is_equal(&valu, &tempo);
        }

        if (!are_the_same) {
            printf("Problem during the codomain computation");
            return 1;
        }
    }

    printf("Gluing compute tests .................................................. PASSED\n");
    return 0;
}

int
test_gluing_Eval(unsigned int Ntest)
{

    for (unsigned int repetition = 0; repetition < Ntest; repetition++) {

        jac_point_t P, Q, TP, TQ, P16, Q16;
        ec_curve_t EC_domain[2];
        ec_curve_init_from_A(&EC_domain[0], &A0);
        ec_twist(&EC_domain[1], &EC_domain[0]);

        uint64_t iter_to_16 = TORSION_EVEN_POWER - 5;

        jac_rational_basis(&P, &Q, &TP, &TQ, &EC_domain[1], &EC_domain[0]);

        jac_dbl_iter(&P16, &P, iter_to_16, &EC_domain[0]);
        jac_dbl_iter(&Q16, &Q, iter_to_16, &EC_domain[1]);

        uint8_t vv[2];
        randombytes(vv, 2);
        uint16_t case_action = ((uint16_t)(vv[0] << 8) | vv[1]) & 8191;
        basis_components_16_t compo_16;
        compo_16.N = all_valid_gluing_16_cases[case_action][0];
        compo_16.s1 = all_valid_gluing_16_cases[case_action][1];
        compo_16.s2 = all_valid_gluing_16_cases[case_action][2];
        compo_16.s3 = all_valid_gluing_16_cases[case_action][3];
        compo_16.s4 = all_valid_gluing_16_cases[case_action][4];

        quad_point_t T16[4];

        construct_symplectic_basis_16(T16, &P16, &Q16, &EC_domain[0], &EC_domain[1], compo_16);

        gluing_isog_dim4_t gluing_iso;

        uint8_t position_theta = compute_theta_position(compo_16);
        uint8_t position_codomain = compute_codomain_position(compo_16);

        gluing_compute(&gluing_iso, T16, &P16, &Q16, EC_domain, position_codomain, position_theta);

        theta_point_dim4_t domain;
        theta_copy(&domain, &gluing_iso.theta_product_domain);
        apply_change_coords_theta(&domain, gluing_iso.position_theta_structure);

        bool twitt = (position_codomain > 1);

        for (int k = 0; k < 2; k++) {
            quad_point_t TK;
            theta_point_dim4_t TK_square_theta;
            quad_point_DBL(&TK, &T16[2 * (!twitt) + k], &EC_domain[!twitt]);
            quad_point_to_dim4_vec(&TK_square_theta, &TK, !twitt);
            dim4_vector_to_theta_point(&TK_square_theta,
                                       &gluing_iso.theta_product_domain,
                                       gluing_iso.basis_swap,
                                       gluing_iso.position_theta_structure);

            dot_prod(&TK_square_theta, &TK_square_theta, &domain);
            hadamard(&TK_square_theta, &TK_square_theta);

            theta_point_dim4_t TK_theta;

            gluing_eval(&TK_theta, &T16[2 * (!twitt) + k], &gluing_iso, EC_domain, !twitt);

            uint32_t are_the_same = fp_is_zero(&ZERO);
            fp_t valu, tempo;
            fp_sqr(&valu, &TK_theta[0]);
            fp_inv(&valu);
            fp_mul(&valu, &valu, &TK_square_theta[0]);

            for (int i = 0; i < 16; i++) {
                fp_sqr(&tempo, &TK_theta[i]);
                fp_inv(&tempo);
                fp_mul(&tempo, &tempo, &TK_square_theta[i]);
                // fp_print("tempo is ",&tempo);
                are_the_same &= (fp_is_equal(&valu, &tempo) || fp_is_zero(&tempo));
                // printf(" Are they equal %d\n", are_the_same );
            }

            if (!are_the_same) {
                printf("Problem with generic gluing evaluation \n");
                return 1;
            }
        }

        theta_point_dim4_t TT[2];
        gluing_eval_special_basis(&TT[0], &TT[1], &T16[2 * twitt], &T16[2 * twitt + 1], &gluing_iso, &EC_domain[twitt]);

        for (int k = 0; k < 2; k++) {
            quad_point_t TK;
            theta_point_dim4_t TK_square_theta;
            quad_point_DBL(&TK, &T16[2 * (twitt) + k], &EC_domain[twitt]);
            quad_point_to_dim4_vec(&TK_square_theta, &TK, twitt);
            dim4_vector_to_theta_point(&TK_square_theta,
                                       &gluing_iso.theta_product_domain,
                                       gluing_iso.basis_swap,
                                       gluing_iso.position_theta_structure);

            dot_prod(&TK_square_theta, &TK_square_theta, &domain);
            hadamard(&TK_square_theta, &TK_square_theta);

            uint32_t are_the_same = fp_is_zero(&ZERO);
            fp_t valu, tempo;
            fp_sqr(&valu, &TT[k][0]);
            fp_inv(&valu);
            fp_mul(&valu, &valu, &TK_square_theta[0]);

            for (int i = 0; i < 16; i++) {
                fp_sqr(&tempo, &TT[k][i]);
                fp_inv(&tempo);
                fp_mul(&tempo, &tempo, &TK_square_theta[i]);
                // fp_print("tempo is ",&tempo);
                are_the_same &= (fp_is_equal(&valu, &tempo) || fp_is_zero(&tempo));
                // printf(" Are they equal %d\n", are_the_same );
            }

            if (!are_the_same) {
                printf("Problem with special gluing evaluation \n");
                return 1;
            }
        }
    }

    printf("Gluing evaluation tests ............................................... PASSED\n");
    return 0;
}

int
test_hypercube_automorphisms()
{
    // Test gray-transform.
    aut_t aut;
    init_aut(&aut);

    uint8_t bin, pos1, pos2;

    // testing gray-code.
    for (uint8_t i = 1; i < 16; i++) {
        pos1 = binary_to_gray_pos(i - 1);
        pos2 = binary_to_gray_pos(i);

        if (test_int_to_ker[pos1 ^ pos2] == 255) {
            printf("Problem with grey-code implementation (dif not power 2) \n");
            return 1;
        }

        bin = gray_pos_to_binary(pos1);

        if (bin != i - 1) {
            printf("Problem with grey-code implementation (not inverse) \n");
            return 1;
        }
    }

    // testing Gray-code + standard automorphism.
    for (uint8_t j = 0; j < 24; j++) {
        for (uint8_t i = 0; i < 4; i++) {
            aut.sigma[i] = Sigmas_test[j][i];
            // printf("automorphism %d --> %d\n", i, aut.sigma[i]);
        }

        for (uint8_t i = 1; i < 16; i++) {
            pos1 = binary_to_pos(i - 1, &aut);
            pos2 = binary_to_pos(i, &aut);
            if (test_int_to_ker[pos1 ^ pos2] == 255) {
                printf("Problem with grey + automorphism implementation (dif not power 2) \n");
                return 1;
            }
            bin = pos_to_binary(pos1, &aut);
            if (bin != i - 1) {
                printf("Problem with grey + automorphism implementation (not inverse) \n");
                return 1;
            }
        }
    }

    // testing Gray-code + standard automorphism.
    extra_aut_t extra;

    for (uint8_t mas = 0; mas < 3; mas += 2) {
        init_extra_aut(&extra, (3 << mas));
        uint8_t paths[2] = { ((3 << mas) & 1), 2 };

        for (uint8_t r = 0; r < 2; r++) {
            extra.path = paths[r];
            for (uint8_t v = 0; v < 4; v++) {
                for (uint8_t sgm1 = 0; sgm1 < 2; sgm1++) {
                    extra.sigma1[0] = sgm1;
                    extra.sigma1[1] = (1 - sgm1);

                    for (uint8_t j = 0; j < 24; j++) {
                        for (uint8_t i = 0; i < 4; i++) {
                            extra.sigma2[i] = Sigmas_test[j][i];
                        }

                        for (uint8_t i = 1; i < 16; i++) {
                            pos1 = binary_to_pos_extra(i - 1, &extra);
                            pos2 = binary_to_pos_extra(i, &extra);
                            if (test_int_to_ker[pos1 ^ pos2] == 255 && ((pos1 ^ pos2) != (3 << mas))) {
                                printf("Problem with grey + automorphism + extra implementation (dif non trivial) \n");
                                return 1;
                            }
                            bin = pos_to_binary_extra(pos1, &extra);
                            if (bin != i - 1) {
                                printf("Problem with grey + automorphism + extra implementation (not inverse) \n");
                                return 1;
                            }
                        }
                    }
                }
            }
        }
    }

    printf("Tesseract automorphism tests .......................................... PASSED\n");
    return 0;
}

bool
test_generic_codomain(int iterations)
{
    // float dt_mus;
    // clock_t t;
    theta_point_dim4_t ker[5];
    theta_point_dim4_t inv_codomain, codomain;
    aut_t aut;
    init_aut(&aut);

    // dt_mus = 0;
    for (int n = 0; n < iterations; n++) {

        for (uint8_t j = 0; j < 24; j++) {
            for (uint8_t i = 0; i < 4; i++) {
                aut.sigma[i] = Sigmas_test[j][i];
            }

            for (int l = 0; l < 4; l++) {
                for (int i = 0; i < 16; i++) {
                    fp_random_test(&ker[l][i]);
                }
            }
            // t = clock();

            // test compute double
            compute_codomain_and_inverse(&inv_codomain, &codomain, ker, &aut);

            // dt_mus = dt_mus + (1000000. * (float)(clock()-t) / CLOCKS_PER_SEC);
            dot_prod(&inv_codomain, &inv_codomain, &codomain);
            uint32_t all_equal = fp_is_zero(&ZERO);
            fp_t mem;
            fp_copy(&mem, &inv_codomain[0]);

            for (int i = 0; i < 16; i++) {
                all_equal &= fp_is_equal(&mem, &inv_codomain[i]);
            }

            if (!all_equal) {
                printf("Problem in computing codomain generic isogeny (double)\n");
                return 1;
            }

            // test compute_simple
            compute_codomain_simple(&inv_codomain, ker, &aut);

            dot_prod(&inv_codomain, &inv_codomain, &codomain);
            all_equal = fp_is_zero(&ZERO);
            fp_copy(&mem, &inv_codomain[0]);

            for (int i = 0; i < 16; i++) {
                all_equal &= fp_is_equal(&mem, &inv_codomain[i]);
            }

            if (!all_equal) {
                printf("Problem in computing codomain generic isogeny (single)\n");
                return 1;
            }
            // test compute_simple_fast
            compute_codomain_simple_fast(&inv_codomain, ker, &aut);

            dot_prod(&inv_codomain, &inv_codomain, &codomain);
            all_equal = fp_is_zero(&ZERO);
            fp_copy(&mem, &inv_codomain[0]);

            for (int i = 0; i < 16; i++) {
                all_equal &= fp_is_equal(&mem, &inv_codomain[i]);
            }

            if (!all_equal) {
                printf("Problem in computing codomain generic isogeny (fast)\n");
                return 1;
            }
        }
    }
    // printf("Combined codomain and inverse computation tests: PASSED [%.3f µs].\n",dt_mus/(24*iterations));

    printf("Generic codomain tests ................................................ PASSED\n");
    return 0;
}

int
test_first_isogenies(unsigned int Ntest)
{

    uint8_t automorphism_list[24][4];
    uint8_t index_auto = 0;

    for (unsigned int repetition = 0; repetition < Ntest; repetition++) {

        // printf("repetition = %d\n",repetition );

        jac_point_t P, Q, TP, TQ, P16, Q16;
        ec_curve_t domain_EC[2];
        ec_curve_init_from_A(&domain_EC[0], &A0);
        ec_twist(&domain_EC[1], &domain_EC[0]);

        uint64_t iter_to_16 = TORSION_EVEN_POWER - 5;

        jac_rational_basis(&P, &Q, &TP, &TQ, &domain_EC[1], &domain_EC[0]);

        jac_dbl_iter(&P16, &P, iter_to_16, &domain_EC[0]);
        jac_dbl_iter(&Q16, &Q, iter_to_16, &domain_EC[1]);

        uint8_t vv[2];
        randombytes(vv, 2);
        uint16_t case_action = ((uint16_t)(vv[0] << 8) | vv[1]) & 8191;
        basis_components_16_t compo_16;
        compo_16.N = all_valid_gluing_16_cases[case_action][0];
        compo_16.s1 = all_valid_gluing_16_cases[case_action][1];
        compo_16.s2 = all_valid_gluing_16_cases[case_action][2];
        compo_16.s3 = all_valid_gluing_16_cases[case_action][3];
        compo_16.s4 = all_valid_gluing_16_cases[case_action][4];

        // Problem cases
        //  N1 = 5, s1 = 12, s2 = 13, s3 = 11, s4 = 1 // path 0 (with mask 12)
        //  N1 = 5, s1 = 7, s2 = 13, s3 = 13, s4 = 12 // path 1 (with mask 3).

        // quick method to test a case.
        /*
        compo_16.N = 3;
        compo_16.s1 = 7;
        compo_16.s2 = 1;
        compo_16.s3 = 0;
        compo_16.s4 = 5;
        */

        quad_point_t T16[4];

        construct_symplectic_basis_16(T16, &P16, &Q16, &domain_EC[0], &domain_EC[1], compo_16);

        first_isogenies_t first_isogenies;
        uint8_t position_theta = compute_theta_position(compo_16);
        uint8_t position_codomain = compute_codomain_position(compo_16);

        theta_point_dim4_t f1_T16[5];
        bool twitt = (position_codomain > 1);
        uint8_t mask = 3 << (2 * twitt);

        // Copy domain_EC
        copy_curve(&(first_isogenies.domain_EC[0]), &domain_EC[0]);
        copy_curve(&(first_isogenies.domain_EC[1]), &domain_EC[1]);
        // printf("After copy curve\n");

        // Computing the gluing f1: E^4 --> A1
        gluing_compute(&(first_isogenies.f1),
                       T16, // 4 elements
                       &P16,
                       &Q16,
                       domain_EC, // 2 elements
                       position_codomain,
                       position_theta);
        // printf("End gluing compute\n");

        // Evaluating T16
        // First, points not on the same curve as the auxiliary point.
        gluing_eval(&f1_T16[2 * (!twitt)], &T16[2 * (!twitt)], &(first_isogenies.f1), domain_EC, !twitt);
        gluing_eval(&f1_T16[2 * (!twitt) + 1], &T16[2 * (!twitt) + 1], &(first_isogenies.f1), domain_EC, !twitt);

        // Second, points on the same curve as the auxiliary point.
        gluing_eval_special_basis(&f1_T16[2 * twitt],
                                  &f1_T16[2 * twitt + 1],
                                  &T16[2 * twitt],
                                  &T16[2 * twitt + 1],
                                  &(first_isogenies.f1),
                                  &domain_EC[twitt]);

        // the auxiliary point
        /*
        quad_point_t fourth_point;
        quad_point_ADD(&fourth_point, T16[2 * (!twitt) ], T16[2 * (!twitt) +1], &domain_EC[!twitt]);
        gluing_eval(&f1_T16[4],&fourth_point, &(first_isogenies.f1), &domain_EC[0], !twitt );
        mask = 15 ^ mask;
        */
        invert_theta_point(&f1_T16[4], &(first_isogenies.f1).inv_codomain_auxiliary_point);

        for (uint8_t i = 0; i < 5; i++) {
            inline_hadamard(&f1_T16[i]);
        }

        // printf("End gluing eval\n");
        //  Computing f2: A1 --> A2

        theta_point_dim4_t HSK_8[5];
        for (uint8_t i = 0; i < 5; i++) {
            squared(&HSK_8[i], &f1_T16[i]);
            inline_hadamard(&HSK_8[i]);
        }

        extra_aut_t extra;
        init_extra_aut(&extra, mask);

        // printf("Testing fast version\n");
        find_hamilton_path_extra_fast(&extra, HSK_8, mask);

        bool in_the_list = false;
        for (uint8_t i = 0; i < index_auto; i++) {
            if ((automorphism_list[i][0] == extra.sigma2[0]) && (automorphism_list[i][1] == extra.sigma2[1]) &&
                (automorphism_list[i][2] == extra.sigma2[2]) && (automorphism_list[i][3] == extra.sigma2[3])) {
                in_the_list = true;
            }
        }

        if (!in_the_list) {
            automorphism_list[index_auto][0] = extra.sigma2[0];
            automorphism_list[index_auto][1] = extra.sigma2[1];
            automorphism_list[index_auto][2] = extra.sigma2[2];
            automorphism_list[index_auto][3] = extra.sigma2[3];
            index_auto += 1;
        }

        compute_codomain_simple_extra(&(first_isogenies.codomain_f2.inv_dual_null_point), HSK_8, &extra);
        first_isogenies.codomain_f2.arith_precomp = false;

        uint32_t test_T, test_zeros;
        test_T = fp_is_zero(&ZERO);
        test_zeros = 0;

        for (int g = 0; g < 5; g++) {

            for (int i = 0; i < 16; i++) {
                fp_mul(&HSK_8[g][i], &HSK_8[g][i], &first_isogenies.codomain_f2.inv_dual_null_point[i]);
            }
        }

        for (int g = 0; g < 4; g++) {
            for (int i = 0; i < 16; i++) {
                test_T &= fp_is_equal(&HSK_8[g][i], &HSK_8[g][i ^ (1 << g)]);
            }
        }
        // Compute last point
        for (int i = 0; i < 16; i++) {
            test_T &= fp_is_equal(&HSK_8[4][i], &HSK_8[4][i ^ mask]);
            test_zeros |= ~fp_is_zero(&HSK_8[4][i]);
        }

        if (!test_T) {
            printf("Problem with the T points of the product theta structure\n");
            printf("N1 = %d, s1 = %d, s2 = %d, s3 = %d, s4 = %d\n",
                   compo_16.N,
                   compo_16.s1,
                   compo_16.s2,
                   compo_16.s3,
                   compo_16.s4);

            return 1;
        }
        if (!test_zeros) {
            printf("codomain is all zeros\n");
            printf("N1 = %d, s1 = %d, s2 = %d, s3 = %d, s4 = %d\n",
                   compo_16.N,
                   compo_16.s1,
                   compo_16.s2,
                   compo_16.s3,
                   compo_16.s4);

            return 1;
        }

        for (uint8_t i = 0; i < 5; i++) {
            inline_hadamard(&HSK_8[i]);
        }

        if (!is_correct_gen_isog_test(&first_isogenies.codomain_f2, HSK_8)) {
            printf("2 torsion points are not compatible with theta action \n");
            printf("N1 = %d, s1 = %d, s2 = %d, s3 = %d, s4 = %d\n",
                   compo_16.N,
                   compo_16.s1,
                   compo_16.s2,
                   compo_16.s3,
                   compo_16.s4);

            return 1;
        }
    }
    /*
    printf("{");
    for(uint8_t i = 0; i < index_auto; i++)
    {
        printf("{%d,%d,%d,%d}, ",
    automorphism_list[i][0],automorphism_list[i][1],automorphism_list[i][2],automorphism_list[i][3]);
    }
    printf("}");
    */

    printf("First isogenies tests ................................................. PASSED\n");
    return 0;
}

bool
isog_test(int iterations)
{
    bool OK = true;
    float dt_mus;
    clock_t t;
    theta_point_dim4_t ker[5], P, Q;
    theta_struct_dim4_t codomain;

    printf("Testing a bogus generic isogeny over GF(p): \n");

    dt_mus = 0;
    for (int n = 0; n < iterations; n++) {
        for (int l = 0; l < 4; l++) {
            for (int i = 0; i < 16; i++) {
                fp_random_test(&ker[l][i]);
            }
        }
        t = clock();
        gen_isog_compute_codomain(&codomain, ker, false);
        dt_mus = dt_mus + (1000000. * (float)(clock() - t) / CLOCKS_PER_SEC);
    }
    printf("Generic codomain computation tests: PASSED [%.3f µs].\n", dt_mus / iterations);

    dt_mus = 0;
    for (int n = 0; n < iterations; n++) {
        for (int l = 0; l < 4; l++) {
            for (int i = 0; i < 16; i++) {
                fp_random_test(&ker[l][i]);
            }
        }
        t = clock();
        gen_isog_compute_codomain(&codomain, ker, true);
        dt_mus = dt_mus + (1000000. * (float)(clock() - t) / CLOCKS_PER_SEC);
    }
    printf("Generic robust codomain computation tests: PASSED [%.3f µs].\n", dt_mus / iterations);

    uint8_t mask = 3;

    dt_mus = 0;
    for (int n = 0; n < iterations; n++) {
        for (int l = 0; l < 5; l++) {
            for (int i = 0; i < 16; i++) {
                fp_random_test(&ker[l][i]);
            }
        }
        t = clock();
        gen_isog_compute_codomain_extra(&codomain, ker, mask, false);
        dt_mus = dt_mus + (1000000. * (float)(clock() - t) / CLOCKS_PER_SEC);
    }
    printf("Special generic codomain computation tests: PASSED [%.3f µs].\n", dt_mus / iterations);

    dt_mus = 0;
    for (int n = 0; n < iterations; n++) {
        for (int l = 0; l < 5; l++) {
            for (int i = 0; i < 16; i++) {
                fp_random_test(&ker[l][i]);
            }
        }
        t = clock();
        gen_isog_compute_codomain_extra(&codomain, ker, mask, true);
        dt_mus = dt_mus + (1000000. * (float)(clock() - t) / CLOCKS_PER_SEC);
    }
    printf("Special generic slow codomain computation tests: PASSED [%.3f µs].\n", dt_mus / iterations);

    dt_mus = 0;
    for (int n = 0; n < iterations; n++) {
        for (int i = 0; i < 16; i++) {
            fp_random_test(&P[i]);
        }
        t = clock();
        gen_isog_eval(&Q, &P, &codomain);
        dt_mus = dt_mus + (1000000. * (float)(clock() - t) / CLOCKS_PER_SEC);
    }
    printf("Evaluation tests: PASSED [%.3f µs].\n", dt_mus / iterations);

    dt_mus = 0;
    for (int n = 0; n < iterations; n++) {
        for (int i = 0; i < 16; i++) {
            fp_random_test(&P[i]);
        }
        t = clock();
        theta_double(&Q, &P, &codomain);
        dt_mus = dt_mus + (1000000. * (float)(clock() - t) / CLOCKS_PER_SEC);
    }
    printf("Point duplication tests: PASSED [%.3f µs].\n", dt_mus / iterations);

    return OK;
}

int
main(int argc, char *argv[])
{
    // int iterations = 100;
    // int null;
    // null = test_gluing(1);

    // test_codomain(iterations);
    // return !isog_test(iterations);

    uint32_t seed[12] = { 0 };
    int iterations = 20; // max 8192
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
    printf("Testing gluing isogenies over GF(p):\n\n");
    print_seed(seed);

    res = test_precomputation();
    res |= test_gluing_domain_shape(iterations);
    res |= test_gluing_product_theta_domain(iterations);
    res |= test_gluing_isotropic_basis(iterations);
    res |= test_gluing_compute(iterations);
    res |= test_gluing_Eval(iterations);
    res |= test_generic_codomain(iterations);
    res |= test_hypercube_automorphisms();
    res |= test_first_isogenies(iterations);
    res |= test_gluing_barycentric_lifting(iterations);

    // isog_test(iterations);

    // test_gluing_16_torsion_points(4352); //Any number bellow 8192

    //! isog_test(iterations);

    
    return res;
}
