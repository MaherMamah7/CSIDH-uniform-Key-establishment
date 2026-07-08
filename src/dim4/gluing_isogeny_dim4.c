/** @file
 *
 * @author(s) Max Duparc
 *
 * @brief dim 4 Gluing functions
 */

#include <isog_dim4.h>
#include <stdio.h>

const uint8_t gluing_int_to_ker[16] = { 255, 0, 1, 4, 2, 255, 255, 255, 3, 255, 255, 255, 4, 255, 255, 255 };
const uint8_t gluing_bin_to_pos[4][16] = {
    {15, 7, 3, 11, 9, 5, 4, 6, 2, 10, 8, 0, 1, 13, 12, 14}, 
    {15, 11, 3, 7, 6, 10, 8, 9, 1, 5, 4, 0, 2, 14, 12, 13}, 
    {15, 13, 12, 14, 6, 5, 1, 9, 8, 10, 2, 0, 4, 7, 3, 11}, 
    {15, 14, 12, 13, 9, 10, 2, 6, 4, 5, 1, 0, 8, 11, 3, 7}
}; 

const uint8_t gluing_theta_change_variable_positions[4][16] = {
    { 0, 7, 10, 13, 2, 5, 8, 15, 3, 4,  9, 14, 1, 6, 11, 12 },
    { 0, 5, 11, 14, 3, 6, 8, 13, 1, 4, 10, 15, 2, 7,  9, 12 },
    { 0, 7,  9, 14, 1, 6, 8, 15, 3, 4, 10, 13, 2, 5, 11, 12 },
    { 0, 6, 11, 13, 3, 5, 8, 14, 2, 4,  9, 15, 1, 7, 10, 12 }
};


uint8_t compute_codomain_position(const basis_components_16_t compo)
{
uint16_t unique_id = ((uint16_t)(compo.N & 3)<< 8) | ((compo.s1 & 3)<< 6)| ((compo.s2 & 3)<< 4)| ((compo.s3 & 3) << 2)| (compo.s4 & 3); 
uint8_t ans; 

    switch (unique_id) {
        case 279: case 317: case 371: case 465: case 791: case 829: case 883: case 977: ans = 0; break;
        case 285: case 311: case 369: case 467: case 797: case 823: case 881: case 979: ans = 3; break;
        case 327: case 380: case 461: case 468: case 839: case 892: case 973: case 980: ans = 1; break;
        case 333: case 372: case 455: case 476: case 845: case 884: case 967: case 988: ans = 2; break;
        default: ans = 255; break; 

    }

    return ans; 
}

uint8_t compute_theta_position(const basis_components_16_t compo)
{
uint16_t unique_id = ((uint16_t)(compo.N & 3)<< 8) | ((compo.s1 & 3)<< 6)| ((compo.s2 & 3)<< 4)| ((compo.s3 & 3) << 2)| (compo.s4 & 3); 
uint8_t ans; 
    switch (unique_id) {
        case 279: case 285: case 311: case 317: case 791: case 797: case 823: case 829: ans = 0; break;
        case 327: case 333: case 455: case 461: case 839: case 845: case 967: case 973: ans = 1; break;
        case 369: case 371: case 465: case 467: case 881: case 883: case 977: case 979: ans = 2; break; 
        case 372: case 380: case 468: case 476: case 884: case 892: case 980: case 988: ans = 3; break; 
        default: ans = 255; break; 
    }

    return ans; 
}

uint8_t
from_ec_to_theta_basis(theta_point_dim1_t *theta_dim1, const ec_point_t *P, const ec_point_t *Q)
// Cost: 3a
{

    fp_add(&(*theta_dim1)[0], &P->x, &P->z);
    fp_sub(&(*theta_dim1)[1], &P->x, &P->z);

    #ifndef NDEBUG
    if (fp_is_zero(&(*theta_dim1)[0]) || fp_is_zero(&(*theta_dim1)[1])) {
        printf("Error: Canonical 4 torsion basis not well defined, P = (1,±1)\n");
    }
    #endif

    fp_t r;
    fp_add(&r, &Q->x, &Q->z);
    uint32_t signfliped_Q = fp_is_zero(&r);
    // Swap if Q =(1:-1).
    fp_cswap(&(*theta_dim1)[0], &(*theta_dim1)[1], signfliped_Q);

    return (uint8_t)(~signfliped_Q & 15);
}

uint8_t
compute_product_theta_coord(theta_point_dim4_t *domain_product_theta_point, const ec_point_t *P, const ec_point_t *Q)
// Cost: 3M + 5S + 3a

{

    fp_t a2, ab, b2;
    theta_point_dim1_t theta_dim1;
    uint8_t basis_switch;

    basis_switch = from_ec_to_theta_basis(&theta_dim1, P, Q);

    fp_sqr(&a2, &theta_dim1[0]);
    fp_sqr(&b2, &theta_dim1[1]);
    fp_mul(&ab, &theta_dim1[0], &theta_dim1[1]);

    // coordinates of Hamming weight 0
    fp_sqr(&(*domain_product_theta_point)[0], &a2);
    // coordinates of Hamming weight 1
    fp_mul(&(*domain_product_theta_point)[1], &a2, &ab);
    fp_copy(&(*domain_product_theta_point)[2], &(*domain_product_theta_point)[1]);
    fp_copy(&(*domain_product_theta_point)[4], &(*domain_product_theta_point)[1]);
    fp_copy(&(*domain_product_theta_point)[8], &(*domain_product_theta_point)[1]);
    // coordinates of Hamming weight 2
    fp_sqr(&(*domain_product_theta_point)[3], &ab);
    fp_copy(&(*domain_product_theta_point)[5], &(*domain_product_theta_point)[3]);
    fp_copy(&(*domain_product_theta_point)[6], &(*domain_product_theta_point)[3]);
    fp_copy(&(*domain_product_theta_point)[9], &(*domain_product_theta_point)[3]);
    fp_copy(&(*domain_product_theta_point)[10], &(*domain_product_theta_point)[3]);
    fp_copy(&(*domain_product_theta_point)[12], &(*domain_product_theta_point)[3]);
    // coordinates of Hamming weight 3
    fp_mul(&(*domain_product_theta_point)[14], &b2, &ab);
    fp_copy(&(*domain_product_theta_point)[7], &(*domain_product_theta_point)[14]);
    fp_copy(&(*domain_product_theta_point)[11], &(*domain_product_theta_point)[14]);
    fp_copy(&(*domain_product_theta_point)[13], &(*domain_product_theta_point)[14]);
    // coordinates of Hamming weight 4
    fp_sqr(&(*domain_product_theta_point)[15], &b2);

    return basis_switch;
}

void
apply_change_coords_theta(theta_point_dim4_t *out, const uint8_t position_theta)
// Cost: 32a
{
    uint8_t i;
    theta_point_dim4_t tmp;

    fp_t tmp1[2], tmp2[2];

    for (i = 0; i < 4; i++) {

        fp_add(&tmp1[0],
               &(*out)[gluing_theta_change_variable_positions[position_theta][4 * i]],
               &(*out)[gluing_theta_change_variable_positions[position_theta][4 * i + 1]]);
        fp_sub(&tmp1[1],
               &(*out)[gluing_theta_change_variable_positions[position_theta][4 * i]],
               &(*out)[gluing_theta_change_variable_positions[position_theta][4 * i + 1]]);
        fp_add(&tmp2[0],
               &(*out)[gluing_theta_change_variable_positions[position_theta][4 * i + 2]],
               &(*out)[gluing_theta_change_variable_positions[position_theta][4 * i + 3]]);
        fp_sub(&tmp2[1],
               &(*out)[gluing_theta_change_variable_positions[position_theta][4 * i + 2]],
               &(*out)[gluing_theta_change_variable_positions[position_theta][4 * i + 3]]);
        fp_add(&tmp[4 * i], &tmp1[0], &tmp2[0]);
        fp_sub(&tmp[4 * i + 1], &tmp1[0], &tmp2[0]);
        fp_add(&tmp[4 * i + 2], &tmp1[1], &tmp2[1]);
        fp_sub(&tmp[4 * i + 3], &tmp1[1], &tmp2[1]);
    }

    for (i = 0; i < 16; i++) {
        fp_copy(&(*out)[i], &tmp[i]);
    }
}

void
quad_point_to_dim4_vec(theta_point_dim4_t *out, const quad_point_t *T, const bool twist)
{ // Given a point T, compute its tensor product.

    theta_point_dim2_t u12, u34;
    ec_point_t T_mont[4], tmp;
    uint8_t i, j;
    uint32_t test_zero;
    fp_t set_zero, set_one, twist_tmp;

    // TODO: Do this better.
    fp_set_zero(&set_zero);
    fp_set_one(&set_one);

    // Deals with zero points.
    for (i = 0; i < 4; i++) {
        jac_to_ec(&tmp, &(*T)[i]);
        test_zero = ec_is_zero(&tmp);

        fp_select(&T_mont[i].x, &tmp.x, &set_one, test_zero);
        fp_select(&T_mont[i].z, &tmp.z, &set_zero, test_zero);
    }

    fp_mul(&u12[0], &T_mont[0].x, &T_mont[1].x);
    fp_mul(&u12[1], &T_mont[0].x, &T_mont[1].z);
    fp_mul(&u12[2], &T_mont[0].z, &T_mont[1].x);
    fp_mul(&u12[3], &T_mont[0].z, &T_mont[1].z);

    fp_mul(&u34[0], &T_mont[2].x, &T_mont[3].x);
    fp_mul(&u34[1], &T_mont[2].x, &T_mont[3].z);
    fp_mul(&u34[2], &T_mont[2].z, &T_mont[3].x);
    fp_mul(&u34[3], &T_mont[2].z, &T_mont[3].z);

    // Handles the twist in constant time.
    test_zero = -(int32_t)(twist);
    fp_neg(&twist_tmp, &u12[1]);
    fp_select(&u12[1], &u12[1], &twist_tmp, test_zero);
    fp_neg(&twist_tmp, &u12[2]);
    fp_select(&u12[2], &u12[2], &twist_tmp, test_zero);
    fp_neg(&twist_tmp, &u34[1]);
    fp_select(&u34[1], &u34[1], &twist_tmp, test_zero);
    fp_neg(&twist_tmp, &u34[2]);
    fp_select(&u34[2], &u34[2], &twist_tmp, test_zero);

    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            fp_mul(&(*out)[4 * i + j], &u12[i], &u34[j]);
        }
    }
}

void
ec_points_to_dim2_vectors(theta_point_dim2_t *u,
                          theta_point_dim2_t *v,
                          const add_components_t *bary1,
                          const add_components_t *bary2,
                          const uint32_t twist,
                          const uint8_t mask)
{ // Given two barycentric coordinates. compute their respective tensor product
    // Cost: 9 M + 4 a

    // Compute u1 x u2
    fp_mul(&(*u)[0], &bary1->u, &bary2->u);
    fp_mul(&(*u)[1], &bary1->u, &bary2->w);
    fp_mul(&(*u)[2], &bary1->w, &bary2->u);
    fp_mul(&(*u)[3], &bary1->w, &bary2->w);

    // Compute v1 x v2
    fp_mul(&(*v)[0], &bary1->v, &bary2->v);
    fp_neg(&(*v)[1], &(*v)[0]);
    fp_select(&(*v)[0], &(*v)[0], &(*v)[1], twist);
    fp_add(&(*u)[mask], &(*u)[mask], &(*v)[0]);


    // Compute v1 x u2
    fp_mul(&(*v)[0 ^ (mask & 2)], &bary1->v, &bary2->u);
    fp_mul(&(*v)[1 ^ (mask & 2)], &bary1->v, &bary2->w);
    fp_set_zero(&(*v)[2 ^ mask]); 
    //fp_set_zero(&(*v)[3 ^ (mask & 2)]); 

    // Compute u1 x v2
    fp_mul(&(*v)[3 ^ mask], &bary1->u, &bary2->v);
    fp_add(&(*v)[0 ^ (mask & 1)], &(*v)[0 ^ (mask & 1)], &(*v)[3 ^ mask]); 
    fp_mul(&(*v)[3 ^ mask], &bary1->w, &bary2->v);
    fp_add(&(*v)[2 ^ (mask & 1)], &(*v)[2 ^ (mask & 1)], &(*v)[3 ^ mask]); 
    
    fp_set_zero(&(*v)[3 ^ mask]); 

}

void
ec_couplepoints_to_dim4_vectors(theta_point_dim4_t *real_vec,
                                theta_point_dim4_t *im_vec,
                                const quad_point_t *P,
                                const quad_point_t *Q,
                                const ec_curve_t *E,
                                const bool twistP,
                                const bool twistQ)
{
    /* Cost:
    - jac_to_xz_add_components: 11M + 5S + 7a (*4)
    - ec_points_to_dim2_vectors: 9 M + 3 a (*2)
    - u12 X u34: 16 M
    - v12 X u34: 12 M
    - u12 X v34: 12 M + 12 a
    - v12 X v34: 9 M + 9 a
    Total: 111 M + 20 S + 55 a
    */
    add_components_t bary1, bary2;
    theta_point_dim2_t u12, v12, u34, v34;
    uint32_t twist_P_pm_Q = -(int32_t)(twistP ^ twistQ);
    uint8_t i, j, mask12, mask34, masktot;


    i = jac_to_xz_add_components(&bary1, &(*P)[0], &(*Q)[0], E, twistP, twistQ);
    j = jac_to_xz_add_components(&bary2, &(*P)[1], &(*Q)[1], E, twistP, twistQ);

    mask12 = ((i << 1) | j);

    fp_select(&bary1.u, &bary1.u, &ONE, fp_is_zero(&bary1.w)); // Deals with the 0 + 0 case
    fp_select(&bary2.u, &bary2.u, &ONE, fp_is_zero(&bary2.w)); // Deals with the 0 + 0 case

    // Compute u12 and v12
    ec_points_to_dim2_vectors(&u12, &v12, &bary1, &bary2, twist_P_pm_Q, mask12);

    i = jac_to_xz_add_components(&bary1, &(*P)[2], &(*Q)[2], E, twistP, twistQ);
    j = jac_to_xz_add_components(&bary2, &(*P)[3], &(*Q)[3], E, twistP, twistQ);

    fp_select(&bary1.u, &bary1.u, &ONE, fp_is_zero(&bary1.w)); // Deals with the 0 + 0 case
    fp_select(&bary2.u, &bary2.u, &ONE, fp_is_zero(&bary2.w)); // Deals with the 0 + 0 case

    mask34 =  ((i << 1) | j);
    masktot = (mask12 << 2) | mask34; 

    // Compute u34 and v34
    ec_points_to_dim2_vectors(&u34, &v34, &bary1, &bary2, twist_P_pm_Q, mask34);

    // Compute u12 X u34
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            fp_mul(&(*real_vec)[(i<< 2) ^ j], &u12[i], &u34[j]);
        }
    }

    // Compute v12 X u34
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 4; j++) {
            fp_mul(&(*im_vec)[((i ^ mask12) << 2) ^ j], &v12[i ^ mask12], &u34[j]);
        }
    }
    // Security as we have uncontroled add.
    fp_set_zero(&(*im_vec)[((3 ^ mask12) << 2)]);
    fp_set_zero(&(*im_vec)[((3 ^ mask12) << 2) ^ 1]);
    fp_set_zero(&(*im_vec)[((3 ^ mask12) << 2) ^ 2]);
    fp_set_zero(&(*im_vec)[((3 ^ mask12) << 2) ^ 3]);


    // Compute u12 X v34
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 3; j++) {
            fp_mul(&(*im_vec)[masktot ^ 15], &u12[i], &v34[j ^ mask34]);
            fp_add(&(*im_vec)[(i<<2) ^ (j ^ mask34)], &(*im_vec)[(i<<2) ^ (j ^ mask34)], &(*im_vec)[masktot ^ 15]);
        }
    }


    // Twist parameters
    for (i = 0; i < 3; i++) {
        fp_neg(&v34[3 ^ mask34], &v34[i ^ mask34]);
        fp_select(&v34[i ^ mask34], &v34[i ^ mask34], &v34[3 ^ mask34], twist_P_pm_Q);
    }

    // Compute v12 X v34
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            fp_mul(&(*im_vec)[masktot ^ 15], &v12[i ^ mask12], &v34[j ^ mask34]);
            fp_add(&(*real_vec)[((i<< 2) ^ j) ^ masktot], &(*real_vec)[((i<< 2) ^ j) ^ masktot], &(*im_vec)[masktot ^ 15]);
        }
    }

    fp_set_zero(&(*im_vec)[masktot ^ 15]);
}

void
dim4_vector_to_theta_point(theta_point_dim4_t *out,
                           const theta_point_dim4_t *product_theta_domain,
                           const uint8_t basis_twist,
                           const uint8_t position_theta)
{
    // Cost: 16M + 96 a
    theta_point_dim4_t tmp;
    hadamard(&tmp, out);

    for (uint8_t i = 0; i < 16; i++) {
        fp_copy(&(*out)[i], &tmp[(basis_twist ^ i)]);
    }
    dot_prod(out, out, product_theta_domain);
    apply_change_coords_theta(out, position_theta);
}




void
gluing_compute_codomain_fast(theta_point_dim4_t *out, const theta_point_dim4_t *HSK_8, const uint8_t position_codomain)
{

    fp_t red[16], blue[16];
    uint8_t i, i0, i1, i2, i3; 

    // SETUP
    i0 = gluing_bin_to_pos[position_codomain][0];
    i1 = gluing_bin_to_pos[position_codomain][1];
    i  = gluing_bin_to_pos[position_codomain][2];
    i2 = gluing_bin_to_pos[position_codomain][14];
    i3 = gluing_bin_to_pos[position_codomain][15];

    fp_copy(&red[i1],&HSK_8[gluing_int_to_ker[i0 ^ i]][i0] );
    fp_copy(&blue[i2],&HSK_8[gluing_int_to_ker[i2 ^ i3]][i3] );


    for(i = 1; i < 14; i++)
    {
        i0 = i1;
        i3 = i2;
        i1 = gluing_bin_to_pos[position_codomain][i + 1];
        i2 = gluing_bin_to_pos[position_codomain][14 - i];
        fp_mul(&red[i1], &red[i0],&HSK_8[gluing_int_to_ker[i0 ^ i1]][i0]);
        fp_mul(&blue[i2], &blue[i3],&HSK_8[gluing_int_to_ker[i2 ^ i3]][i3]);
    }

    // Write solution
    // last element
    i0 = i1;
    i1 = gluing_bin_to_pos[position_codomain][15];    
    fp_mul(&(*out)[i1], &red[i0],&HSK_8[gluing_int_to_ker[i0 ^ i1]][i0]);

    // first element (Special rabbit)
    i1 = i3;
    i3 = i2;
    i2 = gluing_bin_to_pos[position_codomain][0];
    fp_mul(&blue[i2], &blue[i1],&HSK_8[gluing_int_to_ker[i3 ^ i1]][i3]);
    fp_mul(&(*out)[i2], &blue[i2],&HSK_8[gluing_int_to_ker[i2 ^ i1]][i1]);

    for(i = 1; i < 15; i++)
    {
        i0 = gluing_bin_to_pos[position_codomain][i];
        fp_mul(&(*out)[i0],&red[i0],&blue[i0]);
    }
}



void
couplepoint_to_product_of_alt_sum_theta_point(theta_point_dim4_t *out,
                                              const quad_point_t *P,
                                              const quad_point_t *Q,
                                              const ec_curve_t *E, // Always E, not E^t.
                                              const theta_point_dim4_t *product_theta_domain,
                                              const bool twistP,
                                              const bool twistQ,
                                              const uint8_t basis_swap,
                                              const uint8_t position_theta_change_variable)
{
    theta_point_dim4_t real_vec, im_vec;
    uint8_t i;
    fp_t tmp;
    uint32_t twist = -(int32_t)(twistP ^ twistQ);

    ec_couplepoints_to_dim4_vectors(&real_vec, &im_vec, P, Q, E, twistP, twistQ);

    dim4_vector_to_theta_point(&real_vec, product_theta_domain, basis_swap, position_theta_change_variable);
    dim4_vector_to_theta_point(&im_vec, product_theta_domain, basis_swap, position_theta_change_variable);

    squared(&real_vec, &real_vec);
    squared(&im_vec, &im_vec);

    for (i = 0; i < 16; i++) {
        fp_neg(&tmp, &im_vec[i]);
        fp_select(&im_vec[i], &im_vec[i], &tmp, twist);
        fp_sub(&(*out)[i], &real_vec[i], &im_vec[i]);
    }
}

void
couplepoint_to_product_of_alt_sum_theta_point_special(
    theta_point_dim4_t *out,
    const quad_point_t *P,
    const quad_point_t *Q,
    const ec_curve_t *E, // The curve where we both P and Q are defined.
    const theta_point_dim4_t *product_theta_domain,
    const bool twist,
    const uint8_t basis_swap,
    const uint8_t position_theta_change_variable)
{
    quad_point_t tmp, tmp2;
    theta_point_dim4_t T1_theta, T2_theta;

    quad_point_ADD(&tmp, P, Q, E); // P + Qaux
    quad_point_to_dim4_vec(&T1_theta, &tmp, twist);
    dim4_vector_to_theta_point(&T1_theta, product_theta_domain, basis_swap, position_theta_change_variable);
    quad_point_neg(&tmp2, P);
    quad_point_ADD(&tmp, &tmp2, Q, E); // P - Qaux
    quad_point_to_dim4_vec(&T2_theta, &tmp, twist);
    dim4_vector_to_theta_point(&T2_theta, product_theta_domain, basis_swap, position_theta_change_variable);

    dot_prod(out, &T1_theta, &T2_theta);
}



void
gluing_eval_basis(theta_point_dim4_t *out,
                          const quad_point_t *in, //4 ellements
                          const gluing_isog_dim4_t *gluing_isogeny,
                          const ec_curve_t *domain_EC // Always E
                        )
{
    bool twitt = (gluing_isogeny->position_codomain > 1);


    // First, points not on the same curve as the auxiliary point.
    for(uint8_t i =0; i < 2; i++){
        couplepoint_to_product_of_alt_sum_theta_point(&out[2 * (!twitt) + i],
                                                      &in[2 * (!twitt)+ i],
                                                      &gluing_isogeny->auxiliary_point,
                                                      domain_EC,
                                                      &gluing_isogeny->theta_product_domain,
                                                      !twitt,
                                                      twitt,
                                                      gluing_isogeny->basis_swap,
                                                      gluing_isogeny->position_theta_structure);

        couplepoint_to_product_of_alt_sum_theta_point(&out[2 * (twitt) + i],
                                                      &in[2 * (twitt) + i],
                                                      &gluing_isogeny->auxiliary_point,
                                                      domain_EC,
                                                      &gluing_isogeny->theta_product_domain,
                                                      twitt,
                                                      twitt,
                                                      gluing_isogeny->basis_swap,
                                                      gluing_isogeny->position_theta_structure);
    }
    for(uint8_t i =0; i < 4; i++)
    {
        inline_hadamard(&out[i]);
        for (uint8_t j = 0; j < 16; j++) {
            fp_mul(&out[i][j], &out[i][j], &gluing_isogeny->inv_codomain_auxiliary_point[j]);
        }
    }
}


void
gluing_eval_special_basis(theta_point_dim4_t *out1,
                          theta_point_dim4_t *out2,
                          const quad_point_t *T1,
                          const quad_point_t *T2,
                          const gluing_isog_dim4_t *gluing_isogeny,
                          const ec_curve_t *E // The curve upon P1, P2 and P_aux are defined.
)
{
    uint8_t i;

    bool twist_Qaux = (gluing_isogeny->position_codomain > 1);

    quad_point_t tmp;
    theta_point_dim4_t T1_theta, T2_theta;

    quad_point_to_dim4_vec(&T1_theta, T1, twist_Qaux);
    dim4_vector_to_theta_point(&T1_theta,
                               &gluing_isogeny->theta_product_domain,
                               gluing_isogeny->basis_swap,
                               gluing_isogeny->position_theta_structure);
    quad_point_to_dim4_vec(&T2_theta, T2, twist_Qaux);
    dim4_vector_to_theta_point(&T2_theta,
                               &gluing_isogeny->theta_product_domain,
                               gluing_isogeny->basis_swap,
                               gluing_isogeny->position_theta_structure);

    quad_point_ADD(&tmp, T1, &gluing_isogeny->auxiliary_point, E); // 2T1 + T2
    quad_point_to_dim4_vec(out1, &tmp, twist_Qaux);
    dim4_vector_to_theta_point(out1,
                               &gluing_isogeny->theta_product_domain,
                               gluing_isogeny->basis_swap,
                               gluing_isogeny->position_theta_structure);
    quad_point_ADD(&tmp, T2, &gluing_isogeny->auxiliary_point, E); // T1 + 2T2
    quad_point_to_dim4_vec(out2, &tmp, twist_Qaux);
    dim4_vector_to_theta_point(out2,
                               &gluing_isogeny->theta_product_domain,
                               gluing_isogeny->basis_swap,
                               gluing_isogeny->position_theta_structure);

    dot_prod(out1, out1, &T2_theta);

    dot_prod(out2, out2, &T1_theta);

    inline_hadamard(out1);
    inline_hadamard(out2);

    for (i = 0; i < 16; i++) {
        fp_mul(&(*out1)[i], &(*out1)[i], &gluing_isogeny->inv_codomain_auxiliary_point[i]);
        fp_mul(&(*out2)[i], &(*out2)[i], &gluing_isogeny->inv_codomain_auxiliary_point[i]);
    }
}


void
gluing_compute(gluing_isog_dim4_t *gluing_isogeny,
               const quad_point_t *T16, // 4 elements
               const jac_point_t *P_16,
               const jac_point_t *Q_16,
               const ec_curve_t *domain_EC, // 2 elements
               const uint8_t position_codomain,
               const uint8_t position_theta_structure)
{

    gluing_isogeny->position_codomain = position_codomain;
    gluing_isogeny->position_theta_structure = position_theta_structure;

    quad_point_t twist_torsion_basis_8[3];
    uint8_t i;
    jac_point_t jP_4, jQ_4;
    ec_point_t P4, Q4;

    // Compute a symplectic 4 torsion basis.
    DBL(&jP_4, P_16, &domain_EC[0]);
    DBL(&jP_4, &jP_4, &domain_EC[0]);
    jac_to_ec(&P4, &jP_4);
    DBL(&jQ_4, Q_16, &domain_EC[1]);
    DBL(&jQ_4, &jQ_4, &domain_EC[1]);
    jac_to_ec(&Q4, &jQ_4);

    // PART I: Compute the product theta null point corresponding to our basis.
    gluing_isogeny->basis_swap = compute_product_theta_coord(&gluing_isogeny->theta_product_domain, &P4, &Q4);

    // PART II, compute the codomain.
    // compute the auxiliary point, T3 + T4 if codomain_pos > 1, T1 + T2 otherwise.
    bool twitt = (position_codomain > 1);
    quad_point_t P_m_Q, PPP_p_Q, P_p_QQQ, tmp;

    // Select the right auxiliary point.
    quad_point_ADD(&gluing_isogeny->auxiliary_point, &T16[2 * twitt], &T16[2 * twitt + 1], &domain_EC[twitt]);

    // Compute the 8 torsions points on the twist of the auxiliary point.
    quad_point_ADD(&twist_torsion_basis_8[2], &T16[2 * (!twitt)], &T16[2 * (!twitt) + 1], &domain_EC[!twitt]);
    quad_point_DBL(&twist_torsion_basis_8[0], &T16[2 * (!twitt)], &domain_EC[!twitt]);
    quad_point_DBL(&twist_torsion_basis_8[1], &T16[2 * (!twitt) + 1], &domain_EC[!twitt]);
    quad_point_DBL(&twist_torsion_basis_8[2], &twist_torsion_basis_8[2], &domain_EC[!twitt]);

    // Compute the 8 torsion points on the same curve as the auxiliary point.
    quad_point_neg(&tmp, &T16[2 * twitt + 1]);
    quad_point_ADD(&P_m_Q, &T16[2 * twitt], &tmp, &domain_EC[twitt]);
    quad_point_DBL(&tmp, &T16[2 * twitt], &domain_EC[twitt]);
    quad_point_ADD(&PPP_p_Q, &tmp, &gluing_isogeny->auxiliary_point, &domain_EC[twitt]);
    quad_point_DBL(&tmp, &T16[2 * twitt + 1], &domain_EC[twitt]);
    quad_point_ADD(&P_p_QQQ, &tmp, &gluing_isogeny->auxiliary_point, &domain_EC[twitt]);

    theta_point_dim4_t HSK8[5], P_m_Q_theta, tmp_theta;

    // Points on the twist and the double of the auxiliary points.
    couplepoint_to_product_of_alt_sum_theta_point(&HSK8[2 * (!twitt)],
                                                  &gluing_isogeny->auxiliary_point,
                                                  &twist_torsion_basis_8[0],
                                                  &domain_EC[0],
                                                  &gluing_isogeny->theta_product_domain,
                                                  twitt,
                                                  !twitt,
                                                  gluing_isogeny->basis_swap,
                                                  position_theta_structure);

    couplepoint_to_product_of_alt_sum_theta_point(&HSK8[2 * (!twitt) + 1],
                                                  &gluing_isogeny->auxiliary_point,
                                                  &twist_torsion_basis_8[1],
                                                  &domain_EC[0],
                                                  &gluing_isogeny->theta_product_domain,
                                                  twitt,
                                                  !twitt,
                                                  gluing_isogeny->basis_swap,
                                                  position_theta_structure);

    couplepoint_to_product_of_alt_sum_theta_point(&HSK8[4],
                                                  &gluing_isogeny->auxiliary_point,
                                                  &twist_torsion_basis_8[2],
                                                  &domain_EC[0],
                                                  &gluing_isogeny->theta_product_domain,
                                                  twitt,
                                                  !twitt,
                                                  gluing_isogeny->basis_swap,
                                                  position_theta_structure);

    // Points on the curve
    // Lift P - Q
    quad_point_to_dim4_vec(&P_m_Q_theta, &P_m_Q, twitt);
    dim4_vector_to_theta_point(
        &P_m_Q_theta, &gluing_isogeny->theta_product_domain, gluing_isogeny->basis_swap, position_theta_structure);

    // Compute (3P + Q)(-P + Q)
    quad_point_to_dim4_vec(&tmp_theta, &PPP_p_Q, twitt);
    dim4_vector_to_theta_point(
        &tmp_theta, &gluing_isogeny->theta_product_domain, gluing_isogeny->basis_swap, position_theta_structure);
    dot_prod(&HSK8[2 * twitt], &P_m_Q_theta, &tmp_theta);

    // Compute (P + 3Q)(P - Q)
    quad_point_to_dim4_vec(&tmp_theta, &P_p_QQQ, twitt);
    dim4_vector_to_theta_point(
        &tmp_theta, &gluing_isogeny->theta_product_domain, gluing_isogeny->basis_swap, position_theta_structure);
    dot_prod(&HSK8[2 * twitt + 1], &P_m_Q_theta, &tmp_theta);

    for (i = 0; i < 5; i++) {
        inline_hadamard(&HSK8[i]);
    }

    gluing_compute_codomain_fast(&gluing_isogeny->inv_codomain_auxiliary_point, HSK8, position_codomain);
}

//
// Functions used for testing. Do not use in final implementation.
//

// DEBUG CONST.
const uint8_t inv_mod_16[16] = { 0, 1, 0, 11, 0, 13, 0, 7, 0, 9, 0, 3, 0, 5, 0, 15 };

// DEBUG FUNCTION.
void
jac_small_MUL(jac_point_t *Q, const jac_point_t *P, const uint8_t k, const ec_curve_t *curve)
{

    jac_point_init(Q);

    for (uint8_t i = 0; i < k; i++) {
        ADD(Q, Q, P, curve);
    }
}

// DEBUG FUNCTION
void
construct_symplectic_basis_4(quad_point_t *out,
                             const jac_point_t *P,
                             const jac_point_t *Q,
                             const ec_curve_t *curve,
                             const ec_curve_t *twist,
                             const basis_components_16_t compo,
                             const bool full)
{

    uint8_t N, s1, s2, s3, s4, alpha_s1, alpha_s2, alpha_s3, alpha_s4;
    N = compo.N & 3;
    s1 = compo.s1 & 3;
    s2 = compo.s2 & 3;
    s3 = compo.s3 & 3;
    s4 = compo.s4 & 3;
    alpha_s1 = (N * s1) & 3;
    alpha_s2 = (N * s2) & 3;
    alpha_s3 = (N * s3) & 3;
    alpha_s4 = (N * s4) & 3;

    // Compute the T points
    // T1
    jac_small_MUL(&out[0][0], P, N, curve);
    jac_point_init(&out[0][1]);
    jac_small_MUL(&out[0][2], P, s3, curve);
    jac_small_MUL(&out[0][3], P, (4 - s1) & 3, curve);

    // T2
    jac_point_init(&out[1][0]);
    jac_small_MUL(&out[1][1], P, N, curve);
    jac_small_MUL(&out[1][2], P, s2, curve);
    jac_small_MUL(&out[1][3], P, s4, curve);

    // T3
    copy_jac_point(&out[2][0], Q);
    jac_point_init(&out[2][1]);
    jac_small_MUL(&out[2][2], Q, alpha_s4, twist);
    jac_small_MUL(&out[2][3], Q, (4 - alpha_s2) & 3, twist);

    // T4
    jac_point_init(&out[3][0]);
    copy_jac_point(&out[3][1], Q);
    jac_small_MUL(&out[3][2], Q, alpha_s1, twist);
    jac_small_MUL(&out[3][3], Q, alpha_s3, twist);

    // Compute the S-points
    if (full) {
        // S1
        jac_small_MUL(&out[4][0], Q, (4 - N) & 3, twist);
        jac_point_init(&out[4][1]);
        jac_point_init(&out[4][2]);
        jac_point_init(&out[4][3]);

        // S2
        jac_point_init(&out[5][0]);
        jac_small_MUL(&out[5][1], Q, (4 - N) & 3, twist);
        jac_point_init(&out[4][2]);
        jac_point_init(&out[4][3]);

        // S3
        jac_point_init(&out[6][0]);
        jac_point_init(&out[6][1]);
        jac_small_MUL(&out[6][2], P, (4 - alpha_s3) & 3, curve);
        jac_small_MUL(&out[6][3], P, alpha_s1, curve);

        // S4
        jac_point_init(&out[7][0]);
        jac_point_init(&out[7][1]);
        jac_small_MUL(&out[7][2], P, (4 - alpha_s2) & 3, curve);
        jac_small_MUL(&out[7][3], P, (4 - alpha_s4) & 3, curve);
    }
}

// DEBUG FUNCTION.
void
construct_symplectic_basis_16(quad_point_t *out,
                              const jac_point_t *P,
                              const jac_point_t *Q,
                              const ec_curve_t *curve,
                              const ec_curve_t *twist,
                              const basis_components_16_t compo)
{

    uint8_t N, s1, s2, s3, s4, alpha, alpha_s1, alpha_s2, alpha_s3, alpha_s4;
    N = compo.N & 15;
    alpha = inv_mod_16[N];
    s1 = compo.s1 & 15;
    s2 = compo.s2 & 15;
    s3 = compo.s3 & 15;
    s4 = compo.s4 & 15;
    alpha_s1 = (alpha * s1) & 15;
    alpha_s2 = (alpha * s2) & 15;
    alpha_s3 = (alpha * s3) & 15;
    alpha_s4 = (alpha * s4) & 15;

    // printf("N= %d, s1 = %d,s2 = %d,s3 = %d,s4 = %d\n", N,s1, s2, s3, s4);
    //  T1
    jac_small_MUL(&out[0][0], P, N, curve);
    jac_point_init(&out[0][1]);
    jac_small_MUL(&out[0][2], P, s3, curve);
    jac_small_MUL(&out[0][3], P, (16 - s1) & 15, curve);

    // T2
    jac_point_init(&out[1][0]);
    jac_small_MUL(&out[1][1], P, N, curve);
    jac_small_MUL(&out[1][2], P, s2, curve);
    jac_small_MUL(&out[1][3], P, s4, curve);

    // T3
    copy_jac_point(&out[2][0], Q);
    jac_point_init(&out[2][1]);
    jac_small_MUL(&out[2][2], Q, alpha_s4, twist);
    jac_small_MUL(&out[2][3], Q, (16 - alpha_s2) & 15, twist);

    // T4
    jac_point_init(&out[3][0]);
    copy_jac_point(&out[3][1], Q);
    jac_small_MUL(&out[3][2], Q, alpha_s1, twist);
    jac_small_MUL(&out[3][3], Q, alpha_s3, twist);
}

void
gluing_eval(theta_point_dim4_t *out,
            const quad_point_t *P,
            const gluing_isog_dim4_t *gluing_isogeny,
            const ec_curve_t *domain_EC, // Always E, not E^t.
            bool twist_P)
{
    uint8_t i;
    bool twist_Qaux = (gluing_isogeny->position_codomain > 1);

    /*

    if (twist_P == twist_Qaux) {
        couplepoint_to_product_of_alt_sum_theta_point_special(out,
                                                              P,
                                                              &gluing_isogeny->auxiliary_point,
                                                              &domain_EC[twist_P],
                                                              &gluing_isogeny->theta_product_domain,
                                                              twist_P,
                                                              twist_Qaux,
                                                              gluing_isogeny->basis_swap,
                                                              gluing_isogeny->position_theta_structure);

    } else {*/

        couplepoint_to_product_of_alt_sum_theta_point(out,
                                                      P,
                                                      &gluing_isogeny->auxiliary_point,
                                                      &domain_EC[0],
                                                      &gluing_isogeny->theta_product_domain,
                                                      twist_P,
                                                      twist_Qaux,
                                                      gluing_isogeny->basis_swap,
                                                      gluing_isogeny->position_theta_structure);
    //}

    inline_hadamard(out);
    for (i = 0; i < 16; i++) {
        fp_mul(&(*out)[i], &(*out)[i], &gluing_isogeny->inv_codomain_auxiliary_point[i]);
    }
}




/*
void to_prod_vec(fp_t *vec, const fp_t *a, const fp_t *b){
        // Cost: 4S+5M
        fp_t aa[4], bb[4], ab[5];
        int s, i;

        // aa=[a,a^2,a^3,a^4]
        fp_copy(&aa[0],a);// a
        fp_sqr(&aa[1],&aa[0]);// a^2
        fp_mul(&aa[2],&aa[1],&aa[0]);// a^3
        fp_sqr(&aa[3],&aa[2]);// a^4

        // bb=[b,b^2,b^3,b^4]
        fp_copy(&bb[0],b);// Idem for bb
        fp_sqr(&bb[1],&bb[0]);
        fp_mul(&bb[2],&bb[1],&bb[0]);
        fp_sqr(&bb[3],&bb[2]);

        // ab=[a^4,a^3b,a^2b^2,ab^3,b^4]
        fp_copy(&ab[0],&aa[3]);
        for(int i=1;i<4; i++){
                fp_mul(&ab[i],&aa[3-i],&bb[i-1]);
        }
        fp_copy(&ab[4],&bb[3]);

        for(int i0=0; i0<2; i0++){
                for(int i1=0; i1<2; i1++){
                        for(int i2=0; i2<2; i2++){
                                for(int i3=0; i3<2; i3++){
                                        s = i0+i1+i2+i3;
                                        i = i0+2*i1+4*i2+8*i3;
                                        fp_copy(&vec[i],&ab[s]);
                                        if (s&1){// If s==1 mod 2
                                                fp_neg(&vec[i],&vec[i]);
                                        }
                                }
                        }
                }
        }

}



// Recursive version of the gluing compute no longer used.
// Only keeped because I strugled so much with it that I do not have the heart of putting it of its missery.
void
gluing_compute_codomain_recursive(theta_point_dim4_t *out,
                                  const theta_point_dim4_t *HSK8,
                                  const uint8_t position,
                                  const uint8_t start_pos,
                                  const uint8_t end_pos)
{

    if (end_pos - start_pos <= 4) {
        uint8_t i0 = gluing_isomorphism_codomain(gluing_bin_to_pos[start_pos], position);
        uint8_t i1 = gluing_isomorphism_codomain(gluing_bin_to_pos[start_pos + 1], position);
        uint8_t i2 = gluing_isomorphism_codomain(gluing_bin_to_pos[start_pos + 2], position);
        uint8_t i3 = gluing_isomorphism_codomain(gluing_bin_to_pos[start_pos + 3], position);
        uint8_t l0 = i0 ^ i1;
        uint8_t l1 = i1 ^ i2;
        uint8_t l2 = i2 ^ i3;
        if (start_pos == 0) {
            l0 = i0 ^ i2;
            // Total cost: 7M
            fp_mul(&(*out)[i0], &HSK8[gluing_int_to_ker[l1]][i1], &HSK8[gluing_int_to_ker[l2]][i3]);
            fp_mul(&(*out)[i1], &HSK8[gluing_int_to_ker[l1]][i2], &HSK8[gluing_int_to_ker[l2]][i3]);
            fp_mul(&(*out)[i3], &HSK8[gluing_int_to_ker[l1]][i1], &HSK8[gluing_int_to_ker[l2]][i2]);

            fp_mul(&(*out)[i2], &HSK8[gluing_int_to_ker[l0]][i0], &(*out)[i0]);
            fp_mul(&(*out)[i0], &HSK8[gluing_int_to_ker[l0]][i2], &(*out)[i0]);
            fp_mul(&(*out)[i1], &HSK8[gluing_int_to_ker[l0]][i0], &(*out)[i1]);
            fp_mul(&(*out)[i3], &HSK8[gluing_int_to_ker[l0]][i0], &(*out)[i3]);

        } else {
            // Total cost: 6M
            fp_mul(&(*out)[i0], &HSK8[gluing_int_to_ker[l1]][i2], &HSK8[gluing_int_to_ker[l2]][i3]);
            fp_mul(&(*out)[i3], &HSK8[gluing_int_to_ker[l0]][i0], &HSK8[gluing_int_to_ker[l1]][i1]);

            fp_mul(&(*out)[i1], &HSK8[gluing_int_to_ker[l0]][i0], &(*out)[i0]);
            fp_mul(&(*out)[i0], &HSK8[gluing_int_to_ker[l0]][i1], &(*out)[i0]);
            fp_mul(&(*out)[i2], &HSK8[gluing_int_to_ker[l2]][i3], &(*out)[i3]);
            fp_mul(&(*out)[i3], &HSK8[gluing_int_to_ker[l2]][i2], &(*out)[i3]);
        }

    } else {
        uint8_t mid_pos = start_pos + ((end_pos - start_pos) >> 1);
        gluing_compute_codomain_recursive(out, HSK8, position, start_pos, mid_pos);
        gluing_compute_codomain_recursive(out, HSK8, position, mid_pos, end_pos);

        fp_t top_right, bot_left;
        uint8_t i_mid, i_midm1, l, i, q;

        // Compute top_right factor
        i_mid = gluing_isomorphism_codomain(gluing_bin_to_pos[mid_pos], position);
        fp_copy(&top_right, &(*out)[i_mid]);

        i_midm1 = gluing_isomorphism_codomain(gluing_bin_to_pos[mid_pos - 1], position);

        l = i_mid ^ i_midm1;
        fp_mul(&top_right, &top_right, &HSK8[gluing_int_to_ker[l]][i_mid]);

        // Compute bot_left factor
        fp_copy(&bot_left, &(*out)[i_midm1]);
        fp_mul(&bot_left, &bot_left, &HSK8[gluing_int_to_ker[l]][i_midm1]);

        // Multiplying by top_right factor
        for ( q = start_pos; q < mid_pos; q++) {
            i = gluing_isomorphism_codomain(gluing_bin_to_pos[q], position);
            fp_mul(&(*out)[i], &(*out)[i], &top_right);
        }

        // Multiplying by bot_left factor
        for ( q = mid_pos; q < end_pos; q++) {
            i = gluing_isomorphism_codomain(gluing_bin_to_pos[q], position);
            fp_mul(&(*out)[i], &(*out)[i], &bot_left);
        }
    }
}


uint8_t
gluing_isomorphism_codomain(const uint8_t n, const uint8_t position)
{
    return (((n & 1) << (position)) | (((n >> 1) & 1) << (1 ^ position)) | (((n >> 2) & 1) << (2 ^ position)) |
            (((n >> 3) & 1) << (3 ^ position)));
}


void
gluing_compute_codomain(theta_point_dim4_t *out, const theta_point_dim4_t *HSK8, const uint8_t position_codomain)
{

    fp_t top_right, bot_left;
    uint8_t len, pad, mid_up, mid_low, i0, i1, i2, i3, l0, l1, l2;

    // INIT LOOP:
    // Rabbit head
    i0 = gluing_isomorphism_codomain(gluing_bin_to_pos[0], position_codomain);
    i1 = gluing_isomorphism_codomain(gluing_bin_to_pos[1], position_codomain);
    i2 = gluing_isomorphism_codomain(gluing_bin_to_pos[2], position_codomain);
    i3 = gluing_isomorphism_codomain(gluing_bin_to_pos[3], position_codomain);
    l0 = i0 ^ i2;
    l1 = i1 ^ i2;
    l2 = i2 ^ i3;

    fp_mul(&(*out)[i0], &HSK8[gluing_int_to_ker[l1]][i1], &HSK8[gluing_int_to_ker[l2]][i3]);
    fp_mul(&(*out)[i1], &HSK8[gluing_int_to_ker[l1]][i2], &HSK8[gluing_int_to_ker[l2]][i3]);
    fp_mul(&(*out)[i3], &HSK8[gluing_int_to_ker[l1]][i1], &HSK8[gluing_int_to_ker[l2]][i2]);

    fp_mul(&(*out)[i2], &HSK8[gluing_int_to_ker[l0]][i0], &(*out)[i0]);
    fp_mul(&(*out)[i0], &HSK8[gluing_int_to_ker[l0]][i2], &(*out)[i0]);
    fp_mul(&(*out)[i1], &HSK8[gluing_int_to_ker[l0]][i0], &(*out)[i1]);
    fp_mul(&(*out)[i3], &HSK8[gluing_int_to_ker[l0]][i0], &(*out)[i3]);

    // Rabbit tail:
    for (pad = 4; pad < 16; pad += 4) {
        i0 = gluing_isomorphism_codomain(gluing_bin_to_pos[pad], position_codomain);
        i1 = gluing_isomorphism_codomain(gluing_bin_to_pos[pad ^ 1], position_codomain);
        i2 = gluing_isomorphism_codomain(gluing_bin_to_pos[pad ^ 2], position_codomain);
        i3 = gluing_isomorphism_codomain(gluing_bin_to_pos[pad ^ 3], position_codomain);
        l0 = i0 ^ i1;
        l1 = i1 ^ i2;
        l2 = i2 ^ i3;

        // Total cost: 6M
        fp_mul(&(*out)[i0], &HSK8[gluing_int_to_ker[l1]][i2], &HSK8[gluing_int_to_ker[l2]][i3]);
        fp_mul(&(*out)[i3], &HSK8[gluing_int_to_ker[l0]][i0], &HSK8[gluing_int_to_ker[l1]][i1]);

        fp_mul(&(*out)[i1], &HSK8[gluing_int_to_ker[l0]][i0], &(*out)[i0]);
        fp_mul(&(*out)[i0], &HSK8[gluing_int_to_ker[l0]][i1], &(*out)[i0]);
        fp_mul(&(*out)[i2], &HSK8[gluing_int_to_ker[l2]][i3], &(*out)[i3]);
        fp_mul(&(*out)[i3], &HSK8[gluing_int_to_ker[l2]][i2], &(*out)[i3]);
    }
    
    // MAIN LOOP:
    for (len = 4; len < 16; len <<= 1) {
        for (pad = 0; pad < 16; pad += (len << 1)) {
            mid_up = pad ^ len;
            mid_low = pad ^ (len - 1);
            i0 = gluing_isomorphism_codomain(gluing_bin_to_pos[mid_low], position_codomain);
            i1 = gluing_isomorphism_codomain(gluing_bin_to_pos[mid_up], position_codomain);
            l0 = i0 ^ i1;

            fp_mul(&top_right, &(*out)[i1], &HSK8[gluing_int_to_ker[l0]][i1]);
            fp_mul(&bot_left, &(*out)[i0], &HSK8[gluing_int_to_ker[l0]][i0]);

            for (l0 = pad; l0 < mid_up; l0++) {
                i0 = gluing_isomorphism_codomain(gluing_bin_to_pos[l0], position_codomain);
                i1 = gluing_isomorphism_codomain(gluing_bin_to_pos[l0 ^ len], position_codomain);
                fp_mul(&(*out)[i0], &(*out)[i0], &top_right);
                fp_mul(&(*out)[i1], &(*out)[i1], &bot_left);
            }
        }
    }
}

*/

