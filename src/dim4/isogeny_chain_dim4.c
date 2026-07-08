/** @file
 *
 * @author(s) Pierrick Dartois, Max Duparc
 *
 * @brief dim 4 isogeny chain functions
 */

#include <isog_dim4.h>
#include <stdbool.h>
#include <assert.h>



void
compute_first_isogenies(first_isogenies_t *first_isogenies,
                        const ec_curve_t *domain_EC,
                        const quad_point_t *T16,
                        const jac_point_t *P_16,
                        const jac_point_t *Q_16,
                        const uint8_t position_codomain,
                        const uint8_t position_theta_structure)
{

    theta_point_dim4_t f1_T16[5];
    bool twitt = (position_codomain > 1);
    uint8_t mask = 3 << (2 * twitt);

    // Copy domain_EC
    copy_curve(&(first_isogenies->domain_EC[0]), &domain_EC[0]);
    copy_curve(&(first_isogenies->domain_EC[1]), &domain_EC[1]);

    // Computing the gluing f1: E^4 --> A1
    gluing_compute(&(first_isogenies->f1),
                   T16, // 4 elements
                   P_16,
                   Q_16,
                   domain_EC, // 2 elements
                   position_codomain,
                   position_theta_structure);

    // Evaluating T16
    // First, points not on the same curve as the auxiliary point.
    gluing_eval(&f1_T16[2 * (!twitt)], &T16[2 * (!twitt)], &(first_isogenies->f1), &domain_EC[0], !twitt);
    gluing_eval(&f1_T16[2 * (!twitt) + 1], &T16[2 * (!twitt) + 1], &(first_isogenies->f1), &domain_EC[0], !twitt);

    // Second, points on the same curve as the auxiliary point.
    gluing_eval_special_basis(&f1_T16[2 * twitt],
                              &f1_T16[2 * twitt + 1],
                              &T16[2 * twitt],
                              &T16[2 * twitt + 1],
                              &(first_isogenies->f1),
                              &domain_EC[twitt]);

    // the auxiliary point
    invert_theta_point(&f1_T16[4], &(first_isogenies->f1).inv_codomain_auxiliary_point);    // cost 48M

    // Apply Hadamard
    for (uint8_t i = 0; i < 5; i++) {
        inline_hadamard(&f1_T16[i]);
    }
    // Computing f2: A1 --> A2
    gen_isog_compute_codomain_extra(&(first_isogenies->codomain_f2), f1_T16, mask, false);
}

// Given a point P on E^4, compute f_2(f_1(P)).
void
first_isogenies_eval(theta_point_dim4_t *out,
                     const quad_point_t *P,
                     const bool twist,
                     const first_isogenies_t *first_isogenies)
{
    gluing_eval(out, P, &(first_isogenies->f1), &first_isogenies->domain_EC[0], twist);
    inline_hadamard(out);
    gen_isog_eval(out, out, &(first_isogenies->codomain_f2));
}


// Given a basis B on E^4, compute f_2(f_1(B)).
void
first_isogenies_eval_basis(theta_point_dim4_t *out,
                     const quad_point_t *P,
                     const first_isogenies_t *first_isogenies)
{
    gluing_eval_basis(out, P, &first_isogenies->f1, &first_isogenies->domain_EC[0]);
    for(uint8_t i = 0; i < 4; i++)
    {
    inline_hadamard(&out[i]);
    gen_isog_eval(&out[i], &out[i], &(first_isogenies->codomain_f2));
    }
}


#ifndef NDEBUG

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

bool
is_correct_gen_isog(theta_struct_dim4_t *domain, const theta_point_dim4_t *ker8, const uint8_t hight)
{
    theta_point_dim4_t ker2[4], ker2_bis[4];
    uint8_t s = 0;

    for (uint8_t i = 0; i < 4; i++) {
        theta_double_iter(&ker2[i], &ker8[i], domain, hight);
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
    
    //for (uint8_t j = 0; j < 16; j++) {
    //    printf("U_{%u}(0) is zero? %u\n", j, fp_is_zero(&domain->inv_dual_null_point[j]));
    //}
    //for (uint8_t j = 0; j < 16; j++) {
    //    printf("theta_{%u}(0) is zero? %u\n", j, fp_is_zero(&domain->null_point[j]));
    //}

    for (uint8_t j = 0; j < 16; j++) {
        if (fp_is_zero(&domain->null_point[j])) {
            printf("Zeros everywhere\n");
            return false;
        }
    }
    return true;
}
#endif

void
compute_chain_codomain(theta_struct_dim4_t *codomain, const quad_point_t *ker, first_isogenies_t *first_isogenies)
{
    uint16_t strat_idx = 0, limit_idx = 0, prev = 0, mem_idx = 0, n = TORSION_EVEN_POWER - 4;

    quad_point_t domain_ker_elts[MAX_IDX_STRAT_INIT][4];
    theta_point_dim4_t ker_elts[MAX_IDX_STRAT][4];
    uint16_t memory[MAX_IDX_STRAT];

    copy_theta_struct(codomain, &(first_isogenies->codomain_f2));

    for (uint8_t i = 0; i < 4; i++) {
        quad_point_copy(&domain_ker_elts[0][i], &ker[i]);
    }

    theta_copy(&codomain->inv_dual_null_point, &first_isogenies->codomain_f2.inv_dual_null_point);
    codomain->arith_precomp = false;

    for (uint16_t k = 0; k < n; k++) {
        while (prev != n - 1 - k) {
            prev += STRATEGY[strat_idx];
            memory[mem_idx] = STRATEGY[strat_idx];
            mem_idx += 1;


            // Perform the doublings and update kernel elements
            // Prevent the last unnecessary doublings for first isogeny computation
            if (k == 0) {
                if (prev != n - 1) {

                    for (uint8_t j = 0; j < 2; j++) {
                        quad_point_dbl_iter(&domain_ker_elts[limit_idx + 1][j],
                                            &domain_ker_elts[limit_idx][j],
                                            STRATEGY[strat_idx],
                                            &(first_isogenies->domain_EC[0]));
                        quad_point_dbl_iter(&domain_ker_elts[limit_idx + 1][j ^ 2],
                                            &domain_ker_elts[limit_idx][j ^ 2],
                                            STRATEGY[strat_idx],
                                            &(first_isogenies->domain_EC[1]));
                    }
                    limit_idx += 1;
                }

            } else {

                for (uint8_t i = 0; i < 4; i++) {
                    theta_double_iter(
                        &ker_elts[limit_idx + 1][i], &ker_elts[limit_idx][i], codomain, STRATEGY[strat_idx]);
                }
                limit_idx += 1;
            }
            // Update bookkeeping variable
            strat_idx += 1;
        }

        // Compute the codomain from the 8-torsion
        if (k > 0) {
            #ifndef NDEBUG
                // Check correctness of the codomain
                if (k != n-1){
                    assert(is_correct_gen_isog(codomain, ker_elts[limit_idx], 2)); 
                }
            #endif
            // Caution because zeros can appear (search_path=true)
            if ((k <= 5)) {
                gen_isog_compute_codomain(codomain, ker_elts[limit_idx], true);
            } else {
                gen_isog_compute_codomain(codomain, ker_elts[limit_idx], false);
            }
        }

        // "Remove" elements from "list" i.e. update limit_idx and prev
        if ((k > 0) && (limit_idx > 0)) {
            limit_idx -= 1;
        }
        if (mem_idx > 0) {
            mem_idx -= 1;
            prev -= memory[mem_idx];
        }

        // Push through points for the next step
        if (k == 0) {
            for (uint16_t i = 0; i <= limit_idx; i++) {
                first_isogenies_eval_basis(ker_elts[i], domain_ker_elts[i], first_isogenies);
            }
        } else {
            for (uint16_t i = 0; i <= limit_idx; i++) {
                for (uint8_t j = 0; j < 4; j++) {
                    gen_isog_eval(&ker_elts[i][j], &ker_elts[i][j], codomain);
                }
            }
        }
    }
}



