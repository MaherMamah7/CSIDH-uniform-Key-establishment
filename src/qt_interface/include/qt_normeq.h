#ifndef QT_NORMEQ_H
#define QT_NORMEQ_H

#include "intbig.h"
#include "constants.h"
#include "fp_constants.h"

/** @brief Type for qt_normeq output
 *
 * @typedef qt_output_t
 *
 * @struct qt_output
 */
typedef struct qt_output_t
{
    digit_t N1[NWORDS_ORDER];
    digit_t Nb1[NWORDS_ORDER];
    digit_t A1[NWORDS_ORDER];
    digit_t A2[NWORDS_ORDER];
    digit_t B1[NWORDS_ORDER];
    digit_t B2[NWORDS_ORDER];
    digit_t C1[NWORDS_ORDER];
    digit_t C2[NWORDS_ORDER];
    digit_t D1[NWORDS_ORDER];
    digit_t D2[NWORDS_ORDER];
    digit_t E1[NWORDS_ORDER];
    digit_t E2[NWORDS_ORDER];
    // Norm of the reduced ideal N mod 4
    uint8_t Nmod4;
    // value of a mod 2 where alpha = a+b*omega and N generate the reduced ideal
    //(a+b mod 2) where N and gen=a+b pi generators of shortest equivalent ideal
    bool two_isog_choice;
} qt_output_t;

/**
 * @brief Sample a (near-)reduced ideal (N, alpha = a + b*omega).
 *
 * N is a random prime slightly below sqrt(p) with -p a QR mod N, so the ideal
 * is (close to) reduced. omega = (1 + sqrt(-p))/2.
 *
 * @returns 1 on success, 0 on error.
 */
int qt_sample_ideal(ibz_t *N, ibz_t *alpha_a, ibz_t *alpha_b, gmp_randstate_t rng);

/**
 * @brief Solve the norm equation N1 + N2 = 2^{e_sol} for the ideal (N, alpha)
 * and fill sol (a qt_output_t consumable by compute_action).
 *
 * @returns 1 on success, 0 on failure (no solution / bad input).
 */
int qt_normeq(qt_output_t *sol, const ibz_t *N, const ibz_t *alpha_a, const ibz_t *alpha_b);

/**
 * @brief Convenience: sample a reduced ideal and solve its norm equation.
 *
 * @returns 1 on success, 0 on failure.
 */
int qt_sample_and_normeq(qt_output_t *sol, gmp_randstate_t rng);

#endif