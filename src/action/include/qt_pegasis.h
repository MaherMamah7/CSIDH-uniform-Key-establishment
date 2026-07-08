/** @file
 *
 * @author Pierrick Dartois
 *
 * @brief Main functions for the qt-pegasis algorithm
 */

#ifndef QT_PEGASIS_H
#define QT_PEGASIS_H

#include <fp.h>
#include <mp.h>
#include <intbig.h>
#include <ec.h>
#include <isog_dim4.h>
// TODO: include quadratic
#include <stdio.h>
#include <assert.h>
#include <qt_normeq.h>

/** @defgroup Structures Quadratic structures
 * @{
 */

/** @brief Structure that contains integers from the norm equation solution.
 *
 * @typedef norm_eq_sol_t
 *
 * @struct norm_eq_sol_t
 */
typedef  qt_output_t norm_eq_sol_t;


/** // end Quadratic structures
 * @}
 */

/** @defgroup Functions Main functions
 * @{
 */

/**
 * @brief Computes the result of the qt_pegasis action given a norm equation solution.
 *
 * @param Ea (output) resulting curve [a]*E, where a is the input ideal.
 * @param Eabar (output) resulting curve [bar{a}]*E.
 * @param E (input) starting curve.
 * @param sol (input) norm equation solution associated to the input ideal a.
 */
void
compute_action(ec_curve_t *Ea, ec_curve_t *Eabar, ec_curve_t *E, const norm_eq_sol_t *sol);

/** // end Main functions
 * @}
 */

#endif