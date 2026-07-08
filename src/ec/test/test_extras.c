#include "test_extras.h"
#include "rng.h"

// Make n random-ish field elements (for tests only!).
void
fp_random_test(fp_t *a)
{
    uint8_t tmp[FP_ENCODED_BYTES];

    randombytes(tmp, sizeof(tmp));

    fp_decode(a, tmp);
}

void
ec_random_normalized_test(ec_point_t *P, const ec_curve_t *curve)
{
    fp_set_one(&P->z);
    while (1) {
        fp_random_test(&P->x);
        if (x_is_on_curve(&P->y, &P->x, curve)) {
            fp_sqrt(&P->y);
            break;
        }
    }
}

void
ec_random_test(ec_point_t *P, const ec_curve_t *curve)
{
    ec_random_normalized_test(P, curve);
    fp_random_test(&P->z);
    fp_mul(&P->x, &P->x, &P->z);
    fp_mul(&P->y, &P->y, &P->z);
}

void
xz_random_normalized_test(xz_point_t *P, const ec_curve_t *curve)
{
    fp_set_one(&P->z);
    fp_t y2;
    while (1) {
        fp_random_test(&P->x);
        if (x_is_on_curve(&y2, &P->x, curve)) {
            break;
        }
    }
}

void
xz_random_test(xz_point_t *P, const ec_curve_t *curve)
{
    xz_random_normalized_test(P, curve);
    fp_random_test(&P->z);
    fp_mul(&P->x, &P->x, &P->z);
}