/** @file
 *
 * @author(s) Pierrick Dartois, Max Duparc
 *
 * @brief generic dim 4 isogeny functions
 */

#include <isog_dim4.h>
#include <stdint.h>

const uint8_t second_iso_bin_to_pos[3][16] = {
    { 9, 5, 1, 0, 2, 3, 7,  6,  4,  8, 10, 11, 15, 14, 12, 13 }, //  if mask = 12
    { 6, 5, 4, 0, 1, 3, 2, 10, 11,  9,  8, 12, 13, 14, 15,  7 }, //  if mask = 3
    { 0, 1, 3, 2, 6, 7, 5,  4, 12, 13, 15, 14, 10, 11,  9,  8 }  // gray code
};

const uint8_t int_to_ker[16] = { 255, 0, 1, 255, 2, 255, 255, 255, 3, 255, 255, 255, 255, 255, 255, 255 };

const uint8_t Sigmas[24][4] = {
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

// List of sufficient automorphisms to ensure proper evalaution of any points.
const uint8_t Sigma_fast[3][4] = {
    { 0, 1, 2, 3 },
    { 0, 1, 3, 2 },
    { 3, 2, 0, 1 }
};

// List of sufficient automorphisms to ensure proper evalaution of any points in the second isogeny
const uint8_t Sigma_extra_fast[6][4] = {
    { 2, 1, 3, 0 },
    { 3, 2, 0, 1 },
    { 1, 2, 3, 0 },
    { 3, 2, 1, 0 },
    { 3, 1, 2, 0 },
    { 3, 1, 0, 2 }
};


void
compute_codomain_simple_fast(theta_point_dim4_t *inv_null_pt_dual, const theta_point_dim4_t *HSK_8, const aut_t *aut)
{

    fp_t red[16], blue[16];
    uint8_t i, i0, i1, i2, i3; 

    uint8_t pos[16];
    for(i =0; i < 16; i++) pos[i] = binary_to_pos(i, aut);

    // setup
    i0 = pos[0];
    i1 = pos[1];
    i2 = pos[14];
    i3 = pos[15];

    fp_copy(&red[i1],&HSK_8[int_to_ker[i0 ^ i1]][i0] );
    fp_copy(&blue[i2],&HSK_8[int_to_ker[i2 ^ i3]][i3] );

    for(i = 1; i < 14; i++)
    {
        i0 = i1;
        i3 = i2;
        i1 = pos[i + 1];
        i2 = pos[14 - i];
        fp_mul(&red[i1], &red[i0],&HSK_8[int_to_ker[i0 ^ i1]][i0]);
        fp_mul(&blue[i2], &blue[i3],&HSK_8[int_to_ker[i2 ^ i3]][i3]);
    }

    // Write solution
    i0 = i1;
    i1 = pos[15];
    i3 = i2;
    i2 = pos[0];
    
    fp_mul(&(*inv_null_pt_dual)[i1], &red[i0],&HSK_8[int_to_ker[i0 ^ i1]][i0]);
    fp_mul(&(*inv_null_pt_dual)[i2], &blue[i3],&HSK_8[int_to_ker[i2 ^ i3]][i3]);

    for(i = 1; i < 15; i++)
    {
        i0 = pos[i];
        fp_mul(&(*inv_null_pt_dual)[i0],&red[i0],&blue[i0]);
    }
}

void
compute_codomain_simple_extra_fast(theta_point_dim4_t *inv_null_pt_dual, const theta_point_dim4_t *HSK_8, const extra_aut_t *aut)
{

    fp_t red[16], blue[16];
    uint8_t i, i0, i1, i2, i3; 

    uint8_t pos[16];
    for(i = 0; i < 16; i++) pos[i] = binary_to_pos_extra(i, aut);

    // setup
    i0 = pos[0];
    i1 = pos[1];
    i2 = pos[14];
    i3 = pos[15];

    fp_copy(&red[i1],&HSK_8[gluing_int_to_ker[i0 ^ i1]][i0] );
    fp_copy(&blue[i2],&HSK_8[gluing_int_to_ker[i2 ^ i3]][i3] );

    for(i = 1; i < 14; i++)
    {
        i0 = i1;
        i3 = i2;
        i1 = pos[i + 1];
        i2 = pos[14 - i];
        fp_mul(&red[i1], &red[i0],&HSK_8[gluing_int_to_ker[i0 ^ i1]][i0]);
        fp_mul(&blue[i2], &blue[i3],&HSK_8[gluing_int_to_ker[i2 ^ i3]][i3]);
    }

    // Write solution
    i0 = i1;
    i1 = pos[15];
    i3 = i2;
    i2 = pos[0];
    
    fp_mul(&(*inv_null_pt_dual)[i1], &red[i0],&HSK_8[gluing_int_to_ker[i0 ^ i1]][i0]);
    fp_mul(&(*inv_null_pt_dual)[i2], &blue[i3],&HSK_8[gluing_int_to_ker[i2 ^ i3]][i3]);

    for(i = 1; i < 15; i++)
    {
        i0 = pos[i];
        fp_mul(&(*inv_null_pt_dual)[i0],&red[i0],&blue[i0]);
    }

}


// functions to construct the zero map of a given HIIP problem. 

void
compute_zero_map(theta_point_zero_map *maps, const theta_point_dim4_t *HSK_8)
{
    uint8_t g, pad, m;
    for (g = 0; g < 4; g++) {
        for (pad = 0; pad < 16; pad += (1 << (g + 1))) {
            for (m = pad; m < (pad + (1 << g)); m++) {
                maps[g][m] = (fp_is_zero(&HSK_8[g][m]) == 0);
                maps[g][m ^ (1 << g)] = maps[g][m];
            }
        }
    }
}


void
compute_zero_map_extra(theta_point_zero_map *maps,
                       const theta_point_dim4_t *HSK_8,
                       const uint8_t mask,
                       const uint8_t shift)
{
    uint8_t pad, m;
    compute_zero_map(maps, HSK_8);

    {
        // Compute matrix of zeros for the 5th point.
        for (pad = 0; pad < 16; pad += (1 << (shift + 2))) {
            for (m = pad; m < (pad + (1 << (shift + 1))); m++) {
                maps[4][m] = (fp_is_zero(&HSK_8[4][m]) == 0);
                maps[4][m ^ mask] = maps[4][m];
            }
        }
    }
}

// Functions to find valid hamiltonian path
static inline bool
is_correct_path(aut_t *aut, const theta_point_dim4_t *HSK_8)
{
    uint8_t i, ip1, l;
    bool res = true;
    ip1 = binary_to_pos(0, aut);
    for (uint8_t q = 0; q < 15; q++) {
        i = ip1;
        ip1 = binary_to_pos(q + 1, aut);
        l = i ^ ip1;
        res = res & (!fp_is_zero(&HSK_8[int_to_ker[l]][i]));
    }
    return res;
}

static inline bool
is_correct_path_fast(aut_t *aut, const theta_point_zero_map *HSK_8)
{
    uint8_t i, ip1, l;
    bool res = true;
    ip1 = binary_to_pos(0, aut);

    for (uint8_t q = 0; q < 15; q++) {
        i = ip1;
        ip1 = binary_to_pos(q + 1, aut);
        l = i ^ ip1;
        res = res && HSK_8[int_to_ker[l]][i];
    }
    return -res;
}

void
find_hamilton_path(aut_t *aut, const theta_point_dim4_t *HSK_8)
{
    // Returns the correct Tesseract automorphism from the Gray Hamilton path
    // to a correct one that does not go through zero.
    // Assumes aut is already initialised.
    // Not used in practice because too slow.
    aut_t aut2;
    bool sel;
    for (uint8_t i = 0; i < 24; i++) {
        for (uint8_t j = 0; j < 4; j++) {
            aut2.sigma[j] = Sigmas[i][j];
        }

        sel = is_correct_path(&aut2, HSK_8);
        aut_select(aut, aut, &aut2, sel);
    }
}

void
find_hamilton_path_fast(aut_t *aut, const theta_point_dim4_t *HSK_8)
{
    // Returns the correct Tesseract automorphism from the Gray Hamilton path
    // to a correct one that does not go through zero.
    // Assumes aut is already initialised.
    // Not used in practice because too slow.
    aut_t aut2;

    theta_point_zero_map maps[4];
    compute_zero_map(maps, HSK_8);

    bool sel;
    for (uint8_t i = 0; i < 3; i++) {
        for (uint8_t j = 0; j < 4; j++) {
            aut2.sigma[j] = Sigma_fast[i][j];
        }

        sel = is_correct_path_fast(&aut2, maps);
        aut_select(aut, aut, &aut2, sel);
    }
}

void
gen_isog_compute_codomain(theta_struct_dim4_t *codomain, const theta_point_dim4_t *ker, const bool search_path)
{
    theta_point_dim4_t HSK_8[4];
    uint8_t i;
    for (i = 0; i < 4; i++) {
        squared(&HSK_8[i], &ker[i]);
        inline_hadamard(&HSK_8[i]);
    }

    aut_t aut;
    init_aut(&aut);
    if (search_path) { // Usually not needed, used as a caution
        find_hamilton_path_fast(&aut, HSK_8);
    }

    //compute_codomain_simple(&codomain->inv_dual_null_point, HSK_8, &aut);
    compute_codomain_simple_fast(&codomain->inv_dual_null_point, HSK_8, &aut);
    codomain->arith_precomp = false;
}

void
gen_isog_eval(theta_point_dim4_t *imP, const theta_point_dim4_t *P, const theta_struct_dim4_t *codomain)
{
    // Cost 16 M + 16 S + 128 a
    squared(imP, P);
    inline_hadamard(imP);
    dot_prod(imP, imP, &(codomain->inv_dual_null_point));
    inline_hadamard(imP);
}

// Reliable method to check if a point is zero or not.
static inline uint32_t
is_correct_path_extra(extra_aut_t *extra, const theta_point_dim4_t *HSK_8)
{
    uint8_t i, ip1, l;
    uint32_t res = 0xffffffff;
    ip1 = binary_to_pos_extra(0, extra);
    for (uint8_t q = 0; q < 15; q++) {
        i = ip1;
        ip1 = binary_to_pos_extra(q + 1, extra);
        l = i ^ ip1;
        res = res & (~fp_is_zero(&HSK_8[gluing_int_to_ker[l]][i]));
    }
    return res;
}

// Faster method
static inline bool
is_correct_path_extra_fast(extra_aut_t *extra, const theta_point_zero_map *HSK_8)
{
    uint8_t i, ip1, l;
    bool res = true;
    ip1 = binary_to_pos_extra(0, extra);
    for (uint8_t q = 0; q < 15; q++) {
        i = ip1;
        ip1 = binary_to_pos_extra(q + 1, extra);
        l = i ^ ip1;
        res = res && HSK_8[gluing_int_to_ker[l]][i];
    }
    return -res;
}

void
find_hamilton_path_extra_slow(extra_aut_t *extra, const theta_point_dim4_t *HSK_8, const uint8_t mask)
{
    // Very slow (but robust) method that should not be used.
    // Returns the correct Tesseract automorphism from the Gray Hamilton path
    // to a correct one that does not go through zero.
    // Assumes aut is already initialised.
    extra_aut_t extra2;
    init_extra_aut(&extra2, mask);

    uint32_t sel;
    uint8_t g, m, pad, r;

    uint8_t paths[2] = { (extra2.mask & 1), 2 };
    // test all endomorphisms and valid paths.
    for (r = 0; r < 2; r++) {
        extra2.path = paths[r];
        for (g = 0; g < 24; g++) {
            for (m = 0; m < 4; m++) {
                extra2.sigma2[m] = Sigmas[g][m]; // TODO: reduce to more specific precomputed cases.
            }

            for (extra2.v = 0; extra2.v < 4; extra2.v++) {
                for (pad = 0; pad < 2; pad++) {
                    extra2.sigma1[0] = pad;
                    extra2.sigma1[1] = pad ^ 1;

                    sel = is_correct_path_extra(&extra2, HSK_8);
                    // printf("is_correct_path_extra_slow = %d\n",sel);
                    extra_aut_select(extra, extra, &extra2, sel);
                }
            }
        }
    }
}

// This works but can be improved.
void
find_hamilton_path_extra_fast(extra_aut_t *extra, const theta_point_dim4_t *HSK_8, const uint8_t mask)
{
    // Returns the correct Tesseract automorphism from the Gray Hamilton path
    // to a correct one that does not go through zero.
    // Assumes aut is already initialised.

    theta_point_zero_map maps[5];
    uint8_t g, pad, m, r;
    uint32_t sel;
    extra_aut_t extra2;
    init_extra_aut(&extra2, mask);

    compute_zero_map_extra(maps, HSK_8, mask, extra2.shift);

    uint8_t paths[2] = { (extra2.mask & 1), 2 };
    // test all endomorphisms and valid paths.
    for (r = 0; r < 2; r++) {
        extra2.path = paths[r];
        for (g = 0; g < 6; g++) {
            for (m = 0; m < 4; m++) {
                extra2.sigma2[m] = Sigma_extra_fast[g][m];
            }
            extra2.v = 0;

            for (pad = 0; pad < 2; pad++) {
                extra2.sigma1[0] = pad;
                extra2.sigma1[1] = pad ^ 1;

                sel = is_correct_path_extra_fast(&extra2, maps);
                // printf("is_correct_path_extra_fast = %d\n",sel);
                extra_aut_select(extra, extra, &extra2, sel);
            }
        }
    }
}

void
gen_isog_compute_codomain_extra(theta_struct_dim4_t *codomain,
                                const theta_point_dim4_t *ker,
                                const uint8_t mask,
                                const bool slow)
{
    theta_point_dim4_t HSK_8[5];
    for (uint8_t i = 0; i < 5; i++) {
        squared(&HSK_8[i], &ker[i]);
        inline_hadamard(&HSK_8[i]);
    }

    extra_aut_t extra;
    init_extra_aut(&extra, mask);

    if (slow) { // Usually not needed, used as a caution
        find_hamilton_path_extra_slow(&extra, HSK_8, mask);
    } else {
        find_hamilton_path_extra_fast(&extra, HSK_8, mask);
    }
#ifndef NDEBUG
    printf("Found path\n");
#endif

    compute_codomain_simple_extra_fast(&codomain->inv_dual_null_point, HSK_8, &extra);
    codomain->arith_precomp = false;

#ifndef NDEBUG
    printf("compute_codomain_simple_extra_forloops\n");
#endif
}




// DEBUG FUNCTIONS

void
compute_codomain_simple(theta_point_dim4_t *inv_null_pt_dual, const theta_point_dim4_t *HSK_8, const aut_t *aut)
{
    // Const computations:
    // C(2^g) = (2C(2^{g-1})+2^g+2)M, C(2)=0.
    // C(2^g) = (g2^g - 2) M.
    // C(2^4) = 62 M.
    fp_t top_right, bot_left;
    uint8_t len, pad, mid_up, mid_low, i0, i1, l;
    // INIT LOOP:
    for (pad = 0; pad < 16; pad += 2) {
        i0 = binary_to_pos(pad, aut);
        i1 = binary_to_pos(pad ^ 1, aut);
        l = i0 ^ i1;
        fp_copy(&(*inv_null_pt_dual)[i0], &HSK_8[int_to_ker[l]][i1]);
        fp_copy(&(*inv_null_pt_dual)[i1], &HSK_8[int_to_ker[l]][i0]);
    }

    // MAIN LOOP:
    for (len = 2; len < 16; len <<= 1) {
        for (pad = 0; pad < 16; pad += (len << 1)) {
            mid_up = pad ^ len;
            mid_low = pad ^ (len - 1);
            i0 = binary_to_pos(mid_low, aut);
            i1 = binary_to_pos(mid_up, aut);
            l = i0 ^ i1;

            fp_mul(&top_right, &(*inv_null_pt_dual)[i1], &HSK_8[int_to_ker[l]][i1]);
            fp_mul(&bot_left, &(*inv_null_pt_dual)[i0], &HSK_8[int_to_ker[l]][i0]);

            for (l = pad; l < mid_up; l++) {
                i0 = binary_to_pos(l, aut);
                i1 = binary_to_pos(l ^ len, aut);
                fp_mul(&(*inv_null_pt_dual)[i0], &(*inv_null_pt_dual)[i0], &top_right);
                fp_mul(&(*inv_null_pt_dual)[i1], &(*inv_null_pt_dual)[i1], &bot_left);
            }
        }
    }
}

void
compute_codomain_and_inverse(theta_point_dim4_t *inv_null_pt_dual,
                             theta_point_dim4_t *null_pt_dual,
                             const theta_point_dim4_t *HSK_8,
                             const aut_t *aut)
{

    fp_t inv_top_right, inv_bot_left, top_right, bot_left;
    uint8_t len, pad, mid_up, mid_low, i0, i1, i2, i3, l0, l1, l2;

    // INIT LOOP:
    for (pad = 0; pad < 16; pad += 4) {
        i0 = binary_to_pos(pad, aut);
        i1 = binary_to_pos(pad ^ 1, aut);
        i2 = binary_to_pos(pad ^ 2, aut);
        i3 = binary_to_pos(pad ^ 3, aut);

        l0 = i0 ^ i1;
        l1 = i1 ^ i2;
        l2 = i2 ^ i3;

        // compute inv_null_pt_dual
        fp_mul(&(*inv_null_pt_dual)[i0], &HSK_8[int_to_ker[l1]][i2], &HSK_8[int_to_ker[l2]][i3]);
        fp_mul(&(*inv_null_pt_dual)[i3], &HSK_8[int_to_ker[l0]][i0], &HSK_8[int_to_ker[l1]][i1]);

        fp_mul(&(*inv_null_pt_dual)[i1], &HSK_8[gluing_int_to_ker[l0]][i0], &(*inv_null_pt_dual)[i0]);
        fp_mul(&(*inv_null_pt_dual)[i0], &HSK_8[gluing_int_to_ker[l0]][i1], &(*inv_null_pt_dual)[i0]);
        fp_mul(&(*inv_null_pt_dual)[i2], &HSK_8[gluing_int_to_ker[l2]][i3], &(*inv_null_pt_dual)[i3]);
        fp_mul(&(*inv_null_pt_dual)[i3], &HSK_8[gluing_int_to_ker[l2]][i2], &(*inv_null_pt_dual)[i3]);

        // compute null_pt_dual
        fp_mul(&(*null_pt_dual)[i1], &HSK_8[int_to_ker[l0]][i1], &HSK_8[int_to_ker[l2]][i2]);
        fp_mul(&(*null_pt_dual)[i2], &HSK_8[int_to_ker[l1]][i2], &(*null_pt_dual)[i1]);
        fp_mul(&(*null_pt_dual)[i1], &HSK_8[int_to_ker[l1]][i1], &(*null_pt_dual)[i1]);

        fp_copy(&(*null_pt_dual)[i0], &(*inv_null_pt_dual)[i3]);
        fp_copy(&(*null_pt_dual)[i3], &(*inv_null_pt_dual)[i0]);
    }

    // MAIN LOOP:
    for (len = 4; len < 16; len <<= 1) {
        for (pad = 0; pad < 16; pad += (len << 1)) {
            mid_up = pad ^ len;
            mid_low = pad ^ (len - 1);
            i0 = binary_to_pos(mid_low, aut);
            i1 = binary_to_pos(mid_up, aut);
            l0 = i0 ^ i1;

            fp_mul(&inv_top_right, &(*inv_null_pt_dual)[i1], &HSK_8[int_to_ker[l0]][i1]);
            fp_mul(&top_right, &(*null_pt_dual)[i1], &HSK_8[int_to_ker[l0]][i0]);
            fp_mul(&inv_bot_left, &(*inv_null_pt_dual)[i0], &HSK_8[int_to_ker[l0]][i0]);
            fp_mul(&bot_left, &(*null_pt_dual)[i0], &HSK_8[int_to_ker[l0]][i1]);

            i0 = binary_to_pos(pad, aut);
            i1 = binary_to_pos(mid_low ^ len, aut);
            fp_mul(&(*inv_null_pt_dual)[i0], &(*inv_null_pt_dual)[i0], &inv_top_right);
            fp_mul(&(*inv_null_pt_dual)[i1], &(*inv_null_pt_dual)[i1], &inv_bot_left);
            fp_copy(&(*null_pt_dual)[i1], &(*inv_null_pt_dual)[i0]);
            fp_copy(&(*null_pt_dual)[i0], &(*inv_null_pt_dual)[i1]);

            i0 = binary_to_pos(mid_low, aut);
            i1 = binary_to_pos(mid_up, aut);
            fp_mul(&(*inv_null_pt_dual)[i0], &(*inv_null_pt_dual)[i0], &inv_top_right);
            fp_mul(&(*null_pt_dual)[i0], &(*null_pt_dual)[i0], &top_right);
            fp_mul(&(*inv_null_pt_dual)[i1], &(*inv_null_pt_dual)[i1], &inv_bot_left);
            fp_mul(&(*null_pt_dual)[i1], &(*null_pt_dual)[i1], &bot_left);

            for (l0 = pad + 1; l0 < mid_low; l0++) {
                i0 = binary_to_pos(l0, aut);
                i1 = binary_to_pos(l0 ^ len, aut);
                fp_mul(&(*inv_null_pt_dual)[i0], &(*inv_null_pt_dual)[i0], &inv_top_right);
                fp_mul(&(*null_pt_dual)[i0], &(*null_pt_dual)[i0], &top_right);
                fp_mul(&(*inv_null_pt_dual)[i1], &(*inv_null_pt_dual)[i1], &inv_bot_left);
                fp_mul(&(*null_pt_dual)[i1], &(*null_pt_dual)[i1], &bot_left);
            }
        }
    }
}


void
compute_codomain_simple_extra(theta_point_dim4_t *inv_null_pt_dual,
                                       const theta_point_dim4_t *HSK_8,
                                       const extra_aut_t *extra)
{
    // Const computations:
    // C(2^g) = (2C(2^{g-1})+2^g+2)M, C(2)=0.
    // C(2^g) = (g2^g - 2) M.
    // C(2^4) = 62 M.
    fp_t top_right, bot_left;
    uint8_t len, pad, mid_up, mid_low, i0, i1, l;
    // INIT LOOP:
    for (pad = 0; pad < 16; pad += 2) {
        i0 = binary_to_pos_extra(pad, extra);
        i1 = binary_to_pos_extra(pad ^ 1, extra);
        l = i0 ^ i1;
        fp_copy(&(*inv_null_pt_dual)[i0], &HSK_8[gluing_int_to_ker[l]][i1]);
        fp_copy(&(*inv_null_pt_dual)[i1], &HSK_8[gluing_int_to_ker[l]][i0]);
    }

    // MAIN LOOP:
    for (len = 2; len < 16; len <<= 1) {
        for (pad = 0; pad < 16; pad += (len << 1)) {
            mid_up = pad ^ len;
            mid_low = pad ^ (len - 1);
            i0 = binary_to_pos_extra(mid_low, extra);
            i1 = binary_to_pos_extra(mid_up, extra);
            l = i0 ^ i1;

            fp_mul(&top_right, &(*inv_null_pt_dual)[i1], &HSK_8[gluing_int_to_ker[l]][i1]);
            fp_mul(&bot_left, &(*inv_null_pt_dual)[i0], &HSK_8[gluing_int_to_ker[l]][i0]);

            for (l = pad; l < mid_up; l++) {
                i0 = binary_to_pos_extra(l, extra);
                i1 = binary_to_pos_extra(l ^ len, extra);
                fp_mul(&(*inv_null_pt_dual)[i0], &(*inv_null_pt_dual)[i0], &top_right);
                fp_mul(&(*inv_null_pt_dual)[i1], &(*inv_null_pt_dual)[i1], &bot_left);
            }
        }
    }
}






















/*
void
compute_codomain_and_inverse_recurse(theta_point_dim4_t *inv_null_pt_dual, theta_point_dim4_t *null_pt_dual, const
theta_point_dim4_t *HSK_8, const uint8_t *int_to_ker, const aut_t *aut, const uint8_t p_min, const uint8_t p_max){
        // Const computations:
        // TO DO
        if (p_max-p_min<=4)
        {
                // 9M
        uint8_t i0 = binary_to_pos(p_min, aut);
        uint8_t i1 = binary_to_pos(p_min + 1, aut);
        uint8_t i2 = binary_to_pos(p_min + 2, aut);
        uint8_t i3 = binary_to_pos(p_min + 3, aut);
        uint8_t l0 = i0 ^ i1;
        uint8_t l1 = i1 ^ i2;
        uint8_t l2 = i2 ^ i3;

                fp_mul(&(*inv_null_pt_dual)[i0], &HSK_8[int_to_ker[l1]][i2], &HSK_8[int_to_ker[l2]][i3]);
        fp_mul(&(*inv_null_pt_dual)[i3], &HSK_8[int_to_ker[l0]][i0], &HSK_8[int_to_ker[l1]][i1]);

        fp_mul(&(*inv_null_pt_dual)[i1], &HSK_8[gluing_int_to_ker[l0]][i0], &(*inv_null_pt_dual)[i0]);
        fp_mul(&(*inv_null_pt_dual)[i0], &HSK_8[gluing_int_to_ker[l0]][i1], &(*inv_null_pt_dual)[i0]);
        fp_mul(&(*inv_null_pt_dual)[i2], &HSK_8[gluing_int_to_ker[l2]][i3], &(*inv_null_pt_dual)[i3]);
        fp_mul(&(*inv_null_pt_dual)[i3], &HSK_8[gluing_int_to_ker[l2]][i2], &(*inv_null_pt_dual)[i3]);


                fp_mul(&(*null_pt_dual)[i1], &HSK_8[int_to_ker[l0]][i1], &HSK_8[int_to_ker[l2]][i2]);
                fp_mul(&(*null_pt_dual)[i2], &HSK_8[int_to_ker[l1]][i2], &(*null_pt_dual)[i1]);
                fp_mul(&(*null_pt_dual)[i1], &HSK_8[int_to_ker[l1]][i1], &(*null_pt_dual)[i1]);

                fp_copy(&(*null_pt_dual)[i0],&(*inv_null_pt_dual)[i3]);
                fp_copy(&(*null_pt_dual)[i3],&(*inv_null_pt_dual)[i0]);

    } else {
                uint8_t p_mid = p_min + ((p_max-p_min)>>1);
                compute_codomain_and_inverse_recurse(inv_null_pt_dual,null_pt_dual, HSK_8, int_to_ker, aut, p_min,
p_mid); compute_codomain_and_inverse_recurse(inv_null_pt_dual,null_pt_dual, HSK_8, int_to_ker, aut, p_mid, p_max);


                fp_t inv_top_right, inv_bot_left, top_right, bot_left;
                uint8_t i_mid, i_midm1, l, i, q;

                // Compute top_right factor
                i_mid = binary_to_pos(p_mid, aut);

                fp_copy(&inv_top_right,&(*inv_null_pt_dual)[i_mid]);
                fp_copy(&top_right,&(*null_pt_dual)[i_mid]);

                i_midm1 = binary_to_pos(p_mid-1, aut);
                l = i_mid ^ i_midm1;


                fp_mul(&inv_top_right,&inv_top_right,&HSK_8[int_to_ker[l]][i_mid]);
                fp_mul(&top_right,&top_right,&HSK_8[int_to_ker[l]][i_midm1]);

                // Compute bot_left factor
                fp_copy(&inv_bot_left,&(*inv_null_pt_dual)[i_midm1]);
                fp_copy(&bot_left,&(*null_pt_dual)[i_midm1]);

                fp_mul(&inv_bot_left,&inv_bot_left,&HSK_8[int_to_ker[l]][i_midm1]);
                fp_mul(&bot_left,&bot_left,&HSK_8[int_to_ker[l]][i_mid]);

                // Multiplying by top_right factor
                for(q = p_min + 1 ; q < p_mid ; q++){
                        i = binary_to_pos(q, aut);
                        fp_mul(&(*inv_null_pt_dual)[i],&(*inv_null_pt_dual)[i],&inv_top_right);
                        fp_mul(&(*null_pt_dual)[i],&(*null_pt_dual)[i],&top_right);

                }
                // Multiplying by bot_left factor
                for(q = p_mid ; q < p_max -1 ; q++){
                        i = binary_to_pos(q, aut);
                        fp_mul(&(*inv_null_pt_dual)[i],&(*inv_null_pt_dual)[i],&inv_bot_left);
                        fp_mul(&(*null_pt_dual)[i],&(*null_pt_dual)[i],&bot_left);
                }

                i = binary_to_pos(p_min, aut);
                l = binary_to_pos(p_max-1, aut);

                fp_mul(&(*inv_null_pt_dual)[i],&(*inv_null_pt_dual)[i],&inv_top_right);
                fp_copy(&(*null_pt_dual)[l],&(*inv_null_pt_dual)[i]);
                fp_mul(&(*inv_null_pt_dual)[l],&(*inv_null_pt_dual)[l],&inv_bot_left);
                fp_copy(&(*null_pt_dual)[i],&(*inv_null_pt_dual)[l]);



        }
}


void
compute_codomain_and_inverse(theta_point_dim4_t *inv_null_pt_dual,theta_point_dim4_t *null_pt_dual, const
theta_point_dim4_t *HSK_8, const uint8_t *int_to_ker, const aut_t *aut)
        {
                compute_codomain_and_inverse_recurse(inv_null_pt_dual, null_pt_dual ,HSK_8, int_to_ker, aut, 0, 16);
}



//Fonction no longer used.
void
compute_codomain_simple_recurse(theta_point_dim4_t *inv_null_pt_dual, const theta_point_dim4_t *HSK_8,
        const uint8_t *int_to_ker, const aut_t *aut, const uint8_t p_min, const uint8_t p_max){
        // Const computations:
        // C(2^g) = (2C(2^{g-1})+2^g+2)M, C(2)=0.
        // C(2^g) = (g2^g - 2) M.
        // C(2^4) = 62 M.
        if (p_max-p_min<=2){
                uint8_t i0 = binary_to_pos(p_min, aut);
                uint8_t i1 = binary_to_pos(p_min+1, aut);
                uint8_t l = i0^i1;
                // x1*P0^{-1}(x1)
                fp_copy(&(*inv_null_pt_dual)[i0],&HSK_8[int_to_ker[l]][i1]);
                // x0*P1^{-1}(x0)
                fp_copy(&(*inv_null_pt_dual)[i1],&HSK_8[int_to_ker[l]][i0]);
        }
        else{
                uint8_t p_mid = p_min+((p_max-p_min)>>1);
                compute_codomain_simple_recurse(inv_null_pt_dual, HSK_8, int_to_ker, aut, p_min, p_mid);
                compute_codomain_simple_recurse(inv_null_pt_dual, HSK_8, int_to_ker, aut, p_mid, p_max);

                fp_t top_right, bot_left;
                uint8_t i_mid, i_midm1, l, i;

                // Compute top_right factor
                i_mid = binary_to_pos(p_mid, aut);
                fp_copy(&top_right,&(*inv_null_pt_dual)[i_mid]);

                i_midm1 = binary_to_pos(p_mid-1,aut);
                l = i_mid^i_midm1;
                fp_mul(&top_right,&top_right,&HSK_8[int_to_ker[l]][i_mid]);

                // Compute bot_left factor
                fp_copy(&bot_left,&(*inv_null_pt_dual)[i_midm1]);

                fp_mul(&bot_left,&bot_left,&HSK_8[int_to_ker[l]][i_midm1]);

                // Multiplying by top_right factor
                for(uint8_t q=p_min;q<p_mid;q++){
                        i = binary_to_pos(q, aut);
                        fp_mul(&(*inv_null_pt_dual)[i],&(*inv_null_pt_dual)[i],&top_right);
                }
                // Multiplying by bot_left factor
                for(uint8_t q=p_mid;q<p_max;q++){
                        i = binary_to_pos(q, aut);
                        fp_mul(&(*inv_null_pt_dual)[i],&(*inv_null_pt_dual)[i],&bot_left);
                }
        }
}



void
compute_codomain_simple_old(theta_point_dim4_t *inv_null_pt_dual, const theta_point_dim4_t *HSK_8,
        const uint8_t *int_to_ker, const aut_t *aut){
        compute_codomain_simple_recurse(inv_null_pt_dual, HSK_8, int_to_ker, aut, 0, 16);
}
*/

/* Functions for "special generic" isogenies where one more torsion point is needed
for codomain (e.g. f_2: A_1->A_2 where T_1, T_2, T_3, T_4, T_1 + T2 or T_3+T_4 are needed) */
/*
void
compute_codomain_simple_recurse_extra(theta_point_dim4_t *inv_null_pt_dual, const theta_point_dim4_t *HSK_8,
        const uint8_t *int_to_ker, const extra_aut_t *extra, const uint8_t p_min,
        const uint8_t p_max){
        // Const computations:
        // C(2^g) = (2C(2^{g-1})+2^g+2)M, C(2)=0.
        // C(2^g) = (g2^g - 2) M.
        // C(2^4) = 62 M.
        if (p_max-p_min<=2){
                uint8_t i0 = binary_to_pos_extra(p_min, extra);
                uint8_t i1 = binary_to_pos_extra(p_min+1, extra);
                uint8_t l = i0^i1;
                // x1*P0^{-1}(x1)
                fp_copy(&(*inv_null_pt_dual)[i0],&HSK_8[int_to_ker[l]][i1]);
                // x0*P1^{-1}(x0)
                fp_copy(&(*inv_null_pt_dual)[i1],&HSK_8[int_to_ker[l]][i0]);
        }
        else{
                uint8_t p_mid = p_min+((p_max-p_min)>>1);
                compute_codomain_simple_recurse_extra(inv_null_pt_dual, HSK_8, int_to_ker, extra, p_min, p_mid);
                compute_codomain_simple_recurse_extra(inv_null_pt_dual, HSK_8, int_to_ker, extra, p_mid, p_max);

                fp_t top_right, bot_left;
                uint8_t i_mid, i_midm1, l, i;

                // Compute top_right factor
                i_mid = binary_to_pos_extra(p_mid, extra);
                fp_copy(&top_right,&(*inv_null_pt_dual)[i_mid]);

                i_midm1 = binary_to_pos_extra(p_mid-1, extra);
                l = i_mid^i_midm1;
                fp_mul(&top_right,&top_right,&HSK_8[int_to_ker[l]][i_mid]);

                // Compute bot_left factor
                fp_copy(&bot_left,&(*inv_null_pt_dual)[i_midm1]);

                fp_mul(&bot_left,&bot_left,&HSK_8[int_to_ker[l]][i_midm1]);

                // Multiplying by top_right factor
                for(uint8_t q=p_min;q<p_mid;q++){
                        i = binary_to_pos_extra(q, extra);
                        fp_mul(&(*inv_null_pt_dual)[i],&(*inv_null_pt_dual)[i],&top_right);
                }
                // Multiplying by bot_left factor
                for(uint8_t q=p_mid;q<p_max;q++){
                        i = binary_to_pos_extra(q, extra);
                        fp_mul(&(*inv_null_pt_dual)[i],&(*inv_null_pt_dual)[i],&bot_left);
                }
        }
}

void
compute_codomain_simple_extra(theta_point_dim4_t *inv_null_pt_dual, const theta_point_dim4_t *HSK_8,
        const uint8_t *int_to_ker, const extra_aut_t *extra){
        compute_codomain_simple_recurse_extra(inv_null_pt_dual, HSK_8, int_to_ker, extra, 0, 16);
}
*/
