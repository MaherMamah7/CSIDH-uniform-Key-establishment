/** @file
 *
 * @author(s) Pierrick Dartois, Max Duparc
 *
 * @brief Dim 4 isogenies using theta structure.  
 */

#ifndef ISOG_DIM4_H
#define ISOG_DIM4_H

#include <fp.h>
#include <ec.h>
#include <stdint.h>

extern const uint8_t gluing_int_to_ker[16];
extern const uint8_t gluing_bin_to_pos[4][16];
extern const uint8_t gluing_theta_change_variable_positions[4][16];
extern const uint8_t second_iso_bin_to_pos[3][16]; 
extern const uint8_t ind_to_ker[16];
extern const uint8_t splitting_automorphism[12];

/** @defgroup dim4 Dimension 4 Isogenies
 * @{
 */


/** @defgroup Struct Structures
 * @{
 */

/** @brief Projective theta point in dimension 4
 *
 * @typedef theta_point_dim4_t
 *
 * @struct theta_point_dim4_t
 *
 * Given by 16 coordinates in level 2
 */
typedef fp_t theta_point_dim4_t[16];

/** @brief Projective theta point in dimension 2
 *
 * @typedef theta_point_dim2_t
 *
 * @struct theta_point_dim2_t
 *
 * Given by 4 coordinates in level 2
 */
typedef fp_t theta_point_dim2_t[4];

/** @brief Projective theta point in dimension 1
 *
 * @typedef theta_point_dim1_t
 *
 * @struct theta_point_dim1_t
 *
 * Given by 2 coordinates in level 2
 */
typedef fp_t theta_point_dim1_t[2];

/** @brief Theta structure in dimension 4
 *
 * @typedef theta_struct_dim4_t
 *
 * @struct theta_struct_dim4_t
 */
typedef struct theta_struct_dim4_t
{
    theta_point_dim4_t inv_dual_null_point; // Stores the inverse dual theta null point
                                            // by default (for isogeny evaluations).
                                            // This saves a batch inversion for every
                                            // isogeny computation
    bool arith_precomp;
    theta_point_dim4_t null_point;
    theta_point_dim4_t inv_null_point;
    theta_point_dim4_t inv_codomain_dual_null_point_sq;
} theta_struct_dim4_t;

/** @brief Quadruple of points
 *
 * @typedef quad_point_t
 *
 * @struct quad_point_t
 */
typedef struct jac_point_t quad_point_t[4];

/** @brief Gluing isogeny E^4 --> A, where E is an elliptic curve
 *
 * @typedef gluing_isog_dim4_t
 *
 * @struct gluing_isog_dim4_t
 */
typedef struct gluing_isog_dim4_t
{
    theta_point_dim4_t theta_product_domain;
    quad_point_t auxiliary_point;
    theta_point_dim4_t inv_codomain_auxiliary_point;
    uint8_t basis_swap;
    uint8_t position_codomain;        // Type of the codomain graph.
    uint8_t position_theta_structure; // Change of theta coordinates
} gluing_isog_dim4_t;

/** @brief Composition of the two first isogenies f_1: E^4 --> A_1
 * and f_2: A_1 --> A_2
 *
 * @typedef first_isogenies_t
 *
 * @struct first_isogenies_t
 */
typedef struct first_isogenies_t
{
    ec_curve_t domain_EC[2]; // Stores the starting elliptic curve domain_EC[0] 
    // and its twist domain_EC[1]
    gluing_isog_dim4_t f1;
    theta_struct_dim4_t codomain_f2;
} first_isogenies_t;


/** @brief A 5-tuple of point that describe the symplectic basis of order 16 above the kernel of our gluing isogeny.  
 * (See section 4.1 and Appendix A of https://eprint.iacr.org/2025/1859 for more details. )
 * @typedef basis_components_16_t
 *
 * @struct basis_components_16_t
 */
typedef struct basis_components_16_t
{
    uint8_t N;
    uint8_t s1;
    uint8_t s2;
    uint8_t s3;
    uint8_t s4;
} basis_components_16_t;

/** @brief Encodes automorphisms of the Tesseract
 * (hypercube of dim 4) in order to determine
 * deviations from the canonical Hamiltonian path
 * (Gray spanning tree).
 *
 * @typedef aut_t
 *
 * @struct aut_t
 *
 * aut: v --> (v o sigma)^s. Contains data to encode
 * aut^{-1}: w --> (w^s) o sigma^{-1}
 */
typedef struct aut_t
{
    // Table encoding a permutation of {0,...,3}
    uint8_t sigma[4];
    // integer to xor with. Useless because Hamilton
    // path are invariant by rotation
    // unsigned int s;

    // Indicates if the inverse of sigma is computed
    // bool is_inverted;
    // Inverse of sigma
    // uint8_t sigmam1[4];
} aut_t;


/** @brief Encodes extra automorphisms of the
 * enhanced Tesseract (with e_1+e_2 or e_3+e_4
 * added to the edges).
 * (see https://www.pmf.ni.ac.rs/filomat-content/2020/34-8/34-8-25-11857.pdf,
 * Theorem 3.6)
 *
 * @typedef extra_aut_t
 *
 * @struct extra_aut_t
 */
typedef struct extra_aut_t
{
    // mask for the support
    // (e.g. 3 for support {0,1} i.e 1100)
    uint8_t mask;
    // Starting index of the support
    // (e.g. 0 for {0,1} i.e. 1100)
    uint8_t shift;
    // permutation of indices ({0,1}) outside the support
    uint8_t sigma1[2];
    // vector translate outside the support
    uint8_t v;
    // permutation on the support on (Z/2Z)^2
    // (written as a permutation of {0,...,3}
    // in binary representation)
    uint8_t sigma2[4];
    // integer designing which path to use
    uint8_t path; 

    // Indicates if the inverse of sigma is computed
    // bool is_inverted;
    // Inverse of sigma2 (sigma1 is its own inverse)
    // uint8_t sigma2m1[4];
} extra_aut_t;

/** @brief 0-1 description of a 4 torsion point above the kernel and of the topology of the graph underlying the HIIP problem. 
 *
 * @typedef theta_point_zero_map
 *
 * @struct theta_point_zero_map
 *
 * Incidence maps by 16 coordinates in level 2
 */
typedef bool theta_point_zero_map[16];


/** // end Structures
 * @}
 */

/** @defgroup utility Trees, paths and graphs
 * @{
 */




/**
 * @brief maps an int (in binary) to its respective index in the gray spanning tree.
 *
 * @param bin (input) an integer in binary. 
 * @return its index position in the gray spanning tree
 *
 */
static inline uint8_t
binary_to_gray_pos(const uint8_t bin)
{
    return bin ^ (bin >> 1);
}

/**
 * @brief maps an index in the gray spanning tree to its respective integer in binary.
 *
 * @param pos (input) the position of an index in the gray spanning tree 
 * @return the corresponding number (in binary).  
 *
 */
static inline uint8_t
gray_pos_to_binary(const uint8_t pos)
{
    uint8_t bin = pos, mask = pos;
    for (uint8_t i = 0; i < 4; i++) {
        mask >>= 1;
        bin ^= mask;
    }
    return bin;
}


/**
 * @brief given an automorphism, compute its inverse.
 *
 * @param inverse (output) a list corresponding to the inverse of the automorphism. 
 * @param aut (input) the position of an index in the gray spanning tree 
 */
static inline void
invert_aut(uint8_t *inverse ,aut_t const *aut)
{

    //if (!(aut->is_inverted)) {
    for (uint8_t i = 0; i < 4; i++) {
        inverse[aut->sigma[i]] = i;
    }
    //}
}

/**
 * @brief maps an int (in binary) to its respective index in the gray spanning tree postcomposed by an automorphism of the Tesseract.
 *
 * @param bin (input) an integer in binary. 
 * @param aut (input) an automorphism of the tesseract. 
 * @return its index position in the gray spanning tree after application of the automophism
 */
static inline uint8_t
binary_to_pos(const uint8_t bin, const aut_t *aut)
{
    uint8_t pos1, pos2;
    pos1 = binary_to_gray_pos(bin);
    pos2 = 0;
    for (uint8_t i = 0; i < 4; i++) {
        pos2 ^= ((pos1 & (1 << (aut->sigma[i]))) >> (aut->sigma[i])) << i;
    }

    return pos2;
}

/**
 * @brief DEBUG:inverse function of binary_to_pos.
 *
 * @param pos (input) the position of an index in the gray spanning tree (up to automorphism). 
 * @param aut (input) an automorphism of the tesseract. 
 * @return its index position in the gray spanning tree after application of the automophism
 */
static inline uint8_t
pos_to_binary(const uint8_t pos, const aut_t *aut)
{
    uint8_t sigmam1[4];
    invert_aut(sigmam1 ,aut);

    uint8_t bin = 0;
    // bin1 = bin^(aut->s);
    for (uint8_t i = 0; i < 4; i++) {
        bin ^= ((pos & (1 << (sigmam1[i]))) >> (sigmam1[i])) << i;
    }
    return gray_pos_to_binary(bin);
}

/**
 * @brief init an automorphism to the identity.
 *
 * @param aut (input/output) the automophism to init. 
 */
static inline void
init_aut(aut_t *aut)
{
    // aut->s = 0;
    for (uint8_t i = 0; i < 4; i++) {
        aut->sigma[i] = i;
    }
}

/**
 * @brief select an automorphism between 2.
 *
 * @param aut (output) the selected automophism. 
 * @param aut1 (input) the first automophism. 
 * @param aut2 (input) the second automophism. 
 * @param sel (input) the slection parameter. 
 */
static inline void
aut_select(aut_t *aut, const aut_t *aut1, const aut_t *aut2, const bool sel)
{
    uint8_t mask = -sel;
    for (uint8_t i = 0; i < 4; i++) {
        aut->sigma[i] = (aut1->sigma[i]) ^ (((aut1->sigma[i]) ^ (aut2->sigma[i])) & mask);
    }
}


/**
 * @brief given an extra-automorphism, compute its inverse.
 *
 * @param inverse (output) a list corresponding to the inverse of the automorphism. 
 * @param aut (input) the position of an index in the gray spanning tree 
 */
static inline void
invert_extra_aut(uint8_t *inverse ,extra_aut_t const *aut)
{
    //if (!(aut->is_inverted)) {
    for (uint8_t i = 0; i < 4; i++) {
        inverse[aut->sigma2[i]] = i;
    }
    //aut->is_inverted = true;
    //}
}

/**
 * @brief maps an int (in binary) to its respective index in the gray spanning tree postcomposed by an automorphism of the extended Tesseract.
 *
 * @param bin (input) an integer in binary. 
 * @param extra (input) an automorphism of the extended tesseract. 
 * @return its index position in the gray spanning tree after application of the automophism
 */
static inline uint8_t
binary_to_pos_extra(const uint8_t bin, const extra_aut_t *extra){

    uint8_t pos1, pos2, pos3, pos4, pos5;
    pos2 = second_iso_bin_to_pos[extra->path][bin]; 
    //binary_to_gray_pos(bin);

    // Support
    pos1 = pos2 & (extra->mask);
    pos3 = extra->sigma2[pos1 >> (extra->shift)] << (extra->shift);


    // Out of support
    pos4 = (pos2 & (15 ^ (extra->mask))) >> (2-(extra->shift));
    pos5 = 0;
    for (uint8_t i = 0; i < 2; i++) {
        pos5 ^= ((pos4 & (1 << (extra->sigma1[i]))) >> (extra->sigma1[i])) << i;
    }
    pos5 = (pos5 ^ (extra->v)) << (2-(extra->shift));
    return pos3 ^ pos5;
}

/**
 * @brief DEBUG: inverse function of binary_to_pos_extra.
 *
 * @param pos (input) the position of an index in the gray spanning tree (up to automorphism). 
 * @param extra (input) an automorphism of the extended tesseract. 
 * @return its index position in the gray spanning tree after application of the automophism
 */
static inline uint8_t
pos_to_binary_extra(const uint8_t pos, extra_aut_t *extra)
{

    uint8_t sigma2m1[4];
    invert_extra_aut(sigma2m1, extra);
    uint8_t bin1, bin2, bin3, bin4, bin5;
    uint8_t sel = 0;

    // Support
    bin1 = pos & (extra->mask);
    bin3 = sigma2m1[bin1 >> (extra->shift)] << (extra->shift);

    // Out of support
    bin4 = (pos & (15 ^ (extra->mask))) >> (2-(extra->shift));
    bin5 = 0;
    for (uint8_t i = 0; i < 2; i++) {
        bin5 ^= ((bin4 & (1 << (extra->sigma1[i]))) >> (extra->sigma1[i])) << i;
    }
    bin5 = ( bin5 ^ (extra->v)) << (2-(extra->shift));
    bin2 = bin3 ^ bin5;


    bin3 = 0;
    for(bin1 = 0; bin1 < 16; bin1++)
    {
        sel = -(uint8_t)(second_iso_bin_to_pos[extra->path][bin1] == bin2);
        bin3 ^=  ((bin1) & sel);

    }

    return bin3;

}

/**
 * @brief init an automorphism to the identity.
 *
 * @param aut (input/output) the automophism to init. 
 * @param mask (input) either 3 or 12, defines the type of extended tesseract. 
 */
static inline void
init_extra_aut(extra_aut_t *aut, const uint8_t mask)
{
    aut->mask = mask;
    aut->shift = (!(mask & 1)) + (!((mask & 2) >> 1));
    aut->v = 0; 
    aut->path =0;

    aut -> sigma1[0] = 0;
    aut -> sigma1[1] = 1;
    
    for (uint8_t i = 0; i < 4; i++) {
        aut -> sigma2[i] = i;
    }
}

/**
 * @brief select an automorphism of the extended automorphism between 2.
 *
 * @param aut (output) the selected automophism. 
 * @param aut1 (input) the first automophism. 
 * @param aut2 (input) the second automophism. 
 * @param sel (input) the slection parameter. 
 */
static inline void
extra_aut_select(extra_aut_t *aut, const extra_aut_t *aut1, const extra_aut_t *aut2, const bool sel)
{
    uint8_t mask = -sel;
    aut->mask = (aut1->mask) ^ (((aut1->mask) ^ (aut2->mask)) & mask);
    aut->shift = (aut1->shift) ^ (((aut1->shift) ^ (aut2->shift)) & mask);
    aut->v = (aut1->v) ^ (((aut1->v) ^ (aut2->v)) & mask);
    aut->path = (aut1->path) ^ (((aut1->path) ^ (aut2->path)) & mask);

    aut->sigma1[0] = (aut1->sigma1[0]) ^ (((aut1->sigma1[0]) ^ (aut2->sigma1[0])) & mask);
    aut->sigma1[1] = (aut1->sigma1[1]) ^ (((aut1->sigma1[1]) ^ (aut2->sigma1[1])) & mask);


    for (uint8_t i = 0; i < 4; i++) {
        aut->sigma2[i] = (aut1->sigma2[i]) ^ (((aut1->sigma2[i]) ^ (aut2->sigma2[i])) & mask);
    }
}

/** //  utility Trees, paths and graphs
 * @}
 */

/** @defgroup  Theta-structures Theta structures
 * @{
 */

/** @defgroup  Theta-points Theta points
 * @{
 */

/**
 * @brief Set all coordinates in a theta_point_dim4_t to zero.
 *
 * @param out (input/output) the theta_point_dim4_t point.  
 *
 */
static inline void
init_theta_dim_4(theta_point_dim4_t *out)
{
    for (uint8_t i = 0; i < 16; i++) {
        fp_set_zero(&(*out)[i]);
    }
}

/**
 * @brief Hadamard transform on a vector of 2 elements.
 *
 * @param out (output) The hadamard transformed vector.   
 * @param in (input) The input vector of 2 elements.  
 *
 */
static inline void
hadamard_2(fp_t *out, const fp_t *in)
{
    fp_add(&out[0], &in[0], &in[1]);
    fp_sub(&out[1], &in[0], &in[1]);
}

/**
 * @brief Hadamard transform on a vector of 4 elements.
 *
 * @param out (output) The hadamard transformed vector.   
 * @param in (input) The input vector of 4 elements.  
 *
 */
static inline void
hadamard_4(fp_t *out, const fp_t *in)
{
    fp_t in1[2], in2[2], out1[2], out2[2];

    fp_copy(&in1[0], &in[0]);
    fp_copy(&in1[1], &in[1]);
    fp_copy(&in2[0], &in[2]);
    fp_copy(&in2[1], &in[3]);

    hadamard_2(out1, in1);
    hadamard_2(out2, in2);

    fp_add(&out[0], &out1[0], &out2[0]);
    fp_add(&out[1], &out1[1], &out2[1]);
    fp_sub(&out[2], &out1[0], &out2[0]);
    fp_sub(&out[3], &out1[1], &out2[1]);
}

/**
 * @brief Hadamard transform on a vector of 8 elements.
 *
 * @param out (output) The hadamard transformed vector.   
 * @param in (input) The input vector of 8 elements.  
 *
 */
static inline void
hadamard_8(fp_t *out, const fp_t *in)
{
    fp_t in1[4], in2[4], out1[4], out2[4];

    for (uint8_t i = 0; i < 4; i++) {
        fp_copy(&in1[i], &in[i]);
        fp_copy(&in2[i], &in[i + 4]);
    }

    hadamard_4(out1, in1);
    hadamard_4(out2, in2);

    for (uint8_t i = 0; i < 4; i++) {
        fp_add(&out[i], &out1[i], &out2[i]);
        fp_sub(&out[i + 4], &out1[i], &out2[i]);
    }
}

/**
 * @brief Hadamard transform on thetapoint of dim 4.
 *
 * @param out (output) The hadamard transformed thetapoint.   
 * @param in (input) The input thetapoint of dim 4.  
 *
 */
static inline void
hadamard(theta_point_dim4_t *out, const theta_point_dim4_t *in)
{
    fp_t in1[8], in2[8], out1[8], out2[8];

    for (uint8_t i = 0; i < 8; i++) {
        fp_copy(&in1[i], &(*in)[i]);
        fp_copy(&in2[i], &(*in)[i + 8]);
    }

    hadamard_8(out1, in1);
    hadamard_8(out2, in2);

    for (uint8_t i = 0; i < 8; i++) {
        fp_add(&(*out)[i], &out1[i], &out2[i]);
        fp_sub(&(*out)[i + 8], &out1[i], &out2[i]);
    }
}

/**
 * @brief inline Hadamard transform based on Fast Walsh–Hadamard–Fourier Transform
 *  More memory efficient. 
 *
 * @param out (input/output) The vector to hadamard transform.   
 */
static inline void
inline_hadamard(theta_point_dim4_t *out)
{   
    uint8_t len, pad, l, i;
    fp_t x; 
    for(len = 1; len < 16; len <<= 1)
    {
        for(pad = 0; pad < 16; pad += (len<<1))
        {
            for(l = pad; l < (pad ^ len); l++)
            {
                i = l ^ len;
                fp_copy(&x, &(*out)[l]);
                fp_add(&(*out)[l], &(*out)[l], &(*out)[i]); 
                fp_sub(&(*out)[i], &x, &(*out)[i]); 

            }
        }
    }
}

/**
 * @brief square all coefficient of thetapoint of dim 4.
 *
 * @param out (output) The squared thetapoint   
 * @param in (input) The input thetapoint of dim 4.  
 *
 */
static inline void
squared(theta_point_dim4_t *out, const theta_point_dim4_t *in)
{
    for (uint8_t i = 0; i < 16; i++) {
        fp_sqr(&(*out)[i], &(*in)[i]);
    }
}

/**
 * @brief perform dot product of 2 thetapoints of dim 4.
 *
 * @param out (output) The product thetapoint.   
 * @param in1 (input) The first thetapoint of dim 4.  
 * @param in2 (input) The second thetapoint of dim 4. 
 *
 */
static inline void
dot_prod(theta_point_dim4_t *out, const theta_point_dim4_t *in1, const theta_point_dim4_t *in2)
{
    for (uint8_t i = 0; i < 16; i++) {
        fp_mul(&(*out)[i], &(*in1)[i], &(*in2)[i]);
    }
}

/**
 * @brief Copy a thetapoint to another. 
 *
 * @param out (output) where to copy the thetapoint.   
 * @param in (input) The thetapoint to copy.  
 *
 */
static inline void
theta_copy(theta_point_dim4_t *out, const theta_point_dim4_t *in)
{
    for (uint8_t i = 0; i < 16; i++) {
        fp_copy(&(*out)[i], &(*in)[i]);
    }
}

/**
 * @brief invert a theta point. 
 *
 * @param out (output) where to write the inverted thetapoint.   
 * @param in (input) The thetapoint to invert.  
 *
 */
static inline void
invert_theta_point(theta_point_dim4_t *out, const theta_point_dim4_t *in)
{
    fp_t inv[16];
    for(uint8_t i=0; i<16; i++){
        fp_copy(&inv[i],&(*in)[i]);
    }
    fp_proj_batched_inv(inv, 16);
    for(uint8_t i=0; i<16; i++){
        fp_copy(&(*out)[i],&inv[i]);
    }
}
/**
 * @brief init a theta point to 0. 
 *
 * @param P (input/output) The thetapoint to initialise.  
 *
 */
static inline void
theta_init(theta_point_dim4_t *P){
    for(uint8_t i=0; i<16; i++){
        fp_set_zero(&(*P)[i]);
    }
}

/**
 * @brief Copy a theta-structure. 
 *
 * @param out (output) where to copy the theta-structure.  
 * @param in (input)  the theta-structure to copy  
 *
 */
static inline void
copy_theta_struct(theta_struct_dim4_t *out, const theta_struct_dim4_t *in){
    theta_copy(&out->inv_dual_null_point,&(in->inv_dual_null_point));
    theta_copy(&out->null_point,&in->null_point);
    theta_copy(&out->inv_null_point,&in->inv_null_point);
    theta_copy(&out->inv_codomain_dual_null_point_sq,&in->inv_codomain_dual_null_point_sq);
    out->arith_precomp = in->arith_precomp;
}

/** // end utility functions for theta structures
 * @}
 */

 /** @defgroup Algebra Algebra
 * @{
 */

void theta_double(theta_point_dim4_t *out, const theta_point_dim4_t *in, theta_struct_dim4_t *theta_struct);

void theta_double_iter(theta_point_dim4_t *out, const theta_point_dim4_t *in, theta_struct_dim4_t *theta_struct, int n);

/** // end Algebra
 * @}
 */
 /** // end Theta-structures
 * @}
 */


/** @defgroup quad-points Tuple of 4 points
 * @{
 */

/**
 * @brief Double all coefficients of a quad_point. 
 *
 * @param out (output) the doubled quad_point.   
 * @param P (input) The input quad_point.  
 * @param E (input) The curve on which P is defined.   
 *
 */
static inline void
quad_point_DBL(quad_point_t *out, const quad_point_t *P, const ec_curve_t *E)
{
    for (uint8_t i = 0; i < 4; i++) {
        DBL(&(*out)[i], &(*P)[i], E);
    }
}

/**
 * @brief Iterate doubling of all coefficients of a quad_point. 
 *
 * @param out (output) the doubled quad_point.   
 * @param P (input) The input quad_point. 
 * @param n (input) Number of iteration. 
 * @param E (input) The curve on which P is defined.
 *    
 * out = [2^n]*P.
 */
static inline void quad_point_dbl_iter(quad_point_t *out, const quad_point_t *P, const unsigned int  n, ec_curve_t *E)
{
    for (uint8_t i = 0; i < 4; i++) {
        jac_dbl_iter(&(*out)[i], &(*P)[i], n, E);
    }
}

/**
 * @brief Add all coefficients of two quad_points. 
 *
 * @param out (output) the added quad_point.   
 * @param P (input) The first input quad_point.  
 * @param Q (input) The second input quad_point.  
 * @param E (input) The curve on which P and Q are defined.   
 *
 * @note Make sure out is not P or Q, as this can lead to errors. 
 */
static inline void
quad_point_ADD(quad_point_t *out, const quad_point_t *P, const quad_point_t *Q, const ec_curve_t *E)
{
    for (uint8_t i = 0; i < 4; i++) {
        ADD(&(*out)[i], &(*P)[i], &(*Q)[i], E);
    }
}

/**
 * @brief negate all coefficients of a quad_point
 *
 * @param out (output) the negated quad_point.   
 * @param P (input) The input quad_point.  
 *
 */
static inline void
quad_point_neg(quad_point_t *out, const quad_point_t *P)
{
    for (uint8_t i = 0; i < 4; i++) {
        jac_neg(&(*out)[i], &(*P)[i]);
    }
}

/**
 * @brief Copy a quad_point to another. 
 *
 * @param out (output) where to copy the quad_point.   
 * @param in (input) The quad_point to copy.  
 *
 */
static inline void quad_point_copy(quad_point_t *out, const quad_point_t *in)
{
    for (uint8_t i = 0; i < 4; i++) {
        copy_jac_point(&(*out)[i], &(*in)[i]);
    }
}

/** // quad-Points
 * @}
 */


/** @defgroup isogenies dimension 4 isogenies
 * @{
 */

/**
 * @brief Given B = < T_0,T_1,T_2,T_3 > a maximal isotropic subgroup of A[8] expressed in a consistent theta structure, compute the codomain of the isogeny f: A --> B
 * with ker(f) is generated by  [4]T_0,[4]T_1,[4]T_2,[4]T_3.   
 *
 * @param codomain (output) The codomain of the isogeny f: A --> B of kernel [4]B.
 * @param ker (input) {T_0,...,T_3} a maximal isotropic subgroup of A[8] expressed in a consistent theta structure.
 * @param search_path (input) a bool indicating wether or not we test several possible Hamiltonian path to compute the codomain.
 *
 */
void gen_isog_compute_codomain(theta_struct_dim4_t *codomain,
                               const theta_point_dim4_t *ker,
                               bool search_path);
                               
/**
 * @brief Given P in consistent theta coordinates and the codomain of f, compute f(P).     
 *
 * @param imP (output) The image of f(P). 
 * @param P (input) a point P in theta coordinates consistent with f. 
 * @param codomain (input) the codomain of f.  
 *
 */
void gen_isog_eval(theta_point_dim4_t *imP, const theta_point_dim4_t *P, const theta_struct_dim4_t *codomain);

/**
 * @brief Given B = < T_0,T_1,T_2,T_3> a maximal isotropic subgroup of A[8] expressed in a consistent theta structure, compute the codomain of the isogeny f: A --> B
 * with ker(f) is generated by  [4]T_0,[4]T_1,[4]T_2,[4]T_3. This functions is specific to the 2nd isogeny of qt-Pegasis and takes an additional fith point
 * T_4 = T_0 + T_2 or T_4 = T_2 + T_3. 
 *
 * @param codomain (output) The codomain of the isogeny f: A --> B of kernel [4]B.
 * @param ker (input) {T_0,...,T_4} a maximal isotropic subgroup of A[8] expressed in a consistent theta structure.
 * @param mask (input) an indicator of the nature of T_4. ( mask = 3 --> T_4 = T_0 + T_1/ mask = 12 --> T_4 = T_2 + T_3)
 * @param search_path (input) a bool indicating wether or not we test several possible Hamiltonian path to compute the codomain.
 */
void gen_isog_compute_codomain_extra(theta_struct_dim4_t *codomain,
                                     const theta_point_dim4_t *ker,
                                     const uint8_t mask, 
                                     const bool search_path
                                    );


/**
 * @brief DEBUG: function that solves the HIIP in the generic case and compute both codomain and its inverse. Used for debug.  
 *
 * @param inv_null_pt_dual (output) The inverse null theta point of the codomain of f.
 * @param null_pt_dual (output) The null theta point of the codomain of f.
 * @param HSK_8 (input) the hadamard of the square of the 8 points above the kernel of f. 
 * @param aut (input) an automorphism of the Tesseract.
 */
void
compute_codomain_and_inverse(theta_point_dim4_t *inv_null_pt_dual,theta_point_dim4_t *null_pt_dual, const theta_point_dim4_t *HSK_8,  const aut_t *aut); 

/**
 * @brief DEBUG: function that solves the HIIP in the generic case and compute the inverse of the codomain. slower method no longer used.   
 * 
 * @param inv_null_pt_dual (output) The null theta point of the codomain of f.
 * @param HSK_8 (input) the hadamard of the square of the 8 points above the kernel of f. 
 * @param aut (input) an automorphism of the Tesseract.
 */
void compute_codomain_simple(theta_point_dim4_t *inv_null_pt_dual, const theta_point_dim4_t *HSK_8, const aut_t *aut);

/**
 * @brief DEBUG: function that solves the HIIP in the generic case and compute the inverse of the codomain. Specialised function used in the second isogeny of qt-Pegasis
 *  slower method no longer used.   
 * 
 * @param inv_null_pt_dual (output) The null theta point of the codomain of f.
 * @param HSK_8 (input) the hadamard of the square of the 8 points above the kernel of f. 
 * @param extra (input) an automorphism of the extended Tesseract.
 */
void
compute_codomain_simple_extra(theta_point_dim4_t *inv_null_pt_dual, const theta_point_dim4_t *HSK_8, const extra_aut_t *extra);


/**
 * @brief function that finds an automorphism of the extented Tesseract such that all edges are non-zero. Used to compute the HIIP.  
 * 
 * @param extra (output) an valid automorphism of the extended Tesseract.  
 * @param HSK_8 (input) the hadamard of the square of the 8 points above the kernel of f. 
 * @param mask (input) an indicator of the nature of the extended Tesseract.
 */
void
find_hamilton_path_extra_fast(extra_aut_t *extra, 
	const theta_point_dim4_t *HSK_8, const uint8_t mask);


/**
 * @brief function that solves the HIIP in the generic case and compute the inverse of the codomain.   
 * 
 * @param inv_null_pt_dual (output) The null theta point of the codomain of f.
 * @param HSK_8 (input) the hadamard of the square of the 8 points above the kernel of f. 
 * @param aut (input) an automorphism of the extended Tesseract.
 */
void
compute_codomain_simple_fast(theta_point_dim4_t *inv_null_pt_dual, const theta_point_dim4_t *HSK_8, const aut_t *aut);


/** // end isogeny functions
 * @} 
 */

/** @defgroup gluing gluing isogenies
 *  @{
 */

/**
 * @brief Select the adequate precomputed change of theta coordinates used for gluing.
 *
 * @param compo (input) integers N,s1,s2,s3,s4 specifying the 4-torsion basis used to generate the curve.
 * @return position_theta_structure:  an uint from 0,...,3 detailing which change of theta coordintes we shall apply.
 *
 */
uint8_t compute_theta_position(const basis_components_16_t compo);

/**
 * @brief Select the adequate precomputed spanning tree in the codomain graph used to precomputed domain.
 *
 * @param compo (input) integers N,s1,s2,s3,s4 specifying the 4-torsion basis used to generate the curve.
 * @return position_codomain: an uint from 0,...,3 detailing which shape does the codomain graph has.
 *
 */
uint8_t compute_codomain_position(const basis_components_16_t compo);

/**
 * @brief compute all gluing components useful for gluing
 *
 * @param gluing_isogeny (output) all data necesarry to use the gluing isogeny.
 * @param T16 (input) the 16 torsion symplectic basis <T1, T2, T3, T4> above the kernel of the gluing isogeny
 * @param P_16 (input) the 16 torsion point above the 2 torsion defined over E
 * @param Q_16 (input) the 16 torsion point above the 2 torsion defined over E^t
 * @param domain_EC (input) the curve and its twist (E and E^t).
 * @param position_codomain (input) an uint from 0,...,3 detailing where the zeros are in the codomain graph.
 * @param position_theta_structure (input) an uint from 0,...,3 detailing which change of theta coordintes we shall
 * apply.
 *
 */
void gluing_compute(gluing_isog_dim4_t *gluing_isogeny,
                    const quad_point_t *T16, // 4 elements
                    const jac_point_t *P_16,
                    const jac_point_t *Q_16,
                    const ec_curve_t *domain_EC, // 2 elements
                    const uint8_t position_codomain,
                    const uint8_t position_theta_structure);


/**
 * @brief: Evaluate a basis thought a gluing isogeny
 *
 * @param out (output) computed image of the point in the dual theta structure
 * @param in   (input) a basis of 4 4-tuple of points on E^4 (for 0,1) and its twist (for 2,3).
 * @param gluing_isogeny (input) a gluing isogeny precomputed using gluing_compute.
 * @param domain_EC (input) the curve E0.
 */
void gluing_eval_basis(theta_point_dim4_t *out,
                          const quad_point_t *in, //4 elements
                          const gluing_isog_dim4_t *gluing_isogeny,
                          const ec_curve_t *domain_EC 
                        );

/**
 * @brief DEBUG: Evaluate a point thought a gluing isogeny
 *
 * @param out (output) computed image of the point in the dual theta structure
 * @param P   (input) a 4-tuple of points on E^4 or its twist.
 * @param gluing_isogeny (input) a gluing isogeny precomputed using gluing_compute.
 * @param E (input) the curve E0.
 * @param twistP (input) Wether P is defined on the E (twist_P = 0) or on the twist (twist_P = 1).
 */
void gluing_eval(theta_point_dim4_t *out,
                 const quad_point_t *P,
                 const gluing_isog_dim4_t *gluing_isogeny,
                 const ec_curve_t *E, // Always E0
                 bool twistP);


/**
 * @brief Special function to evaluate T1 and T2 when the auxiliary point used to evaluate gluing isogenies is T1 + T2.
 * This prevents a unavoidable technicality where barycentric coordinates are not well defined, while beeing also slightly faster. 
 *
 * @param out1 computed image of T1.
 * @param out2 computed image of T2.
 * @param T1 a 4-tuple of points.
 * @param T2 a 4-tuple of points.
 * @param gluing_isogeny a gluing isogeny precomputed using gluing_compute.
 * @param E the curve where both T1, T2 and T1+T2 are defined.
 */

void gluing_eval_special_basis(theta_point_dim4_t *out1,
                               theta_point_dim4_t *out2,
                               const quad_point_t *T1,
                               const quad_point_t *T2,
                               const gluing_isog_dim4_t *gluing_isogeny,
                               const ec_curve_t *E // The curve upon P1, P2 and P_aux are defined.
                                );

/**
 * @brief Compute the codomain of the gluing isogeny. More specifically, For a given P, compute Phi(P) using its product
 * with 5 specific 8-torsion points.
 *
 * @param out (output) computed codomain.
 * @param HSK8 (input) List (P+T_i)(P-T_i) where T_i are 8 torsion points above the kernel of the form [T_1,T_2,T_3,T_4,
 * T_(1+2) or T_(3+4)] depending on position.
 * @param position (input) an uint from 0,...,3 detailing where the zeros are in the HSK8.
 *
 */
void gluing_compute_codomain(theta_point_dim4_t *out, const theta_point_dim4_t *HSK8, const uint8_t position);


/**
 * @brief Compute the product of theta-null point of the dim 1 theta structure generated by (P_4, Q_4) in the product of
 * Kummerline on E^4.
 *
 * @param product_theta_point (output) computed product.
 * @param P (output) The 4 torsion point P defined on E. It must not be of the form (1:*:1) or (-1:*:1).
 * @param Q (output) The 4 torsion point Q defined on E^t. It can be of the form (1:*:1) or (-1:*:1).
 * @return basis_swap: An uint_8 either 0 or 15, used to perfom conditional bitswitch depending on the form of Q.
 *
 */
uint8_t compute_product_theta_coord(theta_point_dim4_t *product_theta_point, const ec_point_t *P, const ec_point_t *Q);

/**
 * @brief Apply the precomputed change of theta coordinate to go from a product theta structure to desired theta
 * structure.
 *
 * @param out (output/input) theta point on the product theta structure --> theta point on the desired theta structure
 * @param position_theta_change_variable (input) An integer between 0,...,3 specifying which change of varianble must be
 * applied.
 *
 */
void apply_change_coords_theta(theta_point_dim4_t *out, const uint8_t position_theta_change_variable);

/**
 * @brief Transform a vector P in the product Kummerlines to the desired theta structure on E^4.
 *
 * @param out (output/input) kronecker product of points over E^4 --> theta point on the desired theta structure.
 * @param product_theta_domain (input) kronecker product of theta structure of dim 1.
 * @param basis_swap (input) An uint_8 either 0 or 15, used to perfom conditional bitswitch.
 * @param position_theta_change_variable( input) An integer between 0,...,3 specifying which change of varianble must be
 * applied.
 *
 */
void dim4_vector_to_theta_point(theta_point_dim4_t *out,
                                const theta_point_dim4_t *product_theta_domain,
                                const uint8_t basis_swap,
                                const uint8_t position_theta_change_variable);

/**
 * @brief Transform a quad_point T over E^4 or (E^t)^4 into a dim 4 vector in the product Kummerlines.
 *
 * @param out (output) The vector P in the product Kummerlines.
 * @param T (input) A quad_point.
 * @param twist (input)a bool specifying if it is on the twist or not
 *
 */
void quad_point_to_dim4_vec(theta_point_dim4_t *out, const quad_point_t *T, const bool twist);

/**
 * @brief Given P, Q, compute 2 vectors real and im such that real - i_Pi_Q * im =  P+Q and real + i_Pi_Q * im = P-Q in
 * the product of kummerlines.
 * @param real_vec (output) the real part of the vector.
 * @param im_vec (output) the imaginary part of the vector. 
 * @param P (input) the quad_point P.
 * @param Q (input) the quad_point Q.
 * @param E (input) the curve E (never the twist).
 * @param twistP (input) A bool saying if P is on on twist or not.
 * @param twistQ (input) A bool saying if Q is on on twist or not.
 *
 */
void ec_couplepoints_to_dim4_vectors(theta_point_dim4_t *real_vec,
                                     theta_point_dim4_t *im_vec,
                                     const quad_point_t *P,
                                     const quad_point_t *Q,
                                     const ec_curve_t *E,
                                     const bool twistP,
                                     const bool twistQ);

/**
 * @brief Given P, Q, compute both (P+Q)(P-Q) in the theta structure.
 * @param out (output) (P+Q)(P-Q).
 * @param P (input) the quad_point P.
 * @param Q (input) the quad_point Q.
 * @param E (input) the curve E (never the twist).
 * @param product_theta_domain (input) The product of theta-null points of dimension 1. 
 * @param twistP (input) A bool saying if P is on on twist or not.
 * @param twistQ (input) A bool saying if Q is on on twist or not.
 * @param basis_swap (input) An uint_8 either 0 or 15, used to perfom conditional bitswitch.
 * @param position_theta_change_variable (input) An integer between 0,...,3 specifying which change of varianble must be
 * applied.
 *
 */
void couplepoint_to_product_of_alt_sum_theta_point(theta_point_dim4_t *out,
                                                   const quad_point_t *P,
                                                   const quad_point_t *Q,
                                                   const ec_curve_t *E,
                                                   const theta_point_dim4_t *product_theta_domain,
                                                   const bool twistP,
                                                   const bool twistQ,
                                                   const uint8_t basis_swap,
                                                   const uint8_t position_theta_change_variable);


/**
 * @brief DEBUG: bitpermutations used for the computation of the codomain. Applies the following bit-premutations.
 *      0 --> id,
 *      1 --> (01)(23),
 *      2 --> (02)(13),
 *      3 --> (03)(12)
 * @param n the integer to bitpermute.
 * @param position a int between 0 to 3.
 * @return  the bitpermuted integer.
 *
 */
uint8_t gluing_isomorphism_codomain(const uint8_t n, const uint8_t position);


/**
 * @brief DEBUG: Given coefficients defining a valid symplectic basis, precompute the associate symplectic basis of
 * order 4.
 *
 * @param out (output) The symplectic basis generated from our basis coefficients. It is of the form <T_1, T_2, T_3,
 * T_4> + <S_1,S_2,S_3,S_4> with the S point conditionally computed.
 * @param P (input) the 4 torsion point in P in E jacobain form.
 * @param Q (input) the 4 torsion point in Q in E^t jacobain form.
 * @param curve (input) the curve E.
 * @param twist (input) the curve E^t.
 * @param compo (input) the coefficient defining the basis.
 * @param full (input) a boolean specifying wether we compute the S points or not.
 *
 */
void construct_symplectic_basis_4(quad_point_t *out,
                                  const jac_point_t *P,
                                  const jac_point_t *Q,
                                  const ec_curve_t *curve,
                                  const ec_curve_t *twist,
                                  const basis_components_16_t compo,
                                  const bool full);

/**
 * @brief DEBUG: Given coefficients defining a valid symplectic basis, precompute the associate symplectic basis of
 * order 16.
 *
 * @param out (output) The symplectic basis generated from our basis coefficients. It is of the form <T_1, T_2, T_3,
 * T_4>.
 * @param P (input) the 16 torsion point in P in E jacobain form.
 * @param Q (input) the 16 torsion point in Q in E^t jacobain form.
 * @param curve (input) the curve E.
 * @param twist (input) the curve E^t.
 * @param compo the coefficient defining the basis.
 *
 */
void construct_symplectic_basis_16(quad_point_t *out,
                                   const jac_point_t *P,
                                   const jac_point_t *Q,
                                   const ec_curve_t *curve,
                                   const ec_curve_t *twist,
                                   const basis_components_16_t compo);

/** // end gluing isogenies
 * @}
 */

/** @defgroup isogeny chains
 * @{
 */

/**
 * @brief Computes the two first isogenies of the chain f1: E^4 --> A1 (gluing)
 * and f2: A1 --> A2
 *
 * @param first_isogenies (output) result f2*f1
 * @param domain_EC (input) domain curve E and its twist
 * @param T16 (input) basis of 16-torsion points above the kernel of f2*f1 (4 points).
 * @param P_16 (input) 16-torsion point [2^{e-2}]P from the basis of E[2^f].
 * @param Q_16 (input) 16-torsion point [2^{e-2}]Q from the basis of E[2^f] (on the twist).
 * @param position_codomain (input) an uint from 0,...,3 detailing where the zeros are in 
 * the codomain graph.
 * @param position_theta_structure (input) an uint from 0,...,3 detailing which change of 
 * theta coordintes we shall apply.
 */
void 
compute_first_isogenies(first_isogenies_t *first_isogenies, const ec_curve_t *domain_EC, 
    const quad_point_t *T16, const jac_point_t *P_16, const jac_point_t *Q_16, 
    const uint8_t position_codomain, const uint8_t position_theta_structure);

/**
 * @brief Evaluates the two first isogenies of the chain f1: E^4 --> A1 (gluing)
 * and f2: A1 --> A2
 *
 * @param out (output) result f2*f1(P)
 * @param P (input) point to evaluate
 * @param twist (input) a boolean which is 1 if P is on the twist and 0 otherwise
 * @param first_isogenies (input) the two first isogenies f2*f1
 */
void
first_isogenies_eval(theta_point_dim4_t *out, const quad_point_t *P, const bool twist, 
    const first_isogenies_t *first_isogenies);

/**
 * @brief Computes the codomain of the 4d 2-isogeny chain E^4 --> A1 --> ... --> Ae
 *
 * @param codomain (output) resulting codomain
 * @param ker (input) basis of 2^{e+2}-torsion points above the kernel (4 points)
 * @param first_isogenies (input) the precomputed two first isogenies f2*f1
 */
void
compute_chain_codomain(theta_struct_dim4_t *codomain, const quad_point_t *ker, 
    first_isogenies_t *first_isogenies);

/** // end isogeny chains
 * @}
 */

/** @defgroup splitting
 * @{
 */

/**
 * @brief Finds the position of the splitting change of theta coordinates matrix
 *
 * @param N1 integer mod 4 (from the norm equation solution).
 * @param Nb1 integer mod 4 (from the norm equation solution).
 * @param A1 integer mod 4 (from the norm equation solution).
 * @param A2 integer mod 4 (from the norm equation solution).
 * 
 * @return the index of the splitting change of theta coordinates matrix (between 0 and 2).
 */
uint8_t compute_theta_position_splitting(const uint8_t N1,const uint8_t Nb1,const uint8_t A1,
    const uint8_t A2);

/**
 * @brief Extracts the curves Ea, Eabar from the codomain theta structure Ea x Eabar x A
 *
 * @param Ea elliptic curve to extract
 * @param Eabar elliptic curve to extract
 * @param codomain the codomain theta structure
 * @param position_theta_change_variable the position of the splitting change 
 * of theta coordinates matrix
 */
void codomain_to_curves(ec_curve_t *Ea, ec_curve_t *Eabar, const theta_struct_dim4_t *codomain, 
    const uint8_t position_theta_change_variable);

/** // end splitting
 * @}
 */

/** //end Dimension 4 Isogenies
 * @} 
 */

#endif
