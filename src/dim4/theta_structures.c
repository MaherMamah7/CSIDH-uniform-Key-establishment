#include <isog_dim4.h>

static inline void
theta_struct_arith_precomp(theta_struct_dim4_t *theta_struct)
{
    // Cost: 3*(3*16-4) M + 16 S + 128 a = 132 M + 16 S + 128 a
    invert_theta_point(&theta_struct->null_point, &theta_struct->inv_dual_null_point);
    hadamard(&theta_struct->null_point, &theta_struct->null_point);

    invert_theta_point(&theta_struct->inv_null_point, &theta_struct->null_point);

    squared(&theta_struct->inv_codomain_dual_null_point_sq, &theta_struct->null_point);
    hadamard(&theta_struct->inv_codomain_dual_null_point_sq, &theta_struct->inv_codomain_dual_null_point_sq);

    invert_theta_point(&theta_struct->inv_codomain_dual_null_point_sq, &theta_struct->inv_codomain_dual_null_point_sq);

    theta_struct->arith_precomp = true;
}

void
theta_double(theta_point_dim4_t *out, const theta_point_dim4_t *in, theta_struct_dim4_t *theta_struct)
{
    // Cost (without precomputation): 32 M + 32 S + 128 a
    if (!theta_struct->arith_precomp) {
        theta_struct_arith_precomp(theta_struct);
    }

    squared(out, in);
    hadamard(out, out);
    squared(out, out);
    dot_prod(out, out, &theta_struct->inv_codomain_dual_null_point_sq);
    hadamard(out, out);
    dot_prod(out, out, &theta_struct->inv_null_point);
}

void
theta_double_iter(theta_point_dim4_t *out, const theta_point_dim4_t *in, theta_struct_dim4_t *theta_struct, int n)
{

    theta_copy(out, in);
    for (int i = 0; i < n; i++) {
        theta_double(out, out, theta_struct);
    }
}
