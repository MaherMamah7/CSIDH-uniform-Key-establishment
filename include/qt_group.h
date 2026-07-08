#ifndef QT_GROUP_H
#define QT_GROUP_H
#include <inttypes.h>

typedef struct qt_group_elem qt_group_elem_t;
typedef struct qt_set_elem qt_set_elem_t;

// initializers/finalizers
void qt_group_elem_init(qt_group_elem_t *g);
void qt_group_elem_finalize(qt_group_elem_t *g);
void qt_set_elem_init(qt_set_elem_t *s);
void qt_set_elem_finalize(qt_set_elem_t *s);

// getting instatiations
void qt_group_get_neutral(qt_group_elem_t *g);
int qt_group_sample_uniformly(qt_group_elem_t *g);
void qt_set_starting_point(qt_set_elem_t *s);
int32_t qt_group_near_generator(qt_group_elem_t *gen);

// checking validity of normal form
int qt_group_is_valid(const qt_group_elem_t *g);
int qt_set_is_valid(const qt_set_elem_t *s);

// tests
int qt_group_is_neutral(const qt_group_elem_t *g);
int qt_group_is_equal(const qt_group_elem_t *a, const qt_group_elem_t *b);
int qt_set_is_equal(const qt_set_elem_t *a, const qt_set_elem_t *b);

// transformations
int qt_set_twist(qt_set_elem_t *result, const qt_set_elem_t *s);
void qt_group_act(qt_set_elem_t *result, const qt_set_elem_t *s, const qt_group_elem_t *g);
void qt_group_multiply(qt_group_elem_t *result, const qt_group_elem_t *a, const qt_group_elem_t *b);
void qt_group_inverse(qt_group_elem_t *inverse, const qt_group_elem_t *g);

// optimizations
void qt_group_power(qt_group_elem_t *result, const qt_group_elem_t *a, const uint32_t n);
void qt_sample_and_act(qt_set_elem_t *result, qt_group_elem_t *sampled, const qt_set_elem_t *s);
#endif
