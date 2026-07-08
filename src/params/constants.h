#include <intbig.h> // Where ibz_t is defined

extern const uint64_t TORSION_EVEN_POWER;//e_2
extern const ibz_t CHARACTERISTIC;
extern const ibz_t TORSION_EVEN;// 2^e_2
extern const ibz_t TORSION_ODD;// c*l^e_l
extern const uint32_t MAX_IDX_STRAT_INIT; // For strategies
extern const uint32_t MAX_IDX_STRAT;

#if defined(TEST_KANI)
extern const uint64_t L;// l
extern const uint64_t TORSION_L_POWER;// e_l
extern const ibz_t TORSION_L;// l^e_l
extern const ibz_t TORSION_COMPLEMENT_L;// c*2^e_2
// a1**2+a2**2+l**fl=2**f2
extern const uint64_t CONST_FL;
extern const uint64_t CONST_F2;
extern const ibz_t CONST_A1;
extern const ibz_t CONST_A2;
#elif defined(SIDH_ATTACK)
extern const uint64_t L;// l
extern const uint64_t TORSION_L_POWER;// e_l
extern const ibz_t TORSION_L;// l^e_l
extern const ibz_t TORSION_COMPLEMENT_L;// c*2^e_2
// a1**2+a2**2+l**fl=2**f2
extern const uint64_t CONST_FL;
extern const uint64_t CONST_F2;
extern const ibz_t CONST_A1;
extern const ibz_t CONST_A2;
#endif

/* Do not write below that line */

#ifdef P_5_2_248
#define NWORDS_FIELD 5
#define NWORDS_ORDER 4
#define FP_ENCODED_BYTES 40
#endif
#ifdef P_239_2_194_3_107
#define NWORDS_FIELD 7
#define NWORDS_ORDER 6
#define FP_ENCODED_BYTES 56
#endif
#ifdef P_p500
#define NWORDS_FIELD 9
#define NWORDS_ORDER 8
#define FP_ENCODED_BYTES 72
#endif
#ifdef P_500
#define NWORDS_FIELD 9
#define NWORDS_ORDER 8
#define FP_ENCODED_BYTES 72
extern const uint16_t STRATEGY[495];
#endif
#ifdef P_1000
#define NWORDS_FIELD 17
#define NWORDS_ORDER 16
#define FP_ENCODED_BYTES 136
extern const uint16_t STRATEGY[999];
#endif
#ifdef P_1500
#define NWORDS_FIELD 26
#define NWORDS_ORDER 24
#define FP_ENCODED_BYTES 208
extern const uint16_t STRATEGY[1517];
#endif
#ifdef P_2000
#define NWORDS_FIELD 34
#define NWORDS_ORDER 32
#define FP_ENCODED_BYTES 272
extern const uint16_t STRATEGY[2009];
#endif
#ifdef P_4000
#define NWORDS_FIELD 68
#define NWORDS_ORDER 63
#define FP_ENCODED_BYTES 544
extern const uint16_t STRATEGY[4019];
#endif
