/**
 * @file    avxcommon.h
 * @author  Cagri Balkesen <cagri.balkesen@inf.ethz.ch>
 * @date    Tue Dec 11 18:24:10 2012
 * @version $Id $
 * 
 * @brief   Common AVX code, kernels etc. used by implementations.
 * 
 * 
 */
#ifndef AVXCOMMON_H
#define AVXCOMMON_H

//#include <immintrin.h> /* AVX intrinsics */
#include    <x86intrin.h>

namespace multiAttrSort {

/* just to enable compilation with g++ */
#if defined(__cplusplus)
#undef restrict
#define restrict __restrict__
#endif

typedef struct block4_64  {int64_t val[4]; } block4_64;
typedef struct block4_32  {int32_t val[4]; } block4_32;

typedef struct block8_64  {int64_t val[8]; } block8_64;
typedef struct block8_32  {int32_t oid[8]; } block8_32;

typedef struct block16_64 {int64_t val[16];} block16_64;
typedef struct block16_32 {int32_t oid[16];} block16_32;
typedef struct block16_16 {int16_t val[16];} block16_16;

typedef struct block64_32 {int32_t val[64];} block64_32;

typedef struct block256_16 {int16_t val[256];} block256_16;
typedef struct block256_32 {int32_t oid[256];} block256_32;

/** 
 * There are 2 ways to implement branches: 
 *     1) With conditional move instr.s using inline assembly (IFELSEWITHCMOVE).
 *     2) With software predication (IFELSEWITHPREDICATION).
 *     3) With normal if-else
 */
#define IFELSEWITHCMOVE       0
#define IFELSEWITHPREDICATION 1
#define IFELSEWITHNORMAL      0

/** the given address are in int32_t*, load and convert them to 64-bank-size*/
#define LOAD8_32T64(REGL, REGH, ADDR)                                         \
    do {                                                                \
        REGL = _mm256_cvtepi32_pd(_mm_load_si128((__m128i const*) ADDR));                   \
        REGH = _mm256_cvtepi32_pd(_mm_load_si128((__m128i const*)(((block4_32 *)ADDR) + 1)));  \
    } while(0)

/** Load 2 AVX 256-bit registers from the given address */
#define LOAD8_64(REGL, REGH, ADDR)                                         \
    do {                                                                \
        REGL = _mm256_load_pd((double const *) ADDR);                   \
        REGH = _mm256_load_pd((double const *)(((block4_64 *)ADDR) + 1));  \
    } while(0)

/** Load unaligned 2 AVX 256-bit registers from the given address */
#define LOAD8U_64(REGL, REGH, ADDR)                                        \
    do {                                                                \
        REGL = _mm256_loadu_pd((double const *) ADDR);                  \
        REGH = _mm256_loadu_pd((double const *)(((block4_64 *)ADDR) + 1)); \
    } while(0)

#define LOAD8U_32T64(REGL, REGH, ADDR)                                         \
    do {                                                                \
        REGL = _mm256_cvtepi32_pd(_mm_loadu_si128((__m128i const*) ADDR));                   \
        REGH = _mm256_cvtepi32_pd(_mm_loadu_si128((__m128i const*)(((block4_32 *)ADDR) + 1)));  \
    } while(0)

/** convert REGL, REGH to 32-bank-size, and store to the given address with uint32_t* */
#define STORE8_64T32(ADDR, REGL, REGH)                                    \
    do {                                                            \
        _mm_store_si128((__m128i *) ADDR, _mm256_cvtpd_epi32(REGL));                     \
        _mm_store_si128((__m128i *)(((block4_32 *) ADDR) + 1), _mm256_cvtpd_epi32(REGH));   \
    } while(0)

/** Store 2 AVX 256-bit registers to the given address */
#define STORE8_64(ADDR, REGL, REGH)                                    \
    do {                                                            \
        _mm256_store_pd((double *) ADDR, REGL);                     \
        _mm256_store_pd((double *)(((block4_64 *) ADDR) + 1), REGH);   \
    } while(0)

/** Store unaligned 2 AVX 256-bit registers to the given address */
#define STORE8U_64(ADDR, REGL, REGH)                                   \
    do {                                                            \
        _mm256_storeu_pd((double *) ADDR, REGL);                    \
        _mm256_storeu_pd((double *)(((block4_64 *) ADDR) + 1), REGH);  \
    } while(0)

#define STORE8U_64T32(ADDR, REGL, REGH)                                    \
    do {                                                            \
        _mm_storeu_si128((__m128i *) ADDR, _mm256_cvtpd_epi32(REGL));                     \
        _mm_storeu_si128((__m128i *)(((block4_32 *) ADDR) + 1), _mm256_cvtpd_epi32(REGH));   \
    } while(0)

#if 0
/**
 * @note Reversing 64-bit values in an AVX register. It will be possible with
 * single _mm256_permute4x64_pd() instruction in AVX2.
 */
#define REVERSE(REG)                                    \
    do {                                                \
        /* first reverse each 128-bit lane */           \
        REG = _mm256_permute_pd(REG, 0x5);              \
        /* now shuffle 128-bit lanes */                 \
        REG = _mm256_permute2f128_pd(REG, REG, 0x1);    \
    } while(0)
#endif

#define REVERSE_64(REG)                                    \
    do {                                                \
        REG = _mm256_permute4x64_pd(REG, 0x1b);              \
    } while(0)

/* reverse the whole register with 16-bank-size */
#define REVERSE_16(REG)                                    \
    do {                                                \
		/* first reverse each low 64-bit lane and high 64-bit lane*/           \
		REG = _mm256_shufflelo_epi16(REG, 0x1b);  		\
		REG = _mm256_shufflehi_epi16(REG, 0x1b);  		\
														\
        /* now shuffle 64-bit lanes */                 \
		REG = _mm256_permute4x64_epi64(REG, 0x1b);    \
    } while(0)

/* reverse the whole register with 32-bank-size */
#define REVERSE_32(REG)                                    \
    do {                                                \
    	/* first reverse each 128-bit lane */           \
        REG = _mm256_shuffle_epi32(REG, 0x1b);              \
        /* now shuffle 128-bit lanes */                 \
		REG = _mm256_permute2f128_si256(REG, REG, 0x1);    \
    } while(0)

/** Bitonic merge kernel for 2 x 4 elements after the reversing step. */
#define BITONIC4_64(O1, O2, XO1, XO2, A, B, XA, XB)                     \
    do {                                                                \
        /* Level-1 comparisons */                                       \
		__m256i mask1 = _mm256_cmpgt_epi64(_mm256_castpd_si256(A), _mm256_castpd_si256(B)); \
        __m256d l1 = _mm256_blendv_pd(A, B, _mm256_castsi256_pd(mask1)); \
        __m256d h1 = _mm256_blendv_pd(B, A, _mm256_castsi256_pd(mask1)); \
        																\
		__m256d l1_oid = _mm256_blendv_pd(XA, XB, _mm256_castsi256_pd(mask1)); \
		__m256d h1_oid = _mm256_blendv_pd(XB, XA, _mm256_castsi256_pd(mask1));  \
																		\
        /* Level-1 shuffles */                                          \
        __m256d l1p = _mm256_permute2f128_pd(l1, h1, 0x20);             \
        __m256d h1p = _mm256_permute2f128_pd(l1, h1, 0x31);             \
        																\
        __m256d l1p_oid = _mm256_permute2f128_pd(l1_oid, h1_oid, 0x20);             \
        __m256d h1p_oid = _mm256_permute2f128_pd(l1_oid, h1_oid, 0x31);             \
                                                                        \
        /* Level-2 comparisons */                                       \
		__m256i mask2 = _mm256_cmpgt_epi64(_mm256_castpd_si256(l1p), _mm256_castpd_si256(h1p)); \
        __m256d l2 = _mm256_blendv_pd(l1p, h1p, _mm256_castsi256_pd(mask2)); \
        __m256d h2 = _mm256_blendv_pd(h1p, l1p, _mm256_castsi256_pd(mask2)); \
        																\
		__m256d l2_oid = _mm256_blendv_pd(l1p_oid, h1p_oid, _mm256_castsi256_pd(mask2)); \
		__m256d h2_oid = _mm256_blendv_pd(h1p_oid, l1p_oid, _mm256_castsi256_pd(mask2));  \
                                                                        \
        /* Level-2 shuffles */                                          \
        __m256d l2p = _mm256_shuffle_pd(l2, h2, 0x0);                   \
        __m256d h2p = _mm256_shuffle_pd(l2, h2, 0xF);                   \
        																\
        __m256d l2p_oid = _mm256_shuffle_pd(l2_oid, h2_oid, 0x0);                   \
        __m256d h2p_oid = _mm256_shuffle_pd(l2_oid, h2_oid, 0xF);                   \
                                                                        \
        /* Level-3 comparisons */                                       \
		__m256i mask3 = _mm256_cmpgt_epi64(_mm256_castpd_si256(l2p), _mm256_castpd_si256(h2p)); \
        __m256d l3 = _mm256_blendv_pd(l2p, h2p, _mm256_castsi256_pd(mask3)); \
        __m256d h3 = _mm256_blendv_pd(h2p, l2p, _mm256_castsi256_pd(mask3)); \
        																\
		__m256d l3_oid = _mm256_blendv_pd(l2p_oid, h2p_oid, _mm256_castsi256_pd(mask3)); \
		__m256d h3_oid = _mm256_blendv_pd(h2p_oid, l2p_oid, _mm256_castsi256_pd(mask3));  \
                                                                        \
        /* Level-3 shuffles implemented with unpcklps unpckhps */       \
        /* AVX cannot shuffle both inputs from same 128-bit lane */     \
        /* so we need 2 more instructions for this operation. */        \
        __m256d l4 = _mm256_unpacklo_pd(l3, h3);                        \
        __m256d h4 = _mm256_unpackhi_pd(l3, h3);                        \
        O1 = _mm256_permute2f128_pd(l4, h4, 0x20);                      \
        O2 = _mm256_permute2f128_pd(l4, h4, 0x31);                      \
        																\
        __m256d l4_oid = _mm256_unpacklo_pd(l3_oid, h3_oid);                        \
        __m256d h4_oid = _mm256_unpackhi_pd(l3_oid, h3_oid);                        \
        XO1 = _mm256_permute2f128_pd(l4_oid, h4_oid, 0x20);                      \
        XO2 = _mm256_permute2f128_pd(l4_oid, h4_oid, 0x31);                      \
    } while(0)


/** Bitonic merge network for 2 x 8 elements without reversing B */
#define BITONIC8(O1, O2, O3, O4, A1, A2, B1, B2)                        \
    do {                                                                \
        /* Level-0 comparisons */                                       \
        __m256d l11 = _mm256_min_pd(A1, B1);                            \
        __m256d l12 = _mm256_min_pd(A2, B2);                            \
        __m256d h11 = _mm256_max_pd(A1, B1);                            \
        __m256d h12 = _mm256_max_pd(A2, B2);                            \
                                                                        \
        BITONIC4(O1, O2, l11, l12);                                     \
        BITONIC4(O3, O4, h11, h12);                                     \
    } while(0)


/** Bitonic merge kernel for 2 x 4 elements of int64_t */
#define BITONIC_MERGE4(O1, O2, XO1, XO2, A, B, XA, XB)                                    \
    do {                                                                \
        /* reverse the order of input register B */                     \
        REVERSE_64(B);                                                     \
        REVERSE_64(XB);                                                     \
        BITONIC4_64(O1, O2, XO1, XO2, A, B, XA, XB);                                         \
    } while(0)

/** Bitonic merge kernel for 2 x 16 elements of int16_t, X-- refer to those for oids */
/** Note that here XB2 and XB1 is already reversed from input; so just reverse each of them is fine */
#define BITONIC_MERGE16_16(OL, OH, XOL1, XOL2, XOH1, XOH2, A, B, XA1, XA2, XB1, XB2)  \
        do {                                            \
            /* reverse the order of input B */          \
            REVERSE_16(B);                             \
            REVERSE_32(XB1);                             \
            REVERSE_32(XB2);                             \
                                                        \
			/* Level-0 comparisons */                   \
            __m256i l0 = _mm256_min_epi16(A, B);        \
            __m256i h0 = _mm256_max_epi16(A, B);        \
            											\
			__m256i mask0 = _mm256_cmpeq_epi16(A, l0);			\
			__m256 low = reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask0, 0))); \
			__m256 high = reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask0, 1))); \
														\
			__m256 ol0_1 = _mm256_blendv_ps(reinterpret_cast<__m256>(XB1), reinterpret_cast<__m256>(XA1), low);											\
			__m256 ol0_2 = _mm256_blendv_ps(reinterpret_cast<__m256>(XB2), reinterpret_cast<__m256>(XA2), high);										\
														\
			__m256 oh0_1 = _mm256_blendv_ps(reinterpret_cast<__m256>(XA1), reinterpret_cast<__m256>(XB1), low);											\
			__m256 oh0_2 = _mm256_blendv_ps(reinterpret_cast<__m256>(XA2), reinterpret_cast<__m256>(XB2), high);										\
														\
			/* Level-0 shuffle */                  											\
			__m256i l0_sf = _mm256_permute2f128_si256(l0, h0, 0x20);					\
			__m256i h0_sf = _mm256_permute2f128_si256(l0, h0, 0x31);					\
														\
			__m256 ol0_1_sf = ol0_1;										\
			__m256 ol0_2_sf = oh0_1;										\
			__m256 oh0_1_sf = ol0_2;										\
			__m256 oh0_2_sf = oh0_2;										\
														\
														\
            /* Level-1 comparisons */                   \
            __m256i l1 = _mm256_min_epi16(l0_sf, h0_sf);        \
            __m256i h1 = _mm256_max_epi16(l0_sf, h0_sf);        \
            											\
			__m256i mask1 = _mm256_cmpeq_epi16(l0_sf, l1);			\
			low = reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask1, 0))); \
			high = reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask1, 1))); \
														\
			__m256 ol1_1 = _mm256_blendv_ps(oh0_1_sf, ol0_1_sf, low);											\
			__m256 ol1_2 = _mm256_blendv_ps(oh0_2_sf, ol0_2_sf, high);										\
														\
			__m256 oh1_1 = _mm256_blendv_ps(ol0_1_sf, oh0_1_sf, low);											\
			__m256 oh1_2 = _mm256_blendv_ps(ol0_2_sf, oh0_2_sf, high);										\
																						\
			/* Level-1 shuffle */                   \
			__m256i l1_sf = _mm256_unpacklo_epi64(l1, h1);					\
			__m256i h1_sf = _mm256_unpackhi_epi64(l1, h1);					\
															\
			__m256 ol1_1_sf = _mm256_permute2f128_ps(ol1_1, oh1_1, 0x20);					\
			__m256 ol1_2_sf	= _mm256_permute2f128_ps(ol1_2, oh1_2, 0x20);\
																			\
			__m256 oh1_1_sf = _mm256_permute2f128_ps(ol1_1, oh1_1, 0x31);					\
			__m256 oh1_2_sf	= _mm256_permute2f128_ps(ol1_2, oh1_2, 0x31);	\
														\
			/* Level-2 comparisons */                  \
        	__m256i l2 = _mm256_min_epi16(l1_sf, h1_sf);        \
        	__m256i h2 = _mm256_max_epi16(l1_sf, h1_sf);        \
        														\
			__m256i mask2 = _mm256_cmpeq_epi16(l1_sf, l2); \
			low = reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask2, 0))); \
			high = reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask2, 1))); \
														\
			__m256 ol2_1 = _mm256_blendv_ps(oh1_1_sf, ol1_1_sf, low);											\
			__m256 ol2_2 = _mm256_blendv_ps(oh1_2_sf, ol1_2_sf, high);										\
														\
			__m256 oh2_1 = _mm256_blendv_ps(ol1_1_sf, oh1_1_sf, low);											\
			__m256 oh2_2 = _mm256_blendv_ps(ol1_2_sf, oh1_2_sf, high);										\
																						\
			/* Level-2 shuffle */                   \
			l2 = _mm256_shuffle_epi32(l2, 0xd8);					\
			h2 = _mm256_shuffle_epi32(h2, 0xd8);					\
			__m256i l2_sf = _mm256_unpacklo_epi32(l2, h2);					\
			__m256i h2_sf = _mm256_unpackhi_epi32(l2, h2);					\
																			\
			__m256 ol2_1_sf = reinterpret_cast<__m256>(_mm256_unpacklo_epi64(reinterpret_cast<__m256i>(ol2_1), reinterpret_cast<__m256i>(oh2_1)));					\
			__m256 ol2_2_sf = reinterpret_cast<__m256>(_mm256_unpacklo_epi64(reinterpret_cast<__m256i>(ol2_2), reinterpret_cast<__m256i>(oh2_2)));					\
			\
			__m256 oh2_1_sf = reinterpret_cast<__m256>(_mm256_unpackhi_epi64(reinterpret_cast<__m256i>(ol2_1), reinterpret_cast<__m256i>(oh2_1)));					\
			__m256 oh2_2_sf = reinterpret_cast<__m256>(_mm256_unpackhi_epi64(reinterpret_cast<__m256i>(ol2_2), reinterpret_cast<__m256i>(oh2_2)));					\
																		\
			/* Level-3 comparisons */                  \
			__m256i l3 = _mm256_min_epi16(l2_sf, h2_sf);        \
			__m256i h3 = _mm256_max_epi16(l2_sf, h2_sf);        \
				\
			__m256i mask3 = _mm256_cmpeq_epi16(l2_sf, l3); \
			low = reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask3, 0))); \
			high = reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask3, 1))); \
														\
			__m256 ol3_1 = _mm256_blendv_ps(oh2_1_sf, ol2_1_sf, low);											\
			__m256 ol3_2 = _mm256_blendv_ps(oh2_2_sf, ol2_2_sf, high);										\
														\
			__m256 oh3_1 = _mm256_blendv_ps(ol2_1_sf, oh2_1_sf, low);											\
			__m256 oh3_2 = _mm256_blendv_ps(ol2_2_sf, oh2_2_sf, high);										\
			\
			/* Level-3 shuffle */                   				\
			/* step 1: shuffle each 64 lane of high register*/					\
			__m256i h3_tbb = _mm256_shufflelo_epi16(h3, 0xb1);					\
			h3_tbb = _mm256_shufflehi_epi16(h3_tbb, 0xb1);					\
			/* step 2: shuffle each 64 lane of low register*/					\
			__m256i l3_tbb = _mm256_shufflelo_epi16(l3, 0xb1);					\
			l3_tbb = _mm256_shufflehi_epi16(l3_tbb, 0xb1);					\
			/* step 3: do the blend */									\
			__m256i l3_sf = _mm256_blend_epi16(l3, h3_tbb, 0xaa);					\
			__m256i h3_sf = _mm256_blend_epi16(l3_tbb, h3, 0xaa);					\
			\
			ol3_1 = reinterpret_cast<__m256>(_mm256_shuffle_epi32(reinterpret_cast<__m256i>(ol3_1), 0xd8));					\
			ol3_2 = reinterpret_cast<__m256>(_mm256_shuffle_epi32(reinterpret_cast<__m256i>(ol3_2), 0xd8));					\
			oh3_1 = reinterpret_cast<__m256>(_mm256_shuffle_epi32(reinterpret_cast<__m256i>(oh3_1), 0xd8));					\
			oh3_2 = reinterpret_cast<__m256>(_mm256_shuffle_epi32(reinterpret_cast<__m256i>(oh3_2), 0xd8));					\
			__m256 ol3_1_sf = _mm256_unpacklo_ps(ol3_1, oh3_1);					\
			__m256 ol3_2_sf = _mm256_unpacklo_ps(ol3_2, oh3_2);					\
			__m256 oh3_1_sf = _mm256_unpackhi_ps(ol3_1, oh3_1);					\
			__m256 oh3_2_sf = _mm256_unpackhi_ps(ol3_2, oh3_2);					\
			\
			/* Level-4 comparisons */                  \
			__m256i l4 = _mm256_min_epi16(l3_sf, h3_sf);        \
			__m256i h4 = _mm256_max_epi16(l3_sf, h3_sf);        \
			\
			__m256i mask4 = _mm256_cmpeq_epi16(l3_sf, l4); \
			low = reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask4, 0))); \
			high = reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask4, 1))); \
														\
			__m256 ol4_1 = _mm256_blendv_ps(oh3_1_sf, ol3_1_sf, low);											\
			__m256 ol4_2 = _mm256_blendv_ps(oh3_2_sf, ol3_2_sf, high);										\
														\
			__m256 oh4_1 = _mm256_blendv_ps(ol3_1_sf, oh3_1_sf, low);											\
			__m256 oh4_2 = _mm256_blendv_ps(ol3_2_sf, oh3_2_sf, high);										\
			\
			/* Level-4 shuffle */                  \
			__m256i l4_sf = _mm256_permute4x64_epi64(l4, 0xd8);					\
			__m256i h4_sf = _mm256_permute4x64_epi64(h4, 0xd8);					\
			OL = _mm256_unpacklo_epi16(l4_sf, h4_sf);	\
			OH = _mm256_unpackhi_epi16(l4_sf, h4_sf);	\
			                  \
			__m256i ol4_1_sf = _mm256_permute4x64_epi64(reinterpret_cast<__m256i>(ol4_1), 0xd8);					\
			__m256i ol4_2_sf = _mm256_permute4x64_epi64(reinterpret_cast<__m256i>(ol4_2), 0xd8);					\
			__m256i oh4_1_sf = _mm256_permute4x64_epi64(reinterpret_cast<__m256i>(oh4_1), 0xd8);					\
			__m256i oh4_2_sf = _mm256_permute4x64_epi64(reinterpret_cast<__m256i>(oh4_2), 0xd8);					\
			\
			XOL1 = _mm256_unpacklo_epi32(ol4_1_sf, oh4_1_sf);	\
			XOL2 = _mm256_unpackhi_epi32(ol4_1_sf, oh4_1_sf);	\
			XOH1 = _mm256_unpacklo_epi32(ol4_2_sf, oh4_2_sf);	\
			XOH2 = _mm256_unpackhi_epi32(ol4_2_sf, oh4_2_sf);	\
        } while(0)

/** Bitonic merge kernel for 2 x 8 elements of int32_t, X-- refer to those for oids */
#define BITONIC_MERGE8_32(O1, O2, XO1, XO2, A, B, XA, XB)  \
        do {                                            \
            /* reverse the order of input B */          \
            REVERSE_32(B);                             \
            REVERSE_32(XB);                             \
                                                        \
            /* Level-1 comparisons */                   \
            __m256i l1 = _mm256_min_epi32(A, B);        \
            __m256i h1 = _mm256_max_epi32(A, B);        \
            											\
			__m256i mask1 = _mm256_cmpeq_epi32(A, l1); \
            __m256i l1_oid = reinterpret_cast<__m256i>(_mm256_blendv_ps(reinterpret_cast<__m256>(XB), reinterpret_cast<__m256>(XA), _mm256_castsi256_ps(mask1))); \
            __m256i h1_oid = reinterpret_cast<__m256i>(_mm256_blendv_ps(reinterpret_cast<__m256>(XA), reinterpret_cast<__m256>(XB), _mm256_castsi256_ps(mask1))); \
            																				\
			/* Level-1 shuffle */                   \
			__m256i l1_sf = _mm256_permute2f128_si256(l1, h1, 0x20);					\
			__m256i h1_sf = _mm256_permute2f128_si256(l1, h1, 0x31);					\
															\
			__m256i l1_oid_sf = _mm256_permute2f128_si256(l1_oid, h1_oid, 0x20);					\
			__m256i h1_oid_sf = _mm256_permute2f128_si256(l1_oid, h1_oid, 0x31);					\
														\
			/* Level-2 comparisons */                  \
        	__m256i l2 = _mm256_min_epi32(l1_sf, h1_sf);        \
        	__m256i h2 = _mm256_max_epi32(l1_sf, h1_sf);        \
        														\
			__m256i mask2 = _mm256_cmpeq_epi32(l1_sf, l2); \
			__m256i l2_oid = reinterpret_cast<__m256i>(_mm256_blendv_ps(reinterpret_cast<__m256>(h1_oid_sf), reinterpret_cast<__m256>(l1_oid_sf), _mm256_castsi256_ps(mask2))); \
			__m256i h2_oid = reinterpret_cast<__m256i>(_mm256_blendv_ps(reinterpret_cast<__m256>(l1_oid_sf), reinterpret_cast<__m256>(h1_oid_sf), _mm256_castsi256_ps(mask2)));        \
													\
			/* Level-2 shuffle */                   \
			__m256i l2_sf = _mm256_unpacklo_epi64(l2, h2);					\
			__m256i h2_sf = _mm256_unpackhi_epi64(l2, h2);					\
																			\
			__m256i l2_oid_sf = _mm256_unpacklo_epi64(l2_oid, h2_oid);					\
			__m256i h2_oid_sf = _mm256_unpackhi_epi64(l2_oid, h2_oid);					\
																		\
			/* Level-3 comparisons */                  \
			__m256i l3 = _mm256_min_epi32(l2_sf, h2_sf);        \
			__m256i h3 = _mm256_max_epi32(l2_sf, h2_sf);        \
				\
			__m256i mask3 = _mm256_cmpeq_epi32(l2_sf, l3); \
			__m256i l3_oid = reinterpret_cast<__m256i>(_mm256_blendv_ps(reinterpret_cast<__m256>(h2_oid_sf), reinterpret_cast<__m256>(l2_oid_sf), _mm256_castsi256_ps(mask3))); \
			__m256i h3_oid = reinterpret_cast<__m256i>(_mm256_blendv_ps(reinterpret_cast<__m256>(l2_oid_sf), reinterpret_cast<__m256>(h2_oid_sf), _mm256_castsi256_ps(mask3)));        \
				\
			/* Level-3 shuffle */                   \
			l3 = _mm256_shuffle_epi32(l3, 0xd8);					\
			h3 = _mm256_shuffle_epi32(h3, 0xd8);					\
			__m256i l3_sf = _mm256_unpacklo_epi32(l3, h3);					\
			__m256i h3_sf = _mm256_unpackhi_epi32(l3, h3);					\
			\
			l3_oid = _mm256_shuffle_epi32(l3_oid, 0xd8);					\
			h3_oid = _mm256_shuffle_epi32(h3_oid, 0xd8);					\
			__m256i l3_oid_sf = _mm256_unpacklo_epi32(l3_oid, h3_oid);					\
			__m256i h3_oid_sf = _mm256_unpackhi_epi32(l3_oid, h3_oid);					\
			\
			/* Level-4 comparisons */                  \
			__m256i l4 = _mm256_min_epi32(l3_sf, h3_sf);        \
			__m256i h4 = _mm256_max_epi32(l3_sf, h3_sf);        \
			\
			__m256i mask4 = _mm256_cmpeq_epi32(l3_sf, l4); \
			__m256i l4_oid = reinterpret_cast<__m256i>(_mm256_blendv_ps(reinterpret_cast<__m256>(h3_oid_sf), reinterpret_cast<__m256>(l3_oid_sf), _mm256_castsi256_ps(mask4))); \
			__m256i h4_oid = reinterpret_cast<__m256i>(_mm256_blendv_ps(reinterpret_cast<__m256>(l3_oid_sf), reinterpret_cast<__m256>(h3_oid_sf), _mm256_castsi256_ps(mask4)));        \
				\
			/* Level-4 shuffle */                  \
			__m256i l4_sf = _mm256_permute4x64_epi64(l4, 0xd8);					\
			__m256i h4_sf = _mm256_permute4x64_epi64(h4, 0xd8);					\
			O1 = _mm256_unpacklo_epi32(l4_sf, h4_sf);	\
			O2 = _mm256_unpackhi_epi32(l4_sf, h4_sf);	\
			                  \
			__m256i l4_oid_sf = _mm256_permute4x64_epi64(l4_oid, 0xd8);					\
			__m256i h4_oid_sf = _mm256_permute4x64_epi64(h4_oid, 0xd8);					\
			XO1 = _mm256_unpacklo_epi32(l4_oid_sf, h4_oid_sf);	\
			XO2 = _mm256_unpackhi_epi32(l4_oid_sf, h4_oid_sf);	\
        } while(0)

#if 0	//original version with _m256_min_pd
/** Bitonic merge kernel for 2 x 8 elements int64_t, X-- refer to those for oids */
#define BITONIC_MERGE8_64(O1, O2, O3, O4, XO1, XO2, XO3, XO4, A1, A2, B1, B2, XA1, XA2, XB1, XB2)  \
        do {                                            \
            /* reverse the order of input B */          \
            REVERSE_64(B1);                             \
            REVERSE_64(B2);                             \
            REVERSE_64(XB1);                             \
            REVERSE_64(XB2);                             \
                                                        \
            /* Level-0 comparisons */                   \
            __m256d l11_val = _mm256_min_pd(A1, B2);        \
            __m256d h11_val = _mm256_max_pd(A1, B2);        \
            __m256d l12_val = _mm256_min_pd(A2, B1);        \
            __m256d h12_val = _mm256_max_pd(A2, B1);        \
            											\
			__m256i mask1 = _mm256_cmpeq_epi64(_mm256_castpd_si256(A1), _mm256_castpd_si256(l11_val)); \
            __m256d l11_oid = _mm256_blendv_pd(XB2, XA1, _mm256_castsi256_pd(mask1)); \
            __m256d h11_oid = _mm256_blendv_pd(XA1, XB2, _mm256_castsi256_pd(mask1));        \
            																				\
            __m256i mask2 = _mm256_cmpeq_epi64(_mm256_castpd_si256(A2), _mm256_castpd_si256(l12_val)); \
            __m256d l12_oid = _mm256_blendv_pd(XB1, XA2, _mm256_castsi256_pd(mask2));        \
            __m256d h12_oid = _mm256_blendv_pd(XA2, XB1, _mm256_castsi256_pd(mask2));        \
                                                        \
            BITONIC4_64(O1, O2, XO1, XO2, l11_val, l12_val, l11_oid, l12_oid);                 \
            BITONIC4_64(O3, O4, XO3, XO4, h11_val, h12_val, h11_oid, h12_oid);                 \
        } while(0)
#endif

#if 1	//version with _m256_cmpgt and _blendv
/** Bitonic merge kernel for 2 x 8 elements int64_t, X-- refer to those for oids */
#define BITONIC_MERGE8_64(O1, O2, O3, O4, XO1, XO2, XO3, XO4, A1, A2, B1, B2, XA1, XA2, XB1, XB2)  \
        do {                                            \
            /* reverse the order of input B */          \
            REVERSE_64(B1);                             \
            REVERSE_64(B2);                             \
            REVERSE_64(XB1);                             \
            REVERSE_64(XB2);                             \
                                                        \
            /* Level-0 comparisons */                   \
			__m256i mask1 = _mm256_cmpgt_epi64(_mm256_castpd_si256(A1), _mm256_castpd_si256(B2)); \
            __m256d l11_val = _mm256_blendv_pd(A1, B2, _mm256_castsi256_pd(mask1)); \
            __m256d h11_val = _mm256_blendv_pd(B2, A1, _mm256_castsi256_pd(mask1)); \
            __m256i mask2 = _mm256_cmpgt_epi64(_mm256_castpd_si256(A2), _mm256_castpd_si256(B1)); \
            __m256d l12_val = _mm256_blendv_pd(A2, B1, _mm256_castsi256_pd(mask2)); \
            __m256d h12_val = _mm256_blendv_pd(B1, A2, _mm256_castsi256_pd(mask2)); \
            											\
            __m256d l11_oid = _mm256_blendv_pd(XA1, XB2, _mm256_castsi256_pd(mask1)); \
            __m256d h11_oid = _mm256_blendv_pd(XB2, XA1, _mm256_castsi256_pd(mask1));        \
            																				\
            __m256d l12_oid = _mm256_blendv_pd(XA2, XB1, _mm256_castsi256_pd(mask2));        \
            __m256d h12_oid = _mm256_blendv_pd(XB1, XA2, _mm256_castsi256_pd(mask2));        \
                                                        \
            BITONIC4_64(O1, O2, XO1, XO2, l11_val, l12_val, l11_oid, l12_oid);                 \
            BITONIC4_64(O3, O4, XO3, XO4, h11_val, h12_val, h11_oid, h12_oid);                 \
        } while(0)
#endif

/** Bitonic merge kernel for 2 x 16 elements */
#define BITONIC_MERGE16(O1, O2, O3, O4, O5, O6, O7, O8,         \
                        A1, A2, A3, A4, B1, B2, B3, B4)         \
        do {                                                    \
            /** Bitonic merge kernel for 2 x 16 elemenets */    \
            /* reverse the order of input B */                  \
            REVERSE(B1);                                        \
            REVERSE(B2);                                        \
            REVERSE(B3);                                        \
            REVERSE(B4);                                        \
                                                                \
            /* Level-0 comparisons */                           \
            __m256d l01 = _mm256_min_pd(A1, B4);                \
            __m256d l02 = _mm256_min_pd(A2, B3);                \
            __m256d l03 = _mm256_min_pd(A3, B2);                \
            __m256d l04 = _mm256_min_pd(A4, B1);                \
            __m256d h01 = _mm256_max_pd(A1, B4);                \
            __m256d h02 = _mm256_max_pd(A2, B3);                \
            __m256d h03 = _mm256_max_pd(A3, B2);                \
            __m256d h04 = _mm256_max_pd(A4, B1);                \
                                                                \
            BITONIC8(O1, O2, O3, O4, l01, l02, l03, l04);       \
            BITONIC8(O5, O6, O7, O8, h01, h02, h03, h04);       \
        } while(0)


/** 
 * There are 3 ways to implement branches:
 *     1) With conditional move instr.s using inline assembly (IFELSEWITHCMOVE).
 *     2) With software predication (IFELSEWITHPREDICATION).
 *     3) With normal if-else
 */
#if IFELSEWITHCMOVE
#define IFELSECONDMOVE(NXT, INA, INB, INCR)                             \
    do {                                                                \
        register block4 * tmpA, * tmpB;                                 \
        register int64_t tmpKey;                                        \
                                                                        \
        __asm__ ( "mov %[A], %[tmpA]\n"         /* tmpA <-- inA      */ \
                  "add %[INC], %[A]\n"          /* inA += 4          */ \
                  "mov %[B], %[tmpB]\n"         /* tmpB <-- inB      */ \
                  "mov (%[tmpA]), %[tmpKey]\n"  /* tmpKey <-- *inA   */ \
                  "add %[INC], %[B]\n"          /* inB += 4          */ \
                  "mov %[tmpA], %[NEXT]\n"      /* next <-- A        */ \
                  "cmp (%[tmpB]), %[tmpKey]\n"  /* cmp(tmpKey,*inB ) */ \
                  "cmovnc %[tmpB], %[NEXT]\n"   /* if(A>=B) next<--B */ \
                  "cmovnc %[tmpA], %[A]\n"      /* if(A>=B) A<--oldA */ \
                  "cmovc %[tmpB], %[B]\n"       /* if(A<B)  B<--oldB */ \
                  : [A] "=r" (INA), [B] "=r" (INB), [NEXT] "=r" (NXT),  \
                    [tmpA] "=r" (tmpA), [tmpB] "=r" (tmpB),             \
                    [tmpKey] "=r" (tmpKey)                              \
                  : "0" (INA), "1" (INB), [INC] "i" (INCR)              \
                  :                                                     \
                  );                                                    \
    } while(0)

#elif IFELSEWITHPREDICATION
#define IFELSECONDMOVE(NXT, NXT_OID, INA, INB, INA_OID, INB_OID, T)                 \
    do {                                                    \
        int8_t cmp = *((T *)INA) < *((T *)INB); \
        NXT  = cmp ? INA : INB;                             \
        NXT_OID = cmp ? INA_OID : INB_OID;					\
        INA += cmp;                                         \
        INB += !cmp;                                        \
        INA_OID += cmp;										\
        INB_OID += !cmp;									\
    } while(0)

#elif IFELSEWITHNORMAL
#define IFELSECONDMOVE(NXT, INA, INB, INCR)                 \
            do {                                            \
                if(*((int64_t *)INA) < *((int64_t *)INB)) { \
                    NXT = INA;                              \
                    INA ++;                                 \
                }                                           \
                else {                                      \
                    NXT = INB;                              \
                    INB ++;                                 \
                }                                           \
            } while(0)                                      \

#endif

}

#endif /* AVXCOMMON_H */
