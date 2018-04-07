#ifndef MERGESORT_H
#define MERGESORT_H

//#include "sorter.h"
#include "types.h"
#include "barrier.h"
#include "common.h"
#include <cmath>
#include "avxcommon.h"
#include "params.h"
#include <cstring>
#include <algorithm>
#include <pthread.h>
#include <iostream>
#include <vector>

namespace multiAttrSort {

#if 0
struct merge_node_t {
	element_t *buffer;
	volatile uint32_t count;
	volatile uint32_t head;
	volatile uint32_t tail;
} __attribute__((packed));
#endif
	
//class Mergesort: public Sorter {
class Mergesort {
public:
	Mergesort() {}

	//virtual ~Mergesort(){}

	//void do_sort_pack(BAT_pack_t* column, setting_t* setting) override;
	template <class T>
	void do_sort(BAT_t<T>* column, setting_t* setting);

//private:
	template <class T>
	static void* mergesort_thread(void* param);

	/********************************************************** 
	 *                phase 1: partitioning phase             *
	 **********************************************************/
	template <class T>
	static void partitioning_phase(BAT_t<T>*** partitions, mergesort_params_t<T>* setting);
	template <class T>
	static void radix_partitioning(BAT_t<T>* restrict outCol, BAT_t<T>* restrict inCol,
			uint64_t* restrict hist, uint32_t bitshift, uint32_t nradixbits);

	/********************************************************** 
	 *                phase 2: sorting phase                  *
	 **********************************************************/
	template <class T>
	static void sorting_phase(BAT_t<T>** partitions, mergesort_params_t<T>* setting);

	template <class T>
	inline static void __attribute__((always_inline))
	flip(T** inptr_val, uint64_t nitems, uint32_t bitwidth);

	template <class T>
	static void avxmergesort_anytype(T** inputptr_val, surrogate_t** inputptr_oid,
			T** outputptr_val, surrogate_t** outputptr_oid, uint64_t nitems);
	
	/**
	 * Important: the memory of inputptr and outputptr may be swapped
	 * @param block_size MUST be power of two, equal or larger than 128 (hard-coded)
	 */
#if 0
	template <class T>
	inline static void __attribute__((always_inline))
	avxmergesort_block_aligned(T** inputptr_val, surrogate_t** inputptr_oid,
			T** outputptr_val, surrogate_t** outputptr_oid, uint32_t block_size);
#endif
	inline static void __attribute__((always_inline))
	avxmergesort_block_uint64_aligned(uint64_t** inputptr_val, surrogate_t** inputptr_oid,
			uint64_t** outputptr_val, surrogate_t** outputptr_oid, uint32_t block_size);

	inline static void __attribute__((always_inline))
	avxmergesort_block_uint32_aligned(uint32_t** inputptr_val, surrogate_t** inputptr_oid,
			uint32_t** outputptr_val, surrogate_t** outputptr_oid, uint32_t block_size);

	inline static void __attribute__((always_inline))
	avxmergesort_block_uint16_aligned(uint16_t** inputptr_val, surrogate_t** inputptr_oid,
			uint16_t** outputptr_val, surrogate_t** outputptr_oid, uint32_t block_size);


	//this function cannot be made inline because it is called *outside* the class?
	inline static void __attribute__((always_inline))
	//static void
	avxmergesort_rem_uint16_aligned(uint16_t** inputptr_val, surrogate_t** inputptr_oid,
			uint16_t** outputptr_val, surrogate_t** outputptr_oid, uint32_t nitems);

	inline static void __attribute__((always_inline))
	//static void
	avxmergesort_rem_uint32_aligned(uint32_t** inputptr_val, surrogate_t** inputptr_oid,
			uint32_t** outputptr_val, surrogate_t** outputptr_oid, uint32_t nitems);

	inline static void __attribute__((always_inline))
	//static void
	avxmergesort_rem_uint64_aligned(uint64_t** inputptr_val, surrogate_t** inputptr_oid,
			uint64_t** outputptr_val, surrogate_t** outputptr_oid, uint32_t nitems);

	inline static void __attribute__((always_inline))
	inregister_sort_int64_aligned(int64_t * in_val, int64_t * out_val,
			int32_t * in_oid, int32_t * out_oid);

	inline static void __attribute__((always_inline))
	inregister_sort_int32_aligned(int32_t * in_val, int32_t * out_val,
			int32_t * in_oid, int32_t * out_oid);

	inline static void __attribute__((always_inline))
	inregister_sort_int16_aligned(int16_t * in_val, int16_t * out_val,
			int32_t * in_oid, int32_t * out_oid);

	inline static const int __attribute__((always_inline))
	movemask_epi16_lo(const __m256i &a);

	inline static const int __attribute__((always_inline))
	movemask_epi16_hi(const __m256i &a);
	/*************************************************************
	 *                phase 3: multi-way merge phase             *
	 *************************************************************/
#if 0
	template <class T>
	static void multiwaymerge_phase(mergesort_params_t<T>* setting);

	static int64_t avx_multiway_merge(element_t *output,
							BAT_pack_t **parts,
							uint32_t nparts,
							element_t *fifobuffer,
							uint32_t bufnelements);

	/**
	  * Read from <inA> and <inB> to the fifo buffer managed by <node>
	  *	 
	  *
	  */
	static uint32_t readmerge_parallel_decomposed(merge_node_t *node,
									element_t **inA,
									element_t **inB,
									uint32_t lenA,
									uint32_t lenB,
									uint32_t fifosize);

	inline static void __attribute__((always_inline))
	parallel_read(element_t **A, element_t **B, element_t **Out,
			uint32_t *ri, uint32_t *li, uint32_t *oi, uint32_t *outnslots,
			uint32_t lenA, uint32_t lenB);
	
	static void direct_copy_avx(merge_node_t *dest, merge_node_t *src, uint32_t fifosize);
	static uint32_t direct_copy_to_output_avx(element_t *dest, merge_node_t *src, uint32_t fifosize);
	static void simd_memcpy(void *dst, void *src, size_t sz);

	static void merge_parallel_decomposed(merge_node_t * node,
	                          merge_node_t * right,
	                          merge_node_t * left,
	                          uint32_t fifosize,
	                          uint8_t rightdone, uint8_t leftdone);

	static uint64_t mergestore_parallel_decomposed(merge_node_t *right,
									merge_node_t *left,
									element_t **output,
									uint32_t fifosize,
									uint8_t rightnode,
									uint8_t leftdone);
#endif

	/***********************************************************
	 *                    sort/merge kernels                   *
	 ***********************************************************/
#if 1
	inline static void __attribute__((always_inline))
	merge8_int64_varlen_aligned(int64_t * restrict inpA_val,
	                      			int64_t * restrict inpB_val,
	                      			int64_t * restrict Out_val,
									surrogate_t * restrict inpA_oid,
									surrogate_t * restrict inpB_oid,
									surrogate_t * restrict Out_oid,
	                      			const uint32_t lenA,
	                      			const uint32_t lenB);

	inline static void __attribute__((always_inline))
	merge8_int32_varlen_aligned(int32_t * restrict inpA_val,
	                      			int32_t * restrict inpB_val,
	                      			int32_t * restrict Out_val,
									surrogate_t * restrict inpA_oid,
									surrogate_t * restrict inpB_oid,
									surrogate_t * restrict Out_oid,
	                      			const uint32_t lenA,
	                      			const uint32_t lenB);

	inline static void __attribute__((always_inline))
	merge8_int32_varlen_unaligned(int32_t * restrict inpA_val,
	                      			int32_t * restrict inpB_val,
	                      			int32_t * restrict Out_val,
									surrogate_t * restrict inpA_oid,
									surrogate_t * restrict inpB_oid,
									surrogate_t * restrict Out_oid,
	                      			const uint32_t lenA,
	                      			const uint32_t lenB);

	inline static void __attribute__((always_inline))
	merge16_int16_varlen_aligned(int16_t * restrict inpA_val,
	                      			int16_t * restrict inpB_val,
	                      			int16_t * restrict Out_val,
									surrogate_t * restrict inpA_oid,
									surrogate_t * restrict inpB_oid,
									surrogate_t * restrict Out_oid,
	                      			const uint32_t lenA,
	                      			const uint32_t lenB);
#endif
	
	inline static void __attribute__((always_inline))
	merge16_varlen_aligned(int64_t * restrict inpA, 
                       int64_t * restrict inpB,
                       int64_t * restrict Out, 
                       const uint32_t lenA, 
                       const uint32_t lenB);
#if 0
	inline static void __attribute__((always_inline))
	inregister_sort_int64_aligned(int64_t * in_val, int64_t * out_val,
			int32_t * in_oid, int32_t * out_oid);
#endif

	inline static void __attribute__((always_inline))
	merge4_int64_eqlen_aligned(int64_t * const inpA, int64_t * const inpB,
							surrogate_t * const inpA_oid, surrogate_t * const inpB_oid,
							int64_t * const out, surrogate_t * const out_oid, const uint32_t len);
	
	inline static void __attribute__((always_inline))
	merge8_int32_eqlen_aligned(int32_t * const inpA, int32_t * const inpB,
							surrogate_t * const inpA_oid, surrogate_t * const inpB_oid,
							int32_t * const out, surrogate_t * const out_oid, const uint32_t len);

	inline static void __attribute__((always_inline))
	merge16_int16_eqlen_aligned(int16_t * const inpA, int16_t * const inpB,
			surrogate_t * const inpA_oid, surrogate_t * const inpB_oid,
			int16_t * const out, surrogate_t * const out_oid, const uint32_t len);

	/**
	 * key issue: what is the different with merge16_varlen_aligned?
	 *
	 * 
	 */
#if 0
	inline static void __attribute__((always_inline))
	merge16kernel(element_t * restrict A, element_t * restrict B, element_t * restrict Out,
			uint32_t *ri, uint32_t *li, uint32_t *oi, uint32_t *outnslots,
			uint32_t rend, uint32_t lend);	
	
	inline static void __attribute__((always_inline))
	mergestore16kernel(element_t * restrict A, element_t * restrict B,
				element_t **Out, uint32_t * ri, uint32_t * li,
				uint32_t rend, uint32_t lend);

#endif
};

#if 0
template <class T>
void Mergesort::do_sort_single_threaded(BAT_t<T>* column, setting_t* setting) {
	/* directly call avxmergesort_anytype to do the sorting
	 * NOTE: do the value flipping if necessary
	 */
	uint32_t bitwidth = setting->bitwidth;
	uint64_t nitem = column->num_elements;
	//flip the values if the bitwidth is 8/16/32/64 because SIMD-bank uses signed integer
	if ((bitwidth % 8) == 0) {
		flip<T>(&(column->values), nitem, bitwidth);
	}

	T* outptr_val = (T *) malloc_aligned(sizeof(T) * nitem);
	surrogate_t *outptr_oid = (surrogate_t *) malloc_aligned(sizeof(surrogate_t) * nitem);

	switch(setting->intrinType) {
		case IntrinsicsType::AVX:
			avxmergesort_anytype<T>(&(column->values), &(column->oids), &outptr_val, &outptr_oid, nitem);
			break;
		case IntrinsicsType::SSE:
			break;
		case IntrinsicsType::SCALAR:
			break;
	}

}
#endif

template <class T>
void Mergesort::do_sort(BAT_t<T>* column, setting_t* setting) {
	uint32_t nthreads = setting->nthreads;
	uint32_t part_fanout = setting->partition_fanout;
	int rv;
	uint32_t i;
	pthread_barrier_t barrier;

	/** # of elements assigned for each thread **/
	uint32_t numperthr = column->num_elements / nthreads;

	pthread_t* thread_handles = (pthread_t*) malloc_aligned(nthreads*sizeof(pthread_t));
	mergesort_params_t<T>* params = (mergesort_params_t<T> *) malloc_aligned(nthreads*sizeof(mergesort_params_t<T>));

	/** allocate temporary space for partitioning: revised-->both for <values> and <oids>**/
	T* tmpColPartition = NULL;
	tmpColPartition = (T *) malloc_aligned(column->num_elements*sizeof(T) +
			column_padding<T>(nthreads, part_fanout));
	surrogate_t* tmpOIDPartition = NULL;
	tmpOIDPartition = (surrogate_t *) malloc_aligned(column->num_elements*sizeof(surrogate_t) +
			column_padding<surrogate_t>(nthreads, part_fanout));

	/** allocate temporary space for sorting: revised-->both for <values> and <oids> **/
	T* tmpColSort = NULL;
	tmpColSort = (T *) malloc_aligned(column->num_elements*sizeof(T) +
			column_padding<T>(nthreads, part_fanout));
	surrogate_t* tmpOIDSort = NULL;
	tmpOIDSort = (surrogate_t *) malloc_aligned(column->num_elements*sizeof(surrogate_t) +
			column_padding<surrogate_t>(nthreads, part_fanout));

#if 0
	/** allocate histograms arrays, actual allocation is local to threads **/
	uint32_t** histogram = (uint32_t**) malloc_aligned(nthreads*sizeof(uint32_t*));
#endif

	/** initialize the barrier for threads **/
	rv = pthread_barrier_init(&barrier, NULL, nthreads);
	if (rv != 0) {
		printf("[ERROR ] Couldn't create the barrier\n");
		exit(EXIT_FAILURE);
	}

	/** initialize threalcolchunk,
	  * used to pass the per-thread sorted partitions to multi-way merge phase**/
	BAT_t<T>** threadcolchunks = (BAT_t<T>**) malloc_aligned(nthreads * sizeof(BAT_t<T>*));

	T* ptrs_sharedmergebuffer = (T *) malloc_aligned(L3_CACHE_SIZE);//TODO: can make it as runtime parameter

	BAT_t<T>* result = (BAT_t<T> *) malloc_aligned(nthreads * sizeof(BAT_t<T>));

	/** prepare and fork all threads **/
	for (i = 0; i < nthreads; ++i) {

		/** assign values for each thread **/
		params[i].nthreads 			= setting->nthreads;
		params[i].intrinType 		= setting->intrinType;
		params[i].bitwidth 			= setting->bitwidth;
		params[i].partition_fanout 	= setting->partition_fanout;
		params[i].my_tid 			= i;
		params[i].barrier 			= &barrier;

		params[i].columnchunk 		= column->values + i * numperthr;
		params[i].oidchunk			= column->oids + i * numperthr;
		params[i].tmp_partition 	= tmpColPartition + i * (numperthr +
				cache_line_padding<T>(part_fanout));
		params[i].tmp_partition_oid = tmpOIDPartition + i * (numperthr +
				cache_line_padding<surrogate_t>(part_fanout));
		params[i].tmp_sort 			= tmpColSort + i * (numperthr +
				cache_line_padding<T>(part_fanout));
		params[i].tmp_sort_oid		= tmpOIDSort + i * (numperthr +
				cache_line_padding<surrogate_t>(part_fanout));

		params[i].chunksize = (i == (nthreads-1)) ?
			(column->num_elements - i * numperthr) : numperthr;

		params[i].threadcolchunks 	= threadcolchunks;
		params[i].sharedmergebuffer = ptrs_sharedmergebuffer;
		params[i].result			= result;


		/** run the merge sort thread  **/
		rv = pthread_create(&thread_handles[i], NULL, mergesort_thread<T>, (void*)&params[i]);

		if (rv) {
			printf("[ERROR ] return code from pthread_create() is %d\n", rv);
			exit(EXIT_FAILURE);
		}
	}

	/* wait for threads to finish*/
	for (i = 0; i < nthreads; ++i) {
		pthread_join(thread_handles[i], NULL);
	}


	/** extract the sorted columns from <threadcolchunks> or <result>,
	 *  depend on single-thread or multi-thread */
	uint64_t rowIdx = 0;
	uint32_t fid = 0;
	uint64_t part_nelements;
	uint32_t part_rowIdx;
	uint64_t thread_nelements;
	uint32_t thread_rowIdx;
	if (nthreads == 1) {
		/* get the sorted elements from threadcolchunks[0][i] where <i> is the fanout */
		for (fid = 0; fid < part_fanout; fid++) {
			part_nelements = threadcolchunks[0][fid].num_elements;
			for (part_rowIdx = 0; part_rowIdx < part_nelements; ++part_rowIdx) {
				column->values[rowIdx] = threadcolchunks[0][fid].values[part_rowIdx];
				column->oids[rowIdx] = threadcolchunks[0][fid].oids[part_rowIdx];
				//assert(column->oids[rowIdx] < 1638400);	//not a correct insert when we are sorting a subrange
				++rowIdx;
			}
		}
		assert(rowIdx == column->num_elements);	//# of elements should not be changed
	} else {
		/* get the sorted elements from <result> of each thread */
		for (i = 0; i < nthreads; ++i) {
			thread_nelements = result[i].num_elements;
			for (thread_rowIdx = 0; thread_rowIdx < thread_nelements; ++thread_rowIdx) {
				column->values[rowIdx] = result[i].values[thread_rowIdx];
				column->oids[rowIdx] = result[i].oids[thread_rowIdx];
				++rowIdx;
			}
		}
		assert(rowIdx == column->num_elements);	//# of elements should not be changed
	}

	/** destroy **/
	for (i = 0; i < nthreads; ++i) {
		free(threadcolchunks[i]);
	}
	if (nthreads > 1) {
		for (i = 0; i < nthreads; ++i) {
			free(result[i].values);
			free(result[i].oids);
		}
	}

	//free(histogram);
	free(ptrs_sharedmergebuffer);
	free(threadcolchunks);
	free(tmpColSort);
	free(tmpOIDSort);
	free(tmpColPartition);
	free(tmpOIDPartition);
	free(params);
	free(thread_handles);
	free(result);

	pthread_barrier_destroy(&barrier);
}

#if 0
inline void __attribute((always_inline))
Mergesort::inregister_sort_int64_aligned(int64_t * in_val, int64_t * out_val,
		int32_t * in_oid, int32_t * out_oid)
{
/* IACA_START */
    __m256d ra = _mm256_load_si256 (reinterpret_cast<__m256i const *>(in_val));
    __m256d rb = _mm256_load_si256 (reinterpret_cast<__m256i const *>(in_val + 4));
    __m256d rc = _mm256_load_si256 (reinterpret_cast<__m256i const *>(in_val + 8));
    __m256d rd = _mm256_load_si256 (reinterpret_cast<__m256i const *>(in_val + 12));

    __m256d oa = _mm256_cvtps_pd(_mm_load_ps((float const *)(in_oid)));
    __m256d ob = _mm256_cvtps_pd(_mm_load_ps((float const *)(in_oid + 4)));
    __m256d oc = _mm256_cvtps_pd(_mm_load_ps((float const *)(in_oid + 8)));
    __m256d od = _mm256_cvtps_pd(_mm_load_ps((float const *)(in_oid + 12)));

    /* odd-even sorting network begins */
    /* 1st level of comparisons */
    __m256d ra1 = _mm256_min_pd(ra, rb);
    __m256d rb1 = _mm256_max_pd(ra, rb);

    __m256d rc1 = _mm256_min_pd(rc, rd);
    __m256d rd1 = _mm256_max_pd(rc, rd);

    //re-ordering the oids
    __m256i mask = _mm256_cmpeq_epi64(_mm256_castpd_si256(ra), _mm256_castpd_si256(ra1));
    __m256d oa1 = _mm256_blendv_pd(ob, oa, _mm256_castsi256_pd(mask));
    __m256d ob1 = _mm256_blendv_pd(oa, ob, _mm256_castsi256_pd(mask));

    mask = _mm256_cmpeq_epi64(_mm256_castpd_si256(rc), _mm256_castpd_si256(rc1));
    __m256d oc1 = _mm256_blendv_pd(od, oc, _mm256_castsi256_pd(mask));
    __m256d od1 = _mm256_blendv_pd(oc, od, _mm256_castsi256_pd(mask));

    /* 2nd level of comparisons */
    rb = _mm256_min_pd(rb1, rd1);
    rd = _mm256_max_pd(rb1, rd1);

    mask = _mm256_cmpeq_epi64(_mm256_castpd_si256(rb1), _mm256_castpd_si256(rb));
    ob = _mm256_blendv_pd(od1, ob1, _mm256_castsi256_pd(mask));
    od = _mm256_blendv_pd(ob1, od1, _mm256_castsi256_pd(mask));

    /* 3rd level of comparisons */
    __m256d ra2 = _mm256_min_pd(ra1, rc1);
    __m256d rc2 = _mm256_max_pd(ra1, rc1);

    mask = _mm256_cmpeq_epi64(_mm256_castpd_si256(ra1), _mm256_castpd_si256(ra2));
    __m256d oa2 = _mm256_blendv_pd(oc1, oa1, _mm256_castsi256_pd(mask));
    __m256d oc2 = _mm256_blendv_pd(oa1, oc1, _mm256_castsi256_pd(mask));

    /* 4th level of comparisons */
    __m256d rb3 = _mm256_min_pd(rb, rc2);
    __m256d rc3 = _mm256_max_pd(rb, rc2);

    mask = _mm256_cmpeq_epi64(_mm256_castpd_si256(rb), _mm256_castpd_si256(rb3));
    __m256d ob3 = _mm256_blendv_pd(oc2, ob, _mm256_castsi256_pd(mask));
    __m256d oc3 = _mm256_blendv_pd(ob, oc2, _mm256_castsi256_pd(mask));

    /* results are in ra2, rb3, rc3, rd */
    /* re-ordered oids are in oa2, ob3, oc3, od */
    /**
     * Initial data and transposed data looks like following:
     *  a2={ x1  x2  x3  x4  }                      a4={ x1 x5 x9  x13 }
     *  b3={ x5  x6  x7  x8  }  === Transpose ===>  b5={ x2 x6 x10 x14 }
     *  c3={ x9  x10 x11 x12 }                      c5={ x3 x7 x11 x15 }
     *  d={ x13 x14 x15 x16 }                       d4={ x4 x8 x12 x16 }
     */
    /* shuffle x2 and x5 - shuffle x4 and x7 */
    __m256d ra3 = _mm256_unpacklo_pd(ra2, rb3);
    __m256d rb4 = _mm256_unpackhi_pd(ra2, rb3);

    __m256d oa3 = _mm256_unpacklo_pd(oa2, ob3);
    __m256d ob4 = _mm256_unpackhi_pd(oa2, ob3);

    /* shuffle x10 and x13 - shuffle x12 and x15 */
    __m256d rc4 = _mm256_unpacklo_pd(rc3, rd);
    __m256d rd3 = _mm256_unpackhi_pd(rc3, rd);

    __m256d oc4 = _mm256_unpacklo_pd(oc3, od);
    __m256d od3 = _mm256_unpackhi_pd(oc3, od);

    /* shuffle (x3,x7) and (x9,x13) pairs */
    __m256d ra4 = _mm256_permute2f128_pd(ra3, rc4, 0x20);
    __m256d rc5 = _mm256_permute2f128_pd(ra3, rc4, 0x31);

    __m256d oa4 = _mm256_permute2f128_pd(oa3, oc4, 0x20);
    __m256d oc5 = _mm256_permute2f128_pd(oa3, oc4, 0x31);

    /* shuffle (x4,x8) and (x10,x14) pairs */
    __m256d rb5 = _mm256_permute2f128_pd(rb4, rd3, 0x20);
    __m256d rd4 = _mm256_permute2f128_pd(rb4, rd3, 0x31);

    __m256d ob5 = _mm256_permute2f128_pd(ob4, od3, 0x20);
    __m256d od4 = _mm256_permute2f128_pd(ob4, od3, 0x31);

    /* after this, results are in ra4, rb5, rc5, rd4 */
    /* oids are in oa4, ob5, oc5, od4 */
/* IACA_END */
    /* store values */
    _mm256_store_pd((double *) out_val, ra4);
    _mm256_store_pd((double *) (out_val + 4), rb5);
    _mm256_store_pd((double *) (out_val + 8), rc5);
    _mm256_store_pd((double *) (out_val + 12), rd4);

    //convert the oids from 64 to 32
    __m128 oa4_s = _mm256_cvtpd_ps(oa4);
    __m128 ob5_s = _mm256_cvtpd_ps(ob5);
    __m128 oc5_s = _mm256_cvtpd_ps(oc5);
    __m128 od4_s = _mm256_cvtpd_ps(od4);

    __m128i mask_store = _mm_set1_epi32(1);
    _mm_maskstore_ps((float *) out_oid, mask_store, oa4_s);
    _mm_maskstore_ps((float *) (out_oid + 4), mask_store, ob5_s);
    _mm_maskstore_ps((float *) (out_oid + 8), mask_store, oc5_s);
    _mm_maskstore_ps((float *) (out_oid + 12), mask_store, od4_s);
}
#endif

inline void __attribute__((always_inline))
Mergesort::inregister_sort_int32_aligned(int32_t * in_val, int32_t * out_val,
		int32_t * in_oid, int32_t * out_oid)
{
	/* IACA_START */
	    __m256i ra = _mm256_load_si256 (reinterpret_cast<__m256i const *>(in_val));
	    __m256i rb = _mm256_load_si256 (reinterpret_cast<__m256i const *>(in_val + 8));
	    __m256i rc = _mm256_load_si256 (reinterpret_cast<__m256i const *>(in_val + 16));
	    __m256i rd = _mm256_load_si256 (reinterpret_cast<__m256i const *>(in_val + 24));
	    __m256i re = _mm256_load_si256 (reinterpret_cast<__m256i const *>(in_val + 32));
	    __m256i rf = _mm256_load_si256 (reinterpret_cast<__m256i const *>(in_val + 40));
	    __m256i rg = _mm256_load_si256 (reinterpret_cast<__m256i const *>(in_val + 48));
	    __m256i rh = _mm256_load_si256 (reinterpret_cast<__m256i const *>(in_val + 56));

	    __m256 oa = _mm256_load_ps (reinterpret_cast<float const *>(in_oid));
	    __m256 ob = _mm256_load_ps (reinterpret_cast<float const *>(in_oid + 8));
	    __m256 oc = _mm256_load_ps (reinterpret_cast<float const *>(in_oid + 16));
	    __m256 od = _mm256_load_ps (reinterpret_cast<float const *>(in_oid + 24));
	    __m256 oe = _mm256_load_ps (reinterpret_cast<float const *>(in_oid + 32));
	    __m256 of = _mm256_load_ps (reinterpret_cast<float const *>(in_oid + 40));
	    __m256 og = _mm256_load_ps (reinterpret_cast<float const *>(in_oid + 48));
	    __m256 oh = _mm256_load_ps (reinterpret_cast<float const *>(in_oid + 56));

	    /* odd-even sorting network begins: using Hibbard's Algorithm.
	     * [[0,1],[2,3],[4,5],[6,7]]
		 *	[[0,2],[1,3],[4,6],[5,7]]
		 *	[[1,2],[5,6],[0,4],[3,7]]
		 *	[[1,5],[2,6]]
		 *	[[1,4],[3,6]]
		 *	[[2,4],[3,5]]
		 *	[[3,4]]
		 *	grouped into 7 parallel operations.
		 *	important: may be optimized, e.g., re-use the variables as in in_register_sort_int64()
	     */
	    /* 1st level of comparisons: [[0,1],[2,3],[4,5],[6,7]]*/
	    __m256i ra1 = _mm256_min_epi32(ra, rb);
	    __m256i rb1 = _mm256_max_epi32(ra, rb);

	    __m256i rc1 = _mm256_min_epi32(rc, rd);
	    __m256i rd1 = _mm256_max_epi32(rc, rd);

	    __m256i re1 = _mm256_min_epi32(re, rf);
	    __m256i rf1 = _mm256_max_epi32(re, rf);

	    __m256i rg1 = _mm256_min_epi32(rg, rh);
	    __m256i rh1 = _mm256_max_epi32(rg, rh);

	    //re-ordering the oids
	    __m256i mask1 = _mm256_cmpeq_epi32(ra, ra1);
	    __m256 oa1 = _mm256_blendv_ps(ob, oa, _mm256_castsi256_ps(mask1));
	    __m256 ob1 = _mm256_blendv_ps(oa, ob, _mm256_castsi256_ps(mask1));

	    __m256i mask2 = _mm256_cmpeq_epi32(rc, rc1);
	    __m256 oc1 = _mm256_blendv_ps(od, oc, _mm256_castsi256_ps(mask2));
	    __m256 od1 = _mm256_blendv_ps(oc, od, _mm256_castsi256_ps(mask2));

	    __m256i mask3 = _mm256_cmpeq_epi32(re, re1);
	    __m256 oe1 = _mm256_blendv_ps(of, oe, _mm256_castsi256_ps(mask3));
	    __m256 of1 = _mm256_blendv_ps(oe, of, _mm256_castsi256_ps(mask3));

	    __m256i mask4 = _mm256_cmpeq_epi32(rg, rg1);
	    __m256 og1 = _mm256_blendv_ps(oh, og, _mm256_castsi256_ps(mask4));
	    __m256 oh1 = _mm256_blendv_ps(og, oh, _mm256_castsi256_ps(mask4));

	    /* 2nd level of comparisons: [[0,2],[1,3],[4,6],[5,7]] */
	    __m256i ra2 = _mm256_min_epi32(ra1, rc1);
	    __m256i rc2 = _mm256_max_epi32(ra1, rc1);

	    __m256i rb2 = _mm256_min_epi32(rb1, rd1);
	    __m256i rd2 = _mm256_max_epi32(rb1, rd1);

	    __m256i re2 = _mm256_min_epi32(re1, rg1);
	    __m256i rg2 = _mm256_max_epi32(re1, rg1);

	    __m256i rf2 = _mm256_min_epi32(rf1, rh1);
	    __m256i rh2 = _mm256_max_epi32(rf1, rh1);

	    //re-ordering the oids
	    mask1 = _mm256_cmpeq_epi32(ra1, ra2);
	    __m256 oa2 = _mm256_blendv_ps(oc1, oa1, _mm256_castsi256_ps(mask1));
	    __m256 oc2 = _mm256_blendv_ps(oa1, oc1, _mm256_castsi256_ps(mask1));

	    mask2 = _mm256_cmpeq_epi32(rb1, rb2);
	    __m256 ob2 = _mm256_blendv_ps(od1, ob1, _mm256_castsi256_ps(mask2));
	    __m256 od2 = _mm256_blendv_ps(ob1, od1, _mm256_castsi256_ps(mask2));

	    mask3 = _mm256_cmpeq_epi32(re1, re2);
	    __m256 oe2 = _mm256_blendv_ps(og1, oe1, _mm256_castsi256_ps(mask3));
	    __m256 og2 = _mm256_blendv_ps(oe1, og1, _mm256_castsi256_ps(mask3));

	    mask4 = _mm256_cmpeq_epi32(rf1, rf2);
	    __m256 of2 = _mm256_blendv_ps(oh1, of1, _mm256_castsi256_ps(mask4));
	    __m256 oh2 = _mm256_blendv_ps(of1, oh1, _mm256_castsi256_ps(mask4));

	    /* 3rd level of comparisons: [[1,2],[5,6],[0,4],[3,7]]*/
	    __m256i rb3 = _mm256_min_epi32(rb2, rc2);
	    __m256i rc3 = _mm256_max_epi32(rb2, rc2);

	    __m256i ra3 = _mm256_min_epi32(ra2, re2);
	    __m256i re3 = _mm256_max_epi32(ra2, re2);

	    __m256i rf3 = _mm256_min_epi32(rf2, rg2);
	    __m256i rg3 = _mm256_max_epi32(rf2, rg2);

	    __m256i rd3 = _mm256_min_epi32(rd2, rh2);
	    __m256i rh3 = _mm256_max_epi32(rd2, rh2);

	    //re-ordering the oids
	    mask1 = _mm256_cmpeq_epi32(rb2, rb3);
	    __m256 ob3 = _mm256_blendv_ps(oc2, ob2, _mm256_castsi256_ps(mask1));
	    __m256 oc3 = _mm256_blendv_ps(ob2, oc2, _mm256_castsi256_ps(mask1));

	    mask2 = _mm256_cmpeq_epi32(ra2, ra3);
	    __m256 oa3 = _mm256_blendv_ps(oe2, oa2, _mm256_castsi256_ps(mask2));
	    __m256 oe3 = _mm256_blendv_ps(oa2, oe2, _mm256_castsi256_ps(mask2));

	    mask3 = _mm256_cmpeq_epi32(rf2, rf3);
	    __m256 of3 = _mm256_blendv_ps(og2, of2, _mm256_castsi256_ps(mask3));
	    __m256 og3 = _mm256_blendv_ps(of2, og2, _mm256_castsi256_ps(mask3));

	    mask4 = _mm256_cmpeq_epi32(rd2, rd3);
	    __m256 od3 = _mm256_blendv_ps(oh2, od2, _mm256_castsi256_ps(mask4));
	    __m256 oh3 = _mm256_blendv_ps(od2, oh2, _mm256_castsi256_ps(mask4));

	    /* 4th level of comparisons : [[1,5],[2,6]]*/
	    __m256i rb4 = _mm256_min_epi32(rb3, rf3);
	    __m256i rf4 = _mm256_max_epi32(rb3, rf3);

	    __m256i rc4 = _mm256_min_epi32(rc3, rg3);
	    __m256i rg4 = _mm256_max_epi32(rc3, rg3);

	    //re-ordering the oids
	    mask1 = _mm256_cmpeq_epi32(rb3, rb4);
	    __m256 ob4 = _mm256_blendv_ps(of3, ob3, _mm256_castsi256_ps(mask1));
	    __m256 of4= _mm256_blendv_ps(ob3, of3, _mm256_castsi256_ps(mask1));

	    mask2 = _mm256_cmpeq_epi32(rc3, rc4);
	    __m256 oc4 = _mm256_blendv_ps(og3, oc3, _mm256_castsi256_ps(mask2));
	    __m256 og4 = _mm256_blendv_ps(oc3, og3, _mm256_castsi256_ps(mask2));

	    /* 5th level of comparisons : [[1,4],[3,6]]*/
	    __m256i rb5 = _mm256_min_epi32(rb4, re3);
	    __m256i re5 = _mm256_max_epi32(rb4, re3);

	    __m256i rd5 = _mm256_min_epi32(rd3, rg4);
	    __m256i rg5 = _mm256_max_epi32(rd3, rg4);

	    //re-ordering the oids
	    mask3 = _mm256_cmpeq_epi32(rb4, rb5);
	    __m256 ob5 = _mm256_blendv_ps(oe3, ob4, _mm256_castsi256_ps(mask3));
	    __m256 oe5= _mm256_blendv_ps(ob4, oe3, _mm256_castsi256_ps(mask3));

	    mask4 = _mm256_cmpeq_epi32(rd3, rd5);
	    __m256 od5 = _mm256_blendv_ps(og4, od3, _mm256_castsi256_ps(mask4));
	    __m256 og5 = _mm256_blendv_ps(od3, og4, _mm256_castsi256_ps(mask4));

	    /* 6th level of comparisons : [[2,4],[3,5]]*/
	    __m256i rc6 = _mm256_min_epi32(rc4, re5);
	    __m256i re6 = _mm256_max_epi32(rc4, re5);

	    __m256i rd6 = _mm256_min_epi32(rd5, rf4);
	    __m256i rf6 = _mm256_max_epi32(rd5, rf4);

	    //re-ordering the oids
	    mask1 = _mm256_cmpeq_epi32(rc4, rc6);
	    __m256 oc6 = _mm256_blendv_ps(oe5, oc4, _mm256_castsi256_ps(mask1));
	    __m256 oe6= _mm256_blendv_ps(oc4, oe5, _mm256_castsi256_ps(mask1));

	    mask2 = _mm256_cmpeq_epi32(rd5, rd6);
	    __m256 od6 = _mm256_blendv_ps(of4, od5, _mm256_castsi256_ps(mask2));
	    __m256 of6 = _mm256_blendv_ps(od5, of4, _mm256_castsi256_ps(mask2));

	    /* 7th level of comparisons : [[3,4]]*/
	    __m256i rd7 = _mm256_min_epi32(rd6, re6);
	   	__m256i re7 = _mm256_max_epi32(rd6, re6);

	    mask3 = _mm256_cmpeq_epi32(rd6, rd7);
	    __m256 od7 = _mm256_blendv_ps(oe6, od6, _mm256_castsi256_ps(mask3));
	    __m256 oe7= _mm256_blendv_ps(od6, oe6, _mm256_castsi256_ps(mask3));

	    /* results are in ra3, rb5, rc6, rd7, re7, rf6, rg5, rh3 */
	    /* re-ordered oids are in oa3, ob5, oc6, od7, oe7, of6, og5, oh3 */
	    /**
	     * Initial data and transposed data looks like following:
	     *  a3={ x1  x2  x3  x4  x5  x6  x7  x8 }
	     *  b5={ x9  x10 x11 x12 x13 x14 x15 x16}
	     *  c6={ x17 x18 x19 x20 x21 x22 x23 x24}
	     *  d7={ x25 x26 x27 x28 x29 x30 x31 x32}  === Transpose ===>
	     *  e7={ x33 x34 x35 x36 x37 x38 x39 x40}
	     *  f6={ x41 x42 x43 x44 x45 x46 x47 x48}
	     *  g5={ x49 x50 x51 x52 x53 x54 x55 x56}
	     *  h3={ x57 x58 x59 x60 x61 x62 x63 x64}
	     *
	     *	After Step 1 and Step 2:
	     *	a4={ x1  x9  x17 x25| x5  x13 x21 x29}
	     *  b4={ x2  x10 x18 x26| x6  x14 x22 x30}
	     *  a5={ x3  x11 x19 x27| x7  x15 x23 x31}
	     *  b5={ x4  x12 x20 x28| x8  x16 x24 x32}  === Transpose ===>
	     *  _______________________________________
	     *  e4={ x33 x34 x35 x36| x37 x45 x53 x61}
	     *  f4={ x41 x42 x43 x44| x38 x46 x54 x62}
	     *  e5={ x49 x50 x51 x52| x39 x47 x55 x63}
	     *  f5={ x57 x58 x59 x60| x40 x48 x56 x64}
	     *
	     *	After Step 3:
	     *	a4={ x1  x9  x17 x25| x5  x13 x21 x29}
	     *  b4={ x2  x10 x18 x26| x6  x14 x22 x30}
	     *  a5={ x3  x11 x19 x27| x7  x15 x23 x31}
	     *  b5={ x4  x12 x20 x28| x8  x16 x24 x32}  === Transpose ===>
	     *  _______________________________________
	     *  e4={ x33 x41 x49 x57| x37 x45 x53 x61}
	     *  f4={ x34 x42 x50 x58| x38 x46 x54 x62}
	     *  e5={ x35 x43 x51 x59| x39 x47 x55 x63}
	     *  f5={ x36 x44 x52 x60| x40 x48 x56 x64}
	     *
	     * final state:
	     *	c ={ x1  x9  x17 x25 x33 x41 x49 x57}
	     *  g ={ x2  x10 x18 x26 x34 x42 x50 x58}
	     *  c1={ x3  x11 x19 x27 x35 x43 x51 x59}
	     *  g1={ x4  x12 x20 x28 x36 x44 x52 x60}
	     *  d ={ x5  x13 x21 x29 x37 x45 x53 x61}
	     *  h ={ x6  x14 x22 x30 x38 x46 x54 x62}
	     *  d1={ x7  x15 x23 x31 x39 x47 x55 x63}
	     *  h1={ x8  x16 x24 x32 x40 x48 x56 x64}
	     *
	     */
	    /* Step (1): shuffle left-up corner*/
	    /* Step 1.1 ==> shuffle x2 and x9 - shuffle x4 and x11 */
	    ra 	= _mm256_unpacklo_epi32(ra3, rb5);
	    rb 	= _mm256_unpackhi_epi32(ra3, rb5);
	    ra1 = _mm256_unpacklo_epi64(ra, rb);
	    rb1 = _mm256_unpackhi_epi64(ra, rb);

	    oa 	= _mm256_unpacklo_ps(oa3, ob5);
	    ob 	= _mm256_unpackhi_ps(oa3, ob5);
	    oa1 = _mm256_castsi256_ps(_mm256_unpacklo_epi64(_mm256_castps_si256(oa), _mm256_castps_si256(ob)));
	    ob1 = _mm256_castsi256_ps(_mm256_unpackhi_epi64(_mm256_castps_si256(oa), _mm256_castps_si256(ob)));

	    /* Step 1.2 ==> shuffle x18 and x25 - shuffle x20 and x27 */
	    ra2	= _mm256_unpacklo_epi32(rc6, rd7);
	    rb2	= _mm256_unpackhi_epi32(rc6, rd7);
	    ra3 = _mm256_unpacklo_epi64(ra2, rb2);
	    rb3 = _mm256_unpackhi_epi64(ra2, rb2);

	    oa2	= _mm256_unpacklo_ps(oc6, od7);
	    ob2	= _mm256_unpackhi_ps(oc6, od7);
	    oa3 = _mm256_castpd_ps(_mm256_unpacklo_pd(_mm256_castps_pd(oa2), _mm256_castps_pd(ob2)));
	    ob3 = _mm256_castpd_ps(_mm256_unpackhi_pd(_mm256_castps_pd(oa2), _mm256_castps_pd(ob2)));

	    /* Step 1.3 shuffle (x3,x11) and (x17,x25) pairs */
	    /* Step 1.4 shuffle (x4,x12) and (x18,x26) pairs */
	    __m256i ra4	= _mm256_unpacklo_epi64(ra1, ra3);
	    __m256i ra5	= _mm256_unpackhi_epi64(ra1, ra3);
	    rb4	= _mm256_unpacklo_epi64(rb1, rb3);
	    rb5	= _mm256_unpackhi_epi64(rb1, rb3);

	    __m256d oa4	= _mm256_unpacklo_pd(_mm256_castps_pd(oa1), _mm256_castps_pd(oa3));
	    __m256d oa5	= _mm256_unpackhi_pd(_mm256_castps_pd(oa1), _mm256_castps_pd(oa3));
	    ob4	= _mm256_castpd_ps(_mm256_unpacklo_pd(_mm256_castps_pd(ob1), _mm256_castps_pd(ob3)));
	    ob5	= _mm256_castpd_ps(_mm256_unpackhi_pd(_mm256_castps_pd(ob1), _mm256_castps_pd(ob3)));


	    /* Step (2): shuffle right-down corner*/
	    /* Step 2.1 */
	    re 	= _mm256_unpacklo_epi32(re7, rf6);
	    rf 	= _mm256_unpackhi_epi32(re7, rf6);
	    re1 = _mm256_unpacklo_epi64(re, rf);
	    rf1 = _mm256_unpackhi_epi64(re, rf);

	    __m256 oe_tmp 	= _mm256_unpacklo_ps(oe7, of6);
	    __m256 of_tmp 	= _mm256_unpackhi_ps(oe7, of6);
	    oe1 = _mm256_castpd_ps(_mm256_unpacklo_pd(_mm256_castps_pd(oe_tmp), _mm256_castps_pd(of_tmp)));
	    of1 = _mm256_castpd_ps(_mm256_unpackhi_pd(_mm256_castps_pd(oe_tmp), _mm256_castps_pd(of_tmp)));

	    /* Step 2.2 */
	    re2	= _mm256_unpacklo_epi32(rg5, rh3);
	    rf2	= _mm256_unpackhi_epi32(rg5, rh3);
	    re3 = _mm256_unpacklo_epi64(re2, rf2);
	    rf3 = _mm256_unpackhi_epi64(re2, rf2);

	    oe2	= _mm256_unpacklo_ps(og5, oh3);
	    of2	= _mm256_unpackhi_ps(og5, oh3);
	    oe3 = _mm256_castpd_ps(_mm256_unpacklo_pd(_mm256_castps_pd(oe2), _mm256_castps_pd(of2)));
	    of3 = _mm256_castpd_ps(_mm256_unpackhi_pd(_mm256_castps_pd(oe2), _mm256_castps_pd(of2)));

	    /* Step 2.3 */
	    /* Step 2.4 */
	    __m256i re4	= _mm256_unpacklo_epi64(re1, re3);
	    re5	= _mm256_unpackhi_epi64(re1, re3);
	    rf4	= _mm256_unpacklo_epi64(rf1, rf3);
	    __m256i rf5	= _mm256_unpackhi_epi64(rf1, rf3);

	    __m256d oe4	= _mm256_unpacklo_pd(_mm256_castps_pd(oe1), _mm256_castps_pd(oe3));
	    oe5	= _mm256_castpd_ps(_mm256_unpackhi_pd(_mm256_castps_pd(oe1), _mm256_castps_pd(oe3)));
	    of4	= _mm256_castpd_ps(_mm256_unpacklo_pd(_mm256_castps_pd(of1), _mm256_castps_pd(of3)));
	    __m256d of5	= _mm256_unpackhi_pd(_mm256_castps_pd(of1), _mm256_castps_pd(of3));


	    /* Step (3): switch right-top and left-down corner*/
	    rc = _mm256_permute2f128_si256(ra4, re4, 0x20);
	    rd = _mm256_permute2f128_si256(ra4, re4, 0x31);

	    rg = _mm256_permute2f128_si256(rb4, rf4, 0x20);
	    rh = _mm256_permute2f128_si256(rb4, rf4, 0x31);

	    rc1 = _mm256_permute2f128_si256(ra5, re5, 0x20);
	    rd1 = _mm256_permute2f128_si256(ra5, re5, 0x31);

	    rg1 = _mm256_permute2f128_si256(rb5, rf5, 0x20);
	    rh1 = _mm256_permute2f128_si256(rb5, rf5, 0x31);

	    oc = _mm256_permute2f128_ps(_mm256_castpd_ps(oa4), _mm256_castpd_ps(oe4), 0x20);
	    od = _mm256_permute2f128_ps(_mm256_castpd_ps(oa4), _mm256_castpd_ps(oe4), 0x31);

	    og = _mm256_permute2f128_ps(ob4, of4, 0x20);
	    oh = _mm256_permute2f128_ps(ob4, of4, 0x31);

	    oc1 = _mm256_permute2f128_ps(_mm256_castpd_ps(oa5), oe5, 0x20);
	    od1 = _mm256_permute2f128_ps(_mm256_castpd_ps(oa5), oe5, 0x31);

	    og1 = _mm256_permute2f128_ps(ob5, _mm256_castpd_ps(of5), 0x20);
	    oh1 = _mm256_permute2f128_ps(ob5, _mm256_castpd_ps(of5), 0x31);

	    /* after this, results are in rc, rg, rc1, rg1, rd, rh, rd1, rh1 */
	/* IACA_END */
	    /* store values */
	    _mm256_store_si256(reinterpret_cast<__m256i *> (out_val), rc);
	    _mm256_store_si256(reinterpret_cast<__m256i *> (out_val + 8), rg);
	    _mm256_store_si256(reinterpret_cast<__m256i *> (out_val + 16), rc1);
	    _mm256_store_si256(reinterpret_cast<__m256i *> (out_val + 24), rg1);
	    _mm256_store_si256(reinterpret_cast<__m256i *> (out_val + 32), rd);
	    _mm256_store_si256(reinterpret_cast<__m256i *> (out_val + 40), rh);
	    _mm256_store_si256(reinterpret_cast<__m256i *> (out_val + 48), rd1);
	    _mm256_store_si256(reinterpret_cast<__m256i *> (out_val + 56), rh1);

	    _mm256_store_ps(reinterpret_cast<float *> (out_oid), oc);
	    _mm256_store_ps(reinterpret_cast<float *> (out_oid + 8), og);
	    _mm256_store_ps(reinterpret_cast<float *> (out_oid + 16), oc1);
	    _mm256_store_ps(reinterpret_cast<float *> (out_oid + 24), og1);
	    _mm256_store_ps(reinterpret_cast<float *> (out_oid + 32), od);
	    _mm256_store_ps(reinterpret_cast<float *> (out_oid + 40), oh);
	    _mm256_store_ps(reinterpret_cast<float *> (out_oid + 48), od1);
	    _mm256_store_ps(reinterpret_cast<float *> (out_oid + 56), oh1);
}



#if 0
inline void __attribute__((always_inline))
Mergesort::inregister_sort_int16_aligned(int16_t * in_val, int16_t * out_val,
		int32_t * in_oid, int32_t * out_oid)
{
	/* IACA_START */
	    __m256i ra = _mm256_load_si256 (reinterpret_cast<__m256i const *>(in_val));
	    __m256i rb = _mm256_load_si256 (reinterpret_cast<__m256i const *>(in_val + 16));

	    __m256i oa_1 = _mm256_load_si256 (reinterpret_cast<__m256i const *>(in_oid));
	    __m256i oa_2 = _mm256_load_si256 (reinterpret_cast<__m256i const *>(in_oid + 8));
	    __m256i ob_1 = _mm256_load_si256 (reinterpret_cast<__m256i const *>(in_oid + 16));
	    __m256i ob_2 = _mm256_load_si256 (reinterpret_cast<__m256i const *>(in_oid + 24));

	    /* odd-even sorting network begins: using Hibbard's Algorithm.
There are 65 comparators in this network,
grouped into 15 parallel operations.

[[0,1],[2,3],[4,5],[6,7],[8,9],[10,11],[12,13],[14,15]]   =>level 1
[[0,2],[1,3],[4,6],[5,7],[8,10],[9,11],[12,14],[13,15]] =>level 2
[[1,2],[5,6],[0,4],[3,7],[9,10],[13,14],[8,12],[11,15]] =>level 3
[[1,5],[2,6],[9,13],[10,14],[0,8],[7,15]] =>level 4
[[1,4],[3,6],[9,12],[11,14]] =>level 5
[[2,4],[3,5],[10,12],[11,13],[1,9],[6,14]] =>level 6
[[3,4],[11,12],[1,8],[2,10],[5,13],[7,14]] =>level 7
[[3,11],[2,8],[4,12],[7,13]] =>level 8
[[3,10],[5,12]] =>level 9
[[3,9],[6,12]] =>level 10
[[3,8],[7,12],[5,9],[6,10]] =>level 11
[[4,8],[7,11]] =>level 12
[[5,8],[7,10]] =>level 13
[[6,8],[7,9]] =>level 14
[[7,8]] =>level 15
		 *	important: may be optimized, e.g., re-use the variables as in in_register_sort_int64()
	     */
	    /* 1st level of comparisons: [[0,1],[2,3],[4,5],[6,7],[8,9],[10,11],[12,13],[14,15]]*/
	    __m256i ra1 = _mm256_min_epi16(ra, rb);

	    //re-ordering the oids
	    __m256i mask1 = _mm256_cmpeq_epi16(ra, ra1);
	    __m256 low = reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask1, 0)));
		__m256 high = reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask1, 1)));
	    __m256 oa1_1 = _mm256_blendv_ps(reinterpret_cast<__m256>(ob_1), reinterpret_cast<__m256>(oa_1), low);
	    __m256 oa1_2 = _mm256_blendv_ps(reinterpret_cast<__m256>(ob_2), reinterpret_cast<__m256>(oa_2), high);


	    _mm256_store_ps(reinterpret_cast<float *>(out_oid), oa1_1);
	    _mm256_store_ps(reinterpret_cast<float *>(out_oid+8), oa1_2);

	    //reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(a, 0)))
}
#endif

#if 1
inline void __attribute__((always_inline))
Mergesort::inregister_sort_int16_aligned(int16_t * in_val, int16_t * out_val,
		int32_t * in_oid, int32_t * out_oid)
{
	//for debug: output matrix and oid
#if 0
	for (int i = 0; i < 16; i++) {
		for (int j = 0; j < 16; ++j) {
			std::cout << in_val[i*16+j] << "[" << in_oid[i*16+j] << "]|";
		}
		std::cout << std::endl;
	}
#endif

	/* IACA_START */
	    __m256i ra = _mm256_load_si256 (reinterpret_cast<__m256i const *>(in_val));
	    __m256i rb = _mm256_load_si256 (reinterpret_cast<__m256i const *>(in_val + 16));
	    __m256i rc = _mm256_load_si256 (reinterpret_cast<__m256i const *>(in_val + 32));
	    __m256i rd = _mm256_load_si256 (reinterpret_cast<__m256i const *>(in_val + 48));
	    __m256i re = _mm256_load_si256 (reinterpret_cast<__m256i const *>(in_val + 64));
	    __m256i rf = _mm256_load_si256 (reinterpret_cast<__m256i const *>(in_val + 80));
	    __m256i rg = _mm256_load_si256 (reinterpret_cast<__m256i const *>(in_val + 96));
	    __m256i rh = _mm256_load_si256 (reinterpret_cast<__m256i const *>(in_val + 112));
	    __m256i ri = _mm256_load_si256 (reinterpret_cast<__m256i const *>(in_val + 128));
	    __m256i rj = _mm256_load_si256 (reinterpret_cast<__m256i const *>(in_val + 144));
	    __m256i rk = _mm256_load_si256 (reinterpret_cast<__m256i const *>(in_val + 160));
	    __m256i rl = _mm256_load_si256 (reinterpret_cast<__m256i const *>(in_val + 176));
	    __m256i rm = _mm256_load_si256 (reinterpret_cast<__m256i const *>(in_val + 192));
	    __m256i rn = _mm256_load_si256 (reinterpret_cast<__m256i const *>(in_val + 208));
	    __m256i ro = _mm256_load_si256 (reinterpret_cast<__m256i const *>(in_val + 224));
	    __m256i rp = _mm256_load_si256 (reinterpret_cast<__m256i const *>(in_val + 240));

	    __m256 oa_1 = _mm256_load_ps (reinterpret_cast<float const *>(in_oid));
	    __m256 oa_2 = _mm256_load_ps (reinterpret_cast<float const *>(in_oid + 8));
	    __m256 ob_1 = _mm256_load_ps (reinterpret_cast<float const *>(in_oid + 16));
	    __m256 ob_2 = _mm256_load_ps (reinterpret_cast<float const *>(in_oid + 24));
	    __m256 oc_1 = _mm256_load_ps (reinterpret_cast<float const *>(in_oid + 32));
	    __m256 oc_2 = _mm256_load_ps (reinterpret_cast<float const *>(in_oid + 40));
	    __m256 od_1 = _mm256_load_ps (reinterpret_cast<float const *>(in_oid + 48));
	    __m256 od_2 = _mm256_load_ps (reinterpret_cast<float const *>(in_oid + 56));
	    __m256 oe_1 = _mm256_load_ps (reinterpret_cast<float const *>(in_oid + 64));
	    __m256 oe_2 = _mm256_load_ps (reinterpret_cast<float const *>(in_oid + 72));
	    __m256 of_1 = _mm256_load_ps (reinterpret_cast<float const *>(in_oid + 80));
	    __m256 of_2 = _mm256_load_ps (reinterpret_cast<float const *>(in_oid + 88));
	    __m256 og_1 = _mm256_load_ps (reinterpret_cast<float const *>(in_oid + 96));
	    __m256 og_2 = _mm256_load_ps (reinterpret_cast<float const *>(in_oid + 104));
	    __m256 oh_1 = _mm256_load_ps (reinterpret_cast<float const *>(in_oid + 112));
	    __m256 oh_2 = _mm256_load_ps (reinterpret_cast<float const *>(in_oid + 120));
	    __m256 oi_1 = _mm256_load_ps (reinterpret_cast<float const *>(in_oid + 128));
	    __m256 oi_2 = _mm256_load_ps (reinterpret_cast<float const *>(in_oid + 136));
	    __m256 oj_1 = _mm256_load_ps (reinterpret_cast<float const *>(in_oid + 144));
	    __m256 oj_2 = _mm256_load_ps (reinterpret_cast<float const *>(in_oid + 152));
	    __m256 ok_1 = _mm256_load_ps (reinterpret_cast<float const *>(in_oid + 160));
	    __m256 ok_2 = _mm256_load_ps (reinterpret_cast<float const *>(in_oid + 168));
	    __m256 ol_1 = _mm256_load_ps (reinterpret_cast<float const *>(in_oid + 176));
	    __m256 ol_2 = _mm256_load_ps (reinterpret_cast<float const *>(in_oid + 184));
	    __m256 om_1 = _mm256_load_ps (reinterpret_cast<float const *>(in_oid + 192));
	    __m256 om_2 = _mm256_load_ps (reinterpret_cast<float const *>(in_oid + 200));
	    __m256 on_1 = _mm256_load_ps (reinterpret_cast<float const *>(in_oid + 208));
	    __m256 on_2 = _mm256_load_ps (reinterpret_cast<float const *>(in_oid + 216));
	    __m256 oo_1 = _mm256_load_ps (reinterpret_cast<float const *>(in_oid + 224));
	    __m256 oo_2 = _mm256_load_ps (reinterpret_cast<float const *>(in_oid + 232));
	    __m256 op_1 = _mm256_load_ps (reinterpret_cast<float const *>(in_oid + 240));
	    __m256 op_2 = _mm256_load_ps (reinterpret_cast<float const *>(in_oid + 248));

	    /* odd-even sorting network begins: using Hibbard's Algorithm.
There are 65 comparators in this network,
grouped into 15 parallel operations.

[[0,1],[2,3],[4,5],[6,7],[8,9],[10,11],[12,13],[14,15]]   =>level 1
[[0,2],[1,3],[4,6],[5,7],[8,10],[9,11],[12,14],[13,15]] =>level 2
[[1,2],[5,6],[0,4],[3,7],[9,10],[13,14],[8,12],[11,15]] =>level 3
[[1,5],[2,6],[9,13],[10,14],[0,8],[7,15]] =>level 4
[[1,4],[3,6],[9,12],[11,14]] =>level 5
[[2,4],[3,5],[10,12],[11,13],[1,9],[6,14]] =>level 6
[[3,4],[11,12],[1,8],[2,10],[5,13],[7,14]] =>level 7
[[3,11],[2,8],[4,12],[7,13]] =>level 8
[[3,10],[5,12]] =>level 9
[[3,9],[6,12]] =>level 10
[[3,8],[7,12],[5,9],[6,10]] =>level 11
[[4,8],[7,11]] =>level 12
[[5,8],[7,10]] =>level 13
[[6,8],[7,9]] =>level 14
[[7,8]] =>level 15
		 *	important: may be optimized, e.g., re-use the variables as in in_register_sort_int64()
	     */
	    /* 1st level of comparisons: [[0,1],[2,3],[4,5],[6,7],[8,9],[10,11],[12,13],[14,15]]*/
	    __m256i ra1 = _mm256_min_epi16(ra, rb);
	    __m256i rb1 = _mm256_max_epi16(ra, rb);

	    __m256i rc1 = _mm256_min_epi16(rc, rd);
	    __m256i rd1 = _mm256_max_epi16(rc, rd);

	    __m256i re1 = _mm256_min_epi16(re, rf);
	    __m256i rf1 = _mm256_max_epi16(re, rf);

	    __m256i rg1 = _mm256_min_epi16(rg, rh);
	    __m256i rh1 = _mm256_max_epi16(rg, rh);

	    __m256i ri1 = _mm256_min_epi16(ri, rj);
	    __m256i rj1 = _mm256_max_epi16(ri, rj);

	    __m256i rk1 = _mm256_min_epi16(rk, rl);
	    __m256i rl1 = _mm256_max_epi16(rk, rl);

	    __m256i rm1 = _mm256_min_epi16(rm, rn);
	    __m256i rn1 = _mm256_max_epi16(rm, rn);

	    __m256i ro1 = _mm256_min_epi16(ro, rp);
	    __m256i rp1 = _mm256_max_epi16(ro, rp);

	    //re-ordering the oids
	    __m256i mask1 = _mm256_cmpeq_epi16(ra, ra1);
	    __m256 low = reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask1, 0)));
	    __m256 high = reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask1, 1)));
	    __m256 oa1_1 = _mm256_blendv_ps(ob_1, oa_1, low);
	    __m256 oa1_2 = _mm256_blendv_ps(ob_2, oa_2, high);

	    __m256 ob1_1 = _mm256_blendv_ps(oa_1, ob_1, low);
	    __m256 ob1_2 = _mm256_blendv_ps(oa_2, ob_2, high);

	    __m256i mask2 = _mm256_cmpeq_epi16(rc, rc1);
	    low	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask2, 0)));
		high = reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask2, 1)));
	    __m256 oc1_1 = _mm256_blendv_ps(od_1, oc_1, low);
	    __m256 oc1_2 = _mm256_blendv_ps(od_2, oc_2, high);

	    __m256 od1_1 = _mm256_blendv_ps(oc_1, od_1, low);
	    __m256 od1_2 = _mm256_blendv_ps(oc_2, od_2, high);

	    __m256i mask3 = _mm256_cmpeq_epi16(re, re1);
	    low	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask3, 0)));
		high = reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask3, 1)));
	    __m256 oe1_1 = _mm256_blendv_ps(of_1, oe_1, low);
	    __m256 oe1_2 = _mm256_blendv_ps(of_2, oe_2, high);

	    __m256 of1_1 = _mm256_blendv_ps(oe_1, of_1, low);
	    __m256 of1_2 = _mm256_blendv_ps(oe_2, of_2, high);

	    __m256i mask4 = _mm256_cmpeq_epi16(rg, rg1);
	    low	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask4, 0)));
		high = reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask4, 1)));
	    __m256 og1_1 = _mm256_blendv_ps(oh_1, og_1, low);
	    __m256 og1_2 = _mm256_blendv_ps(oh_2, og_2, high);

	    __m256 oh1_1 = _mm256_blendv_ps(og_1, oh_1, low);
	    __m256 oh1_2 = _mm256_blendv_ps(og_2, oh_2, high);

	    __m256i mask5 = _mm256_cmpeq_epi16(ri, ri1);
	    low	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask5, 0)));
		high = reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask5, 1)));
	    __m256 oi1_1 = _mm256_blendv_ps(oj_1, oi_1, low);
	    __m256 oi1_2 = _mm256_blendv_ps(oj_2, oi_2, high);

	    __m256 oj1_1 = _mm256_blendv_ps(oi_1, oj_1, low);
	    __m256 oj1_2 = _mm256_blendv_ps(oi_2, oj_2, high);

	    __m256i mask6 = _mm256_cmpeq_epi16(rk, rk1);
	    low	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask6, 0)));
		high = reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask6, 1)));
	    __m256 ok1_1 = _mm256_blendv_ps(ol_1, ok_1, low);
	    __m256 ok1_2 = _mm256_blendv_ps(ol_2, ok_2, high);

	    __m256 ol1_1 = _mm256_blendv_ps(ok_1, ol_1, low);
	    __m256 ol1_2 = _mm256_blendv_ps(ok_2, ol_2, high);

	    __m256i mask7 = _mm256_cmpeq_epi16(rm, rm1);
	    low	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask7, 0)));
		high = reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask7, 1)));
	    __m256 om1_1 = _mm256_blendv_ps(on_1, om_1, low);
	    __m256 om1_2 = _mm256_blendv_ps(on_2, om_2, high);

	    __m256 on1_1 = _mm256_blendv_ps(om_1, on_1, low);
	    __m256 on1_2 = _mm256_blendv_ps(om_2, on_2, high);

	    __m256i mask8 = _mm256_cmpeq_epi16(ro, ro1);
	    low	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask8, 0)));
		high = reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask8, 1)));
	    __m256 oo1_1 = _mm256_blendv_ps(op_1, oo_1, low);
	    __m256 oo1_2 = _mm256_blendv_ps(op_2, oo_2, high);

	    __m256 op1_1 = _mm256_blendv_ps(oo_1, op_1, low);
	    __m256 op1_2 = _mm256_blendv_ps(oo_2, op_2, high);

	    /* 2nd level of comparisons: [[0,2],[1,3],[4,6],[5,7],[8,10],[9,11],[12,14],[13,15]] */
	    /* 2nd level of comparisons: [[a,c],[b,d],[e,g],[f,h],[i,k],[j,l],[m,o],[n,p]] */
	    __m256i ra2 = _mm256_min_epi16(ra1, rc1);
	    __m256i rc2 = _mm256_max_epi16(ra1, rc1);

	    __m256i rb2 = _mm256_min_epi16(rb1, rd1);
	    __m256i rd2 = _mm256_max_epi16(rb1, rd1);

	    __m256i re2 = _mm256_min_epi16(re1, rg1);
	    __m256i rg2 = _mm256_max_epi16(re1, rg1);

	    __m256i rf2 = _mm256_min_epi16(rf1, rh1);
	    __m256i rh2 = _mm256_max_epi16(rf1, rh1);

	    __m256i ri2 = _mm256_min_epi16(ri1, rk1);
	    __m256i rk2 = _mm256_max_epi16(ri1, rk1);

	    __m256i rj2 = _mm256_min_epi16(rj1, rl1);
	    __m256i rl2 = _mm256_max_epi16(rj1, rl1);

	    __m256i rm2 = _mm256_min_epi16(rm1, ro1);
	    __m256i ro2 = _mm256_max_epi16(rm1, ro1);

	    __m256i rn2 = _mm256_min_epi16(rn1, rp1);
	    __m256i rp2 = _mm256_max_epi16(rn1, rp1);

	    //re-ordering the oids
	    mask1 = _mm256_cmpeq_epi16(ra1, ra2);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask1, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask1, 1)));
	    __m256 oa2_1 = _mm256_blendv_ps(oc1_1, oa1_1,  low);
	    __m256 oa2_2 = _mm256_blendv_ps(oc1_2, oa1_2,  high);

	    __m256 oc2_1 = _mm256_blendv_ps(oa1_1, oc1_1,  low);
	    __m256 oc2_2 = _mm256_blendv_ps(oa1_2, oc1_2,  high);

	    mask2 = _mm256_cmpeq_epi16(rb1, rb2);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask2, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask2, 1)));
	    __m256 ob2_1 = _mm256_blendv_ps(od1_1, ob1_1,  low);
	    __m256 ob2_2 = _mm256_blendv_ps(od1_2, ob1_2,  high);

	    __m256 od2_1 = _mm256_blendv_ps(ob1_1, od1_1,  low);
	    __m256 od2_2 = _mm256_blendv_ps(ob1_2, od1_2,  high);

	    mask3 = _mm256_cmpeq_epi16(re1, re2);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask3, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask3, 1)));
	    __m256 oe2_1 = _mm256_blendv_ps(og1_1, oe1_1,  low);
	    __m256 oe2_2 = _mm256_blendv_ps(og1_2, oe1_2,  high);

	    __m256 og2_1 = _mm256_blendv_ps(oe1_1, og1_1,  low);
	    __m256 og2_2 = _mm256_blendv_ps(oe1_2, og1_2,  high);

	    mask4 = _mm256_cmpeq_epi16(rf1, rf2);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask4, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask4, 1)));
	    __m256 of2_1 = _mm256_blendv_ps(oh1_1, of1_1,  low);
	    __m256 of2_2 = _mm256_blendv_ps(oh1_2, of1_2,  high);

	    __m256 oh2_1 = _mm256_blendv_ps(of1_1, oh1_1,  low);
	    __m256 oh2_2 = _mm256_blendv_ps(of1_2, oh1_2,  high);

	    mask5 = _mm256_cmpeq_epi16(ri1, ri2);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask5, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask5, 1)));
	    __m256 oi2_1 = _mm256_blendv_ps(ok1_1, oi1_1,  low);
	    __m256 oi2_2 = _mm256_blendv_ps(ok1_2, oi1_2,  high);

	    __m256 ok2_1 = _mm256_blendv_ps(oi1_1, ok1_1,  low);
	    __m256 ok2_2 = _mm256_blendv_ps(oi1_2, ok1_2,  high);

	    mask6 = _mm256_cmpeq_epi16(rj1, rj2);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask6, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask6, 1)));
	    __m256 oj2_1 = _mm256_blendv_ps(ol1_1, oj1_1,  low);
	    __m256 oj2_2 = _mm256_blendv_ps(ol1_2, oj1_2,  high);

	    __m256 ol2_1 = _mm256_blendv_ps(oj1_1, ol1_1,  low);
	    __m256 ol2_2 = _mm256_blendv_ps(oj1_2, ol1_2,  high);

	    mask7 = _mm256_cmpeq_epi16(rm1, rm2);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask7, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask7, 1)));
	    __m256 om2_1 = _mm256_blendv_ps(oo1_1, om1_1,  low);
	    __m256 om2_2 = _mm256_blendv_ps(oo1_2, om1_2,  high);

	    __m256 oo2_1 = _mm256_blendv_ps(om1_1, oo1_1,  low);
	    __m256 oo2_2 = _mm256_blendv_ps(om1_2, oo1_2,  high);

	    mask8 = _mm256_cmpeq_epi16(rn1, rn2);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask8, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask8, 1)));
	    __m256 on2_1 = _mm256_blendv_ps(op1_1, on1_1,  low);
	    __m256 on2_2 = _mm256_blendv_ps(op1_2, on1_2,  high);

	    __m256 op2_1 = _mm256_blendv_ps(on1_1, op1_1,  low);
	    __m256 op2_2 = _mm256_blendv_ps(on1_2, op1_2,  high);

	    /* 3rd level of comparisons: [[1,2],[5,6],[0,4],[3,7],[9,10],[13,14],[8,12],[11,15]]*/
	    /* 3rd level of comparisons: [[b,c],[f,g],[a,e],[d,h],[j,k],[n,o],[i,m],[l,p]]*/
	    __m256i rb3 = _mm256_min_epi16(rb2, rc2);
	    __m256i rc3 = _mm256_max_epi16(rb2, rc2);

	    __m256i rf3 = _mm256_min_epi16(rf2, rg2);
	    __m256i rg3 = _mm256_max_epi16(rf2, rg2);

	    __m256i ra3 = _mm256_min_epi16(ra2, re2);
	    __m256i re3 = _mm256_max_epi16(ra2, re2);

	    __m256i rd3 = _mm256_min_epi16(rd2, rh2);
	    __m256i rh3 = _mm256_max_epi16(rd2, rh2);

	    __m256i rj3 = _mm256_min_epi16(rj2, rk2);
	    __m256i rk3 = _mm256_max_epi16(rj2, rk2);

	    __m256i rn3 = _mm256_min_epi16(rn2, ro2);
	    __m256i ro3 = _mm256_max_epi16(rn2, ro2);

	    __m256i ri3 = _mm256_min_epi16(ri2, rm2);
	    __m256i rm3 = _mm256_max_epi16(ri2, rm2);

	    __m256i rl3 = _mm256_min_epi16(rl2, rp2);
	    __m256i rp3 = _mm256_max_epi16(rl2, rp2);

	    //re-ordering the oids
	    mask1 = _mm256_cmpeq_epi16(rb2, rb3);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask1, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask1, 1)));
	    __m256 ob3_1 = _mm256_blendv_ps(oc2_1, ob2_1,  low);
	    __m256 ob3_2 = _mm256_blendv_ps(oc2_2, ob2_2,  high);

	    __m256 oc3_1 = _mm256_blendv_ps(ob2_1, oc2_1,  low);
	    __m256 oc3_2 = _mm256_blendv_ps(ob2_2, oc2_2,  high);

	    mask2 = _mm256_cmpeq_epi16(rf2, rf3);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask2, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask2, 1)));
	    __m256 of3_1 = _mm256_blendv_ps(og2_1, of2_1,  low);
	    __m256 of3_2 = _mm256_blendv_ps(og2_2, of2_2,  high);

	    __m256 og3_1 = _mm256_blendv_ps(of2_1, og2_1,  low);
	    __m256 og3_2 = _mm256_blendv_ps(of2_2, og2_2,  high);

	    mask3 = _mm256_cmpeq_epi16(ra2, ra3);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask3, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask3, 1)));
	    __m256 oa3_1 = _mm256_blendv_ps(oe2_1, oa2_1,  low);
	    __m256 oa3_2 = _mm256_blendv_ps(oe2_2, oa2_2,  high);

	    __m256 oe3_1 = _mm256_blendv_ps(oa2_1, oe2_1,  low);
	    __m256 oe3_2 = _mm256_blendv_ps(oa2_2, oe2_2,  high);

	    mask4 = _mm256_cmpeq_epi16(rd2, rd3);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask4, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask4, 1)));
	    __m256 od3_1 = _mm256_blendv_ps(oh2_1, od2_1,  low);
	    __m256 od3_2 = _mm256_blendv_ps(oh2_2, od2_2,  high);

	    __m256 oh3_1 = _mm256_blendv_ps(od2_1, oh2_1,  low);
	    __m256 oh3_2 = _mm256_blendv_ps(od2_2, oh2_2,  high);

	    mask5 = _mm256_cmpeq_epi16(rj2, rj3);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask5, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask5, 1)));
	    __m256 oj3_1 = _mm256_blendv_ps(ok2_1, oj2_1,  low);
	    __m256 oj3_2 = _mm256_blendv_ps(ok2_2, oj2_2,  high);

	    __m256 ok3_1 = _mm256_blendv_ps(oj2_1, ok2_1,  low);
	    __m256 ok3_2 = _mm256_blendv_ps(oj2_2, ok2_2,  high);

	    mask6 = _mm256_cmpeq_epi16(rn2, rn3);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask6, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask6, 1)));
	    __m256 on3_1 = _mm256_blendv_ps(oo2_1, on2_1,  low);
	    __m256 on3_2 = _mm256_blendv_ps(oo2_2, on2_2,  high);

	    __m256 oo3_1 = _mm256_blendv_ps(on2_1, oo2_1,  low);
	    __m256 oo3_2 = _mm256_blendv_ps(on2_2, oo2_2,  high);

	    mask7 = _mm256_cmpeq_epi16(ri2, ri3);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask7, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask7, 1)));
	    __m256 oi3_1 = _mm256_blendv_ps(om2_1, oi2_1,  low);
	    __m256 oi3_2 = _mm256_blendv_ps(om2_2, oi2_2,  high);

	    __m256 om3_1 = _mm256_blendv_ps(oi2_1, om2_1,  low);
	    __m256 om3_2 = _mm256_blendv_ps(oi2_2, om2_2,  high);

	    mask8 = _mm256_cmpeq_epi16(rl2, rl3);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask8, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask8, 1)));
	    __m256 ol3_1 = _mm256_blendv_ps(op2_1, ol2_1,  low);
	    __m256 ol3_2 = _mm256_blendv_ps(op2_2, ol2_2,  high);

	    __m256 op3_1 = _mm256_blendv_ps(ol2_1, op2_1,  low);
	    __m256 op3_2 = _mm256_blendv_ps(ol2_2, op2_2,  high);

	    /* 4th level of comparisons : [[1,5],[2,6],[9,13],[10,14],[0,8],[7,15]]*/
	    /* 4th level of comparisons : [[b,f],[c,g],[j,n],[k,o],[a,i],[h,p]]*/
	    __m256i rb4 = _mm256_min_epi16(rb3, rf3);
	    __m256i rf4 = _mm256_max_epi16(rb3, rf3);

	    __m256i rc4 = _mm256_min_epi16(rc3, rg3);
	    __m256i rg4 = _mm256_max_epi16(rc3, rg3);

	    __m256i rj4 = _mm256_min_epi16(rj3, rn3);
	    __m256i rn4 = _mm256_max_epi16(rj3, rn3);

	    __m256i rk4 = _mm256_min_epi16(rk3, ro3);
	    __m256i ro4 = _mm256_max_epi16(rk3, ro3);

	    __m256i ra4 = _mm256_min_epi16(ra3, ri3);
	    __m256i ri4 = _mm256_max_epi16(ra3, ri3);

	    __m256i rh4 = _mm256_min_epi16(rh3, rp3);
	    __m256i rp4 = _mm256_max_epi16(rh3, rp3);

	    //re-ordering the oids
	    mask1 = _mm256_cmpeq_epi16(rb3, rb4);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask1, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask1, 1)));
	    __m256 ob4_1 = _mm256_blendv_ps(of3_1, ob3_1,  low);
	    __m256 ob4_2 = _mm256_blendv_ps(of3_2, ob3_2,  high);

	    __m256 of4_1 = _mm256_blendv_ps(ob3_1, of3_1,  low);
	    __m256 of4_2 = _mm256_blendv_ps(ob3_2, of3_2,  high);

	    mask2 = _mm256_cmpeq_epi16(rc3, rc4);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask2, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask2, 1)));
	    __m256 oc4_1 = _mm256_blendv_ps(og3_1, oc3_1,  low);
	    __m256 oc4_2 = _mm256_blendv_ps(og3_2, oc3_2,  high);

	    __m256 og4_1 = _mm256_blendv_ps(oc3_1, og3_1,  low);
	    __m256 og4_2 = _mm256_blendv_ps(oc3_2, og3_2,  high);

	    mask3 = _mm256_cmpeq_epi16(rj3, rj4);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask3, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask3, 1)));
	    __m256 oj4_1 = _mm256_blendv_ps(on3_1, oj3_1,  low);
	    __m256 oj4_2 = _mm256_blendv_ps(on3_2, oj3_2,  high);

	    __m256 on4_1 = _mm256_blendv_ps(oj3_1, on3_1,  low);
	    __m256 on4_2 = _mm256_blendv_ps(oj3_2, on3_2,  high);

	    mask4 = _mm256_cmpeq_epi16(rk3, rk4);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask4, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask4, 1)));
	    __m256 ok4_1 = _mm256_blendv_ps(oo3_1, ok3_1,  low);
	    __m256 ok4_2 = _mm256_blendv_ps(oo3_2, ok3_2,  high);

	    __m256 oo4_1 = _mm256_blendv_ps(ok3_1, oo3_1,  low);
	    __m256 oo4_2 = _mm256_blendv_ps(ok3_2, oo3_2,  high);

	    mask5 = _mm256_cmpeq_epi16(ra3, ra4);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask5, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask5, 1)));
	    __m256 oa4_1 = _mm256_blendv_ps(oi3_1, oa3_1,  low);
	    __m256 oa4_2 = _mm256_blendv_ps(oi3_2, oa3_2,  high);

	    __m256 oi4_1 = _mm256_blendv_ps(oa3_1, oi3_1,  low);
	    __m256 oi4_2 = _mm256_blendv_ps(oa3_2, oi3_2,  high);

	    mask6 = _mm256_cmpeq_epi16(rh3, rh4);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask6, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask6, 1)));
	    __m256 oh4_1 = _mm256_blendv_ps(op3_1, oh3_1,  low);
	    __m256 oh4_2 = _mm256_blendv_ps(op3_2, oh3_2,  high);

	    __m256 op4_1 = _mm256_blendv_ps(oh3_1, op3_1,  low);
	    __m256 op4_2 = _mm256_blendv_ps(oh3_2, op3_2,  high);

	    /* 5th level of comparisons : [[1,4],[3,6],[9,12],[11,14]]*/
	    /* 5th level of comparisons : [[b,e],[d,g],[j,m],[l,o]]*/
	    __m256i rb5 = _mm256_min_epi16(rb4, re3);
	    __m256i re5 = _mm256_max_epi16(rb4, re3);

	    __m256i rd5 = _mm256_min_epi16(rd3, rg4);
	    __m256i rg5 = _mm256_max_epi16(rd3, rg4);

	    __m256i rj5 = _mm256_min_epi16(rj4, rm3);
	    __m256i rm5 = _mm256_max_epi16(rj4, rm3);

	    __m256i rl5 = _mm256_min_epi16(rl3, ro4);
	    __m256i ro5 = _mm256_max_epi16(rl3, ro4);

	    //re-ordering the oids
	    mask1 = _mm256_cmpeq_epi16(rb4, rb5);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask1, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask1, 1)));
	    __m256 ob5_1 = _mm256_blendv_ps(oe3_1, ob4_1,  low);
	    __m256 ob5_2 = _mm256_blendv_ps(oe3_2, ob4_2,  high);

	    __m256 oe5_1 = _mm256_blendv_ps(ob4_1, oe3_1,  low);
	    __m256 oe5_2 = _mm256_blendv_ps(ob4_2, oe3_2,  high);

	    mask2 = _mm256_cmpeq_epi16(rd3, rd5);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask2, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask2, 1)));
	    __m256 od5_1 = _mm256_blendv_ps(og4_1, od3_1,  low);
	    __m256 od5_2 = _mm256_blendv_ps(og4_2, od3_2,  high);

	    __m256 og5_1 = _mm256_blendv_ps(od3_1, og4_1,  low);
	    __m256 og5_2 = _mm256_blendv_ps(od3_2, og4_2,  high);

	    mask3 = _mm256_cmpeq_epi16(rj4, rj5);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask3, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask3, 1)));
	    __m256 oj5_1 = _mm256_blendv_ps(om3_1, oj4_1,  low);
	    __m256 oj5_2 = _mm256_blendv_ps(om3_2, oj4_2,  high);

	    __m256 om5_1 = _mm256_blendv_ps(oj4_1, om3_1,  low);
	    __m256 om5_2 = _mm256_blendv_ps(oj4_2, om3_2,  high);

	    mask4 = _mm256_cmpeq_epi16(rl3, rl5);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask4, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask4, 1)));
	    __m256 ol5_1 = _mm256_blendv_ps(oo4_1, ol3_1,  low);
	    __m256 ol5_2 = _mm256_blendv_ps(oo4_2, ol3_2,  high);

	    __m256 oo5_1 = _mm256_blendv_ps(ol3_1, oo4_1,  low);
	    __m256 oo5_2 = _mm256_blendv_ps(ol3_2, oo4_2,  high);

	    /* 6th level of comparisons : [[2,4],[3,5],[10,12],[11,13],[1,9],[6,14]]*/
	    /* 6th level of comparisons : [[c,e],[d,f],[k,m],[l,n],[b,j],[g,o]]*/
	    __m256i rc6 = _mm256_min_epi16(rc4, re5);
	    __m256i re6 = _mm256_max_epi16(rc4, re5);

	    __m256i rd6 = _mm256_min_epi16(rd5, rf4);
	    __m256i rf6 = _mm256_max_epi16(rd5, rf4);

	    __m256i rk6 = _mm256_min_epi16(rk4, rm5);
	    __m256i rm6 = _mm256_max_epi16(rk4, rm5);

	    __m256i rl6 = _mm256_min_epi16(rl5, rn4);
	    __m256i rn6 = _mm256_max_epi16(rl5, rn4);

	    __m256i rb6 = _mm256_min_epi16(rb5, rj5);
	    __m256i rj6 = _mm256_max_epi16(rb5, rj5);

	    __m256i rg6 = _mm256_min_epi16(rg5, ro5);
	    __m256i ro6 = _mm256_max_epi16(rg5, ro5);

	    //re-ordering the oids
	    mask1 = _mm256_cmpeq_epi16(rc4, rc6);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask1, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask1, 1)));
	    __m256 oc6_1 = _mm256_blendv_ps(oe5_1, oc4_1,  low);
	    __m256 oc6_2 = _mm256_blendv_ps(oe5_2, oc4_2,  high);

	    __m256 oe6_1 = _mm256_blendv_ps(oc4_1, oe5_1,  low);
	    __m256 oe6_2 = _mm256_blendv_ps(oc4_2, oe5_2,  high);

	    mask2 = _mm256_cmpeq_epi16(rd5, rd6);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask2, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask2, 1)));
	    __m256 od6_1 = _mm256_blendv_ps(of4_1, od5_1,  low);
	    __m256 od6_2 = _mm256_blendv_ps(of4_2, od5_2,  high);

	    __m256 of6_1 = _mm256_blendv_ps(od5_1, of4_1,  low);
	    __m256 of6_2 = _mm256_blendv_ps(od5_2, of4_2,  high);

	    mask3 = _mm256_cmpeq_epi16(rk4, rk6);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask3, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask3, 1)));
	    __m256 ok6_1 = _mm256_blendv_ps(om5_1, ok4_1,  low);
	    __m256 ok6_2 = _mm256_blendv_ps(om5_2, ok4_2,  high);

	    __m256 om6_1 = _mm256_blendv_ps(ok4_1, om5_1,  low);
	    __m256 om6_2 = _mm256_blendv_ps(ok4_2, om5_2,  high);

	    mask4 = _mm256_cmpeq_epi16(rl5, rl6);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask4, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask4, 1)));
	    __m256 ol6_1 = _mm256_blendv_ps(on4_1, ol5_1,  low);
	    __m256 ol6_2 = _mm256_blendv_ps(on4_2, ol5_2,  high);

	    __m256 on6_1 = _mm256_blendv_ps(ol5_1, on4_1,  low);
	    __m256 on6_2 = _mm256_blendv_ps(ol5_2, on4_2,  high);

	    mask5 = _mm256_cmpeq_epi16(rb5, rb6);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask5, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask5, 1)));
	    __m256 ob6_1 = _mm256_blendv_ps(oj5_1, ob5_1,  low);
	    __m256 ob6_2 = _mm256_blendv_ps(oj5_2, ob5_2,  high);

	    __m256 oj6_1 = _mm256_blendv_ps(ob5_1, oj5_1,  low);
	    __m256 oj6_2 = _mm256_blendv_ps(ob5_2, oj5_2,  high);

	    mask6 = _mm256_cmpeq_epi16(rg5, rg6);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask6, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask6, 1)));
	    __m256 og6_1 = _mm256_blendv_ps(oo5_1, og5_1,  low);
	    __m256 og6_2 = _mm256_blendv_ps(oo5_2, og5_2,  high);

	    __m256 oo6_1 = _mm256_blendv_ps(og5_1, oo5_1,  low);
	    __m256 oo6_2 = _mm256_blendv_ps(og5_2, oo5_2,  high);

	    /* 7th level of comparisons : [[3,4],[11,12],[1,8],[2,10],[5,13],[7,14]]*/
	    /* 7th level of comparisons : [[d,e],[l,m],[b,i],[c,k],[f,n],[h,o]]*/
	    __m256i rd7 = _mm256_min_epi16(rd6, re6);
	    __m256i re7 = _mm256_max_epi16(rd6, re6);

	    __m256i rl7 = _mm256_min_epi16(rl6, rm6);
	    __m256i rm7 = _mm256_max_epi16(rl6, rm6);

	    __m256i rb7 = _mm256_min_epi16(rb6, ri4);
	    __m256i ri7 = _mm256_max_epi16(rb6, ri4);

	    __m256i rc7 = _mm256_min_epi16(rc6, rk6);
	    __m256i rk7 = _mm256_max_epi16(rc6, rk6);

	    __m256i rf7 = _mm256_min_epi16(rf6, rn6);
	    __m256i rn7 = _mm256_max_epi16(rf6, rn6);

	    __m256i rh7 = _mm256_min_epi16(rh4, ro6);
	    __m256i ro7 = _mm256_max_epi16(rh4, ro6);

	    //re-ordering the oids
	    mask1 = _mm256_cmpeq_epi16(rd6, rd7);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask1, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask1, 1)));
	    __m256 od7_1 = _mm256_blendv_ps(oe6_1, od6_1,  low);
	    __m256 od7_2 = _mm256_blendv_ps(oe6_2, od6_2,  high);

	    __m256 oe7_1 = _mm256_blendv_ps(od6_1, oe6_1,  low);
	    __m256 oe7_2 = _mm256_blendv_ps(od6_2, oe6_2,  high);

	    mask2 = _mm256_cmpeq_epi16(rl6, rl7);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask2, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask2, 1)));
	    __m256 ol7_1 = _mm256_blendv_ps(om6_1, ol6_1,  low);
	    __m256 ol7_2 = _mm256_blendv_ps(om6_2, ol6_2,  high);

	    __m256 om7_1 = _mm256_blendv_ps(ol6_1, om6_1,  low);
	    __m256 om7_2 = _mm256_blendv_ps(ol6_2, om6_2,  high);

	    mask3 = _mm256_cmpeq_epi16(rb6, rb7);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask3, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask3, 1)));
	    __m256 ob7_1 = _mm256_blendv_ps(oi4_1, ob6_1,  low);
	    __m256 ob7_2 = _mm256_blendv_ps(oi4_2, ob6_2,  high);

	    __m256 oi7_1 = _mm256_blendv_ps(ob6_1, oi4_1,  low);
	    __m256 oi7_2 = _mm256_blendv_ps(ob6_2, oi4_2,  high);

	    mask4 = _mm256_cmpeq_epi16(rc6, rc7);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask4, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask4, 1)));
	    __m256 oc7_1 = _mm256_blendv_ps(ok6_1, oc6_1,  low);
	    __m256 oc7_2 = _mm256_blendv_ps(ok6_2, oc6_2,  high);

	    __m256 ok7_1 = _mm256_blendv_ps(oc6_1, ok6_1,  low);
	    __m256 ok7_2 = _mm256_blendv_ps(oc6_2, ok6_2,  high);

	    mask5 = _mm256_cmpeq_epi16(rf6, rf7);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask5, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask5, 1)));
	    __m256 of7_1 = _mm256_blendv_ps(on6_1, of6_1,  low);
	    __m256 of7_2 = _mm256_blendv_ps(on6_2, of6_2,  high);

	    __m256 on7_1 = _mm256_blendv_ps(of6_1, on6_1,  low);
	    __m256 on7_2 = _mm256_blendv_ps(of6_2, on6_2,  high);

	    mask6 = _mm256_cmpeq_epi16(rh4, rh7);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask6, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask6, 1)));
	    __m256 oh7_1 = _mm256_blendv_ps(oo6_1, oh4_1,  low);
	    __m256 oh7_2 = _mm256_blendv_ps(oo6_2, oh4_2,  high);

	    __m256 oo7_1 = _mm256_blendv_ps(oh4_1, oo6_1,  low);
	    __m256 oo7_2 = _mm256_blendv_ps(oh4_2, oo6_2,  high);

	    /* 8th level of comparisons : [[3,11],[2,8],[4,12],[7,13]]*/
	    /* 8th level of comparisons : [[d,l],[c,i],[e,m],[h,n]]*/
	    __m256i rd8 = _mm256_min_epi16(rd7, rl7);
	    __m256i rl8 = _mm256_max_epi16(rd7, rl7);

	    __m256i rc8 = _mm256_min_epi16(rc7, ri7);
	    __m256i ri8 = _mm256_max_epi16(rc7, ri7);

	    __m256i re8 = _mm256_min_epi16(re7, rm7);
	    __m256i rm8 = _mm256_max_epi16(re7, rm7);

	    __m256i rh8 = _mm256_min_epi16(rh7, rn7);
	    __m256i rn8 = _mm256_max_epi16(rh7, rn7);

	    //re-ordering the oids
	    mask1 = _mm256_cmpeq_epi16(rd7, rd8);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask1, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask1, 1)));
	    __m256 od8_1 = _mm256_blendv_ps(ol7_1, od7_1,  low);
	    __m256 od8_2 = _mm256_blendv_ps(ol7_2, od7_2,  high);

	    __m256 ol8_1 = _mm256_blendv_ps(od7_1, ol7_1,  low);
	    __m256 ol8_2 = _mm256_blendv_ps(od7_2, ol7_2,  high);

	    mask2 = _mm256_cmpeq_epi16(rc7, rc8);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask2, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask2, 1)));
	    __m256 oc8_1 = _mm256_blendv_ps(oi7_1, oc7_1,  low);
	    __m256 oc8_2 = _mm256_blendv_ps(oi7_2, oc7_2,  high);

	    __m256 oi8_1 = _mm256_blendv_ps(oc7_1, oi7_1,  low);
	    __m256 oi8_2 = _mm256_blendv_ps(oc7_2, oi7_2,  high);

	    mask3 = _mm256_cmpeq_epi16(re7, re8);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask3, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask3, 1)));
	    __m256 oe8_1 = _mm256_blendv_ps(om7_1, oe7_1,  low);
	    __m256 oe8_2 = _mm256_blendv_ps(om7_2, oe7_2,  high);

	    __m256 om8_1 = _mm256_blendv_ps(oe7_1, om7_1,  low);
	    __m256 om8_2 = _mm256_blendv_ps(oe7_2, om7_2,  high);

	    mask4 = _mm256_cmpeq_epi16(rh7, rh8);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask4, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask4, 1)));
	    __m256 oh8_1 = _mm256_blendv_ps(on7_1, oh7_1,  low);
	    __m256 oh8_2 = _mm256_blendv_ps(on7_2, oh7_2,  high);

	    __m256 on8_1 = _mm256_blendv_ps(oh7_1, on7_1,  low);
	    __m256 on8_2 = _mm256_blendv_ps(oh7_2, on7_2,  high);

	    /* 9th level of comparisons : [[3,10],[5,12]]*/
	    /* 9th level of comparisons : [[d,k],[f,m]]*/
	    __m256i rd9 = _mm256_min_epi16(rd8, rk7);
	    __m256i rk9 = _mm256_max_epi16(rd8, rk7);

	    __m256i rf9 = _mm256_min_epi16(rf7, rm8);
	    __m256i rm9 = _mm256_max_epi16(rf7, rm8);

	    //re-ordering the oids
	    mask1 = _mm256_cmpeq_epi16(rd8, rd9);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask1, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask1, 1)));
	    __m256 od9_1 = _mm256_blendv_ps(ok7_1, od8_1,  low);
	    __m256 od9_2 = _mm256_blendv_ps(ok7_2, od8_2,  high);

	    __m256 ok9_1 = _mm256_blendv_ps(od8_1, ok7_1,  low);
	    __m256 ok9_2 = _mm256_blendv_ps(od8_2, ok7_2,  high);

	    mask2 = _mm256_cmpeq_epi16(rf7, rf9);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask2, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask2, 1)));
	    __m256 of9_1 = _mm256_blendv_ps(om8_1, of7_1,  low);
	    __m256 of9_2 = _mm256_blendv_ps(om8_2, of7_2,  high);

	    __m256 om9_1 = _mm256_blendv_ps(of7_1, om8_1,  low);
	    __m256 om9_2 = _mm256_blendv_ps(of7_2, om8_2,  high);

	    /* 10th level of comparisons : [[3,9],[6,12]]*/
	    /* 10th level of comparisons : [[d,j],[g,m]]*/
	    __m256i rd10 = _mm256_min_epi16(rd9, rj6);
	    __m256i rj10 = _mm256_max_epi16(rd9, rj6);

	    __m256i rg10 = _mm256_min_epi16(rg6, rm9);
	    __m256i rm10 = _mm256_max_epi16(rg6, rm9);

	    //re-ordering the oids
	    mask1 = _mm256_cmpeq_epi16(rd9, rd10);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask1, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask1, 1)));
	    __m256 od10_1 = _mm256_blendv_ps(oj6_1, od9_1,  low);
	    __m256 od10_2 = _mm256_blendv_ps(oj6_2, od9_2,  high);

	    __m256 oj10_1 = _mm256_blendv_ps(od9_1, oj6_1,  low);
	    __m256 oj10_2 = _mm256_blendv_ps(od9_2, oj6_2,  high);

	    mask2 = _mm256_cmpeq_epi16(rg6, rg10);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask2, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask2, 1)));
	    __m256 og10_1 = _mm256_blendv_ps(om9_1, og6_1,  low);
	    __m256 og10_2 = _mm256_blendv_ps(om9_2, og6_2,  high);

	    __m256 om10_1 = _mm256_blendv_ps(og6_1, om9_1,  low);
	    __m256 om10_2 = _mm256_blendv_ps(og6_2, om9_2,  high);

	    /* 11th level of comparisons : [[3,8],[7,12],[5,9],[6,10]] */
	    /* 11th level of comparisons : [[d,i],[h,m],[f,j],[g,k]] */
	    __m256i rd11 = _mm256_min_epi16(rd10, ri8);
	    __m256i ri11 = _mm256_max_epi16(rd10, ri8);

	    __m256i rh11 = _mm256_min_epi16(rh8, rm10);
	    __m256i rm11 = _mm256_max_epi16(rh8, rm10);

	    __m256i rf11 = _mm256_min_epi16(rf9, rj10);
	    __m256i rj11 = _mm256_max_epi16(rf9, rj10);

	    __m256i rg11 = _mm256_min_epi16(rg10, rk9);
	    __m256i rk11 = _mm256_max_epi16(rg10, rk9);

	    //re-ordering the oids
	    mask1 = _mm256_cmpeq_epi16(rd10, rd11);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask1, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask1, 1)));
	    __m256 od11_1 = _mm256_blendv_ps(oi8_1, od10_1,  low);
	    __m256 od11_2 = _mm256_blendv_ps(oi8_2, od10_2,  high);

	    __m256 oi11_1 = _mm256_blendv_ps(od10_1, oi8_1,  low);
	    __m256 oi11_2 = _mm256_blendv_ps(od10_2, oi8_2,  high);

	    mask2 = _mm256_cmpeq_epi16(rh8, rh11);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask2, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask2, 1)));
	    __m256 oh11_1 = _mm256_blendv_ps(om10_1, oh8_1,  low);
	    __m256 oh11_2 = _mm256_blendv_ps(om10_2, oh8_2,  high);

	    __m256 om11_1 = _mm256_blendv_ps(oh8_1, om10_1,  low);
	    __m256 om11_2 = _mm256_blendv_ps(oh8_2, om10_2,  high);

	    mask3 = _mm256_cmpeq_epi16(rf9, rf11);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask3, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask3, 1)));
	    __m256 of11_1 = _mm256_blendv_ps(oj10_1, of9_1,  low);
	    __m256 of11_2 = _mm256_blendv_ps(oj10_2, of9_2,  high);

	    __m256 oj11_1 = _mm256_blendv_ps(of9_1, oj10_1,  low);
	    __m256 oj11_2 = _mm256_blendv_ps(of9_2, oj10_2,  high);

	    mask4 = _mm256_cmpeq_epi16(rg10, rg11);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask4, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask4, 1)));
	    __m256 og11_1 = _mm256_blendv_ps(ok9_1, og10_1,  low);
	    __m256 og11_2 = _mm256_blendv_ps(ok9_2, og10_2,  high);

	    __m256 ok11_1 = _mm256_blendv_ps(og10_1, ok9_1,  low);
	    __m256 ok11_2 = _mm256_blendv_ps(og10_2, ok9_2,  high);

	    /* 12th level of comparisons : [[4,8],[7,11]] */
	    /* 12th level of comparisons : [[e,i],[h,l]] */
	    __m256i re12 = _mm256_min_epi16(re8, ri11);
	    __m256i ri12 = _mm256_max_epi16(re8, ri11);

	    __m256i rh12 = _mm256_min_epi16(rh11, rl8);
	    __m256i rl12 = _mm256_max_epi16(rh11, rl8);

	    //re-ordering the oids
	    mask1 = _mm256_cmpeq_epi16(re8, re12);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask1, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask1, 1)));
	    __m256 oe12_1 = _mm256_blendv_ps(oi11_1, oe8_1,  low);
	    __m256 oe12_2 = _mm256_blendv_ps(oi11_2, oe8_2,  high);

	    __m256 oi12_1 = _mm256_blendv_ps(oe8_1, oi11_1,  low);
	    __m256 oi12_2 = _mm256_blendv_ps(oe8_2, oi11_2,  high);

	    mask2 = _mm256_cmpeq_epi16(rh11, rh12);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask2, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask2, 1)));
	    __m256 oh12_1 = _mm256_blendv_ps(ol8_1, oh11_1,  low);
	    __m256 oh12_2 = _mm256_blendv_ps(ol8_2, oh11_2,  high);

	    __m256 ol12_1 = _mm256_blendv_ps(oh11_1, ol8_1,  low);
	    __m256 ol12_2 = _mm256_blendv_ps(oh11_2, ol8_2,  high);

	    /* 13th level of comparisons : [[5,8],[7,10]] */
	    /* 13th level of comparisons : [[f,i],[h,k]] */
	    __m256i rf13 = _mm256_min_epi16(rf11, ri12);
	    __m256i ri13 = _mm256_max_epi16(rf11, ri12);

	    __m256i rh13 = _mm256_min_epi16(rh12, rk11);
	    __m256i rk13 = _mm256_max_epi16(rh12, rk11);

	    //re-ordering the oids
	    mask1 = _mm256_cmpeq_epi16(rf11, rf13);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask1, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask1, 1)));
	    __m256 of13_1 = _mm256_blendv_ps(oi12_1, of11_1,  low);
	    __m256 of13_2 = _mm256_blendv_ps(oi12_2, of11_2,  high);

	    __m256 oi13_1 = _mm256_blendv_ps(of11_1, oi12_1,  low);
	    __m256 oi13_2 = _mm256_blendv_ps(of11_2, oi12_2,  high);

	    mask2 = _mm256_cmpeq_epi16(rh12, rh13);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask2, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask2, 1)));
	    __m256 oh13_1 = _mm256_blendv_ps(ok11_1, oh12_1,  low);
	    __m256 oh13_2 = _mm256_blendv_ps(ok11_2, oh12_2,  high);

	    __m256 ok13_1 = _mm256_blendv_ps(oh12_1, ok11_1,  low);
	    __m256 ok13_2 = _mm256_blendv_ps(oh12_2, ok11_2,  high);

	    /* 14th level of comparisons : [[6,8],[7,9]] */
	    /* 14th level of comparisons : [[g,i],[h,j]] */
	    __m256i rg14 = _mm256_min_epi16(rg11, ri13);
	    __m256i ri14 = _mm256_max_epi16(rg11, ri13);

	    __m256i rh14 = _mm256_min_epi16(rh13, rj11);
	    __m256i rj14 = _mm256_max_epi16(rh13, rj11);

	    //re-ordering the oids
	    mask1 = _mm256_cmpeq_epi16(rg11, rg14);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask1, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask1, 1)));
	    __m256 og14_1 = _mm256_blendv_ps(oi13_1, og11_1,  low);
	    __m256 og14_2 = _mm256_blendv_ps(oi13_2, og11_2,  high);

	    __m256 oi14_1 = _mm256_blendv_ps(og11_1, oi13_1,  low);
	    __m256 oi14_2 = _mm256_blendv_ps(og11_2, oi13_2,  high);

	    mask2 = _mm256_cmpeq_epi16(rh13, rh14);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask2, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask2, 1)));
	    __m256 oh14_1 = _mm256_blendv_ps(oj11_1, oh13_1,  low);
	    __m256 oh14_2 = _mm256_blendv_ps(oj11_2, oh13_2,  high);

	    __m256 oj14_1 = _mm256_blendv_ps(oh13_1, oj11_1,  low);
	    __m256 oj14_2 = _mm256_blendv_ps(oh13_2, oj11_2,  high);


	    /* 15th level of comparisons : [[7,8]] */
	    /* 15th level of comparisons : [[h,i]] */
	    __m256i rh15 = _mm256_min_epi16(rh14, ri14);
	    __m256i ri15 = _mm256_max_epi16(rh14, ri14);

	    mask1 = _mm256_cmpeq_epi16(rh14, rh15);
	    low 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask1, 0)));
	    high 	= reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(mask1, 1)));
	    __m256 oh15_1 = _mm256_blendv_ps(oi14_1, oh14_1,  low);
	    __m256 oh15_2 = _mm256_blendv_ps(oi14_2, oh14_2,  high);

	    __m256 oi15_1 = _mm256_blendv_ps(oh14_1, oi14_1,  low);
	    __m256 oi15_2 = _mm256_blendv_ps(oh14_2, oi14_2,  high);

	    /* results are in ra4, rb7, rc8, rd11, re12, rf13, rg14, rh15, ri15, rj14, rk13, rl12, rm11, rn8, ro7, rp4
	     * oids are in (ra4_1, ra4_2), (rb7_1, rb7_2), ...
	     * then, do the transpose...
	     */
	    /* Step 1: 16-bank-size unpacking */
	    __m256i ra_shuffle 	= _mm256_unpacklo_epi16(ra4, rb7);
	    __m256i rb_shuffle 	= _mm256_unpackhi_epi16(ra4, rb7);

	    __m256i rc_shuffle 	= _mm256_unpacklo_epi16(rc8, rd11);
	    __m256i rd_shuffle 	= _mm256_unpackhi_epi16(rc8, rd11);

	    __m256i re_shuffle 	= _mm256_unpacklo_epi16(re12, rf13);
	    __m256i rf_shuffle 	= _mm256_unpackhi_epi16(re12, rf13);

	    __m256i rg_shuffle 	= _mm256_unpacklo_epi16(rg14, rh15);
	    __m256i rh_shuffle 	= _mm256_unpackhi_epi16(rg14, rh15);

	    __m256i ri_shuffle 	= _mm256_unpacklo_epi16(ri15, rj14);
	    __m256i rj_shuffle 	= _mm256_unpackhi_epi16(ri15, rj14);

	    __m256i rk_shuffle 	= _mm256_unpacklo_epi16(rk13, rl12);
	    __m256i rl_shuffle 	= _mm256_unpackhi_epi16(rk13, rl12);

	    __m256i rm_shuffle 	= _mm256_unpacklo_epi16(rm11, rn8);
	    __m256i rn_shuffle 	= _mm256_unpackhi_epi16(rm11, rn8);

	    __m256i ro_shuffle 	= _mm256_unpacklo_epi16(ro7, rp4);
	    __m256i rp_shuffle 	= _mm256_unpackhi_epi16(ro7, rp4);

	    /* Corresponding operations on oids */
	    /* Note: shuffling oids do not exactly follow the pattern */
	    __m256i oa_1_shuffle 	= _mm256_castps_si256(_mm256_unpacklo_ps(oa4_1, ob7_1));
	    __m256i ob_1_shuffle 	= _mm256_castps_si256(_mm256_unpackhi_ps(oa4_1, ob7_1));
	    __m256i oa_2_shuffle 	= _mm256_castps_si256(_mm256_unpacklo_ps(oa4_2, ob7_2));
	    __m256i ob_2_shuffle 	= _mm256_castps_si256(_mm256_unpackhi_ps(oa4_2, ob7_2));

	    __m256i oc_1_shuffle 	= _mm256_castps_si256(_mm256_unpacklo_ps(oc8_1, od11_1));
	    __m256i od_1_shuffle 	= _mm256_castps_si256(_mm256_unpackhi_ps(oc8_1, od11_1));
	    __m256i oc_2_shuffle 	= _mm256_castps_si256(_mm256_unpacklo_ps(oc8_2, od11_2));
	    __m256i od_2_shuffle 	= _mm256_castps_si256(_mm256_unpackhi_ps(oc8_2, od11_2));

	    __m256i oe_1_shuffle 	= _mm256_castps_si256(_mm256_unpacklo_ps(oe12_1, of13_1));
	    __m256i of_1_shuffle 	= _mm256_castps_si256(_mm256_unpackhi_ps(oe12_1, of13_1));
	    __m256i oe_2_shuffle 	= _mm256_castps_si256(_mm256_unpacklo_ps(oe12_2, of13_2));
	    __m256i of_2_shuffle 	= _mm256_castps_si256(_mm256_unpackhi_ps(oe12_2, of13_2));

	    __m256i og_1_shuffle 	= _mm256_castps_si256(_mm256_unpacklo_ps(og14_1, oh15_1));
	    __m256i oh_1_shuffle 	= _mm256_castps_si256(_mm256_unpackhi_ps(og14_1, oh15_1));
	    __m256i og_2_shuffle 	= _mm256_castps_si256(_mm256_unpacklo_ps(og14_2, oh15_2));
	    __m256i oh_2_shuffle 	= _mm256_castps_si256(_mm256_unpackhi_ps(og14_2, oh15_2));

	    __m256i oi_1_shuffle 	= _mm256_castps_si256(_mm256_unpacklo_ps(oi15_1, oj14_1));
	    __m256i oj_1_shuffle 	= _mm256_castps_si256(_mm256_unpackhi_ps(oi15_1, oj14_1));
	    __m256i oi_2_shuffle 	= _mm256_castps_si256(_mm256_unpacklo_ps(oi15_2, oj14_2));
	    __m256i oj_2_shuffle 	= _mm256_castps_si256(_mm256_unpackhi_ps(oi15_2, oj14_2));

	    __m256i ok_1_shuffle 	= _mm256_castps_si256(_mm256_unpacklo_ps(ok13_1, ol12_1));
	    __m256i ol_1_shuffle 	= _mm256_castps_si256(_mm256_unpackhi_ps(ok13_1, ol12_1));
	    __m256i ok_2_shuffle 	= _mm256_castps_si256(_mm256_unpacklo_ps(ok13_2, ol12_2));
	    __m256i ol_2_shuffle 	= _mm256_castps_si256(_mm256_unpackhi_ps(ok13_2, ol12_2));

	    __m256i om_1_shuffle 	= _mm256_castps_si256(_mm256_unpacklo_ps(om11_1, on8_1));
	    __m256i on_1_shuffle 	= _mm256_castps_si256(_mm256_unpackhi_ps(om11_1, on8_1));
	    __m256i om_2_shuffle 	= _mm256_castps_si256(_mm256_unpacklo_ps(om11_2, on8_2));
	    __m256i on_2_shuffle 	= _mm256_castps_si256(_mm256_unpackhi_ps(om11_2, on8_2));

	    __m256i oo_1_shuffle 	= _mm256_castps_si256(_mm256_unpacklo_ps(oo7_1, op4_1));
	    __m256i op_1_shuffle 	= _mm256_castps_si256(_mm256_unpackhi_ps(oo7_1, op4_1));
	    __m256i oo_2_shuffle 	= _mm256_castps_si256(_mm256_unpacklo_ps(oo7_2, op4_2));
	    __m256i op_2_shuffle 	= _mm256_castps_si256(_mm256_unpackhi_ps(oo7_2, op4_2));

	    /* Step 2: 32-bank-size unpacking*/
	    /* (2.1) first-round */
	    __m256i ra1_shuffle 	= _mm256_unpacklo_epi32(ra_shuffle, rb_shuffle);
	    __m256i rb1_shuffle 	= _mm256_unpackhi_epi32(ra_shuffle, rb_shuffle);

	    __m256i rc1_shuffle 	= _mm256_unpacklo_epi32(rc_shuffle, rd_shuffle);
	    __m256i rd1_shuffle 	= _mm256_unpackhi_epi32(rc_shuffle, rd_shuffle);

	    __m256i re1_shuffle 	= _mm256_unpacklo_epi32(re_shuffle, rf_shuffle);
	    __m256i rf1_shuffle 	= _mm256_unpackhi_epi32(re_shuffle, rf_shuffle);

	    __m256i rg1_shuffle 	= _mm256_unpacklo_epi32(rg_shuffle, rh_shuffle);
	    __m256i rh1_shuffle 	= _mm256_unpackhi_epi32(rg_shuffle, rh_shuffle);

	    __m256i ri1_shuffle 	= _mm256_unpacklo_epi32(ri_shuffle, rj_shuffle);
	    __m256i rj1_shuffle 	= _mm256_unpackhi_epi32(ri_shuffle, rj_shuffle);

	    __m256i rk1_shuffle 	= _mm256_unpacklo_epi32(rk_shuffle, rl_shuffle);
	    __m256i rl1_shuffle 	= _mm256_unpackhi_epi32(rk_shuffle, rl_shuffle);

	    __m256i rm1_shuffle 	= _mm256_unpacklo_epi32(rm_shuffle, rn_shuffle);
	    __m256i rn1_shuffle 	= _mm256_unpackhi_epi32(rm_shuffle, rn_shuffle);

	    __m256i ro1_shuffle 	= _mm256_unpacklo_epi32(ro_shuffle, rp_shuffle);
	    __m256i rp1_shuffle 	= _mm256_unpackhi_epi32(ro_shuffle, rp_shuffle);

	    /* (2.2) second-round */
	    __m256i ra2_shuffle 	= _mm256_unpacklo_epi32(ra1_shuffle, rc1_shuffle);
	    __m256i rb2_shuffle 	= _mm256_unpackhi_epi32(ra1_shuffle, rc1_shuffle);

	    __m256i rc2_shuffle 	= _mm256_unpacklo_epi32(rb1_shuffle, rd1_shuffle);
	    __m256i rd2_shuffle 	= _mm256_unpackhi_epi32(rb1_shuffle, rd1_shuffle);

	    __m256i re2_shuffle 	= _mm256_unpacklo_epi32(re1_shuffle, rg1_shuffle);
	    __m256i rf2_shuffle 	= _mm256_unpackhi_epi32(re1_shuffle, rg1_shuffle);

	    __m256i rg2_shuffle 	= _mm256_unpacklo_epi32(rf1_shuffle, rh1_shuffle);
	    __m256i rh2_shuffle 	= _mm256_unpackhi_epi32(rf1_shuffle, rh1_shuffle);

	    __m256i ri2_shuffle 	= _mm256_unpacklo_epi32(ri1_shuffle, rk1_shuffle);
	    __m256i rj2_shuffle 	= _mm256_unpackhi_epi32(ri1_shuffle, rk1_shuffle);

	    __m256i rk2_shuffle 	= _mm256_unpacklo_epi32(rj1_shuffle, rl1_shuffle);
	    __m256i rl2_shuffle 	= _mm256_unpackhi_epi32(rj1_shuffle, rl1_shuffle);

	    __m256i rm2_shuffle 	= _mm256_unpacklo_epi32(rm1_shuffle, ro1_shuffle);
	    __m256i rn2_shuffle 	= _mm256_unpackhi_epi32(rm1_shuffle, ro1_shuffle);

	    __m256i ro2_shuffle 	= _mm256_unpacklo_epi32(rn1_shuffle, rp1_shuffle);
	    __m256i rp2_shuffle 	= _mm256_unpackhi_epi32(rn1_shuffle, rp1_shuffle);

	    /* Corresponding operations on oids */ /*can use _mm256_blend_pd() to replace*/
	    /* Do not follow the same pattern as oids */
	    __m256d oa2_1_shuffle 	= _mm256_shuffle_pd(_mm256_castsi256_pd(oa_1_shuffle), _mm256_castsi256_pd(oc_1_shuffle), 0);
	    __m256d ob2_1_shuffle 	= _mm256_shuffle_pd(_mm256_castsi256_pd(oa_1_shuffle), _mm256_castsi256_pd(oc_1_shuffle), 0x0f);

	    __m256d oa2_2_shuffle 	= _mm256_shuffle_pd(_mm256_castsi256_pd(oa_2_shuffle), _mm256_castsi256_pd(oc_2_shuffle), 0);
	    __m256d ob2_2_shuffle 	= _mm256_shuffle_pd(_mm256_castsi256_pd(oa_2_shuffle), _mm256_castsi256_pd(oc_2_shuffle), 0x0f);

	    __m256d oc2_1_shuffle 	= _mm256_shuffle_pd(_mm256_castsi256_pd(ob_1_shuffle), _mm256_castsi256_pd(od_1_shuffle), 0);
	    __m256d od2_1_shuffle 	= _mm256_shuffle_pd(_mm256_castsi256_pd(ob_1_shuffle), _mm256_castsi256_pd(od_1_shuffle), 0x0f);

	    __m256d oc2_2_shuffle 	= _mm256_shuffle_pd(_mm256_castsi256_pd(ob_2_shuffle), _mm256_castsi256_pd(od_2_shuffle), 0);
	    __m256d od2_2_shuffle 	= _mm256_shuffle_pd(_mm256_castsi256_pd(ob_2_shuffle), _mm256_castsi256_pd(od_2_shuffle), 0x0f);

	    __m256d oe2_1_shuffle 	= _mm256_shuffle_pd(_mm256_castsi256_pd(oe_1_shuffle), _mm256_castsi256_pd(og_1_shuffle), 0);
	    __m256d of2_1_shuffle 	= _mm256_shuffle_pd(_mm256_castsi256_pd(oe_1_shuffle), _mm256_castsi256_pd(og_1_shuffle), 0x0f);

	    __m256d oe2_2_shuffle 	= _mm256_shuffle_pd(_mm256_castsi256_pd(oe_2_shuffle), _mm256_castsi256_pd(og_2_shuffle), 0);
	    __m256d of2_2_shuffle 	= _mm256_shuffle_pd(_mm256_castsi256_pd(oe_2_shuffle), _mm256_castsi256_pd(og_2_shuffle), 0x0f);

	    __m256d og2_1_shuffle 	= _mm256_shuffle_pd(_mm256_castsi256_pd(of_1_shuffle), _mm256_castsi256_pd(oh_1_shuffle), 0);
	    __m256d oh2_1_shuffle 	= _mm256_shuffle_pd(_mm256_castsi256_pd(of_1_shuffle), _mm256_castsi256_pd(oh_1_shuffle), 0x0f);

	    __m256d og2_2_shuffle 	= _mm256_shuffle_pd(_mm256_castsi256_pd(of_2_shuffle), _mm256_castsi256_pd(oh_2_shuffle), 0);
	    __m256d oh2_2_shuffle 	= _mm256_shuffle_pd(_mm256_castsi256_pd(of_2_shuffle), _mm256_castsi256_pd(oh_2_shuffle), 0x0f);

	    __m256d oi2_1_shuffle 	= _mm256_shuffle_pd(_mm256_castsi256_pd(oi_1_shuffle), _mm256_castsi256_pd(ok_1_shuffle), 0);
	    __m256d oj2_1_shuffle 	= _mm256_shuffle_pd(_mm256_castsi256_pd(oi_1_shuffle), _mm256_castsi256_pd(ok_1_shuffle), 0x0f);

	    __m256d oi2_2_shuffle 	= _mm256_shuffle_pd(_mm256_castsi256_pd(oi_2_shuffle), _mm256_castsi256_pd(ok_2_shuffle), 0);
	    __m256d oj2_2_shuffle 	= _mm256_shuffle_pd(_mm256_castsi256_pd(oi_2_shuffle), _mm256_castsi256_pd(ok_2_shuffle), 0x0f);

	    __m256d ok2_1_shuffle 	= _mm256_shuffle_pd(_mm256_castsi256_pd(oj_1_shuffle), _mm256_castsi256_pd(ol_1_shuffle), 0);
	    __m256d ol2_1_shuffle 	= _mm256_shuffle_pd(_mm256_castsi256_pd(oj_1_shuffle), _mm256_castsi256_pd(ol_1_shuffle), 0x0f);

	    __m256d ok2_2_shuffle 	= _mm256_shuffle_pd(_mm256_castsi256_pd(oj_2_shuffle), _mm256_castsi256_pd(ol_2_shuffle), 0);
	    __m256d ol2_2_shuffle 	= _mm256_shuffle_pd(_mm256_castsi256_pd(oj_2_shuffle), _mm256_castsi256_pd(ol_2_shuffle), 0x0f);

	    __m256d om2_1_shuffle 	= _mm256_shuffle_pd(_mm256_castsi256_pd(om_1_shuffle), _mm256_castsi256_pd(oo_1_shuffle), 0);
	    __m256d on2_1_shuffle 	= _mm256_shuffle_pd(_mm256_castsi256_pd(om_1_shuffle), _mm256_castsi256_pd(oo_1_shuffle), 0x0f);

	    __m256d om2_2_shuffle 	= _mm256_shuffle_pd(_mm256_castsi256_pd(om_2_shuffle), _mm256_castsi256_pd(oo_2_shuffle), 0);
	    __m256d on2_2_shuffle 	= _mm256_shuffle_pd(_mm256_castsi256_pd(om_2_shuffle), _mm256_castsi256_pd(oo_2_shuffle), 0x0f);

	    __m256d oo2_1_shuffle 	= _mm256_shuffle_pd(_mm256_castsi256_pd(on_1_shuffle), _mm256_castsi256_pd(op_1_shuffle), 0);
	    __m256d op2_1_shuffle 	= _mm256_shuffle_pd(_mm256_castsi256_pd(on_1_shuffle), _mm256_castsi256_pd(op_1_shuffle), 0x0f);

	    __m256d oo2_2_shuffle 	= _mm256_shuffle_pd(_mm256_castsi256_pd(on_2_shuffle), _mm256_castsi256_pd(op_2_shuffle), 0);
	    __m256d op2_2_shuffle 	= _mm256_shuffle_pd(_mm256_castsi256_pd(on_2_shuffle), _mm256_castsi256_pd(op_2_shuffle), 0x0f);


	    /* Step 3: 64-bank-size unpacking*/
	    __m256i ra3_shuffle 	= _mm256_unpacklo_epi64(ra2_shuffle, re2_shuffle);
	    __m256i re3_shuffle 	= _mm256_unpackhi_epi64(ra2_shuffle, re2_shuffle);

	    __m256i rb3_shuffle 	= _mm256_unpacklo_epi64(rb2_shuffle, rf2_shuffle);
	    __m256i rf3_shuffle 	= _mm256_unpackhi_epi64(rb2_shuffle, rf2_shuffle);

	    __m256i rc3_shuffle 	= _mm256_unpacklo_epi64(rc2_shuffle, rg2_shuffle);
	    __m256i rg3_shuffle 	= _mm256_unpackhi_epi64(rc2_shuffle, rg2_shuffle);

	    __m256i rd3_shuffle 	= _mm256_unpacklo_epi64(rd2_shuffle, rh2_shuffle);
	    __m256i rh3_shuffle 	= _mm256_unpackhi_epi64(rd2_shuffle, rh2_shuffle);

	    __m256i ri3_shuffle 	= _mm256_unpacklo_epi64(ri2_shuffle, rm2_shuffle);
	    __m256i rm3_shuffle 	= _mm256_unpackhi_epi64(ri2_shuffle, rm2_shuffle);

	    __m256i rj3_shuffle 	= _mm256_unpacklo_epi64(rj2_shuffle, rn2_shuffle);
	    __m256i rn3_shuffle 	= _mm256_unpackhi_epi64(rj2_shuffle, rn2_shuffle);

	    __m256i rk3_shuffle 	= _mm256_unpacklo_epi64(rk2_shuffle, ro2_shuffle);
	    __m256i ro3_shuffle 	= _mm256_unpackhi_epi64(rk2_shuffle, ro2_shuffle);

	    __m256i rl3_shuffle 	= _mm256_unpacklo_epi64(rl2_shuffle, rp2_shuffle);
	    __m256i rp3_shuffle 	= _mm256_unpackhi_epi64(rl2_shuffle, rp2_shuffle);

	    /* corresponding operations on oid*/
	    __m256i oa3_1_shuffle 	= _mm256_permute2f128_si256(_mm256_castpd_si256(oa2_1_shuffle), _mm256_castpd_si256(oe2_1_shuffle), 0x20);
	    __m256i oa3_2_shuffle 	= _mm256_permute2f128_si256(_mm256_castpd_si256(oa2_2_shuffle), _mm256_castpd_si256(oe2_2_shuffle), 0x20);
	    __m256i oe3_1_shuffle 	= _mm256_permute2f128_si256(_mm256_castpd_si256(oa2_1_shuffle), _mm256_castpd_si256(oe2_1_shuffle), 0x31);
	    __m256i oe3_2_shuffle 	= _mm256_permute2f128_si256(_mm256_castpd_si256(oa2_2_shuffle), _mm256_castpd_si256(oe2_2_shuffle), 0x31);

	    __m256i ob3_1_shuffle 	= _mm256_permute2f128_si256(_mm256_castpd_si256(ob2_1_shuffle), _mm256_castpd_si256(of2_1_shuffle), 0x20);
	    __m256i ob3_2_shuffle 	= _mm256_permute2f128_si256(_mm256_castpd_si256(ob2_2_shuffle), _mm256_castpd_si256(of2_2_shuffle), 0x20);
	    __m256i of3_1_shuffle 	= _mm256_permute2f128_si256(_mm256_castpd_si256(ob2_1_shuffle), _mm256_castpd_si256(of2_1_shuffle), 0x31);
	    __m256i of3_2_shuffle 	= _mm256_permute2f128_si256(_mm256_castpd_si256(ob2_2_shuffle), _mm256_castpd_si256(of2_2_shuffle), 0x31);

	    __m256i oc3_1_shuffle 	= _mm256_permute2f128_si256(_mm256_castpd_si256(oc2_1_shuffle), _mm256_castpd_si256(og2_1_shuffle), 0x20);
	    __m256i oc3_2_shuffle 	= _mm256_permute2f128_si256(_mm256_castpd_si256(oc2_2_shuffle), _mm256_castpd_si256(og2_2_shuffle), 0x20);
	    __m256i og3_1_shuffle 	= _mm256_permute2f128_si256(_mm256_castpd_si256(oc2_1_shuffle), _mm256_castpd_si256(og2_1_shuffle), 0x31);
	    __m256i og3_2_shuffle 	= _mm256_permute2f128_si256(_mm256_castpd_si256(oc2_2_shuffle), _mm256_castpd_si256(og2_2_shuffle), 0x31);

	    __m256i od3_1_shuffle 	= _mm256_permute2f128_si256(_mm256_castpd_si256(od2_1_shuffle), _mm256_castpd_si256(oh2_1_shuffle), 0x20);
	    __m256i od3_2_shuffle 	= _mm256_permute2f128_si256(_mm256_castpd_si256(od2_2_shuffle), _mm256_castpd_si256(oh2_2_shuffle), 0x20);
	    __m256i oh3_1_shuffle 	= _mm256_permute2f128_si256(_mm256_castpd_si256(od2_1_shuffle), _mm256_castpd_si256(oh2_1_shuffle), 0x31);
	    __m256i oh3_2_shuffle 	= _mm256_permute2f128_si256(_mm256_castpd_si256(od2_2_shuffle), _mm256_castpd_si256(oh2_2_shuffle), 0x31);

	    __m256i oi3_1_shuffle 	= _mm256_permute2f128_si256(_mm256_castpd_si256(oi2_1_shuffle), _mm256_castpd_si256(om2_1_shuffle), 0x20);
	    __m256i oi3_2_shuffle 	= _mm256_permute2f128_si256(_mm256_castpd_si256(oi2_2_shuffle), _mm256_castpd_si256(om2_2_shuffle), 0x20);
	    __m256i om3_1_shuffle 	= _mm256_permute2f128_si256(_mm256_castpd_si256(oi2_1_shuffle), _mm256_castpd_si256(om2_1_shuffle), 0x31);
	    __m256i om3_2_shuffle 	= _mm256_permute2f128_si256(_mm256_castpd_si256(oi2_2_shuffle), _mm256_castpd_si256(om2_2_shuffle), 0x31);

	    __m256i oj3_1_shuffle 	= _mm256_permute2f128_si256(_mm256_castpd_si256(oj2_1_shuffle), _mm256_castpd_si256(on2_1_shuffle), 0x20);
	    __m256i oj3_2_shuffle 	= _mm256_permute2f128_si256(_mm256_castpd_si256(oj2_2_shuffle), _mm256_castpd_si256(on2_2_shuffle), 0x20);
	    __m256i on3_1_shuffle 	= _mm256_permute2f128_si256(_mm256_castpd_si256(oj2_1_shuffle), _mm256_castpd_si256(on2_1_shuffle), 0x31);
	    __m256i on3_2_shuffle 	= _mm256_permute2f128_si256(_mm256_castpd_si256(oj2_2_shuffle), _mm256_castpd_si256(on2_2_shuffle), 0x31);

	    __m256i ok3_1_shuffle 	= _mm256_permute2f128_si256(_mm256_castpd_si256(ok2_1_shuffle), _mm256_castpd_si256(oo2_1_shuffle), 0x20);
	    __m256i ok3_2_shuffle 	= _mm256_permute2f128_si256(_mm256_castpd_si256(ok2_2_shuffle), _mm256_castpd_si256(oo2_2_shuffle), 0x20);
	    __m256i oo3_1_shuffle 	= _mm256_permute2f128_si256(_mm256_castpd_si256(ok2_1_shuffle), _mm256_castpd_si256(oo2_1_shuffle), 0x31);
	    __m256i oo3_2_shuffle 	= _mm256_permute2f128_si256(_mm256_castpd_si256(ok2_2_shuffle), _mm256_castpd_si256(oo2_2_shuffle), 0x31);

	    __m256i ol3_1_shuffle 	= _mm256_permute2f128_si256(_mm256_castpd_si256(ol2_1_shuffle), _mm256_castpd_si256(op2_1_shuffle), 0x20);
	    __m256i ol3_2_shuffle 	= _mm256_permute2f128_si256(_mm256_castpd_si256(ol2_2_shuffle), _mm256_castpd_si256(op2_2_shuffle), 0x20);
	    __m256i op3_1_shuffle 	= _mm256_permute2f128_si256(_mm256_castpd_si256(ol2_1_shuffle), _mm256_castpd_si256(op2_1_shuffle), 0x31);
	    __m256i op3_2_shuffle 	= _mm256_permute2f128_si256(_mm256_castpd_si256(ol2_2_shuffle), _mm256_castpd_si256(op2_2_shuffle), 0x31);

	    /* Step 4: switch right-top and left-down corner*/
	    __m256i ra4_shuffle = _mm256_permute2f128_si256(ra3_shuffle, ri3_shuffle, 0x20);
	    __m256i ri4_shuffle = _mm256_permute2f128_si256(ra3_shuffle, ri3_shuffle,  0x31);

	    __m256i rb4_shuffle = _mm256_permute2f128_si256(rb3_shuffle, rj3_shuffle, 0x20);
	    __m256i rj4_shuffle = _mm256_permute2f128_si256(rb3_shuffle, rj3_shuffle,  0x31);

	    __m256i rc4_shuffle = _mm256_permute2f128_si256(rc3_shuffle, rk3_shuffle, 0x20);
	    __m256i rk4_shuffle = _mm256_permute2f128_si256(rc3_shuffle, rk3_shuffle,  0x31);

	    __m256i rd4_shuffle = _mm256_permute2f128_si256(rd3_shuffle, rl3_shuffle, 0x20);
	    __m256i rl4_shuffle = _mm256_permute2f128_si256(rd3_shuffle, rl3_shuffle,  0x31);

	    __m256i re4_shuffle = _mm256_permute2f128_si256(re3_shuffle, rm3_shuffle, 0x20);
	    __m256i rm4_shuffle = _mm256_permute2f128_si256(re3_shuffle, rm3_shuffle,  0x31);

	    __m256i rf4_shuffle = _mm256_permute2f128_si256(rf3_shuffle, rn3_shuffle, 0x20);
	    __m256i rn4_shuffle = _mm256_permute2f128_si256(rf3_shuffle, rn3_shuffle,  0x31);

	    __m256i rg4_shuffle = _mm256_permute2f128_si256(rg3_shuffle, ro3_shuffle, 0x20);
	    __m256i ro4_shuffle = _mm256_permute2f128_si256(rg3_shuffle, ro3_shuffle,  0x31);

	    __m256i rh4_shuffle = _mm256_permute2f128_si256(rh3_shuffle, rp3_shuffle, 0x20);
	    __m256i rp4_shuffle = _mm256_permute2f128_si256(rh3_shuffle, rp3_shuffle,  0x31);

	    /* after this, results are in ra4_shuffle, rb4_shuffle, rc4_shuffle, rd4_shuffle, ..., rp4_shuffle*/
	    /* after this, oids are in oa3_1_shuffle, oi3_1_shuffle,
	     	 	 	 	 	 	   ob3_1_shuffle, oj3_1_shuffle,
	     	 	 	 	 	 	   oc3_1_shuffle, ok3_1_shuffle,
	     	 	 	 	 	 	   od3_1_shuffle, ol3_1_shuffle,
	     	 	 	 	 	 	   oe3_1_shuffle, om3_1_shuffle,
	     	 	 	 	 	 	   of3_1_shuffle, on3_1_shuffle,
	     	 	 	 	 	 	   og3_1_shuffle, oo3_1_shuffle,
	     	 	 	 	 	 	   oh3_1_shuffle, op3_1_shuffle,
	     	 	 	 	 	 	   oa3_2_shuffle, oi3_2_shuffle,
	     	 	 	 	 	 	   ob3_2_shuffle, oj3_2_shuffle,
	     	 	 	 	 	 	   oc3_2_shuffle, ok3_2_shuffle,
	     	 	 	 	 	 	   od3_2_shuffle, ol3_2_shuffle,
	     	 	 	 	 	 	   oe3_2_shuffle, om3_2_shuffle,
	     	 	 	 	 	 	   of3_2_shuffle, on3_2_shuffle,
	     	 	 	 	 	 	   og3_2_shuffle, oo3_2_shuffle,
	     	 	 	 	 	 	   oh3_2_shuffle, op3_2_shuffle,
	     	 	 	 	 	 	   */
	/* IACA_END */
	    /* store values */
	    _mm256_store_si256(reinterpret_cast<__m256i *> (out_val), ra4_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *> (out_val + 16), rb4_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *> (out_val + 32), rc4_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *> (out_val + 48), rd4_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *> (out_val + 64), re4_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *> (out_val + 80), rf4_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *> (out_val + 96), rg4_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *> (out_val + 112), rh4_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *> (out_val + 128), ri4_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *> (out_val + 144), rj4_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *> (out_val + 160), rk4_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *> (out_val + 176), rl4_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *> (out_val + 192), rm4_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *> (out_val + 208), rn4_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *> (out_val + 224), ro4_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *> (out_val + 240), rp4_shuffle);
	    /* after this, oids are in oa3_1_shuffle, oi3_1_shuffle,
	     	 	 	 	 	 	   ob3_1_shuffle, oj3_1_shuffle,
	     	 	 	 	 	 	   oc3_1_shuffle, ok3_1_shuffle,
	     	 	 	 	 	 	   od3_1_shuffle, ol3_1_shuffle,
	     	 	 	 	 	 	   oe3_1_shuffle, om3_1_shuffle,
	     	 	 	 	 	 	   of3_1_shuffle, on3_1_shuffle,
	     	 	 	 	 	 	   og3_1_shuffle, oo3_1_shuffle,
	     	 	 	 	 	 	   oh3_1_shuffle, op3_1_shuffle,
	     	 	 	 	 	 	   oa3_2_shuffle, oi3_2_shuffle,
	     	 	 	 	 	 	   ob3_2_shuffle, oj3_2_shuffle,
	     	 	 	 	 	 	   oc3_2_shuffle, ok3_2_shuffle,
	     	 	 	 	 	 	   od3_2_shuffle, ol3_2_shuffle,
	     	 	 	 	 	 	   oe3_2_shuffle, om3_2_shuffle,
	     	 	 	 	 	 	   of3_2_shuffle, on3_2_shuffle,
	     	 	 	 	 	 	   og3_2_shuffle, oo3_2_shuffle,
	     	 	 	 	 	 	   oh3_2_shuffle, op3_2_shuffle,*/

	    _mm256_store_si256(reinterpret_cast<__m256i *>(out_oid), oa3_1_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *>(out_oid + 8), oi3_1_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *>(out_oid + 16), ob3_1_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *>(out_oid + 24), oj3_1_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *>(out_oid + 32), oc3_1_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *>(out_oid + 40), ok3_1_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *>(out_oid + 48), od3_1_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *>(out_oid + 56), ol3_1_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *>(out_oid + 64), oe3_1_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *>(out_oid + 72), om3_1_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *>(out_oid + 80), of3_1_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *>(out_oid + 88), on3_1_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *>(out_oid + 96), og3_1_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *>(out_oid + 104), oo3_1_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *>(out_oid + 112), oh3_1_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *>(out_oid + 120), op3_1_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *>(out_oid + 128), oa3_2_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *>(out_oid + 136), oi3_2_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *>(out_oid + 144), ob3_2_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *>(out_oid + 152), oj3_2_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *>(out_oid + 160), oc3_2_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *>(out_oid + 168), ok3_2_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *>(out_oid + 176), od3_2_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *>(out_oid + 184), ol3_2_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *>(out_oid + 192), oe3_2_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *>(out_oid + 200), om3_2_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *>(out_oid + 208), of3_2_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *>(out_oid + 216), on3_2_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *>(out_oid + 224), og3_2_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *>(out_oid + 232), oo3_2_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *>(out_oid + 240), oh3_2_shuffle);
	    _mm256_store_si256(reinterpret_cast<__m256i *>(out_oid + 248), op3_2_shuffle);
}
#endif

#if 0
inline int __attribute((always_inline))
Mergesort::movemask_epi16_lo(const __m256i &a){

   int ret = 0;
   ret |= _mm256_movemask_ps(reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(a, 0))));
   ret |= _mm256_movemask_ps(reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(a, 1)))) << 8;
   return ret;
}
#endif

inline const int __attribute__((always_inline))
Mergesort::movemask_epi16_lo(const __m256i &a){

   return _mm256_movemask_ps(reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(a, 0))));

}

inline const int __attribute__((always_inline))
Mergesort::movemask_epi16_hi(const __m256i &a){

   return _mm256_movemask_ps(reinterpret_cast<__m256>(_mm256_cvtepi16_epi32(_mm256_extracti128_si256(a, 1))));
}

#if 0	//original one use _mm256_max_pd/_mm256_min_pd
inline void __attribute__((always_inline))
Mergesort::inregister_sort_int64_aligned(int64_t * in_val, int64_t * out_val,
		int32_t * in_oid, int32_t * out_oid)
{
/* IACA_START */
    __m256d ra = _mm256_load_pd ((double const *)(in_val));
    __m256d rb = _mm256_load_pd ((double const *)(in_val + 4));
    __m256d rc = _mm256_load_pd ((double const *)(in_val + 8));
    __m256d rd = _mm256_load_pd ((double const *)(in_val + 12));

    __m256d oa = _mm256_cvtps_pd(_mm_load_ps((float const *)(in_oid)));
    __m256d ob = _mm256_cvtps_pd(_mm_load_ps((float const *)(in_oid + 4)));
    __m256d oc = _mm256_cvtps_pd(_mm_load_ps((float const *)(in_oid + 8)));
    __m256d od = _mm256_cvtps_pd(_mm_load_ps((float const *)(in_oid + 12)));

    /* odd-even sorting network begins */
    /* 1st level of comparisons */
    __m256d ra1 = _mm256_min_pd(ra, rb);
#if 0
    _mm256_store_pd((double *) out_val, oa);
    std::cout << out_val[0] << "\t" << out_val[1]
							<< "\t" << out_val[2]
							<< "\t" << out_val[3] << std::endl;
    std::cout << ((double*)out_val)[0] << "\t" << ((double*)out_val)[1]
							<< "\t" << ((double*)out_val)[2]
							<< "\t" << ((double*)out_val)[3] << std::endl;
#endif
    __m256d rb1 = _mm256_max_pd(ra, rb);

    __m256d rc1 = _mm256_min_pd(rc, rd);
    __m256d rd1 = _mm256_max_pd(rc, rd);

    //re-ordering the oids
    __m256i mask1 = _mm256_cmpeq_epi64(_mm256_castpd_si256(ra), _mm256_castpd_si256(ra1));
    __m256d oa1 = _mm256_blendv_pd(ob, oa, _mm256_castsi256_pd(mask1));
    __m256d ob1 = _mm256_blendv_pd(oa, ob, _mm256_castsi256_pd(mask1));

    __m256i mask2 = _mm256_cmpeq_epi64(_mm256_castpd_si256(rc), _mm256_castpd_si256(rc1));
    __m256d oc1 = _mm256_blendv_pd(od, oc, _mm256_castsi256_pd(mask2));
    __m256d od1 = _mm256_blendv_pd(oc, od, _mm256_castsi256_pd(mask2));

    /* 2nd level of comparisons */
    rb = _mm256_min_pd(rb1, rd1);
    rd = _mm256_max_pd(rb1, rd1);

    mask1 = _mm256_cmpeq_epi64(_mm256_castpd_si256(rb1), _mm256_castpd_si256(rb));
    ob = _mm256_blendv_pd(od1, ob1, _mm256_castsi256_pd(mask1));
    od = _mm256_blendv_pd(ob1, od1, _mm256_castsi256_pd(mask1));

    /* 3rd level of comparisons */
    __m256d ra2 = _mm256_min_pd(ra1, rc1);
    __m256d rc2 = _mm256_max_pd(ra1, rc1);

    mask2 = _mm256_cmpeq_epi64(_mm256_castpd_si256(ra1), _mm256_castpd_si256(ra2));
    __m256d oa2 = _mm256_blendv_pd(oc1, oa1, _mm256_castsi256_pd(mask2));
    __m256d oc2 = _mm256_blendv_pd(oa1, oc1, _mm256_castsi256_pd(mask2));

    /* 4th level of comparisons */
    __m256d rb3 = _mm256_min_pd(rb, rc2);
    __m256d rc3 = _mm256_max_pd(rb, rc2);

    mask1 = _mm256_cmpeq_epi64(_mm256_castpd_si256(rb), _mm256_castpd_si256(rb3));
    __m256d ob3 = _mm256_blendv_pd(oc2, ob, _mm256_castsi256_pd(mask1));
    __m256d oc3 = _mm256_blendv_pd(ob, oc2, _mm256_castsi256_pd(mask1));

    /* results are in ra2, rb3, rc3, rd */
    /* re-ordered oids are in oa2, ob3, oc3, od */
    /**
     * Initial data and transposed data looks like following:
     *  a2={ x1  x2  x3  x4  }                      a4={ x1 x5 x9  x13 }
     *  b3={ x5  x6  x7  x8  }  === Transpose ===>  b5={ x2 x6 x10 x14 }
     *  c3={ x9  x10 x11 x12 }                      c5={ x3 x7 x11 x15 }
     *  d={ x13 x14 x15 x16 }                       d4={ x4 x8 x12 x16 }
     */
    /* shuffle x2 and x5 - shuffle x4 and x7 */
    __m256d ra3 = _mm256_unpacklo_pd(ra2, rb3);
    __m256d rb4 = _mm256_unpackhi_pd(ra2, rb3);

    __m256d oa3 = _mm256_unpacklo_pd(oa2, ob3);
    __m256d ob4 = _mm256_unpackhi_pd(oa2, ob3);

    /* shuffle x10 and x13 - shuffle x12 and x15 */
    __m256d rc4 = _mm256_unpacklo_pd(rc3, rd);
    __m256d rd3 = _mm256_unpackhi_pd(rc3, rd);

    __m256d oc4 = _mm256_unpacklo_pd(oc3, od);
    __m256d od3 = _mm256_unpackhi_pd(oc3, od);

    /* shuffle (x3,x7) and (x9,x13) pairs */
    __m256d ra4 = _mm256_permute2f128_pd(ra3, rc4, 0x20);
    __m256d rc5 = _mm256_permute2f128_pd(ra3, rc4, 0x31);

    __m256d oa4 = _mm256_permute2f128_pd(oa3, oc4, 0x20);
    __m256d oc5 = _mm256_permute2f128_pd(oa3, oc4, 0x31);

    /* shuffle (x4,x8) and (x10,x14) pairs */
    __m256d rb5 = _mm256_permute2f128_pd(rb4, rd3, 0x20);
    __m256d rd4 = _mm256_permute2f128_pd(rb4, rd3, 0x31);

    __m256d ob5 = _mm256_permute2f128_pd(ob4, od3, 0x20);
    __m256d od4 = _mm256_permute2f128_pd(ob4, od3, 0x31);

    /* after this, results are in ra4, rb5, rc5, rd4 */
    /* oids are in oa4, ob5, oc5, od4 */
/* IACA_END */
    /* store values */
    _mm256_store_pd((double *) out_val, ra4);
    _mm256_store_pd((double *) (out_val + 4), rb5);
    _mm256_store_pd((double *) (out_val + 8), rc5);
    _mm256_store_pd((double *) (out_val + 12), rd4);

    //convert the oids from 64 to 32
    __m128 oa4_s = _mm256_cvtpd_ps(oa4);
    __m128 ob5_s = _mm256_cvtpd_ps(ob5);
    __m128 oc5_s = _mm256_cvtpd_ps(oc5);
    __m128 od4_s = _mm256_cvtpd_ps(od4);

    __m128i mask_store = _mm_set1_epi32(1<<31);
    _mm_maskstore_ps((float *) out_oid, mask_store, oa4_s);
    _mm_maskstore_ps((float *) (out_oid + 4), mask_store, ob5_s);
    _mm_maskstore_ps((float *) (out_oid + 8), mask_store, oc5_s);
    _mm_maskstore_ps((float *) (out_oid + 12), mask_store, od4_s);
}
#endif

#if 1	//use _mm256_cmpgt_epi256 to replace _mm256_max_pd/_mm256_min_pd
inline void __attribute__((always_inline))
Mergesort::inregister_sort_int64_aligned(int64_t * in_val, int64_t * out_val,
		int32_t * in_oid, int32_t * out_oid)
{
/* IACA_START */
    __m256i ra = _mm256_load_si256 ((__m256i const *)(in_val));
    __m256i rb = _mm256_load_si256 ((__m256i const *)(in_val + 4));
    __m256i rc = _mm256_load_si256 ((__m256i const *)(in_val + 8));
    __m256i rd = _mm256_load_si256 ((__m256i const *)(in_val + 12));

    __m256d oa = _mm256_cvtps_pd(_mm_load_ps((float const *)(in_oid)));
    __m256d ob = _mm256_cvtps_pd(_mm_load_ps((float const *)(in_oid + 4)));
    __m256d oc = _mm256_cvtps_pd(_mm_load_ps((float const *)(in_oid + 8)));
    __m256d od = _mm256_cvtps_pd(_mm_load_ps((float const *)(in_oid + 12)));

    /* odd-even sorting network begins */
    /* 1st level of comparisons */
    __m256i mask1 = _mm256_cmpgt_epi64(ra, rb);
    __m256i ra1 = reinterpret_cast<__m256i>(_mm256_blendv_pd(reinterpret_cast<__m256d>(ra),
    									reinterpret_cast<__m256d>(rb),
										reinterpret_cast<__m256d>(mask1)));
#if 0
    _mm256_store_pd((double *) out_val, oa);
    std::cout << out_val[0] << "\t" << out_val[1]
							<< "\t" << out_val[2]
							<< "\t" << out_val[3] << std::endl;
    std::cout << ((double*)out_val)[0] << "\t" << ((double*)out_val)[1]
							<< "\t" << ((double*)out_val)[2]
							<< "\t" << ((double*)out_val)[3] << std::endl;
#endif
    __m256i rb1 = reinterpret_cast<__m256i>(_mm256_blendv_pd(reinterpret_cast<__m256d>(rb),
										reinterpret_cast<__m256d>(ra),
										reinterpret_cast<__m256d>(mask1)));

    __m256i mask2 = _mm256_cmpgt_epi64(rc, rd);
    __m256i rc1 = reinterpret_cast<__m256i>(_mm256_blendv_pd(reinterpret_cast<__m256d>(rc),
										reinterpret_cast<__m256d>(rd),
										reinterpret_cast<__m256d>(mask2)));
    __m256i rd1 = reinterpret_cast<__m256i>(_mm256_blendv_pd(reinterpret_cast<__m256d>(rd),
			reinterpret_cast<__m256d>(rc),
			reinterpret_cast<__m256d>(mask2)));

    //re-ordering the oids
    __m256d oa1 = _mm256_blendv_pd(oa, ob, reinterpret_cast<__m256d>(mask1));
    __m256d ob1 = _mm256_blendv_pd(ob, oa, reinterpret_cast<__m256d>(mask1));

    __m256d oc1 = _mm256_blendv_pd(oc, od, reinterpret_cast<__m256d>(mask2));
    __m256d od1 = _mm256_blendv_pd(od, oc, reinterpret_cast<__m256d>(mask2));

    /* 2nd level of comparisons */
    __m256i mask3 = _mm256_cmpgt_epi64(rb1, rd1);
    rb = reinterpret_cast<__m256i>(_mm256_blendv_pd(reinterpret_cast<__m256d>(rb1),
			reinterpret_cast<__m256d>(rd1),
			reinterpret_cast<__m256d>(mask3)));
    rd = reinterpret_cast<__m256i>(_mm256_blendv_pd(reinterpret_cast<__m256d>(rd1),
			reinterpret_cast<__m256d>(rb1),
			reinterpret_cast<__m256d>(mask3)));

    ob = _mm256_blendv_pd(ob1, od1, reinterpret_cast<__m256d>(mask3));
    od = _mm256_blendv_pd(od1, ob1, reinterpret_cast<__m256d>(mask3));

    /* 3rd level of comparisons */
    __m256i mask4 = _mm256_cmpgt_epi64(ra1, rc1);
    __m256i ra2 = reinterpret_cast<__m256i>(_mm256_blendv_pd(reinterpret_cast<__m256d>(ra1),
			reinterpret_cast<__m256d>(rc1),
			reinterpret_cast<__m256d>(mask4)));
    __m256i rc2 = reinterpret_cast<__m256i>(_mm256_blendv_pd(reinterpret_cast<__m256d>(rc1),
			reinterpret_cast<__m256d>(ra1),
			reinterpret_cast<__m256d>(mask4)));

    __m256d oa2 = _mm256_blendv_pd(oa1, oc1, reinterpret_cast<__m256d>(mask4));
    __m256d oc2 = _mm256_blendv_pd(oc1, oa1, reinterpret_cast<__m256d>(mask4));

    /* 4th level of comparisons */
    __m256i mask5 = _mm256_cmpgt_epi64(rb, rc2);
    __m256i rb3 = reinterpret_cast<__m256i>(_mm256_blendv_pd(reinterpret_cast<__m256d>(rb),
			reinterpret_cast<__m256d>(rc2),
			reinterpret_cast<__m256d>(mask5)));
    __m256i rc3 = reinterpret_cast<__m256i>(_mm256_blendv_pd(reinterpret_cast<__m256d>(rc2),
			reinterpret_cast<__m256d>(rb),
			reinterpret_cast<__m256d>(mask5)));

    __m256d ob3 = _mm256_blendv_pd(ob, oc2, reinterpret_cast<__m256d>(mask5));
    __m256d oc3 = _mm256_blendv_pd(oc2, ob, reinterpret_cast<__m256d>(mask5));

#if 0
    _mm256_store_si256((__m256i *) out_val, oc3);
    std::cout << "intermediate result: " << std::endl;
    std::cout << out_val[0] << "\t" << out_val[1]
							<< "\t" << out_val[2]
							<< "\t" << out_val[3] << std::endl;

#endif
    /* results are in ra2, rb3, rc3, rd */
    /* re-ordered oids are in oa2, ob3, oc3, od */
    /**
     * Initial data and transposed data looks like following:
     *  a2={ x1  x2  x3  x4  }                      a4={ x1 x5 x9  x13 }
     *  b3={ x5  x6  x7  x8  }  === Transpose ===>  b5={ x2 x6 x10 x14 }
     *  c3={ x9  x10 x11 x12 }                      c5={ x3 x7 x11 x15 }
     *  d={ x13 x14 x15 x16 }                       d4={ x4 x8 x12 x16 }
     */
    /* shuffle x2 and x5 - shuffle x4 and x7 */
    __m256i ra3 = _mm256_unpacklo_epi64(ra2, rb3);
    __m256i rb4 = _mm256_unpackhi_epi64(ra2, rb3);

    __m256d oa3 = _mm256_unpacklo_pd(oa2, ob3);
    __m256d ob4 = _mm256_unpackhi_pd(oa2, ob3);

    /* shuffle x10 and x13 - shuffle x12 and x15 */
    __m256i rc4 = _mm256_unpacklo_epi64(rc3, rd);
    __m256i rd3 = _mm256_unpackhi_epi64(rc3, rd);

    __m256d oc4 = _mm256_unpacklo_pd(oc3, od);
    __m256d od3 = _mm256_unpackhi_pd(oc3, od);

    /* shuffle (x3,x7) and (x9,x13) pairs */
    __m256i ra4 = _mm256_permute2f128_si256(ra3, rc4, 0x20);
    __m256i rc5 = _mm256_permute2f128_si256(ra3, rc4, 0x31);

    __m256d oa4 = _mm256_permute2f128_pd(oa3, oc4, 0x20);
    __m256d oc5 = _mm256_permute2f128_pd(oa3, oc4, 0x31);

    /* shuffle (x4,x8) and (x10,x14) pairs */
    __m256i rb5 = _mm256_permute2f128_si256(rb4, rd3, 0x20);
    __m256i rd4 = _mm256_permute2f128_si256(rb4, rd3, 0x31);

    __m256d ob5 = _mm256_permute2f128_pd(ob4, od3, 0x20);
    __m256d od4 = _mm256_permute2f128_pd(ob4, od3, 0x31);

#if 0
    _mm256_store_si256((__m256i *) out_val, od4);
    std::cout << "intermediate result: " << std::endl;
    std::cout << out_val[0] << "\t" << out_val[1]
							<< "\t" << out_val[2]
							<< "\t" << out_val[3] << std::endl;

#endif

    /* after this, results are in ra4, rb5, rc5, rd4 */
    /* oids are in oa4, ob5, oc5, od4 */
/* IACA_END */
    /* store values */
    _mm256_store_si256((__m256i *) out_val, ra4);
    _mm256_store_si256((__m256i *) (out_val + 4), rb5);
    _mm256_store_si256((__m256i *) (out_val + 8), rc5);
    _mm256_store_si256((__m256i *) (out_val + 12), rd4);

    //convert the oids from 64 to 32
    __m128 oa4_s = _mm256_cvtpd_ps(oa4);
    __m128 ob5_s = _mm256_cvtpd_ps(ob5);
    __m128 oc5_s = _mm256_cvtpd_ps(oc5);
    __m128 od4_s = _mm256_cvtpd_ps(od4);

#if 0
    _mm256_store_si256((__m256i *) out_val, od4);
    std::cout << "intermediate result: " << std::endl;
    std::cout << out_val[0] << "\t" << out_val[1]
							<< "\t" << out_val[2]
							<< "\t" << out_val[3] << std::endl;

#endif

    __m128i mask_store = _mm_set1_epi32(1<<31);
    _mm_maskstore_ps((float *) out_oid, mask_store, oa4_s);
    _mm_maskstore_ps((float *) (out_oid + 4), mask_store, ob5_s);
    _mm_maskstore_ps((float *) (out_oid + 8), mask_store, oc5_s);
    _mm_maskstore_ps((float *) (out_oid + 12), mask_store, od4_s);
}
#endif

template <class T>
void* Mergesort::mergesort_thread(void* param) {
	mergesort_params_t<T>* setting 	= (mergesort_params_t<T>*) param;
	//int32_t my_tid 		= setting->my_tid;
	int rv;

	//DEBUGMSG(1, "Thread-%d started running ...\n", my_tid);

	BARRIER_ARRIVE(setting->barrier, rv);


	/*********************************************************************
	  *
	  *   Phase 1: Partitioning
	  *
	  ******************************************************************/

	BAT_t<T>**	partitions = NULL;
	partitioning_phase<T>(&partitions, setting);


	BARRIER_ARRIVE(setting->barrier, rv);


	/*********************************************************************
	  *
	  *   Phase 2: sorting of individual partitions
	  *
	  ******************************************************************/
	sorting_phase<T>(partitions, setting);


	BARRIER_ARRIVE(setting->barrier, rv);

	/*********************************************************************
	  *
	  *   Phase 3: apply multi-way merging with in-cache resident buffers
	  *
	  ******************************************************************/
	//BAT mergedCol;
	/**
	 * if nthread==1, the global merged result is just in <threadcolchunks>,
	 * otherwise, the merged results for current thread is in <result> */

	//multiwaymerge_phase<T>(setting);


	//BARRIER_ARRIVE(setting->barrier, rv);

	//clean up
	free(partitions[0]);	//i.e., memory initialized by <pcols> in partitioning_phase
	free(partitions);

	return NULL;
}


template <class T>
void Mergesort::partitioning_phase(BAT_t<T>*** partitions, mergesort_params_t<T>* setting) {

	//TODO:may need to revise this when T=uint8_t
	if (1 == sizeof(T)) {
		//setting->partition_fanout = 256;
		printf("[Error ] Not supporting T=uint8_t yet...");
		exit(EXIT_SUCCESS);
	}

	uint32_t fanout 	= setting->partition_fanout;

	uint32_t nradixbits = log2(fanout);
	uint32_t bitwidth = setting->bitwidth;

	if (bitwidth < nradixbits) {
		nradixbits = bitwidth;
#ifdef DEBUG
		printf("change the nradixbits from original %d to the bitwdith %d\n",
				(int)(log2(fanout)), nradixbits);
#endif
	}

	/**  **/
	BAT_t<T>** colPartitions = (BAT_t<T>**) malloc_aligned(fanout * sizeof(BAT_t<T>*));

	BAT_t<T>* pcols = (BAT_t<T>*) malloc_aligned(fanout * sizeof(BAT_t<T>));
	uint32_t i;
	for (i = 0; i < fanout; ++i) {
		colPartitions[i] = pcols + i;
	}

	BAT_t<T> columnchunk;	//target column to be partitioned
	BAT_t<T> tmpcol;		//output of radix partitioning

	columnchunk.values	 		= setting->columnchunk;
	columnchunk.oids			= setting->oidchunk;
	columnchunk.num_elements 	= setting->chunksize;
   	tmpcol.values				= setting->tmp_partition;
   	tmpcol.oids					= setting->tmp_partition_oid;
	tmpcol.num_elements			= setting->chunksize;

	uint32_t bitshift = bitwidth - nradixbits;	//# of bit needs to shift in order to get radixbits

	/**
	  * initialize the histogram, aligned to the CACHE_LINE_SIZE
	  */
	uint64_t* hist, * histAligned;
	hist = (uint64_t*) calloc(fanout + CACHE_LINE_SIZE/sizeof(uint64_t), sizeof(uint64_t));
	histAligned = (uint64_t*) ALIGNPTR(hist, CACHE_LINE_SIZE);

	radix_partitioning<T>(&tmpcol, &columnchunk, histAligned, bitshift, nradixbits);

	uint64_t offset_values 	= 0;	//for the alignment of each partition
	uint64_t offset_oids	= 0;
	for (i = 0; i < fanout; ++i) {
		BAT_t<T>* part 		= colPartitions[i];
		part->num_elements 	= histAligned[i];
		part->values 		= tmpcol.values + offset_values;
		part->oids 			= tmpcol.oids + offset_oids;

		offset_values 	+= ALIGN_NUM_VALUES(T, histAligned[i]);
		offset_oids		+= ALIGN_NUM_VALUES(surrogate_t, histAligned[i]);

#if 0
		printf("# of elements in partition %d: %lu\n", i, part->num_elements);
		//debug: output the values in each partition
		printf("values in partition %d:\n", i);
		for  (uint32_t tmpId = 0; tmpId < part->num_elements; ++tmpId) {
			printf("%u\t", part->elements[tmpId].value);
		}
#endif
#if 0
		for(uint32_t j = 0; j < part->num_elements; ++j) {
			assert(part->oids[j] < 1638400);
		}
#endif
	}

	/** return partitions **/
	*partitions = colPartitions;	//here we cannot free colPartitions, why??

	free(hist);
	//free(colPartitions);
}

template <class T>
void Mergesort::radix_partitioning(BAT_t<T>* restrict outCol, BAT_t<T>* restrict inCol,
		uint64_t* restrict hist, uint32_t bitshift, uint32_t nradixbits) {

	const uint32_t fanout = 1 << nradixbits;
	const uint32_t num_elements = inCol->num_elements;
	uint32_t i;

	T* input_val 			= inCol->values;
	surrogate_t* input_oid	= inCol->oids;
	T* output_val		 	= outCol->values;
	surrogate_t* output_oid = outCol->oids;

	uint64_t			dst_val[fanout]		__attribute__((aligned(CACHE_LINE_SIZE)));
	uint64_t			dst_oid[fanout]		__attribute__((aligned(CACHE_LINE_SIZE)));
	cacheline_val_t<T>	buffer_val[fanout]	__attribute__((aligned(CACHE_LINE_SIZE)));
	cacheline_oid_t		buffer_oid[fanout]	__attribute__((aligned(CACHE_LINE_SIZE)));

	uint32_t idx = 0;
	const uint32_t mask = ((1 << nradixbits) - 1) << bitshift;
	/** count elements per partition **/
	for (i = 0; i < num_elements; ++i) {
		idx = RADIX_HASH_FUNCTION(*input_val, mask, bitshift);
		hist[idx]++;
		input_val++;
	}

	for (i = 0; i < fanout; ++i) {
		buffer_val[i].data.slot = 0;
	}

	uint64_t offset_val = 0;
	uint64_t offset_oid = 0;
	/** determine the start and end of each partition depending on the counts **/
	for (i = 0; i < fanout; ++i) {
		dst_val[i] = offset_val;
		dst_oid[i] = offset_oid;
		/* for aligning partition-outputs to cacheline
		 * (related to last portion of partitioning_phase())  */
		offset_val += ALIGN_NUM_VALUES(T, hist[i]);
		offset_oid += ALIGN_NUM_VALUES(surrogate_t, hist[i]);
	}

	input_val = inCol->values;
	//input_oid = inCol->oids;
	/** copy elements to their corresponding partitions at appropriate offset **/
	uint32_t slot = 0;
	T* bucket_val = NULL;
	surrogate_t* bucket_oid = NULL;
	for (i = 0; i < num_elements; ++i) {
		idx = RADIX_HASH_FUNCTION(*input_val, mask, bitshift);

		/* store in the cache-resident buffer first*/
		slot = buffer_val[idx].data.slot;
		bucket_val = (T*)(buffer_val + idx);
		bucket_oid = (surrogate_t*)(buffer_oid + idx);
		bucket_val[slot] = *input_val;
		bucket_oid[slot] = *input_oid;
		input_val++;
		input_oid++;
		slot++;

		/* flush the bucket when full */
		if (slot == VALUESSPERCACHELINE(surrogate_t)) {
			//store_nontemp_64B((output + dst[idx]), (buffer + idx));
			//memcpy(output + dst[idx], buffer + idx, CACHE_LINE_SIZE);
			*(cacheline_val_t<T> *)(output_val + dst_val[idx]) = *(cacheline_val_t<T> *)(buffer_val + idx);
			*(cacheline_oid_t *)(output_oid + dst_oid[idx]) = *(cacheline_oid_t *)(buffer_oid + idx);
			slot = 0;
			dst_val[idx] += VALUESSPERCACHELINE(surrogate_t);
			dst_oid[idx] += VALUESSPERCACHELINE(surrogate_t);
		}
		buffer_val[idx].data.slot = slot;
	}

	/* flush the remainder elements in the buffer */
	uint32_t num;
	uint32_t j;
	T* destination_val;
	surrogate_t* destination_oid;
	for (i = 0; i < fanout; ++i) {
		num = buffer_val[i].data.slot;
		if (num > 0) {
			destination_val = output_val + dst_val[i];
			destination_oid = output_oid + dst_oid[i];
			for (j = 0; j < num; ++j) {
				destination_val[j] = buffer_val[i].data.values[j];
				destination_oid[j] = buffer_oid[i].oids[j];
			}
		}
	}
}

template <class T>
inline void __attribute__((always_inline))
Mergesort::flip(T** inptr_val, uint64_t nitems, uint32_t bitwidth) {
	T* values = *inptr_val;
	for (uint64_t idx = 0; idx < nitems; ++idx) {
		values[idx] ^= (1 << (bitwidth-1));	//can use SIMD to speed up???
	}
}

template <class T>
void Mergesort::sorting_phase(BAT_t<T>** partitions, mergesort_params_t<T>* setting) {
	const int fanout = setting->partition_fanout;
	int32_t my_tid = setting->my_tid;
	const IntrinsicsType intrin_type = setting->intrinType;
	int i;
	uint32_t bitwidth = setting->bitwidth;

	setting->threadcolchunks[my_tid] = (BAT_t<T>*) malloc_aligned(fanout * sizeof(BAT_t<T>));

	uint64_t nelements_per_part;
	uint64_t offset_val = 0;
	uint64_t offset_oid = 0;
	T* start_location_val 			= setting->tmp_sort;
	surrogate_t* start_location_oid = setting->tmp_sort_oid;

	T* inptr_val = NULL;			//start location of each partition for values
	surrogate_t* inptr_oid = NULL;	//start location of each partition for oids
	T* outptr_val = NULL;
	surrogate_t* outptr_oid = NULL;
	for (i = 0; i < fanout; ++i) {
		inptr_val 	= partitions[i]->values;
		inptr_oid	= partitions[i]->oids;
		outptr_val 	= start_location_val + offset_val;
		outptr_oid	= start_location_oid + offset_oid;
		nelements_per_part 	= partitions[i]->num_elements;
		offset_val			+= ALIGN_NUM_VALUES(T, nelements_per_part);
		offset_oid			+= ALIGN_NUM_VALUES(surrogate_t, nelements_per_part);

#if 0
		//debug: output the values in each partition （before sorting）
		printf("values in partition %d before sorting (%lu)...(output first 100 and last 100)\n", i, nelements_per_part);
		uint32_t tmpId;
		for  (tmpId = 0; tmpId < 100; ++tmpId) {
			//printf("%u\n", inptr[tmpId].value);
			std::cout << inptr_val[tmpId] << "[" << inptr_oid[tmpId] << "]" << std::endl;
			//printf("%lx\n", ((int64_t *)inptr_val)[tmpId]);
			//printf("value: %x \t oid: %x\n", inptr[tmpId].value, inptr[tmpId].oid);
		}
		for  (tmpId = nelements_per_part-100; tmpId < nelements_per_part; ++tmpId) {
			//printf("%u\n", inptr[tmpId].value);
			std::cout << inptr_val[tmpId] << "[" << inptr_oid[tmpId] << "]" << std::endl;
			//printf("%lx\n", ((int64_t *)inptr_val)[tmpId]);
			//printf("value: %x \t oid: %x\n", inptr[tmpId].value, inptr[tmpId].oid);
		}
		printf("\n");
#endif

		//TODO: may need revise: when T = uint8_t, just copy the content from input to output because the elements are already sorted
		if (1 == sizeof(T)) {
			memcpy(outptr_val, inptr_val, nelements_per_part*sizeof(T));
			memcpy(outptr_oid, inptr_oid, nelements_per_part*sizeof(surrogate_t));
		} else {	//do the sorting for each partition

			//flip the values if the bitwidth is 8/16/32/64 because SIMD-bank uses signed integer
			if ((bitwidth % 8) == 0) {
				flip<T>(&inptr_val, nelements_per_part, bitwidth);
			}

			switch(intrin_type) {
				case IntrinsicsType::AVX:
					avxmergesort_anytype<T>(&inptr_val, &inptr_oid, &outptr_val, &outptr_oid, nelements_per_part);
					break;
				case IntrinsicsType::SSE:
					break;
				case IntrinsicsType::SCALAR:
					break;
			}
		}

		//store the sorted partitions into the global settings
		setting->threadcolchunks[my_tid][i].values 		= outptr_val;
		setting->threadcolchunks[my_tid][i].oids		= outptr_oid;
		setting->threadcolchunks[my_tid][i].num_elements 	= nelements_per_part;
#if 0
		//debug: output the values in each partition （after sorting）
		printf("values in partition %d after sorting:\n", i);
		for  (uint32_t tmpId = 0; tmpId < nelements_per_part; ++tmpId) {
			printf("V--%u\t\%u\n", outptr[tmpId].value, inptr[tmpId].value);
			printf("O--%u\t\%u\n", outptr[tmpId].oid, inptr[tmpId].oid);
			printf("H--%lx\t\%lx\n", ((int64_t *)outptr)[tmpId], ((int64_t *)inptr)[tmpId]);
			printf("---%lu\t\%lu\n", ((int64_t *)outptr)[tmpId], ((int64_t *)inptr)[tmpId]);
		}
		printf("\n");
#endif
	}
}

template <class T>
void Mergesort::avxmergesort_anytype(T** inputptr_val, surrogate_t** inputptr_oid,
		T** outputptr_val, surrogate_t** outputptr_oid, uint64_t nitems) {

	if (nitems <= 0) return;

	T* input_val 			= *inputptr_val;
	surrogate_t* input_oid	= *inputptr_oid;

	T* output_val 			= *outputptr_val;
	surrogate_t* output_oid = *outputptr_oid;


#if 0
		//debug: output the values in each partition （before sorting）
		std::cout << "values avxmergesort_anytype to be sorted...(output first 100 and last 100), total " << nitems << std::endl;
		uint32_t tmpId;
		for  (tmpId = 0; tmpId < 100; ++tmpId) {
			std::cout << input_val[tmpId] << "[" << input_oid[tmpId] << "]" << std::endl;
		}
		for  (tmpId = nitems-100; tmpId < nitems; ++tmpId) {
			std::cout << input_val[tmpId] << "[" << input_oid[tmpId] << "]" << std::endl;
		}
		printf("\n");
#endif


	uint64_t i;
	uint64_t nblocks = nitems / BLOCKSIZE<T>();
	int rem = nitems % BLOCKSIZE<T>();			//comment temporary
	//assert(rem == 0);

	/** each block keeps track of its temporary memory offset **/
	T* ptrs_val[nblocks+1][2]; //[chunk-in, chunk-out-tmp]
	surrogate_t* ptrs_oid[nblocks+1][2];
	uint32_t sizes[nblocks+1];

	/** Divide the input into blocks fitting into L2 cache. **/
	for (i = 0; i <= nblocks; ++i) {
		ptrs_val[i][0] 	= input_val + i * BLOCKSIZE<T>();
		ptrs_val[i][1] 	= output_val + i * BLOCKSIZE<T>();
		ptrs_oid[i][0] 	= input_oid + i * BLOCKSIZE<T>();
		ptrs_oid[i][1]	= output_oid + i * BLOCKSIZE<T>();
		sizes[i]		= BLOCKSIZE<T>();
	}

	/* sort each block */
	switch (sizeof(T)) {
		case 1://uint8_t, cannot enter here because skip sorting for the type
			{
				printf("[ERROR ] Do not sort values with type uint8_t\n");
				exit(EXIT_SUCCESS);
				break;
			}
		case 2://uint16_t
			for (i = 0; i < nblocks; ++i) {
				avxmergesort_block_uint16_aligned((uint16_t **)&ptrs_val[i][0], &ptrs_oid[i][0],
						(uint16_t **)&ptrs_val[i][1], &ptrs_oid[i][1], BLOCKSIZE<T>());
				swap<T>(&ptrs_val[i][0], &ptrs_val[i][1]);
				swap<surrogate_t>(&ptrs_oid[i][0], &ptrs_oid[i][1]);
			}
			break;
		case 4://uint32_t
			for (i = 0; i < nblocks; ++i) {
				avxmergesort_block_uint32_aligned((uint32_t **)&ptrs_val[i][0], &ptrs_oid[i][0],
						(uint32_t **)&ptrs_val[i][1], &ptrs_oid[i][1], BLOCKSIZE<T>());
#if 0
				for (uint32_t idx = 0; idx < BLOCKSIZE<T>(); ++idx) {
					assert(ptrs_oid[i][1][idx] < 1638400);
				}
#endif
				swap<T>(&ptrs_val[i][0], &ptrs_val[i][1]);
				swap<surrogate_t>(&ptrs_oid[i][0], &ptrs_oid[i][1]);

			}
			break;
		case 8:	//uint64_t
			for (i = 0; i < nblocks; ++i) {
				avxmergesort_block_uint64_aligned((uint64_t **)&ptrs_val[i][0], &ptrs_oid[i][0],
						(uint64_t **)&ptrs_val[i][1], &ptrs_oid[i][1], BLOCKSIZE<T>());
				swap<T>(&ptrs_val[i][0], &ptrs_val[i][1]);
				swap<surrogate_t>(&ptrs_oid[i][0], &ptrs_oid[i][1]);
			}
			break;
	}

#if 0
		//debug: output the values after sorting blocks
		std::cout << "output the values after sorting blocks: " << std::endl;
		uint32_t tmpId;
		std::cout << "ptrs_val[0][0]" << std::endl;
		for  (tmpId = 0; tmpId < 100; ++tmpId) {
			std::cout << ptrs_val[0][0][tmpId] << std::endl;
		}
		std::cout << "ptrs_val[0][1]" << std::endl;
		for  (tmpId = nitems-100; tmpId < nitems; ++tmpId) {
			std::cout << ptrs_val[0][1][tmpId] << std::endl;
		}
		printf("\n");
#endif

#if 1
	if (rem) {	//rem != 0 , sort the last block which is less than BLOCKSIZE
		switch (sizeof(T)) {
			case 1://uint8_t, cannot enter here because skip sorting for the type
				{
					printf("[ERROR ] Do not sort values with type uint8_t\n");
					exit(EXIT_SUCCESS);
					break;
				}
			case 2://uint16_t
				{
					avxmergesort_rem_uint16_aligned((uint16_t **)&ptrs_val[i][0], &ptrs_oid[i][0],
							(uint16_t **)&ptrs_val[i][1], &ptrs_oid[i][1], rem);
					swap<T>(&ptrs_val[i][0], &ptrs_val[i][1]);
					swap<surrogate_t>(&ptrs_oid[i][0], &ptrs_oid[i][1]);
				}
				break;
			case 4://uint32_t
				{
					avxmergesort_rem_uint32_aligned((uint32_t **)&ptrs_val[i][0], &ptrs_oid[i][0],
							(uint32_t **)&ptrs_val[i][1], &ptrs_oid[i][1], rem);
					swap<T>(&ptrs_val[i][0], &ptrs_val[i][1]);
					swap<surrogate_t>(&ptrs_oid[i][0], &ptrs_oid[i][1]);
				}
				break;
			case 8:	//uint64_t
				{
					avxmergesort_rem_uint64_aligned((uint64_t **)&ptrs_val[i][0], &ptrs_oid[i][0],
							(uint64_t **)&ptrs_val[i][1], &ptrs_oid[i][1], rem);
					swap<T>(&ptrs_val[i][0], &ptrs_val[i][1]);
					swap<surrogate_t>(&ptrs_oid[i][0], &ptrs_oid[i][1]);
				}
				break;
		}
		sizes[i] = rem;
	}
#endif

	/**
	  * 2-way merge: TODO: may modify to multi-way merge to alleviate the memory access bottleneck
	  */
	nblocks += (rem > 0);

	const uint64_t logN = ceil(log2(nitems));
	uint64_t k, j;
	T* inpA_val;
	T* inpB_val;
	T* out_val;
	surrogate_t* inpA_oid;
	surrogate_t* inpB_oid;
	surrogate_t* out_oid;
	uint32_t sizeA;
	uint32_t sizeB;
	for (i = log2(BLOCKSIZE<T>()); i < logN; ++i) {

		k = 0;
		for (j = 0; j < (nblocks-1); j += 2) {
			inpA_val 	= ptrs_val[j][0];
			inpB_val 	= ptrs_val[j+1][0];
			out_val		= ptrs_val[j][1];

			inpA_oid	= ptrs_oid[j][0];
			inpB_oid	= ptrs_oid[j+1][0];
			out_oid		= ptrs_oid[j][1];

			sizeA 	= sizes[j];
			sizeB	= sizes[j+1];

			//use switch-case: TODO: can the compiler remove the branches during compiling?
			switch (sizeof(T)) {
				case 1:
				{
					printf("[ERROR ] Do not merge values with type uint8_t\n");
					exit(EXIT_SUCCESS);
					break;
				}
				case 2:
					merge16_int16_varlen_aligned((int16_t *)inpA_val, (int16_t *)inpB_val, (int16_t *)out_val,
							inpA_oid, inpB_oid, out_oid, sizeA, sizeB);
					break;
				case 4:
					merge8_int32_varlen_aligned((int32_t *)inpA_val, (int32_t *)inpB_val, (int32_t *)out_val,
							inpA_oid, inpB_oid, out_oid, sizeA, sizeB);
					break;
				case 8:
					merge8_int64_varlen_aligned((int64_t *)inpA_val, (int64_t *)inpB_val, (int64_t *)out_val,
							inpA_oid, inpB_oid, out_oid, sizeA, sizeB);	//must be various length because nblock can be *odd*
					break;
			}

			/* setup new pointers */
			ptrs_val[k][0] 	= out_val;
			ptrs_val[k][1] 	= inpA_val;
			ptrs_oid[k][0]	= out_oid;
			ptrs_oid[k][1]	= inpA_oid;
			sizes[k]	= sizeA + sizeB;
			k++;
		}

		//odd number of blocks
		if (nblocks % 2) {
			ptrs_val[k][0]	= ptrs_val[nblocks-1][0];
			ptrs_val[k][1]	= ptrs_val[nblocks-1][1];
			ptrs_oid[k][0]	= ptrs_oid[nblocks-1][0];
			ptrs_oid[k][1]	= ptrs_oid[nblocks-1][1];
			sizes[k]	= sizes[nblocks-1];
			k++;
		}

		nblocks = k;
	}

	/** finally swap input/output pointers, where output holds the sorted list **/
	*outputptr_val 	= ptrs_val[0][0];
	*inputptr_val 	= ptrs_val[0][1];
	*outputptr_oid	= ptrs_oid[0][0];
	*inputptr_oid	= ptrs_oid[0][1];
}

#if 1
/**
 * Merge two sorted arrays (of int64_t) to a final output using 8-way AVX bitonic merge.
 * Main principle: convert ONE 8*32 register to TWO 4*64 register for consistent issue
 *
 * @param inpA input array A
 * @param inpB input array B
 * @param Out  output array
 * @param lenA size of A
 * @param lenB size of B
 */
inline void __attribute__((always_inline))
Mergesort::merge8_int64_varlen_aligned(int64_t * restrict inpA_val,
                      			int64_t * restrict inpB_val,
                      			int64_t * restrict Out_val,
								surrogate_t * restrict inpA_oid,
								surrogate_t * restrict inpB_oid,
								surrogate_t * restrict Out_oid,
                      			const uint32_t lenA,
                      			const uint32_t lenB)
{
    uint32_t lenA8 = lenA & ~0x7, lenB8 = lenB & ~0x7;
    uint32_t ai = 0, bi = 0;

    int64_t * out_val = Out_val;
    surrogate_t * out_oid = Out_oid;

    if(lenA8 > 8 && lenB8 > 8) {

        register block8_64 * inA_val  = (block8_64 *) inpA_val;
        register block8_64 * inB_val  = (block8_64 *) inpB_val;
        block8_64 * const    endA_val = (block8_64 *) (inpA_val + lenA) - 1;
        block8_64 * const    endB_val = (block8_64 *) (inpB_val + lenB) - 1;

        register block8_32 * inA_oid  = (block8_32 *) inpA_oid;
        register block8_32 * inB_oid  = (block8_32 *) inpB_oid;

        block8_64 * outp_val = (block8_64 *) out_val;
        block8_32 * outp_oid = (block8_32 *) out_oid;

        register block8_64 * next_val = inB_val;
        register block8_32 * next_oid = inB_oid;

        register __m256d outreg1l_val, outreg1h_val;
        register __m256d outreg2l_val, outreg2h_val;

        register __m256d regAl_val, regAh_val;
        register __m256d regBl_val, regBh_val;

        LOAD8_64(regAl_val, regAh_val, inA_val);
        LOAD8_64(regBl_val, regBh_val, next_val);

        register __m256d outreg1l_oid, outreg1h_oid;
        register __m256d outreg2l_oid, outreg2h_oid;

        register __m256d regAl_oid, regAh_oid;
        register __m256d regBl_oid, regBh_oid;

        LOAD8_32T64(regAl_oid, regAh_oid, inA_oid);
        LOAD8_32T64(regBl_oid, regBh_oid, next_oid);

        inA_val++;
        inB_val++;
        inA_oid++;
        inB_oid++;

        BITONIC_MERGE8_64(outreg1l_val, outreg1h_val, outreg2l_val, outreg2h_val,
        				outreg1l_oid, outreg1h_oid, outreg2l_oid, outreg2h_oid,
                       regAl_val, regAh_val, regBl_val, regBh_val,
					   regAl_oid, regAh_oid, regBl_oid, regBh_oid);

        /* store outreg1 */
        STORE8_64(outp_val, outreg1l_val, outreg1h_val);
        outp_val++;

        /*for oids, convert from 64 back to 32 and store*/
        STORE8_64T32(outp_oid, outreg1l_oid, outreg1h_oid);
        outp_oid++;


        while( inA_val < endA_val && inB_val < endB_val ) {

            /* 3 Options : normal-if, cmove-with-assembly, sw-predication */
            IFELSECONDMOVE(next_val, next_oid, inA_val, inB_val, inA_oid, inB_oid, int64_t);

            regAl_val = outreg2l_val;
            regAh_val = outreg2h_val;

            regAl_oid = outreg2l_oid;
            regAh_oid = outreg2h_oid;

            LOAD8_64(regBl_val, regBh_val, next_val);
            LOAD8_32T64(regBl_oid, regBh_oid, next_oid);

            BITONIC_MERGE8_64(outreg1l_val, outreg1h_val, outreg2l_val, outreg2h_val,
            				outreg1l_oid, outreg1h_oid, outreg2l_oid, outreg2h_oid,
                           regAl_val, regAh_val, regBl_val, regBh_val,
    					   regAl_oid, regAh_oid, regBl_oid, regBh_oid);

            /* store outreg1 */
            STORE8_64(outp_val, outreg1l_val, outreg1h_val);
            outp_val++;

            STORE8_64T32(outp_oid, outreg1l_oid, outreg1h_oid);
            outp_oid++;
        }

        /* flush the register to one of the lists */
        int64_t hireg_val[4] __attribute__((aligned(32)));
        _mm256_store_pd ( (double *)hireg_val, outreg2h_val);

        if(*((int64_t *)inA_val) >= *((int64_t*)(hireg_val+3))) {
            /* store the last remaining register values to A */
            inA_val--;
            inA_oid--;
            STORE8_64(inA_val, outreg2l_val, outreg2h_val);
            STORE8_64T32(inA_oid, outreg2l_oid, outreg2h_oid);
        }
        else {
            /* store the last remaining register values to B */
            inB_val--;
            inB_oid--;
            STORE8_64(inB_val, outreg2l_val, outreg2h_val);
            STORE8_64T32(inB_oid, outreg2l_oid, outreg2h_oid);
        }

        ai = ((int64_t *)inA_val - inpA_val);
        bi = ((int64_t *)inB_val - inpB_val);

        inpA_val = (int64_t *)inA_val;
        inpB_val = (int64_t *)inB_val;
        out_val  = (int64_t *)outp_val;

        inpA_oid = (surrogate_t *)inA_oid;
        inpB_oid = (surrogate_t *)inB_oid;
        out_oid  = (surrogate_t *)outp_oid;
    }

    /* serial-merge */
    while(ai < lenA && bi < lenB){
        int64_t * in_val = inpB_val;
        surrogate_t * in_oid = inpB_oid;
        uint32_t cmp = (*inpA_val < *inpB_val);
        uint32_t notcmp = !cmp;

        ai += cmp;
        bi += notcmp;

        if(cmp) {
            in_val = inpA_val;
            in_oid = inpA_oid;
        }

        *out_val = *in_val;
        *out_oid = *in_oid;
        out_val++;
        out_oid++;
        inpA_val += cmp;
        inpB_val += notcmp;
        inpA_oid += cmp;
        inpB_oid += notcmp;
    }

    if(ai < lenA) {
        /* if A has any more items to be output */
#if 1
        if((lenA - ai) >= 8) {
            /* if A still has some times to be output with AVX */
            uint32_t lenA8_ = ((lenA-ai) & ~0x7);
            register block8_64 * inA_val  = (block8_64 *) inpA_val;
            block8_64 * const    endA_val = (block8_64 *) (inpA_val + lenA8_);
            block8_64 * outp_val = (block8_64 *) out_val;

            register block8_32 *inA_oid = (block8_32 *) inpA_oid;
            block8_32 * outp_oid = (block8_32 *) out_oid;

            while(inA_val < endA_val) {
                __m256d regAl, regAh;
                LOAD8U_64(regAl, regAh, inA_val);
                STORE8U_64(outp_val, regAl, regAh);
                outp_val++;
                inA_val++;

                __m256d regAl_oid, regAh_oid;
                LOAD8U_32T64(regAl_oid, regAh_oid, inA_oid);
                STORE8U_64T32(outp_oid, regAl_oid, regAh_oid);
                outp_oid++;
                inA_oid++;
            }

            ai   += ((int64_t*)inA_val - inpA_val);
            inpA_val  = (int64_t *)inA_val;
            out_val   = (int64_t *)outp_val;

            inpA_oid  = (surrogate_t *)inA_oid;
            out_oid   = (surrogate_t *)outp_oid;
        }
#endif

        while(ai < lenA) {
            *out_val = *inpA_val;
            *out_oid = *inpA_oid;
            ai++;
            out_val++;
            out_oid++;
            inpA_val++;
            inpA_oid++;
        }
    }
    else if(bi < lenB) {
        /* if B has any more items to be output */
#if 1
        if((lenB - bi) >= 8) {
            /* if B still has some times to be output with AVX */
            uint32_t lenB8_ = ((lenB-bi) & ~0x7);
            register block8_64 * inB_val  = (block8_64 *) inpB_val;
            block8_64 * const    endB_val = (block8_64 *) (inpB_val + lenB8_);
            block8_64 * outp_val = (block8_64 *) out_val;

            register block8_32 *inB_oid = (block8_32 *) inpB_oid;
            block8_32 * outp_oid = (block8_32 *) out_oid;

            while(inB_val < endB_val) {
                __m256d regBl, regBh;
                LOAD8U_64(regBl, regBh, inB_val);
                STORE8U_64(outp_val, regBl, regBh);
                outp_val++;
                inB_val++;

                __m256d regBl_oid, regBh_oid;
                LOAD8U_32T64(regBl_oid, regBh_oid, inB_oid);
                STORE8U_64T32(outp_oid, regBl_oid, regBh_oid);
                outp_oid++;
                inB_oid++;
            }

            bi   += ((int64_t*)inB_val - inpB_val);
            inpB_val  = (int64_t *)inB_val;
            out_val   = (int64_t *)outp_val;

            inpB_oid  = (surrogate_t *)inB_oid;
            out_oid   = (surrogate_t *)outp_oid;
        }
#endif
        while(bi < lenB) {
            *out_val = *inpB_val;
            *out_oid = *inpB_oid;
            bi++;
            out_val++;
            out_oid++;
            inpB_val++;
            inpB_oid++;
        }
    }
}
#endif

inline void __attribute__((always_inline))
Mergesort::merge16_int16_varlen_aligned(int16_t * restrict inpA_val,
                      			int16_t * restrict inpB_val,
                      			int16_t * restrict Out_val,
								surrogate_t * restrict inpA_oid,
								surrogate_t * restrict inpB_oid,
								surrogate_t * restrict Out_oid,
                      			const uint32_t lenA,
                      			const uint32_t lenB)
{
#if 1
	assert(0 == ((uint64_t)inpB_val & 31));
	assert(0 == ((uint64_t)inpA_val & 31));
	assert(0 == ((uint64_t)inpA_oid & 31));
	assert(0 == ((uint64_t)inpB_oid & 31));
	assert(0 == ((uint64_t)Out_val & 31));
	assert(0 == ((uint64_t)Out_oid & 31));

    uint32_t lenA16 = lenA & ~0xf, lenB16 = lenB & ~0xf;
    uint32_t ai = 0, bi = 0;

    int16_t * out_val = Out_val;
    surrogate_t * out_oid = Out_oid;

#if 0
		//debug: output the values in each partition （before sorting）
		std::cout << "values merge8_int32_varlen_aligned():" << std::endl;
		uint32_t tmpId;
		std::cout << "inpA:" << lenA << std::endl;
		for(tmpId = 0; tmpId < 100; ++tmpId) {
			std::cout << inpA_val[tmpId] << "[" << inpA_oid[tmpId] << "]" << std::endl;
		}
		std::cout << "inpB:" << 100 << std::endl;
		for(tmpId = 0; tmpId < lenB; ++tmpId) {
			std::cout << inpB_val[tmpId] << "[" << inpB_oid[tmpId] << "]" << std::endl;
		}
		printf("\n");
#endif

    if(lenA16 > 16 && lenB16 > 16) {

        register block16_16 * inA_val  = (block16_16 *) inpA_val;
        register block16_16 * inB_val  = (block16_16 *) inpB_val;
        block16_16 * const    endA_val = (block16_16 *) (inpA_val + lenA) - 1;
        block16_16 * const    endB_val = (block16_16 *) (inpB_val + lenB) - 1;

        register block16_32 * inA_oid  = (block16_32 *) inpA_oid;
        register block16_32 * inB_oid  = (block16_32 *) inpB_oid;

        block16_16 * outp_val = (block16_16 *) out_val;
        block16_32 * outp_oid = (block16_32 *) out_oid;

        register block16_16 * next_val = inB_val;
        register block16_32 * next_oid = inB_oid;

        register __m256i outregl_val, outregh_val;

        register __m256i regA_val, regB_val;

        regA_val = _mm256_load_si256((__m256i const *) inA_val);
        regB_val = _mm256_load_si256((__m256i const *) next_val);

        register __m256i outregl1_oid, outregl2_oid, outregh1_oid, outregh2_oid;

        register __m256i regA1_oid, regA2_oid, regB1_oid, regB2_oid;

        regA1_oid = _mm256_load_si256((__m256i const *) inA_oid);
        regA2_oid = _mm256_load_si256((__m256i const *) (((block8_32 *)inA_oid) + 1));
        regB1_oid = _mm256_load_si256((__m256i const *) next_oid);
        regB2_oid = _mm256_load_si256((__m256i const *) (((block8_32 *)next_oid) + 1));

        inA_val++;
        inB_val++;
        inA_oid++;
        inB_oid++;

        /**Note the final two arguments: they are reversed for the ease of implementation**/
        BITONIC_MERGE16_16(outregl_val, outregh_val,
        		outregl1_oid, outregl2_oid, outregh1_oid, outregh2_oid,
        		regA_val, regB_val, regA1_oid, regA2_oid, regB2_oid, regB1_oid);

        /* store outregl */
        assert(0 == ((uint64_t)outp_val & 31));
        _mm256_store_si256((__m256i *)outp_val, outregl_val);
        outp_val++;

        /*store oids*/
        assert(0 == ((uint64_t)outp_oid & 31));
        _mm256_store_si256((__m256i *)outp_oid, outregl1_oid);
        _mm256_store_si256((__m256i *)((block8_32 *)outp_oid + 1), outregl2_oid);
        outp_oid++;


        while( inA_val < endA_val && inB_val < endB_val ) {

            /* 3 Options : normal-if, cmove-with-assembly, sw-predication */
            IFELSECONDMOVE(next_val, next_oid, inA_val, inB_val, inA_oid, inB_oid, int16_t);

            regA_val = outregh_val;
            regA1_oid = outregh1_oid;
            regA2_oid = outregh2_oid;

            regB_val = _mm256_load_si256((__m256i const *) next_val);
            regB1_oid = _mm256_load_si256((__m256i const *) next_oid);
            regB2_oid = _mm256_load_si256((__m256i const *) (((block8_32 *)next_oid) + 1));


            BITONIC_MERGE16_16(outregl_val, outregh_val,
            		outregl1_oid, outregl2_oid, outregh1_oid, outregh2_oid,
            		regA_val, regB_val, regA1_oid, regA2_oid, regB2_oid, regB1_oid);


            /* store outreg1 */
            assert(0 == ((uint64_t)outp_val & 31));
            _mm256_store_si256((__m256i *)outp_val, outregl_val);
            outp_val++;

            assert(0 == ((uint64_t)outp_oid & 31));
            _mm256_store_si256((__m256i *)outp_oid, outregl1_oid);
            _mm256_store_si256((__m256i *)((block8_32 *)outp_oid + 1), outregl2_oid);
            outp_oid++;
        }

        /* flush the register to one of the lists */
        int16_t hireg_val[16] __attribute__((aligned(32)));
        assert(0 == ((uint64_t)hireg_val & 31));
        _mm256_store_si256( (__m256i *)hireg_val, outregh_val);

        if(*((int16_t *)inA_val) >= *((int16_t*)(hireg_val+15))) {
            /* store the last remaining register values to A */
            inA_val--;
            inA_oid--;
            assert(0 == ((uint64_t)inA_val & 31));
            assert(0 == ((uint64_t)inA_oid & 31));
            _mm256_store_si256((__m256i *)inA_val, outregh_val);
            _mm256_store_si256((__m256i *)inA_oid, outregh1_oid);
            _mm256_store_si256((__m256i *)((block8_32 *)inA_oid + 1), outregh2_oid);
        }
        else {
            /* store the last remaining register values to B */
            inB_val--;
            inB_oid--;
            assert(0 == ((uint64_t)inB_val & 31));
            assert(0 == ((uint64_t)inB_oid & 31));
            _mm256_store_si256((__m256i *)inB_val, outregh_val);
            _mm256_store_si256((__m256i *)inB_oid, outregh1_oid);
            _mm256_store_si256((__m256i *)((block8_32 *)inB_oid + 1), outregh2_oid);
        }

        ai = ((int16_t *)inA_val - inpA_val);
        bi = ((int16_t *)inB_val - inpB_val);

        inpA_val = (int16_t *)inA_val;
        inpB_val = (int16_t *)inB_val;
        out_val  = (int16_t *)outp_val;

        inpA_oid = (surrogate_t *)inA_oid;
        inpB_oid = (surrogate_t *)inB_oid;
        out_oid  = (surrogate_t *)outp_oid;
    }

    /* serial-merge */
    while(ai < lenA && bi < lenB){
        int16_t * in_val = inpB_val;
        surrogate_t * in_oid = inpB_oid;
        uint32_t cmp = (*inpA_val < *inpB_val);
        uint32_t notcmp = !cmp;

        ai += cmp;
        bi += notcmp;

        if(cmp) {
            in_val = inpA_val;
            in_oid = inpA_oid;
        }

        *out_val = *in_val;
        *out_oid = *in_oid;
        out_val++;
        out_oid++;
        inpA_val += cmp;
        inpB_val += notcmp;
        inpA_oid += cmp;
        inpB_oid += notcmp;
    }

    if(ai < lenA) {
        /* if A has any more items to be output */
#if 1
        if((lenA - ai) >= 16) {
            /* if A still has some times to be output with AVX */
            uint32_t lenA16_ = ((lenA-ai) & ~0xf);
            register block16_16 * inA_val  = (block16_16 *) inpA_val;
            block16_16 * const    endA_val = (block16_16 *) (inpA_val + lenA16_);
            block16_16 * outp_val = (block16_16 *) out_val;

            register block16_32 * inA_oid  = (block16_32 *) inpA_oid;
            block16_32 * outp_oid = (block16_32 *) out_oid;

            while(inA_val < endA_val) {
                __m256i regA_val;
                regA_val = _mm256_loadu_si256((__m256i const *) inA_val);
                //assert(0 == ((uint64_t)outp_val & 31));
                _mm256_storeu_si256((__m256i *)outp_val, regA_val);
                outp_val++;
                inA_val++;

                __m256i regA1_oid, regA2_oid;
                regA1_oid = _mm256_loadu_si256((__m256i const *) inA_oid);
                regA2_oid = _mm256_loadu_si256((__m256i const *) ((block8_32 *)inA_oid + 1));
                //assert(0 == ((uint64_t)outp_oid & 31));
                _mm256_storeu_si256((__m256i *)outp_oid, regA1_oid);
                _mm256_storeu_si256((__m256i *)((block8_32 *)outp_oid + 1), regA2_oid);
                outp_oid++;
                inA_oid++;
            }

            ai   += ((int16_t*)inA_val - inpA_val);
            inpA_val  = (int16_t *)inA_val;
            out_val   = (int16_t *)outp_val;

            inpA_oid  = (surrogate_t *)inA_oid;
            out_oid   = (surrogate_t *)outp_oid;
        }
#endif

        while(ai < lenA) {
            *out_val = *inpA_val;
            *out_oid = *inpA_oid;
            ai++;
            out_val++;
            out_oid++;
            inpA_val++;
            inpA_oid++;
        }
    }
    else if(bi < lenB) {
        /* if B has any more items to be output */
#if 1
        if((lenB - bi) >= 16) {
            /* if B still has some times to be output with AVX */
            uint32_t lenB16_ = ((lenB-bi) & ~0xf);
            register block16_16 * inB_val  = (block16_16 *) inpB_val;
            block16_16 * const    endB_val = (block16_16 *) (inpB_val + lenB16_);
            block16_16 * outp_val = (block16_16 *) out_val;

            register block16_32 * inB_oid  = (block16_32 *) inpB_oid;
            block16_32 * outp_oid = (block16_32 *) out_oid;

            while(inB_val < endB_val) {
                __m256i regB_val;
                regB_val = _mm256_loadu_si256((__m256i const *) inB_val);
                //assert(0 == ((uint64_t)outp_val & 31));
                _mm256_storeu_si256((__m256i *)outp_val, regB_val);
                outp_val++;
                inB_val++;

                __m256i regB1_oid, regB2_oid;
                regB1_oid = _mm256_loadu_si256((__m256i const *) inB_oid);
                regB2_oid = _mm256_loadu_si256((__m256i const *) ((block8_32 *)inB_oid + 1));
                //assert(0 == ((uint64_t)outp_oid & 31));
                _mm256_storeu_si256((__m256i *)outp_oid, regB1_oid);
                _mm256_storeu_si256((__m256i *)((block8_32 *)outp_oid + 1), regB2_oid);
                outp_oid++;
                inB_oid++;
            }

            bi   += ((int16_t*)inB_val - inpB_val);
            inpB_val  = (int16_t *)inB_val;
            out_val   = (int16_t *)outp_val;

            inpB_oid  = (surrogate_t *)inB_oid;
            out_oid   = (surrogate_t *)outp_oid;
        }
#endif
        while(bi < lenB) {
            *out_val = *inpB_val;
            *out_oid = *inpB_oid;
            bi++;
            out_val++;
            out_oid++;
            inpB_val++;
            inpB_oid++;
        }
    }
#endif
}

inline void __attribute__((always_inline))
Mergesort::merge8_int32_varlen_aligned(int32_t * restrict inpA_val,
                      			int32_t * restrict inpB_val,
                      			int32_t * restrict Out_val,
								surrogate_t * restrict inpA_oid,
								surrogate_t * restrict inpB_oid,
								surrogate_t * restrict Out_oid,
                      			const uint32_t lenA,
                      			const uint32_t lenB)
{
#if 1
	assert(0 == ((uint64_t)inpB_val & 31));
	assert(0 == ((uint64_t)inpA_val & 31));
	assert(0 == ((uint64_t)inpA_oid & 31));
	assert(0 == ((uint64_t)inpB_oid & 31));

    uint32_t lenA8 = lenA & ~0x7, lenB8 = lenB & ~0x7;
    uint32_t ai = 0, bi = 0;

    int32_t * out_val = Out_val;
    surrogate_t * out_oid = Out_oid;

#if 0
		//debug: output the values in each partition （before sorting）
		std::cout << "values merge8_int32_varlen_aligned():" << std::endl;
		uint32_t tmpId;
		std::cout << "inpA:" << lenA << std::endl;
		for(tmpId = 0; tmpId < 100; ++tmpId) {
			std::cout << inpA_val[tmpId] << "[" << inpA_oid[tmpId] << "]" << std::endl;
		}
		std::cout << "inpB:" << 100 << std::endl;
		for(tmpId = 0; tmpId < lenB; ++tmpId) {
			std::cout << inpB_val[tmpId] << "[" << inpB_oid[tmpId] << "]" << std::endl;
		}
		printf("\n");
#endif

    if(lenA8 > 8 && lenB8 > 8) {

        register block8_32 * inA_val  = (block8_32 *) inpA_val;
        register block8_32 * inB_val  = (block8_32 *) inpB_val;
        block8_32 * const    endA_val = (block8_32 *) (inpA_val + lenA) - 1;
        block8_32 * const    endB_val = (block8_32 *) (inpB_val + lenB) - 1;

        register block8_32 * inA_oid  = (block8_32 *) inpA_oid;
        register block8_32 * inB_oid  = (block8_32 *) inpB_oid;

        block8_32 * outp_val = (block8_32 *) out_val;
        block8_32 * outp_oid = (block8_32 *) out_oid;

        register block8_32 * next_val = inB_val;
        register block8_32 * next_oid = inB_oid;

        register __m256i outregl_val, outregh_val;

        register __m256i regA_val, regB_val;

        regA_val = _mm256_load_si256((__m256i const *) inA_val);
        regB_val = _mm256_load_si256((__m256i const *) next_val);

        register __m256i outregl_oid, outregh_oid;

        register __m256i regA_oid, regB_oid;

        regA_oid = _mm256_load_si256((__m256i const *) inA_oid);
        regB_oid = _mm256_load_si256((__m256i const *) next_oid);

        inA_val++;
        inB_val++;
        inA_oid++;
        inB_oid++;

        BITONIC_MERGE8_32(outregl_val, outregh_val, outregl_oid, outregh_oid,
        		regA_val, regB_val, regA_oid, regB_oid);

        /* store outregl */
        assert(0 == ((uint64_t)outp_val & 31));
        _mm256_store_si256((__m256i *)outp_val, outregl_val);
        outp_val++;

        /*store oids*/
        assert(0 == ((uint64_t)outp_oid & 31));
        _mm256_store_si256((__m256i *)outp_oid, outregl_oid);
        outp_oid++;


        while( inA_val < endA_val && inB_val < endB_val ) {

            /* 3 Options : normal-if, cmove-with-assembly, sw-predication */
            IFELSECONDMOVE(next_val, next_oid, inA_val, inB_val, inA_oid, inB_oid, int32_t);

            regA_val = outregh_val;
            regA_oid = outregh_oid;

            regB_val = _mm256_load_si256((__m256i const *) next_val);
            regB_oid = _mm256_load_si256((__m256i const *) next_oid);


            BITONIC_MERGE8_32(outregl_val, outregh_val, outregl_oid, outregh_oid,
            		regA_val, regB_val, regA_oid, regB_oid);


            /* store outreg1 */
            assert(0 == ((uint64_t)outp_val & 31));
            _mm256_store_si256((__m256i *)outp_val, outregl_val);
            outp_val++;

            assert(0 == ((uint64_t)outp_oid & 31));
            _mm256_store_si256((__m256i *)outp_oid, outregl_oid);
            outp_oid++;
        }

        /* flush the register to one of the lists */
        int32_t hireg_val[8] __attribute__((aligned(32)));
        assert(0 == ((uint64_t)hireg_val & 31));
        _mm256_store_si256( (__m256i *)hireg_val, outregh_val);

        if(*((int32_t *)inA_val) >= *((int32_t*)(hireg_val+7))) {
            /* store the last remaining register values to A */
            inA_val--;
            inA_oid--;
            assert(0 == ((uint64_t)inA_val & 31));
            assert(0 == ((uint64_t)inA_oid & 31));
            _mm256_store_si256((__m256i *)inA_val, outregh_val);
            _mm256_store_si256((__m256i *)inA_oid, outregh_oid);
        }
        else {
            /* store the last remaining register values to B */
            inB_val--;
            inB_oid--;
            assert(0 == ((uint64_t)inB_val & 31));
            assert(0 == ((uint64_t)inB_oid & 31));
            _mm256_store_si256((__m256i *)inB_val, outregh_val);
            _mm256_store_si256((__m256i *)inB_oid, outregh_oid);
        }

        ai = ((int32_t *)inA_val - inpA_val);
        bi = ((int32_t *)inB_val - inpB_val);

        inpA_val = (int32_t *)inA_val;
        inpB_val = (int32_t *)inB_val;
        out_val  = (int32_t *)outp_val;

        inpA_oid = (surrogate_t *)inA_oid;
        inpB_oid = (surrogate_t *)inB_oid;
        out_oid  = (surrogate_t *)outp_oid;
    }

    /* serial-merge */
    while(ai < lenA && bi < lenB){
        int32_t * in_val = inpB_val;
        surrogate_t * in_oid = inpB_oid;
        uint32_t cmp = (*inpA_val < *inpB_val);
        uint32_t notcmp = !cmp;

        ai += cmp;
        bi += notcmp;

        if(cmp) {
            in_val = inpA_val;
            in_oid = inpA_oid;
        }

        *out_val = *in_val;
        *out_oid = *in_oid;
        out_val++;
        out_oid++;
        inpA_val += cmp;
        inpB_val += notcmp;
        inpA_oid += cmp;
        inpB_oid += notcmp;
    }

    if(ai < lenA) {
        /* if A has any more items to be output */
#if 1
        if((lenA - ai) >= 8) {
            /* if A still has some times to be output with AVX */
            uint32_t lenA8_ = ((lenA-ai) & ~0x7);
            register block8_32 * inA_val  = (block8_32 *) inpA_val;
            block8_32 * const    endA_val = (block8_32 *) (inpA_val + lenA8_);
            block8_32 * outp_val = (block8_32 *) out_val;

            register block8_32 * inA_oid  = (block8_32 *) inpA_oid;
            block8_32 * outp_oid = (block8_32 *) out_oid;

            while(inA_val < endA_val) {
                __m256i regA_val;
                regA_val = _mm256_loadu_si256((__m256i const *) inA_val);
                //assert(0 == ((uint64_t)outp_val & 31));
                _mm256_storeu_si256((__m256i *)outp_val, regA_val);
                outp_val++;
                inA_val++;

                __m256i regA_oid;
                regA_oid = _mm256_loadu_si256((__m256i const *) inA_oid);
                //assert(0 == ((uint64_t)outp_oid & 31));
                _mm256_storeu_si256((__m256i *)outp_oid, regA_oid);
                outp_oid++;
                inA_oid++;
            }

            ai   += ((int32_t*)inA_val - inpA_val);
            inpA_val  = (int32_t *)inA_val;
            out_val   = (int32_t *)outp_val;

            inpA_oid  = (surrogate_t *)inA_oid;
            out_oid   = (surrogate_t *)outp_oid;
        }
#endif

        while(ai < lenA) {
            *out_val = *inpA_val;
            *out_oid = *inpA_oid;
            ai++;
            out_val++;
            out_oid++;
            inpA_val++;
            inpA_oid++;
        }
    }
    else if(bi < lenB) {
        /* if B has any more items to be output */
#if 1
        if((lenB - bi) >= 8) {
            /* if B still has some times to be output with AVX */
            uint32_t lenB8_ = ((lenB-bi) & ~0x7);
            register block8_32 * inB_val  = (block8_32 *) inpB_val;
            block8_32 * const    endB_val = (block8_32 *) (inpB_val + lenB8_);
            block8_32 * outp_val = (block8_32 *) out_val;

            register block8_32 * inB_oid  = (block8_32 *) inpB_oid;
            block8_32 * outp_oid = (block8_32 *) out_oid;

            while(inB_val < endB_val) {
                __m256i regB_val;
                regB_val = _mm256_loadu_si256((__m256i const *) inB_val);
                //assert(0 == ((uint64_t)outp_val & 31));
                _mm256_storeu_si256((__m256i *)outp_val, regB_val);
                outp_val++;
                inB_val++;

                __m256i regB_oid;
                regB_oid = _mm256_loadu_si256((__m256i const *) inB_oid);
                //assert(0 == ((uint64_t)outp_oid & 31));
                _mm256_storeu_si256((__m256i *)outp_oid, regB_oid);
                outp_oid++;
                inB_oid++;
            }

            bi   += ((int32_t*)inB_val - inpB_val);
            inpB_val  = (int32_t *)inB_val;
            out_val   = (int32_t *)outp_val;

            inpB_oid  = (surrogate_t *)inB_oid;
            out_oid   = (surrogate_t *)outp_oid;
        }
#endif
        while(bi < lenB) {
            *out_val = *inpB_val;
            *out_oid = *inpB_oid;
            bi++;
            out_val++;
            out_oid++;
            inpB_val++;
            inpB_oid++;
        }
    }
#endif
}

inline void __attribute__((always_inline))
Mergesort::merge8_int32_varlen_unaligned(int32_t * restrict inpA_val,
                      			int32_t * restrict inpB_val,
                      			int32_t * restrict Out_val,
								surrogate_t * restrict inpA_oid,
								surrogate_t * restrict inpB_oid,
								surrogate_t * restrict Out_oid,
                      			const uint32_t lenA,
                      			const uint32_t lenB)
{
#if 1

    uint32_t lenA8 = lenA & ~0x7, lenB8 = lenB & ~0x7;
    uint32_t ai = 0, bi = 0;

    int32_t * out_val = Out_val;
    surrogate_t * out_oid = Out_oid;

#if 0
		//debug: output the values in each partition （before sorting）
		std::cout << "values merge8_int32_varlen_aligned():" << std::endl;
		uint32_t tmpId;
		std::cout << "inpA:" << lenA << std::endl;
		for(tmpId = 0; tmpId < 100; ++tmpId) {
			std::cout << inpA_val[tmpId] << "[" << inpA_oid[tmpId] << "]" << std::endl;
		}
		std::cout << "inpB:" << 100 << std::endl;
		for(tmpId = 0; tmpId < lenB; ++tmpId) {
			std::cout << inpB_val[tmpId] << "[" << inpB_oid[tmpId] << "]" << std::endl;
		}
		printf("\n");
#endif

    if(lenA8 > 8 && lenB8 > 8) {

        register block8_32 * inA_val  = (block8_32 *) inpA_val;
        register block8_32 * inB_val  = (block8_32 *) inpB_val;
        block8_32 * const    endA_val = (block8_32 *) (inpA_val + lenA) - 1;
        block8_32 * const    endB_val = (block8_32 *) (inpB_val + lenB) - 1;

        register block8_32 * inA_oid  = (block8_32 *) inpA_oid;
        register block8_32 * inB_oid  = (block8_32 *) inpB_oid;

        block8_32 * outp_val = (block8_32 *) out_val;
        block8_32 * outp_oid = (block8_32 *) out_oid;

        register block8_32 * next_val = inB_val;
        register block8_32 * next_oid = inB_oid;

        register __m256i outregl_val, outregh_val;

        register __m256i regA_val, regB_val;

        regA_val = _mm256_loadu_si256((__m256i const *) inA_val);
        regB_val = _mm256_loadu_si256((__m256i const *) next_val);

        register __m256i outregl_oid, outregh_oid;

        register __m256i regA_oid, regB_oid;

        regA_oid = _mm256_loadu_si256((__m256i const *) inA_oid);
        regB_oid = _mm256_loadu_si256((__m256i const *) next_oid);

        inA_val++;
        inB_val++;
        inA_oid++;
        inB_oid++;

        BITONIC_MERGE8_32(outregl_val, outregh_val, outregl_oid, outregh_oid,
        		regA_val, regB_val, regA_oid, regB_oid);

        /* store outregl */

        _mm256_storeu_si256((__m256i *)outp_val, outregl_val);
        outp_val++;

        /*store oids*/

        _mm256_storeu_si256((__m256i *)outp_oid, outregl_oid);
        outp_oid++;


        while( inA_val < endA_val && inB_val < endB_val ) {

            /* 3 Options : normal-if, cmove-with-assembly, sw-predication */
            IFELSECONDMOVE(next_val, next_oid, inA_val, inB_val, inA_oid, inB_oid, int32_t);

            regA_val = outregh_val;
            regA_oid = outregh_oid;

            regB_val = _mm256_loadu_si256((__m256i const *) next_val);
            regB_oid = _mm256_loadu_si256((__m256i const *) next_oid);


            BITONIC_MERGE8_32(outregl_val, outregh_val, outregl_oid, outregh_oid,
            		regA_val, regB_val, regA_oid, regB_oid);


            /* store outreg1 */

            _mm256_storeu_si256((__m256i *)outp_val, outregl_val);
            outp_val++;


            _mm256_storeu_si256((__m256i *)outp_oid, outregl_oid);
            outp_oid++;
        }

        /* flush the register to one of the lists */
        int32_t hireg_val[8] __attribute__((aligned(32)));

        _mm256_store_si256( (__m256i *)hireg_val, outregh_val);

        if(*((int32_t *)inA_val) >= *((int32_t*)(hireg_val+7))) {
            /* store the last remaining register values to A */
            inA_val--;
            inA_oid--;
            _mm256_storeu_si256((__m256i *)inA_val, outregh_val);
            _mm256_storeu_si256((__m256i *)inA_oid, outregh_oid);
        }
        else {
            /* store the last remaining register values to B */
            inB_val--;
            inB_oid--;
            _mm256_storeu_si256((__m256i *)inB_val, outregh_val);
            _mm256_storeu_si256((__m256i *)inB_oid, outregh_oid);
        }

        ai = ((int32_t *)inA_val - inpA_val);
        bi = ((int32_t *)inB_val - inpB_val);

        inpA_val = (int32_t *)inA_val;
        inpB_val = (int32_t *)inB_val;
        out_val  = (int32_t *)outp_val;

        inpA_oid = (surrogate_t *)inA_oid;
        inpB_oid = (surrogate_t *)inB_oid;
        out_oid  = (surrogate_t *)outp_oid;
    }

    /* serial-merge */
    while(ai < lenA && bi < lenB){
        int32_t * in_val = inpB_val;
        surrogate_t * in_oid = inpB_oid;
        uint32_t cmp = (*inpA_val < *inpB_val);
        uint32_t notcmp = !cmp;

        ai += cmp;
        bi += notcmp;

        if(cmp) {
            in_val = inpA_val;
            in_oid = inpA_oid;
        }

        *out_val = *in_val;
        *out_oid = *in_oid;
        out_val++;
        out_oid++;
        inpA_val += cmp;
        inpB_val += notcmp;
        inpA_oid += cmp;
        inpB_oid += notcmp;
    }

    if(ai < lenA) {
        /* if A has any more items to be output */
#if 1
        if((lenA - ai) >= 8) {
            /* if A still has some times to be output with AVX */
            uint32_t lenA8_ = ((lenA-ai) & ~0x7);
            register block8_32 * inA_val  = (block8_32 *) inpA_val;
            block8_32 * const    endA_val = (block8_32 *) (inpA_val + lenA8_);
            block8_32 * outp_val = (block8_32 *) out_val;

            register block8_32 * inA_oid  = (block8_32 *) inpA_oid;
            block8_32 * outp_oid = (block8_32 *) out_oid;

            while(inA_val < endA_val) {
                __m256i regA_val;
                regA_val = _mm256_loadu_si256((__m256i const *) inA_val);
                //assert(0 == ((uint64_t)outp_val & 31));
                _mm256_storeu_si256((__m256i *)outp_val, regA_val);
                outp_val++;
                inA_val++;

                __m256i regA_oid;
                regA_oid = _mm256_loadu_si256((__m256i const *) inA_oid);
                //assert(0 == ((uint64_t)outp_oid & 31));
                _mm256_storeu_si256((__m256i *)outp_oid, regA_oid);
                outp_oid++;
                inA_oid++;
            }

            ai   += ((int32_t*)inA_val - inpA_val);
            inpA_val  = (int32_t *)inA_val;
            out_val   = (int32_t *)outp_val;

            inpA_oid  = (surrogate_t *)inA_oid;
            out_oid   = (surrogate_t *)outp_oid;
        }
#endif

        while(ai < lenA) {
            *out_val = *inpA_val;
            *out_oid = *inpA_oid;
            ai++;
            out_val++;
            out_oid++;
            inpA_val++;
            inpA_oid++;
        }
    }
    else if(bi < lenB) {
        /* if B has any more items to be output */
#if 1
        if((lenB - bi) >= 8) {
            /* if B still has some times to be output with AVX */
            uint32_t lenB8_ = ((lenB-bi) & ~0x7);
            register block8_32 * inB_val  = (block8_32 *) inpB_val;
            block8_32 * const    endB_val = (block8_32 *) (inpB_val + lenB8_);
            block8_32 * outp_val = (block8_32 *) out_val;

            register block8_32 * inB_oid  = (block8_32 *) inpB_oid;
            block8_32 * outp_oid = (block8_32 *) out_oid;

            while(inB_val < endB_val) {
                __m256i regB_val;
                regB_val = _mm256_loadu_si256((__m256i const *) inB_val);
                //assert(0 == ((uint64_t)outp_val & 31));
                _mm256_storeu_si256((__m256i *)outp_val, regB_val);
                outp_val++;
                inB_val++;

                __m256i regB_oid;
                regB_oid = _mm256_loadu_si256((__m256i const *) inB_oid);
                //assert(0 == ((uint64_t)outp_oid & 31));
                _mm256_storeu_si256((__m256i *)outp_oid, regB_oid);
                outp_oid++;
                inB_oid++;
            }

            bi   += ((int32_t*)inB_val - inpB_val);
            inpB_val  = (int32_t *)inB_val;
            out_val   = (int32_t *)outp_val;

            inpB_oid  = (surrogate_t *)inB_oid;
            out_oid   = (surrogate_t *)outp_oid;
        }
#endif
        while(bi < lenB) {
            *out_val = *inpB_val;
            *out_oid = *inpB_oid;
            bi++;
            out_val++;
            out_oid++;
            inpB_val++;
            inpB_oid++;
        }
    }
#endif
}


inline void __attribute__((always_inline))
Mergesort::avxmergesort_block_uint16_aligned(uint16_t** inputptr_val, surrogate_t** inputptr_oid,
		uint16_t** outputptr_val, surrogate_t** outputptr_oid, uint32_t block_size)
{
    uint16_t* ptrs_val[2];
    surrogate_t* ptrs_oid[2];
    const uint64_t logBSZ = mylog2(block_size);

    ptrs_val[0] = *inputptr_val;
    ptrs_val[1] = *outputptr_val;

    ptrs_oid[0] = *inputptr_oid;
    ptrs_oid[1] = *outputptr_oid;

    /** 1.a) Perform in-register sort (i.e., sorting network) to get sorted seq of K(K=16)*/
    block256_16 * inptr_val 		= (block256_16 *) ptrs_val[0];
    block256_16 * const end_val 	= (block256_16 *) (ptrs_val[0] + block_size);
    block256_32 * inptr_oid 		= (block256_32 *) ptrs_oid[0];
    //block16_32 * const end_oid 	= (block16_32 *) (ptrs_oid[0] + block_size);

    while(inptr_val < end_val) {
        inregister_sort_int16_aligned((int16_t*)inptr_val, (int16_t*)inptr_val,
        		(int32_t*)inptr_oid, (int32_t*)inptr_oid);
        inptr_val++;
        inptr_oid++;
    }

#if 1
    /**
     * 1.b) for itr <- [(logK) .. (logM - 3)]
     *  - Simultaneously merge 4 sequences (using a K by K
     *  network) of length 2^itr to obtain sorted seq. of 2^{itr+1}
     */
    uint64_t j;
    const uint64_t jend = logBSZ - 2;

    for(j = 4; j < jend; j++) {	//important: j starts from 4 because of merge16
        int ptridx = j & 1;
        int16_t * inp = (int16_t *) ptrs_val[ptridx];
        int16_t * out = (int16_t *) ptrs_val[ptridx ^ 1];
        int16_t * const end = (int16_t*) (inp + block_size);

        surrogate_t * inp_oid = ptrs_oid[ptridx];
        surrogate_t * out_oid = ptrs_oid[ptridx ^ 1];

        /**
         *  merge length 2^j lists beginnig at inp and output a
         *  sorted list of length 2^(j+1) starting at out
         */
        const uint64_t inlen  = (1 << j);
        const uint64_t outlen = (inlen << 1);

        while(inp < end) {

        	merge16_int16_eqlen_aligned(inp, inp + inlen,
        	                    inp_oid, inp_oid + inlen,
        	                    out, out_oid, inlen);

            inp += outlen;
            out += outlen;
            inp_oid += outlen;
            out_oid += outlen;

            /* TODO: Try following. */
            /* simultaneous merge of 4 list pairs */
            /* merge 4 seqs simultaneously (always >= 4) */
            /* merge 2 seqs simultaneously (always >= 2) */
        }
    }

    /**
     * 1.c) for itr = (logM - 2), simultaneously merge 2 sequences
     *  (using a 2K by 2K network) of length M/4 to obtain sorted
     *  sequences of M/2.
     */
    uint64_t inlen  = (1 << j);
    int16_t * inp;
    surrogate_t * inp_oid;
    int16_t * out;
    surrogate_t * out_oid;
    int ptridx = j & 1;

    inp = (int16_t *)ptrs_val[ptridx];
    out = (int16_t *)ptrs_val[ptridx ^ 1];

    inp_oid = ptrs_oid[ptridx];
    out_oid = ptrs_oid[ptridx ^ 1];

	merge16_int16_eqlen_aligned(inp, inp + inlen,
		 inp_oid, inp_oid + inlen,
		 out, out_oid, inlen);
	merge16_int16_eqlen_aligned(inp+2*inlen, inp + 3*inlen,
		 inp_oid+2*inlen, inp_oid + 3*inlen,
		 out+2*inlen, out_oid+2*inlen, inlen);

    /* TODO: simultaneous merge of 2 list pairs */
    /**
     * 1.d) for itr = (logM - 1), merge 2 final sequences (using a
     * 4K by 4K network) of length M/2 to get sorted seq. of M.
     */
    j++; /* j=(LOG2_BLOCK_SIZE-1); inputsize M/2 --> outputsize M*/
    inlen  = (1 << j);
    /* now we know that input is out from the last pass */
	merge16_int16_eqlen_aligned(out, out + inlen,
			out_oid, out_oid + inlen,
			inp, inp_oid, inlen);

    /* finally swap input/output ptrs, output is the sorted list */
    * outputptr_val = (uint16_t *)inp;
    * inputptr_val  = (uint16_t *)out;

    * outputptr_oid = inp_oid;
    * inputptr_oid = out_oid;

#endif
}

inline void __attribute__((always_inline))
Mergesort::avxmergesort_block_uint32_aligned(uint32_t** inputptr_val, surrogate_t** inputptr_oid,
		uint32_t** outputptr_val, surrogate_t** outputptr_oid, uint32_t block_size)
{
    uint32_t* ptrs_val[2];
    surrogate_t* ptrs_oid[2];
    const uint64_t logBSZ = mylog2(block_size);

    ptrs_val[0] = *inputptr_val;
    ptrs_val[1] = *outputptr_val;

    ptrs_oid[0] = *inputptr_oid;
    ptrs_oid[1] = *outputptr_oid;

#if 0
		//debug: output the values in each partition （before sorting）
		std::cout << "values avxmergesort_block_uint32_aligned to be sorted...(output first 100 and last 100), blocksize " << block_size << std::endl;
		uint32_t tmpId;
		for  (tmpId = 0; tmpId < 100; ++tmpId) {
			std::cout << ptrs_val[0][tmpId] << "[" << ptrs_oid[0][tmpId] << "]" << std::endl;
		}
		for  (tmpId = block_size-100; tmpId < block_size; ++tmpId) {
			std::cout << ptrs_val[0][tmpId] << "[" << ptrs_oid[0][tmpId] << "]" << std::endl;
		}
		printf("\n");
#endif

    /** 1.a) Perform in-register sort (i.e., sorting network) to get sorted seq of K(K=8)*/
    block64_32 * inptr_val 		= (block64_32 *) ptrs_val[0];
    block64_32 * const end_val 	= (block64_32 *) (ptrs_val[0] + block_size);
    block64_32 * inptr_oid 		= (block64_32 *) ptrs_oid[0];
    //block16_32 * const end_oid 	= (block16_32 *) (ptrs_oid[0] + block_size);

    while(inptr_val < end_val) {
#if 0
        for (uint32_t idx = 0; idx < 64; ++idx) {
        	assert(((int32_t*)inptr_oid)[idx] < 16384);
        }
#endif
        inregister_sort_int32_aligned((int32_t*)inptr_val, (int32_t*)inptr_val,
        		(int32_t*)inptr_oid, (int32_t*)inptr_oid);
        inptr_val++;
        inptr_oid++;
#if 0
        for (uint32_t idx = 0; idx < 64; ++idx) {
        	assert(((int32_t*)inptr_oid)[idx] < 16384);
        }
#endif
    }


    /**
     * 1.b) for itr <- [(logK) .. (logM - 3)]
     *  - Simultaneously merge 4 sequences (using a K by K
     *  network) of length 2^itr to obtain sorted seq. of 2^{itr+1}
     */
    uint64_t j;
    const uint64_t jend = logBSZ - 2;

    for(j = 3; j < jend; ++j)	//important: j starts from 3 because of merge8
    {
        int ptridx = (j-1) & 1;
        int32_t * inp_val = (int32_t *) ptrs_val[ptridx];
        int32_t * out_val = (int32_t *) ptrs_val[ptridx ^ 1];
        int32_t * const end_val = (int32_t*) (inp_val + block_size);

        surrogate_t * inp_oid = ptrs_oid[ptridx];
        surrogate_t * out_oid = ptrs_oid[ptridx ^ 1];

        /**
         *  merge length 2^j lists beginnig at inp and output a
         *  sorted list of length 2^(j+1) starting at out
         */
        const uint64_t inlen  = (1 << j);
        const uint64_t outlen = (inlen << 1);

        while(inp_val < end_val) {

            merge8_int32_eqlen_aligned(inp_val, inp_val + inlen,
            		inp_oid, inp_oid + inlen,
            		out_val, out_oid, inlen);
            inp_val += outlen;
            out_val += outlen;
            inp_oid += outlen;
            out_oid += outlen;
        }
    }
#if 0
    j = 3;
    {
        int ptridx = j & 1;
        int64_t * inp = (int64_t *) ptrs_val[ptridx];
        int64_t * out = (int64_t *) ptrs_val[ptridx ^ 1];
        int64_t * const end = (int64_t*) (inp + block_size);

        /**
         *  merge length 2^j lists beginnig at inp and output a
         *  sorted list of length 2^(j+1) starting at out
         */
        const uint64_t inlen  = (1 << j);
        const uint64_t outlen = (inlen << 1);

        while(inp < end) {

            merge8_eqlen_aligned(inp, inp + inlen, out, inlen);
            inp += outlen;
            out += outlen;
        }
    }
    for(j = 4; j < jend; j++) {
        int ptridx = j & 1;
        int64_t * inp = (int64_t *) ptrs_val[ptridx];
        int64_t * out = (int64_t *) ptrs_val[ptridx ^ 1];
        int64_t * const end = (int64_t*) (inp + block_size);

        /**
         *  merge length 2^j lists beginnig at inp and output a
         *  sorted list of length 2^(j+1) starting at out
         */
        const uint64_t inlen  = (1 << j);
        const uint64_t outlen = (inlen << 1);

        while(inp < end) {

            merge16_eqlen_aligned(inp, inp + inlen, out, inlen);
            inp += outlen;
            out += outlen;

            /* TODO: Try following. */
            /* simultaneous merge of 4 list pairs */
            /* merge 4 seqs simultaneously (always >= 4) */
            /* merge 2 seqs simultaneously (always >= 2) */
        }
    }
#endif
    /**
     * 1.c) for itr = (logM - 2), simultaneously merge 2 sequences
     *  (using a 2K by 2K network) of length M/4 to obtain sorted
     *  sequences of M/2.
     */
    uint64_t inlen  = (1 << j);
    int32_t * inp_val;
    surrogate_t * inp_oid;
    int32_t * out_val;
    surrogate_t * out_oid;
    int ptridx = (j-1) & 1;

    inp_val = (int32_t *)ptrs_val[ptridx];
    out_val = (int32_t *)ptrs_val[ptridx ^ 1];

    inp_oid = ptrs_oid[ptridx];
    out_oid = ptrs_oid[ptridx ^ 1];

    merge8_int32_eqlen_aligned(inp_val, inp_val + inlen,
    		inp_oid, inp_oid + inlen,
    		out_val, out_oid, inlen);
    merge8_int32_eqlen_aligned(inp_val+2*inlen, inp_val + 3*inlen,
    		inp_oid+2*inlen, inp_oid + 3*inlen,
    		out_val+2*inlen, out_oid+2*inlen, inlen);

    /* TODO: simultaneous merge of 2 list pairs */
    /**
     * 1.d) for itr = (logM - 1), merge 2 final sequences (using a
     * 4K by 4K network) of length M/2 to get sorted seq. of M.
     */
    j++; /* j=(LOG2_BLOCK_SIZE-1); inputsize M/2 --> outputsize M*/
    inlen  = (1 << j);
    /* now we know that input is out from the last pass */
    merge8_int32_eqlen_aligned(out_val, out_val + inlen,
    		out_oid, out_oid + inlen,
    		inp_val, inp_oid, inlen);

    /* finally swap input/output ptrs, output is the sorted list */
    * outputptr_val = (uint32_t *)inp_val;
    * inputptr_val  = (uint32_t *)out_val;

    * outputptr_oid = inp_oid;
    * inputptr_oid = out_oid;

#if 0
        for (uint32_t idx = 0; idx < block_size; ++idx) {
        	assert(((int32_t*)inptr_oid)[idx] < 16384);
        }
#endif

#if 0
		//debug: output the values after sorting blocks
		std::cout << "output the values after sorting blocks: " << std::endl;
		uint32_t tmpId;
		for  (tmpId = 0; tmpId < 100; ++tmpId) {
			std::cout << (*outputptr_val)[tmpId] << std::endl;
		}
		for  (tmpId = block_size-100; tmpId < block_size; ++tmpId) {
			std::cout << (*outputptr_val)[tmpId] << std::endl;
		}
		printf("\n");
#endif
}

inline void __attribute__((always_inline))
Mergesort::avxmergesort_block_uint64_aligned(uint64_t** inputptr_val, surrogate_t** inputptr_oid,
		uint64_t** outputptr_val, surrogate_t** outputptr_oid, uint32_t block_size)
{
    uint64_t* ptrs_val[2];
    surrogate_t* ptrs_oid[2];
    const uint64_t logBSZ = mylog2(block_size);

    ptrs_val[0] = *inputptr_val;
    ptrs_val[1] = *outputptr_val;

    ptrs_oid[0] = *inputptr_oid;
    ptrs_oid[1] = *outputptr_oid;

    /** 1.a) Perform in-register sort (i.e., sorting network) to get sorted seq of K(K=4)*/
    block16_64 * inptr_val 		= (block16_64 *) ptrs_val[0];
    block16_64 * const end_val 	= (block16_64 *) (ptrs_val[0] + block_size);
    block16_32 * inptr_oid 		= (block16_32 *) ptrs_oid[0];
    //block16_32 * const end_oid 	= (block16_32 *) (ptrs_oid[0] + block_size);

    while(inptr_val < end_val) {
        inregister_sort_int64_aligned((int64_t*)inptr_val, (int64_t*)inptr_val,
        		(int32_t*)inptr_oid, (int32_t*)inptr_oid);
        inptr_val++;
        inptr_oid++;
    }

    /**
     * 1.b) for itr <- [(logK) .. (logM - 3)]
     *  - Simultaneously merge 4 sequences (using a K by K
     *  network) of length 2^itr to obtain sorted seq. of 2^{itr+1}
     */
    uint64_t j;
    const uint64_t jend = logBSZ - 2;

    for(j = 2; j < jend; j++)
    {
        int ptridx = j & 1;
        int64_t * inp_val = (int64_t *) ptrs_val[ptridx];
        int64_t * out_val = (int64_t *) ptrs_val[ptridx ^ 1];
        int64_t * const end_val = (int64_t*) (inp_val + block_size);

        surrogate_t * inp_oid = ptrs_oid[ptridx];
        surrogate_t * out_oid = ptrs_oid[ptridx ^ 1];

        /**
         *  merge length 2^j lists beginning at inp and output a
         *  sorted list of length 2^(j+1) starting at out
         */
        const uint64_t inlen  = (1 << j);
        const uint64_t outlen = (inlen << 1);

        while(inp_val < end_val) {

            merge4_int64_eqlen_aligned(inp_val, inp_val + inlen, inp_oid, inp_oid + inlen, out_val, out_oid, inlen);
            inp_val += outlen;
            out_val += outlen;
            inp_oid += outlen;
            out_oid += outlen;
        }
    }
#if 0
    for(j = 3; j < jend; j++)
    {
        int ptridx = j & 1;
        int64_t * inp_val = (int64_t *) ptrs_val[ptridx];
        int64_t * out_val = (int64_t *) ptrs_val[ptridx ^ 1];
        int64_t * const end_val = (int64_t*) (inp_val + block_size);

        surrogate_t * inp_oid = ptrs_oid[ptridx];
        surrogate_t * out_oid = ptrs_oid[ptridx ^ 1];

        /**
         *  merge length 2^j lists beginnig at inp and output a
         *  sorted list of length 2^(j+1) starting at out
         */
        const uint64_t inlen  = (1 << j);
        const uint64_t outlen = (inlen << 1);

        while(inp_val < end) {

            merge8_eqlen_aligned(inp, inp + inlen, out, inlen);
            inp += outlen;
            out += outlen;
        }
    }

    for(j = 4; j < jend; j++) {
        int ptridx = j & 1;
        int64_t * inp = (int64_t *) ptrs_val[ptridx];
        int64_t * out = (int64_t *) ptrs_val[ptridx ^ 1];
        int64_t * const end = (int64_t*) (inp + block_size);

        /**
         *  merge length 2^j lists beginnig at inp and output a
         *  sorted list of length 2^(j+1) starting at out
         */
        const uint64_t inlen  = (1 << j);
        const uint64_t outlen = (inlen << 1);

        while(inp < end) {

            merge16_eqlen_aligned(inp, inp + inlen, out, inlen);
            inp += outlen;
            out += outlen;

            /* TODO: Try following. */
            /* simultaneous merge of 4 list pairs */
            /* merge 4 seqs simultaneously (always >= 4) */
            /* merge 2 seqs simultaneously (always >= 2) */
        }
    }
#endif
    /**
     * 1.c) for itr = (logM - 2), simultaneously merge 2 sequences
     *  (using a 2K by 2K network) of length M/4 to obtain sorted
     *  sequences of M/2.
     */
    uint64_t inlen  = (1 << j);
    int64_t * inp_val;
    surrogate_t * inp_oid;
    int64_t * out_val;
    surrogate_t * out_oid;
    int ptridx = j & 1;

    inp_val = (int64_t *)ptrs_val[ptridx];
    out_val = (int64_t *)ptrs_val[ptridx ^ 1];

    inp_oid = ptrs_oid[ptridx];
    out_oid = ptrs_oid[ptridx ^ 1];

    merge4_int64_eqlen_aligned(inp_val, inp_val + inlen, inp_oid, inp_oid + inlen, out_val, out_oid, inlen);
    merge4_int64_eqlen_aligned(inp_val+2*inlen, inp_val+3*inlen,
    		inp_oid+2*inlen, inp_oid+3*inlen, out_val + 2*inlen, out_oid+2*inlen, inlen);

    /* TODO: simultaneous merge of 2 list pairs */
    /**
     * 1.d) for itr = (logM - 1), merge 2 final sequences (using a
     * 4K by 4K network) of length M/2 to get sorted seq. of M.
     */
    j++; /* j=(LOG2_BLOCK_SIZE-1); inputsize M/2 --> outputsize M*/
    inlen  = (1 << j);
    /* now we know that input is out from the last pass */
    merge4_int64_eqlen_aligned(out_val, out_val + inlen, out_oid, out_oid + inlen, inp_val, inp_oid, inlen);

    /* finally swap input/output ptrs, output is the sorted list */
    * outputptr_val = (uint64_t *)inp_val;
    * inputptr_val  = (uint64_t *)out_val;

    * outputptr_oid = inp_oid;
    * inputptr_oid  = out_oid;
}

inline void __attribute__((always_inline))
Mergesort::merge4_int64_eqlen_aligned(int64_t * const inpA, int64_t * const inpB,
						surrogate_t * const inpA_oid, surrogate_t * const inpB_oid,
						int64_t * const out, surrogate_t * const out_oid, const uint32_t len)
{
    register block4_64 * inA  = (block4_64 *) inpA;
    register block4_64 * inB  = (block4_64 *) inpB;
    block4_64 * const    endA = (block4_64 *) (inpA + len);
    block4_64 * const    endB = (block4_64 *) (inpB + len);

    register block4_32 * inA_oid  = (block4_32 *) inpA_oid;
    register block4_32 * inB_oid  = (block4_32 *) inpB_oid;

    block4_64 * outp = (block4_64 *) out;
    block4_32 * outp_oid = (block4_32 *) out_oid;

    register block4_64 * next = inB;
    register block4_32 * next_oid = inB_oid;

    register __m256d outreg1;
    register __m256d outreg2;

    register __m256d outreg1_oid;
    register __m256d outreg2_oid;

    register __m256d regA = _mm256_load_pd((double const *) inA);
    register __m256d regB = _mm256_load_pd((double const *) next);

    register __m256d regA_oid = _mm256_cvtepi32_pd(_mm_load_si128((__m128i const*) inA_oid));
    register __m256d regB_oid = _mm256_cvtepi32_pd(_mm_load_si128((__m128i const*) next_oid));

    inA++;
    inB++;

    inA_oid++;
    inB_oid++;

    BITONIC_MERGE4(outreg1, outreg2, outreg1_oid, outreg2_oid, regA, regB, regA_oid, regB_oid);

    /* store outreg1 */
    _mm256_store_pd((double *) outp, outreg1);
    _mm_store_si128((__m128i *) outp_oid, _mm256_cvtpd_epi32(outreg1_oid));
    outp++;
    outp_oid++;

    while( inA < endA && inB < endB ) {

        /* 3 Options : normal-if, cmove-with-assembly, sw-predication */
        IFELSECONDMOVE(next, next_oid, inA, inB, inA_oid, inB_oid, int64_t);

        regA = outreg2;
        regB = _mm256_load_pd((double const *) next);

        regA_oid = outreg2_oid;
        regB_oid = _mm256_cvtepi32_pd(_mm_load_si128((__m128i const*) next_oid));

        BITONIC_MERGE4(outreg1, outreg2, outreg1_oid, outreg2_oid, regA, regB, regA_oid, regB_oid);

        /* store outreg1 */
        _mm256_store_pd((double *) outp, outreg1);
        _mm_store_si128((__m128i *) outp_oid, _mm256_cvtpd_epi32(outreg1_oid));
        outp++;
        outp_oid++;
    }

    /* handle remaining items */
    while( inA < endA ) {
        __m256d regA = _mm256_load_pd((double const *) inA);
        __m256d regB = outreg2;

        __m256d regA_oid = _mm256_cvtepi32_pd(_mm_load_si128((__m128i const*) inA_oid));
        __m256d regB_oid = outreg2_oid;

        BITONIC_MERGE4(outreg1, outreg2, outreg1_oid, outreg2_oid, regA, regB, regA_oid, regB_oid);

        _mm256_store_pd((double *) outp, outreg1);
        _mm_store_si128((__m128i *) outp_oid, _mm256_cvtpd_epi32(outreg1_oid));
        inA++;
        inA_oid++;
        outp++;
        outp_oid++;
    }

    while( inB < endB ) {
        __m256d regA = outreg2;
        __m256d regB = _mm256_load_pd((double const *) inB);

        __m256d regA_oid = outreg2_oid;
        __m256d regB_oid = _mm256_cvtepi32_pd(_mm_load_si128((__m128i const*) inB_oid));

        BITONIC_MERGE4(outreg1, outreg2, outreg1_oid, outreg2_oid, regA, regB, regA_oid, regB_oid);

        _mm256_store_pd((double *) outp, outreg1);
        _mm_store_si128((__m128i *) outp_oid, _mm256_cvtpd_epi32(outreg1_oid));
        inB++;
        inB_oid++;
        outp++;
        outp_oid++;
    }

    /* store the last remaining register values */
    _mm256_store_pd((double *) outp, outreg2);
    _mm_store_si128((__m128i *) outp_oid, _mm256_cvtpd_epi32(outreg2_oid));

}

inline void __attribute__((always_inline))
Mergesort::merge16_int16_eqlen_aligned(int16_t * const inpA, int16_t * const inpB,
			surrogate_t * const inpA_oid, surrogate_t * const inpB_oid,
			int16_t * const out, surrogate_t * const out_oid, const uint32_t len)
{
	assert(0 == ((uint64_t)inpA & 31));
	assert(0 == ((uint64_t)inpB & 31));
	assert(0 == ((uint64_t)inpA_oid & 31));
	assert(0 == ((uint64_t)inpB_oid & 31));
	assert(0 == ((uint64_t)out & 31));
	assert(0 == ((uint64_t)out_oid & 31));

    register block16_16 * inA  = (block16_16 *) inpA;
    register block16_16 * inB  = (block16_16 *) inpB;
    block16_16 * const    endA = (block16_16 *) (inpA + len);
    block16_16 * const    endB = (block16_16 *) (inpB + len);

    register block16_32 * inA_oid  = (block16_32 *) inpA_oid;
    register block16_32 * inB_oid  = (block16_32 *) inpB_oid;

    block16_16 * outp = (block16_16 *) out;
    block16_32 * outp_oid = (block16_32 *) out_oid;

    register block16_16 * next = inB;
    register block16_32 * next_oid = inB_oid;

    register __m256i outregl;
    register __m256i outregh;

    register __m256i outregl1_oid, outregl2_oid;
    register __m256i outregh1_oid, outregh2_oid;

    register __m256i regA = _mm256_load_si256((__m256i const *) inA);
    register __m256i regB = _mm256_load_si256((__m256i const *) next);

    register __m256i regA1_oid = _mm256_load_si256((__m256i const *) inA_oid);
    register __m256i regA2_oid = _mm256_load_si256((__m256i const *) ((block8_32 *)inA_oid + 1));
    register __m256i regB1_oid = _mm256_load_si256((__m256i const *) next_oid);
    register __m256i regB2_oid = _mm256_load_si256((__m256i const *) ((block8_32 *)next_oid + 1));

    inA++;
    inB++;

    inA_oid++;
    inB_oid++;

    BITONIC_MERGE16_16(outregl, outregh,
                    outregl1_oid, outregl2_oid, outregh1_oid, outregh2_oid,
                    regA, regB, regA1_oid, regA2_oid, regB2_oid, regB1_oid);

    /* store outreg1 */
    _mm256_store_si256((__m256i *) outp, outregl);
    _mm256_store_si256((__m256i *) outp_oid, outregl1_oid);
    _mm256_store_si256((__m256i *) ((block8_32 *)outp_oid + 1), outregl2_oid);
    outp++;
    outp_oid++;

    while( inA < endA && inB < endB ) {

        /* 3 Options : normal-if, cmove-with-assembly, sw-predication */
        IFELSECONDMOVE(next, next_oid, inA, inB, inA_oid, inB_oid, int16_t);

        regA = outregh;
        regB = _mm256_load_si256((__m256i const *) next);

        regA1_oid = outregh1_oid;
        regA2_oid = outregh2_oid;
        regB1_oid = _mm256_load_si256((__m256i const *) next_oid);
        regB2_oid = _mm256_load_si256((__m256i const *) ((block8_32 *)next_oid + 1));

        BITONIC_MERGE16_16(outregl, outregh,
                        outregl1_oid, outregl2_oid, outregh1_oid, outregh2_oid,
                        regA, regB, regA1_oid, regA2_oid, regB2_oid, regB1_oid);

        /* store outreg1 */
        _mm256_store_si256((__m256i *) outp, outregl);
        _mm256_store_si256((__m256i *) outp_oid, outregl1_oid);
        _mm256_store_si256((__m256i *) ((block8_32 *)outp_oid + 1), outregl2_oid);
        outp++;
        outp_oid++;
    }

    /* handle remaining items */
    while( inA < endA ) {
        __m256i regA = _mm256_load_si256((__m256i const *) inA);
        __m256i regB = outregh;

        __m256i regA1_oid = _mm256_load_si256((__m256i const *) inA_oid);
        __m256i regA2_oid = _mm256_load_si256((__m256i const *) ((block8_32 *)inA_oid + 1));
        __m256i regB1_oid = outregh1_oid;
        __m256i regB2_oid = outregh2_oid;

        BITONIC_MERGE16_16(outregl, outregh,
                        outregl1_oid, outregl2_oid, outregh1_oid, outregh2_oid,
                        regA, regB, regA1_oid, regA2_oid, regB2_oid, regB1_oid);

        _mm256_store_si256((__m256i *) outp, outregl);
        _mm256_store_si256((__m256i *) outp_oid, outregl1_oid);
        _mm256_store_si256((__m256i *) ((block8_32 *)outp_oid + 1), outregl2_oid);
        inA++;
        inA_oid++;
        outp++;
        outp_oid++;
    }

    while( inB < endB ) {
        __m256i regA = outregh;
        __m256i regB = _mm256_load_si256((__m256i const *) inB);

        __m256i regA1_oid = outregh1_oid;
        __m256i regA2_oid = outregh2_oid;
        __m256i regB1_oid = _mm256_load_si256((__m256i const *) inB_oid);
        __m256i regB2_oid = _mm256_load_si256((__m256i const *) ((block8_32 *)inB_oid + 1));

        BITONIC_MERGE16_16(outregl, outregh,
                        outregl1_oid, outregl2_oid, outregh1_oid, outregh2_oid,
                        regA, regB, regA1_oid, regA2_oid, regB2_oid, regB1_oid);

        _mm256_store_si256((__m256i *) outp, outregl);
        _mm256_store_si256((__m256i *) outp_oid, outregl1_oid);
        _mm256_store_si256((__m256i *) ((block8_32 *)outp_oid + 1), outregl2_oid);
        inB++;
        inB_oid++;
        outp++;
        outp_oid++;
    }

    /* store the last remaining register values */
    _mm256_store_si256((__m256i *) outp, outregh);
    _mm256_store_si256((__m256i *) outp_oid, outregh1_oid);
    _mm256_store_si256((__m256i *) ((block8_32 *)outp_oid + 1), outregh2_oid);
}

inline void __attribute__((always_inline))
Mergesort::merge8_int32_eqlen_aligned(int32_t * const inpA, int32_t * const inpB,
						surrogate_t * const inpA_oid, surrogate_t * const inpB_oid,
						int32_t * const out, surrogate_t * const out_oid, const uint32_t len)
{
	assert(0 == ((uint64_t)inpA & 31));
	assert(0 == ((uint64_t)inpB & 31));
	assert(0 == ((uint64_t)inpA_oid & 31));
	assert(0 == ((uint64_t)inpB_oid & 31));
	assert(0 == ((uint64_t)out & 31));
	assert(0 == ((uint64_t)out_oid & 31));

    register block8_32 * inA  = (block8_32 *) inpA;
    register block8_32 * inB  = (block8_32 *) inpB;
    block8_32 * const    endA = (block8_32 *) (inpA + len);
    block8_32 * const    endB = (block8_32 *) (inpB + len);

    register block8_32 * inA_oid  = (block8_32 *) inpA_oid;
    register block8_32 * inB_oid  = (block8_32 *) inpB_oid;

    block8_32 * outp = (block8_32 *) out;
    block8_32 * outp_oid = (block8_32 *) out_oid;

    register block8_32 * next = inB;
    register block8_32 * next_oid = inB_oid;

    register __m256i outreg1;
    register __m256i outreg2;

    register __m256i outreg1_oid;
    register __m256i outreg2_oid;

    register __m256i regA = _mm256_load_si256((__m256i const *) inA);
    register __m256i regB = _mm256_load_si256((__m256i const *) next);

    register __m256i regA_oid = _mm256_load_si256((__m256i const *) inA_oid);
    register __m256i regB_oid = _mm256_load_si256((__m256i const *) next_oid);

    inA++;
    inB++;

    inA_oid++;
    inB_oid++;

    BITONIC_MERGE8_32(outreg1, outreg2, outreg1_oid, outreg2_oid,
    		regA, regB, regA_oid, regB_oid);

    /* store outreg1 */
    _mm256_store_si256((__m256i *) outp, outreg1);
    _mm256_store_si256((__m256i *) outp_oid, outreg1_oid);
    outp++;
    outp_oid++;

    while( inA < endA && inB < endB ) {

        /* 3 Options : normal-if, cmove-with-assembly, sw-predication */
        IFELSECONDMOVE(next, next_oid, inA, inB, inA_oid, inB_oid, int32_t);

        regA = outreg2;
        regB = _mm256_load_si256((__m256i const *) next);

        regA_oid = outreg2_oid;
        regB_oid = _mm256_load_si256((__m256i const *) next_oid);

        BITONIC_MERGE8_32(outreg1, outreg2, outreg1_oid, outreg2_oid,
        		regA, regB, regA_oid, regB_oid);

        /* store outreg1 */
        _mm256_store_si256((__m256i *) outp, outreg1);
        _mm256_store_si256((__m256i *) outp_oid, outreg1_oid);
        outp++;
        outp_oid++;
    }

    /* handle remaining items */
    while( inA < endA ) {
        __m256i regA = _mm256_load_si256((__m256i const *) inA);
        __m256i regB = outreg2;

        __m256i regA_oid = _mm256_load_si256((__m256i const *) inA_oid);
        __m256i regB_oid = outreg2_oid;

        BITONIC_MERGE8_32(outreg1, outreg2, outreg1_oid, outreg2_oid,
        		regA, regB, regA_oid, regB_oid);

        _mm256_store_si256((__m256i *) outp, outreg1);
        _mm256_store_si256((__m256i *) outp_oid, outreg1_oid);
        inA++;
        inA_oid++;
        outp++;
        outp_oid++;
    }

    while( inB < endB ) {
        __m256i regA = outreg2;
        __m256i regB = _mm256_load_si256((__m256i const *) inB);

        __m256i regA_oid = outreg2_oid;
        __m256i regB_oid = _mm256_load_si256((__m256i const *) inB_oid);

        BITONIC_MERGE8_32(outreg1, outreg2, outreg1_oid, outreg2_oid,
        		regA, regB, regA_oid, regB_oid);

        _mm256_store_si256((__m256i *) outp, outreg1);
        _mm256_store_si256((__m256i *) outp_oid, outreg1_oid);
        inB++;
        inB_oid++;
        outp++;
        outp_oid++;
    }

    /* store the last remaining register values */
    _mm256_store_si256((__m256i *) outp, outreg2);
    _mm256_store_si256((__m256i *) outp_oid, outreg2_oid);

}

/**
 * Sorts the last chunk of the input, which is less than BLOCKSIZE tuples.
 * @note This function assumes a hard-coded BLOCKSIZE of 16384 and nitems must
 * be less than 16384.
 *
 * @param inputptr
 * @param outputptr
 * @param nitems
 */
inline void __attribute__((always_inline))
//void
Mergesort::avxmergesort_rem_uint16_aligned(uint16_t** inputptr_val, surrogate_t** inputptr_oid,
		uint16_t** outputptr_val, surrogate_t** outputptr_oid, uint32_t nitems)
{
#if 1
    int16_t * inp_val = (int16_t *) *inputptr_val;
    int16_t * out_val = (int16_t *) *outputptr_val;

    surrogate_t * inp_oid = *inputptr_oid;
    surrogate_t * out_oid = *outputptr_oid;

    /* each chunk keeps track of its temporary memory offset */
    int16_t *ptrs_val[8][2];/* [chunk-in, chunk-out-tmp] */
    surrogate_t *ptrs_oid[8][2];

    uint32_t n = nitems, pos = 0, i = 0;
    uint32_t nxtpow = 8192;/* TODO: infer from nitems, nearest pow2 to nitems */
    uint32_t sizes[8];//XU: originally the constant is 6

    while(n < nxtpow) {
        nxtpow >>= 1;
    }

    while(nxtpow > 128) {
        ptrs_val[i][0] = inp_val + pos;
        ptrs_val[i][1] = out_val + pos;
        ptrs_oid[i][0] = inp_oid + pos;
        ptrs_oid[i][1] = out_oid + pos;
        sizes[i]   = nxtpow;

        avxmergesort_block_uint16_aligned((uint16_t **)&ptrs_val[i][0], &ptrs_oid[i][0],
        				(uint16_t **)&ptrs_val[i][1], &ptrs_oid[i][1], nxtpow);
        assert(std::is_sorted(ptrs_val[i][1], ptrs_val[i][1] + nxtpow));
#if 0
        std::cout << "block " << i << " with output oids:" << std::endl;
        for (uint32_t idx = 0; idx < nxtpow; ++idx) {
        	//std::cout <<(ptrs_oid[i][1][idx]) << std::endl;
        }
		for (uint32_t idx = 0; idx < nxtpow; ++idx) {
			assert(ptrs_oid[i][1][idx] < nitems);
		}
#endif

        pos += nxtpow;
        n   -= nxtpow;
        swap<int16_t>(&ptrs_val[i][0], &ptrs_val[i][1]);
        swap<surrogate_t>(&ptrs_oid[i][0], &ptrs_oid[i][1]);
        i++;

        while(n < nxtpow) {
            nxtpow >>= 1;
        }
    }

    if(n > 0) {
        /* sort last n < 128 items using scalar sort */
        ptrs_val[i][0] = inp_val + pos;
        ptrs_val[i][1] = out_val + pos;
        ptrs_oid[i][0] = inp_oid + pos;
        ptrs_oid[i][1] = out_oid + pos;
        sizes[i]   = n;
        std::vector<element_t<int16_t> > comparedCol;	//here we convert to int16_t because the values are already flipped
    	element_t<int16_t> tmpElement;
    	uint32_t idx;
    	for (idx = 0; idx < n; ++idx) {
    		tmpElement.value 	= ptrs_val[i][0][idx];
    		tmpElement.oid		= ptrs_oid[i][0][idx];
    		comparedCol.push_back(tmpElement);
    	}

        std::sort(comparedCol.begin(), comparedCol.end(), mycomparison<int16_t>());

        for (idx = 0; idx < n; ++idx) {
        	ptrs_val[i][0][idx] = comparedCol[idx].value;
        	ptrs_oid[i][0][idx] = comparedCol[idx].oid;
        }

#if 0
        std::cout << "the block small than 128 with output oids:" << std::endl;
        for (uint32_t tmp = 0; tmp < n; ++tmp) {
        	std::cout <<(ptrs_oid[i][0][tmp]) << std::endl;
        }
		for (uint32_t tmp = 0; tmp < n; ++tmp) {
			assert(ptrs_oid[i][0][tmp] < nitems);
		}
#endif

		/* no need to swap */
		/* XU: why not swap here: try to swap...cannot swap...*/
        //swap<int32_t>(&ptrs_val[i][0], &ptrs_val[i][1]);
        //swap<surrogate_t>(&ptrs_oid[i][0], &ptrs_oid[i][1]);
        i++;
    }

    uint32_t nchunks = i;

    /* merge sorted blocks */
    while(nchunks > 1) {
        uint64_t k = 0;
        for(uint64_t j = 0; j < (nchunks-1); j += 2) {
            int16_t * inpA_val  = ptrs_val[j][0];
            int16_t * inpB_val  = ptrs_val[j+1][0];
            int16_t * out_val   = ptrs_val[j][1];
            surrogate_t * inpA_oid  = ptrs_oid[j][0];
            surrogate_t * inpB_oid  = ptrs_oid[j+1][0];
            surrogate_t * out_oid   = ptrs_oid[j][1];
            uint32_t  sizeA = sizes[j];
            uint32_t  sizeB = sizes[j+1];
#if 0
            std::cout << "merging block " << j << " and " << (j+1) <<
            		" with size " << sizeA << " and " << sizeB << std::endl;

            for (uint32_t idx = 0; idx < sizeA; ++idx) {
            	//std::cout << inpB_oid[idx] << std::endl;
            	assert(inpA_oid[idx] < nitems);
            }
            for (uint32_t idx = 0; idx < sizeB; ++idx) {
            	assert(inpB_oid[idx] < nitems);
            }
#endif
            merge16_int16_varlen_aligned(inpA_val, inpB_val, out_val,
            		inpA_oid, inpB_oid, out_oid, sizeA, sizeB);

#if 0
            std::cout << "the output oids after a merge:" << std::endl;
            for (uint32_t idx = 0; idx < sizeA+sizeB; ++idx) {
            	std::cout << "idx: " << idx << "\t" << out_oid[idx] << std::endl;
            	assert(out_oid[idx] < nitems);
            }
#endif

            /* setup new pointers */
            ptrs_val[k][0] = out_val;
            ptrs_val[k][1] = inpA_val;
            ptrs_oid[k][0] = out_oid;
            ptrs_oid[k][1] = inpA_oid;
            sizes[k]   = sizeA + sizeB;
            k++;
        }

        if((nchunks % 2)) {
            /* just move the pointers */
            ptrs_val[k][0] = ptrs_val[nchunks-1][0];
            ptrs_val[k][1] = ptrs_val[nchunks-1][1];
            ptrs_oid[k][0] = ptrs_oid[nchunks-1][0];
            ptrs_oid[k][1] = ptrs_oid[nchunks-1][1];
            sizes[k]   = sizes[nchunks-1];
            k++;
        }

        nchunks = k;
    }

    /* finally swap input/output pointers, where output holds the sorted list */
    * outputptr_val = (uint16_t *)ptrs_val[0][0];
    * inputptr_val  = (uint16_t *)ptrs_val[0][1];

    * outputptr_oid = ptrs_oid[0][0];
    * inputptr_oid  = ptrs_oid[0][1];

#if 0
            std::cout << "the output oids finally:" << std::endl;
            for (uint32_t idx = 0; idx < nitems; ++idx) {
            	std::cout << "idx: " << idx << "\t" << ptrs_oid[0][0][idx] << std::endl;
            	assert(ptrs_oid[0][0][idx] < nitems);
            }
#endif
#endif
}

/**
 * Sorts the last chunk of the input, which is less than BLOCKSIZE tuples.
 * @note This function assumes a hard-coded BLOCKSIZE of 16384 and nitems must
 * be less than 16384.
 *
 * @param inputptr
 * @param outputptr
 * @param nitems
 */
inline void __attribute__((always_inline))
//void
Mergesort::avxmergesort_rem_uint32_aligned(uint32_t** inputptr_val, surrogate_t** inputptr_oid,
		uint32_t** outputptr_val, surrogate_t** outputptr_oid, uint32_t nitems)
{
#if 0
    int32_t * inp_val = (int32_t *) *inputptr_val;
    int32_t * out_val = (int32_t *) *outputptr_val;

    surrogate_t * inp_oid = *inputptr_oid;
    surrogate_t * out_oid = *outputptr_oid;

    std::vector<element_t<uint32_t> > comparedCol;	//here we convert to int32_t because the values are already flipped
	element_t<uint32_t> tmpElement;
	uint32_t idx;
	for (idx = 0; idx < nitems; ++idx) {
		tmpElement.value 	= inp_val[idx];
		tmpElement.oid		= inp_oid[idx];
		comparedCol.push_back(tmpElement);
	}

    std::sort(comparedCol.begin(), comparedCol.end(), mycomparison<uint32_t>());

    for (idx = 0; idx < nitems; ++idx) {
    	out_val[idx] = comparedCol[idx].value;
    	out_oid[idx] = comparedCol[idx].oid;
    }
#endif

#if 1
    int32_t * inp_val = (int32_t *) *inputptr_val;
    int32_t * out_val = (int32_t *) *outputptr_val;

    surrogate_t * inp_oid = *inputptr_oid;
    surrogate_t * out_oid = *outputptr_oid;

    /* each chunk keeps track of its temporary memory offset */
    int32_t *ptrs_val[8][2];/* [chunk-in, chunk-out-tmp] */
    surrogate_t *ptrs_oid[8][2];

    uint32_t n = nitems, pos = 0, i = 0;
    uint32_t nxtpow = 8192;/* TODO: infer from nitems, nearest pow2 to nitems */
    uint32_t sizes[8];//XU: originally the constant is 6

    while(n < nxtpow) {
        nxtpow >>= 1;
    }

    while(nxtpow > 128) {
        ptrs_val[i][0] = inp_val + pos;
        ptrs_val[i][1] = out_val + pos;
        ptrs_oid[i][0] = inp_oid + pos;
        ptrs_oid[i][1] = out_oid + pos;
        sizes[i]   = nxtpow;

        avxmergesort_block_uint32_aligned((uint32_t **)&ptrs_val[i][0], &ptrs_oid[i][0],
        				(uint32_t **)&ptrs_val[i][1], &ptrs_oid[i][1], nxtpow);
        assert(std::is_sorted(ptrs_val[i][1], ptrs_val[i][1] + nxtpow));
#if 0
        std::cout << "block " << i << " with output oids:" << std::endl;
        for (uint32_t idx = 0; idx < nxtpow; ++idx) {
        	//std::cout <<(ptrs_oid[i][1][idx]) << std::endl;
        }
		for (uint32_t idx = 0; idx < nxtpow; ++idx) {
			assert(ptrs_oid[i][1][idx] < nitems);
		}
#endif

        pos += nxtpow;
        n   -= nxtpow;
        swap<int32_t>(&ptrs_val[i][0], &ptrs_val[i][1]);
        swap<surrogate_t>(&ptrs_oid[i][0], &ptrs_oid[i][1]);
        i++;

        while(n < nxtpow) {
            nxtpow >>= 1;
        }
    }

    if(n > 0) {
        /* sort last n < 128 items using scalar sort */
        ptrs_val[i][0] = inp_val + pos;
        ptrs_val[i][1] = out_val + pos;
        ptrs_oid[i][0] = inp_oid + pos;
        ptrs_oid[i][1] = out_oid + pos;
        sizes[i]   = n;
        std::vector<element_t<int32_t> > comparedCol;	//here we convert to int32_t because the values are already flipped
    	element_t<int32_t> tmpElement;
    	uint32_t idx;
    	for (idx = 0; idx < n; ++idx) {
    		tmpElement.value 	= ptrs_val[i][0][idx];
    		tmpElement.oid		= ptrs_oid[i][0][idx];
    		comparedCol.push_back(tmpElement);
    	}

        std::sort(comparedCol.begin(), comparedCol.end(), mycomparison<int32_t>());

        for (idx = 0; idx < n; ++idx) {
        	ptrs_val[i][0][idx] = comparedCol[idx].value;
        	ptrs_oid[i][0][idx] = comparedCol[idx].oid;
        }

#if 0
        std::cout << "the block small than 128 with output oids:" << std::endl;
        for (uint32_t tmp = 0; tmp < nxtpow; ++tmp) {
        	std::cout <<(ptrs_oid[i][0][tmp]) << std::endl;
        }
		for (uint32_t tmp = 0; tmp < n; ++tmp) {
			assert(ptrs_oid[i][0][tmp] < nitems);
		}
#endif

		/* no need to swap */
		/* XU: why not swap here: try to swap...cannot swap...*/
        //swap<int32_t>(&ptrs_val[i][0], &ptrs_val[i][1]);
        //swap<surrogate_t>(&ptrs_oid[i][0], &ptrs_oid[i][1]);
        i++;
    }

    uint32_t nchunks = i;

    /* merge sorted blocks */
    while(nchunks > 1) {
        uint64_t k = 0;
        for(uint64_t j = 0; j < (nchunks-1); j += 2) {
            int32_t * inpA_val  = ptrs_val[j][0];
            int32_t * inpB_val  = ptrs_val[j+1][0];
            int32_t * out_val   = ptrs_val[j][1];
            surrogate_t * inpA_oid  = ptrs_oid[j][0];
            surrogate_t * inpB_oid  = ptrs_oid[j+1][0];
            surrogate_t * out_oid   = ptrs_oid[j][1];
            uint32_t  sizeA = sizes[j];
            uint32_t  sizeB = sizes[j+1];
#if 0
            std::cout << "merging block " << j << " and " << (j+1) <<
            		" with size " << sizeA << " and " << sizeB << std::endl;

            for (uint32_t idx = 0; idx < sizeA; ++idx) {
            	//std::cout << inpB_oid[idx] << std::endl;
            	assert(inpA_oid[idx] < nitems);
            }
            for (uint32_t idx = 0; idx < sizeB; ++idx) {
            	assert(inpB_oid[idx] < nitems);
            }
#endif
            merge8_int32_varlen_aligned(inpA_val, inpB_val, out_val,
            		inpA_oid, inpB_oid, out_oid, sizeA, sizeB);

#if 0
            for (uint32_t idx = 0; idx < sizeA+sizeB; ++idx) {
            	std::cout << out_oid[idx] << std::endl;
            	assert(out_oid[idx] < nitems);
            }
#endif

            /* setup new pointers */
            ptrs_val[k][0] = out_val;
            ptrs_val[k][1] = inpA_val;
            ptrs_oid[k][0] = out_oid;
            ptrs_oid[k][1] = inpA_oid;
            sizes[k]   = sizeA + sizeB;
            k++;
        }

        if((nchunks % 2)) {
            /* just move the pointers */
            ptrs_val[k][0] = ptrs_val[nchunks-1][0];
            ptrs_val[k][1] = ptrs_val[nchunks-1][1];
            ptrs_oid[k][0] = ptrs_oid[nchunks-1][0];
            ptrs_oid[k][1] = ptrs_oid[nchunks-1][1];
            sizes[k]   = sizes[nchunks-1];
            k++;
        }

        nchunks = k;
    }

    /* finally swap input/output pointers, where output holds the sorted list */
    * outputptr_val = (uint32_t *)ptrs_val[0][0];
    * inputptr_val  = (uint32_t *)ptrs_val[0][1];

    * outputptr_oid = ptrs_oid[0][0];
    * inputptr_oid  = ptrs_oid[0][1];
#endif
}

/**
 *
 * Sorts the last chunk of the input, which is less than BLOCKSIZE tuples.
 * @note This function assumes a hard-coded BLOCKSIZE of 8192 and nitems must
 * be less than 8192. Also, assumes surrogate_t=uint32_t
 *
 * @param inputptr
 * @param outputptr
 * @param nitems
 */

inline void __attribute__((always_inline))
//void
Mergesort::avxmergesort_rem_uint64_aligned(uint64_t** inputptr_val, surrogate_t** inputptr_oid,
		uint64_t** outputptr_val, surrogate_t** outputptr_oid, uint32_t nitems)
{
#if 1
    int64_t * inp_val = (int64_t *) *inputptr_val;
    int64_t * out_val = (int64_t *) *outputptr_val;

    surrogate_t * inp_oid = *inputptr_oid;
    surrogate_t * out_oid = *outputptr_oid;

    /* each chunk keeps track of its temporary memory offset */
    int64_t *ptrs_val[8][2];/* [chunk-in, chunk-out-tmp] */
    surrogate_t *ptrs_oid[8][2];

    uint32_t n = nitems, pos = 0, i = 0;
    uint32_t nxtpow = 4096;/* TODO: infer from nitems, nearest pow2 to nitems */
    uint32_t sizes[8];//XU: originally the constant is 6

    while(n < nxtpow) {
        nxtpow >>= 1;
    }

    while(nxtpow > 128) {
        ptrs_val[i][0] = inp_val + pos;
        ptrs_val[i][1] = out_val + pos;
        ptrs_oid[i][0] = inp_oid + pos;
        ptrs_oid[i][1] = out_oid + pos;
        sizes[i]   = nxtpow;

        avxmergesort_block_uint64_aligned((uint64_t **)&ptrs_val[i][0], &ptrs_oid[i][0],
        				(uint64_t **)&ptrs_val[i][1], &ptrs_oid[i][1], nxtpow);
        assert(std::is_sorted(ptrs_val[i][1], ptrs_val[i][1] + nxtpow));
#if 0
        std::cout << "block " << i << " with output oids:" << std::endl;
        for (uint32_t idx = 0; idx < nxtpow; ++idx) {
        	//std::cout <<(ptrs_oid[i][1][idx]) << std::endl;
        }
		for (uint32_t idx = 0; idx < nxtpow; ++idx) {
			assert(ptrs_oid[i][1][idx] < nitems);
		}
#endif

        pos += nxtpow;
        n   -= nxtpow;
        swap<int64_t>(&ptrs_val[i][0], &ptrs_val[i][1]);
        swap<surrogate_t>(&ptrs_oid[i][0], &ptrs_oid[i][1]);
        i++;

        while(n < nxtpow) {
            nxtpow >>= 1;
        }
    }

    if(n > 0) {
        /* sort last n < 128 items using scalar sort */
        ptrs_val[i][0] = inp_val + pos;
        ptrs_val[i][1] = out_val + pos;
        ptrs_oid[i][0] = inp_oid + pos;
        ptrs_oid[i][1] = out_oid + pos;
        sizes[i]   = n;
        std::vector<element_t<int64_t> > comparedCol;	//here we convert to int64_t because the values are already flipped
    	element_t<int64_t> tmpElement;
    	uint32_t idx;
    	for (idx = 0; idx < n; ++idx) {
    		tmpElement.value 	= ptrs_val[i][0][idx];
    		tmpElement.oid		= ptrs_oid[i][0][idx];
    		comparedCol.push_back(tmpElement);
    	}

        std::sort(comparedCol.begin(), comparedCol.end(), mycomparison<int64_t>());

        for (idx = 0; idx < n; ++idx) {
        	ptrs_val[i][0][idx] = comparedCol[idx].value;
        	ptrs_oid[i][0][idx] = comparedCol[idx].oid;
        }

#if 0
        std::cout << "the block small than 128 with output oids:" << std::endl;
        for (uint32_t tmp = 0; tmp < nxtpow; ++tmp) {
        	std::cout <<(ptrs_oid[i][0][tmp]) << std::endl;
        }
		for (uint32_t tmp = 0; tmp < n; ++tmp) {
			assert(ptrs_oid[i][0][tmp] < nitems);
		}
#endif
		/* no need to swap */
		/* XU: why not swap here: try to swap...cannot swap*/
        //swap<int64_t>(&ptrs_val[i][0], &ptrs_val[i][1]);
        //swap<surrogate_t>(&ptrs_oid[i][0], &ptrs_oid[i][1]);

        i++;
    }

    uint32_t nchunks = i;

    /* merge sorted blocks */
    while(nchunks > 1) {
        uint64_t k = 0;
        for(uint64_t j = 0; j < (nchunks-1); j += 2) {
            int64_t * inpA_val  = ptrs_val[j][0];
            int64_t * inpB_val  = ptrs_val[j+1][0];
            int64_t * out_val   = ptrs_val[j][1];
            surrogate_t * inpA_oid  = ptrs_oid[j][0];
            surrogate_t * inpB_oid  = ptrs_oid[j+1][0];
            surrogate_t * out_oid   = ptrs_oid[j][1];
            uint32_t  sizeA = sizes[j];
            uint32_t  sizeB = sizes[j+1];
#if 0
            std::cout << "merging block " << j << " and " << (j+1) <<
            		" with size " << sizeA << " and " << sizeB << std::endl;

            for (uint32_t idx = 0; idx < sizeA; ++idx) {
            	//std::cout << inpB_oid[idx] << std::endl;
            	assert(inpA_oid[idx] < nitems);
            }
            for (uint32_t idx = 0; idx < sizeB; ++idx) {
            	assert(inpB_oid[idx] < nitems);
            }
#endif
            merge8_int64_varlen_aligned(inpA_val, inpB_val, out_val,
            		inpA_oid, inpB_oid, out_oid, sizeA, sizeB);

#if 0
            for (uint32_t idx = 0; idx < sizeA+sizeB; ++idx) {
            	std::cout << out_oid[idx] << std::endl;
            	assert(out_oid[idx] < nitems);
            }
#endif

            /* setup new pointers */
            ptrs_val[k][0] = out_val;
            ptrs_val[k][1] = inpA_val;
            ptrs_oid[k][0] = out_oid;
            ptrs_oid[k][1] = inpA_oid;
            sizes[k]   = sizeA + sizeB;
            k++;
        }

        if((nchunks % 2)) {
            /* just move the pointers */
            ptrs_val[k][0] = ptrs_val[nchunks-1][0];
            ptrs_val[k][1] = ptrs_val[nchunks-1][1];
            ptrs_oid[k][0] = ptrs_oid[nchunks-1][0];
            ptrs_oid[k][1] = ptrs_oid[nchunks-1][1];
            sizes[k]   = sizes[nchunks-1];
            k++;
        }

        nchunks = k;
    }

    /* finally swap input/output pointers, where output holds the sorted list */
    * outputptr_val = (uint64_t *)ptrs_val[0][0];
    * inputptr_val  = (uint64_t *)ptrs_val[0][1];

    * outputptr_oid = ptrs_oid[0][0];
    * inputptr_oid  = ptrs_oid[0][1];
#endif
}
}
#endif
