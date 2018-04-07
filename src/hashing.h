#ifndef HASHING_H
#define HASHING_H

#include "types.h"
#include "avxcommon.h"

namespace multiAttrSort{

typedef struct {
	//int type;		/* type of index entity */
	//int width;		/* width of hash entries */
	surrogate_t nil;		/* nil representation */
	surrogate_t lim;		/* collision list size, equivalently, tuple # */
	surrogate_t mask;	/* number of hash buckets-1 (power of 2) */
	surrogate_t *Hash;		/* hash table */
	surrogate_t *Link;		/* collision list */
} Hash_t;

//#define ENTRY_NONE ((size_t) 0xFFFFFFFFFFFFFFFF)	//
#define ENTRY_NONE ((surrogate_t) 0xFFFFFFFF)	//
#define mix_bte(X)	((unsigned int) (X))
#define mix_sht(X)	((unsigned int) (X))
#define mix_int(X)	(((X)>>7)^((X)>>13)^((X)>>21)^(X))
#define mix_lng(X)	mix_int((unsigned int) ((X) ^ ((X) >> 32)))

class Hashing{
public:
	Hashing() {}

	/**
	 * find largest power of 2 smaller than or equal to cnt
	 */
	inline static size_t __attribute__((always_inline))
	HASHmask(size_t cnt);

	/**
	 * do the grouping by hashing (with group info from previous round)
	 * NOTE: used by hashtype = 0, i.e., not re-ordering based
	 */
	template <class T>
	inline static void __attribute__((always_inline))
	subgroup_non_first_round(surrogate_t *outGroup, surrogate_t *outGroupNum,
			T *inValues, surrogate_t *inGroup, surrogate_t *inGroupNum, Hash_t* hashPtr);

	/**
	 * do the grouping by hashing (w/o group info from previous round)
	 * NOTE: used by hashtype = 0, i.e., not re-ordering based
	 */
	template <class T>
	inline static void __attribute__((always_inline))
	subgroup_first_round(surrogate_t *outGroup, surrogate_t *outGroupNum, T *inValues, Hash_t* hashPtr);

	/**
	 * NOTE: used by hashtype = 1. i.e., re-ordering based
	 *
	 * INPUT:
	 * columnPtr: its column values to be grouped
	 * start: starting index of values w.r.t. columnPtr->values
	 * nitems: # of items to be grouped
	 * startGrpIdx: initial grouping index
	 *
	 *
	 * OUTPUT:
	 * intermediate_group_idx: 	the group index of target range of values
	 * outGroupNum: # of groups for the range of values
	 *
	 *
	 */
	template <class T>
	inline static void __attribute__((always_inline))
	subgroup_range(surrogate_t *intermediate_group_idx, surrogate_t *outGroupNum,
			BAT_t<T> *columnPtr, surrogate_t start, surrogate_t nitems,
			surrogate_t startGrpIdx, Hash_t* hashPtr, surrogate_t* hashed_per_block);

	template <class T>
	inline static void __attribute__((always_inline))
	subgroup_simd(surrogate_t *outGroup, T *inValues, Hash_t* hashPtr);

	template <class T>
	inline static surrogate_t __attribute__((always_inline))
	hash_any(Hash_t *hashPtr, T *value);

	/**
	 *
	 * OUTPUT:
	 * hashed_per_block: hashed values
	 */
	template <class T>
	inline static void __attribute__((always_inline))
	calculate_hash_in_block(surrogate_t* hashed_per_block, Hash_t *hashPtr,
			const T *inValues, const surrogate_t start_idx, const surrogate_t block_szs);

};


template <class T>
inline void __attribute__((always_inline))
Hashing::subgroup_non_first_round(surrogate_t *outGroup, surrogate_t *outGroupNum,
		T *inValues, surrogate_t *inGroup, surrogate_t *inGroupNum, Hash_t* hashPtr)

{
	surrogate_t value_idx;
	surrogate_t hashed;
	surrogate_t oid;
	surrogate_t group_idx = 0;	//group index (assume the largest group # is 2^32)

	uint32_t inGroup_valid_bits = ceil(log2(*inGroupNum));

	//get bits info
	uint32_t bits = log2(hashPtr->mask + 1) - inGroup_valid_bits;
	/*
	uint32_t bits = 0;
	if (sizeof(T) <= 2) {
		bits = 16;
		//bits = 0;
	} else {
		//bits = 0;
		bits = 8;
	}
	*/

	for (value_idx = 0; value_idx < hashPtr->lim; ++value_idx) {

		/**hashing the value**/
		hashed = hash_any<T>(hashPtr, &inValues[value_idx]);

		hashed = (hashed ^ ((size_t)inGroup[value_idx] << bits)) & hashPtr->mask;

		/**build the hash table on the fly**/
		for (oid = hashPtr->Hash[hashed];
				oid != ENTRY_NONE;
				oid = hashPtr->Link[oid]) {
			if ((inValues[value_idx] == inValues[oid]) &&
					(inGroup[oid] == inGroup[value_idx])) {	//the group exists
				outGroup[value_idx] = outGroup[oid];
				break;
			}
		}
		if (oid == ENTRY_NONE) {
			outGroup[value_idx] = group_idx;
			group_idx++;
			hashPtr->Link[value_idx] = hashPtr->Hash[hashed];
			hashPtr->Hash[hashed] = value_idx;
		}
	}

	*outGroupNum = group_idx;
}

#if 0
template <class T>
inline void __attribute__((always_inline))
subgroup_simd(surrogate_t *outGroup, T *inValues, Hash_t* hashPtr)
{
	size_t value_idx;
	size_t hashed;
	size_t oid;
	surrogate_t group_idx = 0;	//group index (assume the largest group # is 2^32)
	const size_t block_size = L2_CACHE_SIZE / (2 * sizeof(T));

	assert(0 == (hashPtr->lim % block_size));	//no incomplete block

	const size_t block_num = hashPtr->lim / block_size;

	size_t* hashed_per_block = (size_t *) malloc_aligned(block_size * sizeof(size_t));

	for (size_t block_idx = 0; block_idx < block_num; ++block_idx) {

		/**calculate hashing value with SIMD**/
		calculate_hash_in_block(hashed_per_block, hashPtr, inValues, block_idx, block_size);

	}

	free(hashed_per_block);
}
#endif

template <class T>
inline void __attribute__((always_inline))
Hashing::calculate_hash_in_block(surrogate_t* hashed_per_block, Hash_t *hashPtr,
		const T *inValues, const surrogate_t start_idx, const surrogate_t block_sz)
{
	surrogate_t idx = start_idx;
	const surrogate_t end_idx = start_idx + block_sz;
	surrogate_t base_idx = 0;
	surrogate_t bank_num = 0;
	switch(sizeof(T)) {
		case 1:
		case 2:
			for (; idx < end_idx; ++idx, ++base_idx) {
				hashed_per_block[base_idx] = static_cast<surrogate_t>(inValues[idx]) & hashPtr->mask;
			}
			break;
		case 4:	//uint32_t
			bank_num =  8;
			//assert(0 == (block_size % bank_num));
			while ((idx + bank_num - 1) < end_idx) {
				__m256i ra = _mm256_loadu_si256 (reinterpret_cast<__m256i const *>(inValues + idx));
				__m256i rb = _mm256_srli_epi32 (ra, 7);
				__m256i rc = _mm256_srli_epi32 (ra, 13);
				__m256i rd = _mm256_srli_epi32 (ra, 21);

				__m256i re = _mm256_xor_si256 (ra, rb);
				__m256i rf = _mm256_xor_si256 (rc, rd);

				__m256i rg = _mm256_xor_si256 (re, rf);

				//note: we need to do the masking before store back
				//TODO: whether the order to extracting and storing is correct????
				//__m128i low = _mm256_extracti128_si256(rg, 0);
				//__m128i high = _mm256_extracti128_si256(rg, 1);
				//__m256i low_to_store = _mm256_cvtepu32_epi64(low);
				//__m256i high_to_store = _mm256_cvtepu32_epi64(high);

				__m256i mask_vec = _mm256_set1_epi32 (hashPtr->mask);

				_mm256_storeu_si256(reinterpret_cast<__m256i *> (hashed_per_block + base_idx), _mm256_and_si256(rg, mask_vec));

				//_mm256_store_si256(reinterpret_cast<__m256i *> (hashed_per_block + base_idx),
				//		_mm256_and_si256(low_to_store, mask_vec));
				//_mm256_store_si256(reinterpret_cast<__m256i *> (hashed_per_block + base_idx + 4),
				//		_mm256_and_si256(high_to_store, mask_vec));

				idx += bank_num;
				base_idx += bank_num;
			}

			/* Operate on remaining elements */
			while (idx < end_idx) {

				hashed_per_block[base_idx] = (surrogate_t) mix_int(*(const unsigned int *) (inValues + idx)) & hashPtr->mask;

				++idx;
				++base_idx;
			}
			break;
		case 8:	//uint64_t
			bank_num =  4;
			//assert(0 == (block_size % bank_num));
			while ((idx + bank_num - 1) < end_idx) {
				__m256i ra = _mm256_loadu_si256(reinterpret_cast<__m256i const *>(inValues + idx));
				__m256i rb = _mm256_srli_epi64(ra, 32);

				__m256i rc = _mm256_xor_si256(ra, rb);

				__m256i zero = _mm256_setzero_si256();

				__m256i target = _mm256_blend_epi32(rc, zero, 0xAA);

				__m256i target_ra = _mm256_srli_epi32 (target, 7);
				__m256i target_rb = _mm256_srli_epi32 (target, 13);
				__m256i target_rc = _mm256_srli_epi32 (target, 21);

				__m256i target_re = _mm256_xor_si256 (target, target_ra);
				__m256i target_rf = _mm256_xor_si256 (target_rb, target_rc);

				__m256i rg = _mm256_xor_si256 (target_re, target_rf);

				__m256i mask_vec = _mm256_set1_epi64x (hashPtr->mask);

				__m256i result = _mm256_and_si256(rg, mask_vec);

				//store back; currently use sequential store; any better solution???
				//TODO: is the storing order correct?
				hashed_per_block[base_idx] = (surrogate_t)_mm256_extract_epi64 (result, 0);
				hashed_per_block[base_idx + 1] = (surrogate_t)_mm256_extract_epi64 (result, 1);
				hashed_per_block[base_idx + 2] = (surrogate_t)_mm256_extract_epi64 (result, 2);
				hashed_per_block[base_idx + 3] = (surrogate_t)_mm256_extract_epi64 (result, 3);

				//_mm256_store_si256(reinterpret_cast<__m256i *> (hashed_per_block + base_idx), _mm256_and_si256(rg, mask_vec));

				idx += bank_num;
				base_idx += bank_num;
			}

			/* Operate on remaining elements */
			while (idx < end_idx) {

				hashed_per_block[base_idx] = (surrogate_t) mix_lng(
						*(const unsigned long long *) (inValues + idx)) & hashPtr->mask;

				++idx;
				++base_idx;
			}
			break;
	}
}

template <class T>
inline void __attribute__((always_inline))
Hashing::subgroup_range(surrogate_t *intermediate_group_idx, surrogate_t *outGroupNum,
		BAT_t<T> *columnPtr, surrogate_t start, surrogate_t nitems,
		surrogate_t startGrpIdx, Hash_t* hashPtr, surrogate_t* hashed_per_block)
{
	assert(nitems > 0);
	//std::cout << "nitems: " << nitems << std::endl;

	//different from sorting, we still need to do something (i.e., generate group index) when nitem == 1
	if (nitems == 1) {
		//std::cout << "with nitems ==1 " << std::endl;
		*outGroupNum = 1;
		intermediate_group_idx[start] = startGrpIdx;
		return;
	}

	//size_t new_mask = max<size_t>(Hashing::HASHmask(nitems), 1 << 16);
	size_t new_mask = max<size_t>(Hashing::HASHmask(nitems), 1 << 8);

	hashPtr->mask = static_cast<surrogate_t>(new_mask - 1);

	//clear the hash table
	memset((void*)hashPtr->Hash, 0xFF, (hashPtr->mask + 1) * sizeof(surrogate_t));


	//std::cout << "start subgrouping" << std::endl;
	surrogate_t oid;								//oid starts from the range
	surrogate_t group_idx = startGrpIdx;	//group index (assume the largest group # is 2^32)
	T *inValues = columnPtr->values;

	const size_t block_size = L2_CACHE_SIZE / (2 * sizeof(T));
	const size_t block_num = ceil((double)nitems / block_size);	//number of blocks
	const size_t rem_sz = nitems % block_size;		//number of tuples in the final remaining (incomplete) block

	//std::cout << "block_num: " << block_num << std::endl;

	/** Deal with blocks one by one**/
	assert((start + nitems) <= hashPtr->lim);


	surrogate_t hashed_idx; 	//index in *hashed_per_block*
	surrogate_t hashed_val = 0;

	surrogate_t each_block_sz = 0;
	surrogate_t each_block_start = 0;
	surrogate_t each_block_start_within_range = 0;

	for (size_t block_idx = 0; block_idx < block_num; ++block_idx) {

#if 0	//debug
		if (nitems < 500) {
			std::cout << "block_idx: " << block_idx << std::endl;
		}
#endif
		/** calculate hashing value with SIMD for a block**/
		each_block_sz = (rem_sz == 0) ? block_size : rem_sz;
		each_block_start = start + block_idx*block_size;
		assert((start + block_idx*block_size + each_block_sz) <= hashPtr->lim);
		calculate_hash_in_block(hashed_per_block, hashPtr, inValues, each_block_start, each_block_sz);


		each_block_start_within_range = block_idx*block_size;

		/** Calculating the group index for block, by building the hash table on the fly.
		 * TODO: use gathering SIMD intrinsics **/
		for (hashed_idx = 0; hashed_idx < each_block_sz; hashed_idx++) {

#if 0	//debug
		if (nitems < 500) {
			std::cout << "hashed_idx: " << hashed_idx << std::endl;
		}
#endif

			hashed_val = hashed_per_block[hashed_idx];

			for (oid = hashPtr->Hash[hashed_val];
					oid != ENTRY_NONE;
					oid = hashPtr->Link[oid]) {
				if (inValues[each_block_start + hashed_idx] == inValues[start + oid]) {	//the group exists, BE CAREFUL WITH THE IDs
					intermediate_group_idx[each_block_start + hashed_idx] = intermediate_group_idx[start + oid];
					break;
				}
			}
			if (oid == ENTRY_NONE) {
				intermediate_group_idx[each_block_start + hashed_idx] = group_idx;
				group_idx++;
				hashPtr->Link[each_block_start_within_range + hashed_idx] = hashPtr->Hash[hashed_val];
				hashPtr->Hash[hashed_val] = each_block_start_within_range + hashed_idx;	//be careful with the id to be recorded
			}
		}

	}

	//std::cout << "group_idx: " << group_idx << std::endl;
	*outGroupNum = (group_idx - startGrpIdx);
}

template <class T>
inline void __attribute__((always_inline))
Hashing::subgroup_first_round(surrogate_t *outGroup, surrogate_t *outGroupNum, T *inValues, Hash_t* hashPtr)
{
	//std::cout << "start subgrouping" << std::endl;
	surrogate_t value_idx;
	surrogate_t hashed;
	surrogate_t oid;
	surrogate_t group_idx = 0;	//group index (assume the largest group # is 2^32)

	for (value_idx = 0; value_idx < hashPtr->lim; ++value_idx) {

		/**hashing the value**/
		hashed = hash_any<T>(hashPtr, &inValues[value_idx]);

		assert(hashed <= hashPtr->mask);

		/**build the hash table on the fly**/
		for (oid = hashPtr->Hash[hashed];
				oid != ENTRY_NONE;
				oid = hashPtr->Link[oid]) {
			if (inValues[value_idx] == inValues[oid]) {	//the group exists
				outGroup[value_idx] = outGroup[oid];
				break;
			}
		}
		if (oid == ENTRY_NONE) {
			outGroup[value_idx] = group_idx;
			group_idx++;
			hashPtr->Link[value_idx] = hashPtr->Hash[hashed];
			hashPtr->Hash[hashed] = value_idx;
		}
	}

	*outGroupNum = group_idx;
}

template <class T>
inline surrogate_t __attribute__((always_inline))
Hashing::hash_any(Hash_t *hashPtr, T *valuePtr)
{
	switch(sizeof(T)) {
		case 1:	//byte
			assert((hashPtr->mask & 0xFF) == 0xFF);
			return (surrogate_t) mix_bte(*(const unsigned char*) (valuePtr));
		case 2:
			assert((hashPtr->mask & 0xFFFF) == 0xFFFF);
			return (surrogate_t) mix_sht(*(const unsigned short*) (valuePtr));
		case 4:
			return (surrogate_t) mix_int(*(const unsigned int *) (valuePtr)) & hashPtr->mask;
		case 8:
			return (surrogate_t) mix_lng(*(const unsigned long long *) (valuePtr)) & hashPtr->mask;
	}
}

inline size_t __attribute__((always_inline))
Hashing::HASHmask(size_t cnt)
{
	size_t m = cnt;

	/* find largest power of 2 smaller than or equal to cnt */
	m |= m >> 1;
	m |= m >> 2;
	m |= m >> 4;
	m |= m >> 8;
	m |= m >> 16;
#if SIZEOF_BUN == 8
	m |= m >> 32;
#endif
	m -= m >> 1;

	/* if cnt is more than 1/3 into the gap between m and 2*m,
	   double m */
	if (m + m - cnt < 2 * (cnt - m))
		m += m;
	//if (m < BATTINY)
	//	m = BATTINY;
	return m;
}


}

#endif  //HASHING_H
