#ifndef PARTITION_H
#define PARTITION_H

#include "types.h"
#include "avxcommon.h"

namespace multiAttrSort{

class Partition{
public:
	Partition() {}

	/**
	 * NOTE: used by hashtype = 1. i.e., re-ordering based
	 * all operations happen within the range of [start, start + nitems)
	 *
	 * INPUT:
	 * columnPtr: contain the oid list (before reordering). column values are not used in this function.
	 * intermediate_group_idx: the group index to be clustered (re-ordered)
	 * start: starting index of values w.r.t. columnPtr->values
	 * nitems: # of items to be grouped
	 *
	 *
	 *
	 * OUTPUT:
	 * outGroup: the CLUSTERED group information
	 * outOrder: the reordered oids
	 *
	 * AUXILIARY:
	 * histogram: to compute the offset of each partition
	 * startGrpIdx: smallest value of *intermediate_group_idx* within this range
	 * curGrpNum: num of groups within the range
	 *
	 */
	template <class T>
	inline static void __attribute__((always_inline))
	histogram_based_reordering(surrogate_t *outGroup, surrogate_t *outOrder,
			const BAT_t<T> *columnPtr, const surrogate_t * intermediate_group_idx,
			const surrogate_t start, const surrogate_t nitems, surrogate_t *histogram,
			const surrogate_t startGrpIdx, const surrogate_t curGrpNum);

};

template <class T>
inline void __attribute__((always_inline))
Partition::histogram_based_reordering(surrogate_t *outGroup, surrogate_t *outOrder,
		const BAT_t<T> *columnPtr, const surrogate_t * intermediate_group_idx,
		const surrogate_t start, const surrogate_t nitems, surrogate_t *histogram,
		const surrogate_t startGrpIdx, const surrogate_t curGrpNum)
{
	assert(histogram != NULL);
	assert(curGrpNum > 0);


	//TODO: I need the *startGrpIdx* and *curGrpNum* to reduce the size of software-managed buffer
	const uint32_t fanout = curGrpNum;
	//const uint32_t fanout = 100;

	//std::cout << "fanout: " << fanout << std::endl;

	surrogate_t *input_oid	= columnPtr->oids + start;
	const surrogate_t *input_grp_idx = intermediate_group_idx + start;

	surrogate_t* output_oid = outOrder + start;


	/** start position of each partition/group (of oids) **/
	uint32_t			dst_oid[fanout]		__attribute__((aligned(CACHE_LINE_SIZE)));

	/** software-managed buffer **/
	cacheline_oid_enhanced_t		buffer_oid[fanout]	__attribute__((aligned(CACHE_LINE_SIZE)));

	assert(NULL != buffer_oid);

	uint32_t idx = 0;
	//const uint32_t mask = ((1 << nradixbits) - 1) << bitshift;
	/** count elements per partition **/
	for (uint32_t i = 0; i < nitems; ++i) {
		//idx = RADIX_HASH_FUNCTION(*input_val, mask, bitshift);
		idx = (*input_grp_idx) - startGrpIdx;
		assert((idx >= 0) && (idx < curGrpNum));

		histogram[idx]++;
		input_grp_idx++;
	}

	for (uint32_t i = 0; i < fanout; ++i) {
		//std::cout << "i: " << i << std::endl;
		//printf("i: %u\n", i);
		//assert(NULL != (void *)(buffer_oid[i].oids));
		//assert(NULL != *(buffer_oid + i));
		//std::cout << "size of cacheline_oid_enhanced_t: " << sizeof(buffer_oid[i]);
		buffer_oid[i].data.slot = 0;
	}
	//exit(EXIT_SUCCESS);

	uint32_t offset_oid = 0;
	uint32_t out_grp_idx = startGrpIdx;
	uint32_t accum_idx = start;
	/** determine the start and end of each partition depending on the counts **/
	for (uint32_t i = 0; i < fanout; ++i) {
		dst_oid[i] = offset_oid;
		/* TODO:for aligning partition-outputs to cacheline
		 * (related to last portion of partitioning_phase())  */
		//offset_oid += ALIGN_NUM_VALUES(surrogate_t, hist[i]);

		/** Assemble the outGroup **/
		for (uint32_t j = 0; j < histogram[i]; ++j) {
			outGroup[accum_idx] = out_grp_idx;
			accum_idx++;
		}
		out_grp_idx++;

		offset_oid += histogram[i];
	}

	/** copy elements to their corresponding partitions at appropriate offset **/
	uint32_t slot = 0;
	surrogate_t* bucket_oid = NULL;
	input_grp_idx = intermediate_group_idx + start;		//reset the input_grp_idx

	for (uint32_t i = 0; i < nitems; ++i) {
		//idx = RADIX_HASH_FUNCTION(*input_val, mask, bitshift);
		idx = (*input_grp_idx) - startGrpIdx;
		assert((idx >= 0) && (idx < curGrpNum));

		/* store in the cache-resident buffer first*/
		slot = buffer_oid[idx].data.slot;
		bucket_oid = (surrogate_t*)(buffer_oid + idx);
		bucket_oid[slot] = *input_oid;
		input_oid++;
		input_grp_idx++;
		slot++;

		/* flush the bucket when full */
		if (slot == VALUESSPERCACHELINE(surrogate_t)) {
			//store_nontemp_64B((output + dst[idx]), (buffer + idx));
			//memcpy(output + dst[idx], buffer + idx, CACHE_LINE_SIZE);
			*(cacheline_oid_enhanced_t *)(output_oid + dst_oid[idx]) = *(cacheline_oid_enhanced_t *)(buffer_oid + idx);
			slot = 0;
			dst_oid[idx] += VALUESSPERCACHELINE(surrogate_t);
		}
		buffer_oid[idx].data.slot = slot;
	}

	/* flush the remainder elements in the buffer */
	uint32_t num;
	//T* destination_val;
	surrogate_t* destination_oid;
	for (uint32_t i = 0; i < fanout; ++i) {
		num = buffer_oid[i].data.slot;
		if (num > 0) {
			//destination_val = output_val + dst_val[i];
			destination_oid = output_oid + dst_oid[i];
			for (uint32_t j = 0; j < num; ++j) {
				//destination_val[j] = buffer_val[i].data.values[j];
				destination_oid[j] = buffer_oid[i].oids.oids[j];
			}
		}
	}
}

}

#endif  //HASHING_H
