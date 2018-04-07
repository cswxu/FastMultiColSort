/*******************************************************************************
 * Copyright (c) 2016
 * The Hong Kong Polytechnic University, Database Group
 *
 * Author: Wenjian Xu (cswxu AT comp DOT polyu.edu.hk)
 *
 * See file LICENSE.md for details.
 *******************************************************************************/

#include	"composer.h"

namespace multiAttrSort{

//template instantiation
template
void Composer::subsort_mergesort<uint64_t>(BAT_t<uint64_t> *subcolumn, uint64_t start, uint64_t nitems,
		bool asc_desc, uint32_t bitwidth);

template
void Composer::subsort_mergesort<uint32_t>(BAT_t<uint32_t> *subcolumn, uint64_t start, uint64_t nitems,
		bool asc_desc, uint32_t bitwidth);

template
void Composer::subsort_mergesort<uint16_t>(BAT_t<uint16_t> *subcolumn, uint64_t start, uint64_t nitems,
		bool asc_desc, uint32_t bitwidth);

void Composer::SortAllColumns() {
	switch(compose_params_.packtype) {
		case PackOIDType::pack:
		   SortAllColumns_Pack();
		   break;
		case PackOIDType::nonpack:
			SortAllColumns_Nonpack();
		   break;
	}
}

void Composer::createHashTable(Hash_t *hashPtr) {

	/* # of entries for hash table */
	size_t mask = max<size_t>(Hashing::HASHmask(num_rows_), 1 << 16);
	//size_t mask = (1ULL << 30);

	surrogate_t *mem_pool = (surrogate_t *)
			malloc_aligned((mask + num_rows_) * sizeof(surrogate_t));

	assert(NULL != mem_pool);

	//hashPtr->width = 4;	/* regard entries of hash table and collision as uint32_t*/
	hashPtr->mask = static_cast<surrogate_t>(mask - 1);
	hashPtr->nil = ENTRY_NONE;
	hashPtr->lim = static_cast<surrogate_t>(num_rows_);	//tuple number; equivalently, size of collision list
	hashPtr->Link = mem_pool;
	hashPtr->Hash = mem_pool + hashPtr->lim;
}

void Composer::multiRoundSorting() {

	/** do the multiple round sorting as the baseline**/
	surrogate_t *order = (surrogate_t *)malloc_aligned(num_rows_ * sizeof(surrogate_t));
	for (surrogate_t id = 0; id < num_rows_; ++id) {
		order[id] = id;
	}

	surrogate_t *outGroup = (surrogate_t *)malloc_aligned(num_rows_ * sizeof(surrogate_t));
	surrogate_t *inGroup = NULL;
	uint32_t columnid;

	bool isLastRound = false;	//if the last round, no need to record order and group information
	for (columnid = 0; columnid < num_columns_; ++columnid) {

		isLastRound = (columnid == num_columns_ - 1);

		//printf("[INFO ] Sort %u-th column...\n", columnid);
		switch(column_values_[columnid]->GetWidth()) {
			case 1:
				sortByGroup<uint8_t>(outGroup,
						(uint8_t *)(column_values_[columnid]->GetColumn()), order, inGroup,
						compose_params_.column_asc_desc[columnid],
						compose_params_.bitwidth[columnid], isLastRound);
				break;
			case 2:
				sortByGroup<uint16_t>(outGroup,
						(uint16_t *)(column_values_[columnid]->GetColumn()), order, inGroup,
						compose_params_.column_asc_desc[columnid],
						compose_params_.bitwidth[columnid], isLastRound);
				break;
			case 4:
				sortByGroup<uint32_t>(outGroup,
						(uint32_t *)(column_values_[columnid]->GetColumn()), order, inGroup,
						compose_params_.column_asc_desc[columnid],
						compose_params_.bitwidth[columnid], isLastRound);
				break;
			case 8:
				sortByGroup<uint64_t>(outGroup,
						(uint64_t *)(column_values_[columnid]->GetColumn()), order, inGroup,
						compose_params_.column_asc_desc[columnid],
						compose_params_.bitwidth[columnid], isLastRound);
				break;
		}

		HybridTimer timer_ordergroup;
		timer_ordergroup.Start();

		if (NULL == inGroup) {
			inGroup = (surrogate_t *)malloc_aligned(num_rows_ * sizeof(surrogate_t));
		}

		swap<surrogate_t>(&inGroup, &outGroup);

		timer_ordergroup.Stop();
		setting_.time_mmcpyOrderGroup += timer_ordergroup.GetNumCycles();
		setting_.time_orderGroupInfo += timer_ordergroup.GetNumCycles();
	}

	assert(inGroup != NULL);

	//destroy
	free(inGroup);
	free(outGroup);
	free(order);
}

void Composer::multiRoundHashing_reordered() {

	/** we need a re-ordered list now**/
	//printf("Into multiRoundHashing_reordered()\n");
	surrogate_t *inOrder = (surrogate_t *)malloc_aligned(num_rows_ * sizeof(surrogate_t));
	surrogate_t *outOrder = (surrogate_t *)malloc_aligned(num_rows_ * sizeof(surrogate_t));
	for (surrogate_t id = 0; id < num_rows_; ++id) {
		inOrder[id] = id;
		outOrder[id] = id;
	}

	surrogate_t *outGroup = (surrogate_t *)malloc_aligned(num_rows_ * sizeof(surrogate_t));
	surrogate_t *inGroup = NULL;
	uint32_t columnid;

	//histogram
	surrogate_t *histogram = NULL;

	//new
	surrogate_t *intermediate_group_idx = (surrogate_t *)malloc_aligned(num_rows_ * sizeof(surrogate_t));

	Hash_t * hashPtr =  (Hash_t *) malloc_aligned(sizeof(Hash_t));
	createHashTable(hashPtr);
	assert(NULL != hashPtr);

	bool isLastRound = false;	//if the last round, no need to record order and group information
	for (columnid = 0; columnid < num_columns_; ++columnid) {

		isLastRound = (columnid == (num_columns_ - 1));

		//printf("[INFO ] Sort %u-th column...\n", columnid);
		switch(column_values_[columnid]->GetWidth()) {
			case 1:
				hashByGroup_reordered<uint8_t>(outGroup, outOrder,
						(uint8_t *)(column_values_[columnid]->GetColumn()), inOrder,
						inGroup, isLastRound, hashPtr,
						intermediate_group_idx, histogram);
				break;
			case 2:
				hashByGroup_reordered<uint16_t>(outGroup, outOrder,
						(uint16_t *)(column_values_[columnid]->GetColumn()), inOrder,
						inGroup, isLastRound, hashPtr,
						intermediate_group_idx, histogram);
				break;
			case 4:
				hashByGroup_reordered<uint32_t>(outGroup, outOrder,
						(uint32_t *)(column_values_[columnid]->GetColumn()), inOrder,
						inGroup, isLastRound, hashPtr,
						intermediate_group_idx, histogram);
				break;
			case 8:
				hashByGroup_reordered<uint64_t>(outGroup, outOrder,
						(uint64_t *)(column_values_[columnid]->GetColumn()), inOrder,
						inGroup, isLastRound, hashPtr,
						intermediate_group_idx, histogram);
				break;
		}

		if (NULL == inGroup) {
			inGroup = (surrogate_t *)malloc_aligned(num_rows_ * sizeof(surrogate_t));
		}

		if (!isLastRound) {
			swap<surrogate_t>(&inGroup, &outGroup);
			swap<surrogate_t>(&inOrder, &outOrder);
		}
	}

	assert(inGroup != NULL);

#if 0
	//output the outGroup for debugging
	for (uint64_t i = 0; i < num_rows_; i++) {
		std::cout << outGroup[i] << std::endl;
	}
#endif

	//destroy
	free(inGroup);
	free(outGroup);
	free(inOrder);
	free(outOrder);
	free(intermediate_group_idx);

	free(hashPtr->Link);	//link is the start memory position of both hash table and collision list
	free(hashPtr);

	if (histogram != NULL) {	//only when # of rounds = 1, histogram is null
		free(histogram);
	}
}

void Composer::multiRoundHashing_partitioned_reordered() {

	/** we need a re-ordered list now**/
	//printf("Into multiRoundHashing_partitioned_reordered()\n");
	surrogate_t *inOrder = (surrogate_t *)malloc_aligned(num_rows_ * sizeof(surrogate_t));
	surrogate_t *outOrder = (surrogate_t *)malloc_aligned(num_rows_ * sizeof(surrogate_t));
	for (surrogate_t id = 0; id < num_rows_; ++id) {
		inOrder[id] = id;
		outOrder[id] = id;
	}

	surrogate_t *outGroup = (surrogate_t *)malloc_aligned(num_rows_ * sizeof(surrogate_t));
	surrogate_t *inGroup = NULL;
	uint32_t columnid;

	//histogram
	surrogate_t *histogram = NULL;

	//new
	surrogate_t *intermediate_group_idx = (surrogate_t *)malloc_aligned(num_rows_ * sizeof(surrogate_t));

	Hash_t * hashPtr =  (Hash_t *) malloc_aligned(sizeof(Hash_t));
	createHashTable(hashPtr);
	assert(NULL != hashPtr);

	bool isLastRound = false;	//if the last round, no need to record order and group information
	for (columnid = 0; columnid < num_columns_; ++columnid) {

		isLastRound = (columnid == (num_columns_ - 1));

		//printf("[INFO ] Sort %u-th column...\n", columnid);
		switch(column_values_[columnid]->GetWidth()) {
			case 1:
				hashByGroup_partitioned_reordered<uint8_t>(outGroup, outOrder,
						(uint8_t *)(column_values_[columnid]->GetColumn()), inOrder,
						inGroup, isLastRound, hashPtr,
						intermediate_group_idx, histogram);
				break;
			case 2:
				hashByGroup_partitioned_reordered<uint16_t>(outGroup, outOrder,
						(uint16_t *)(column_values_[columnid]->GetColumn()), inOrder,
						inGroup, isLastRound, hashPtr,
						intermediate_group_idx, histogram);
				break;
			case 4:
				hashByGroup_partitioned_reordered<uint32_t>(outGroup, outOrder,
						(uint32_t *)(column_values_[columnid]->GetColumn()), inOrder,
						inGroup, isLastRound, hashPtr,
						intermediate_group_idx, histogram);
				break;
			case 8:
				hashByGroup_partitioned_reordered<uint64_t>(outGroup, outOrder,
						(uint64_t *)(column_values_[columnid]->GetColumn()), inOrder,
						inGroup, isLastRound, hashPtr,
						intermediate_group_idx, histogram);
				break;
		}

		if (NULL == inGroup) {
			inGroup = (surrogate_t *)malloc_aligned(num_rows_ * sizeof(surrogate_t));
		}

		if (!isLastRound) {
			swap<surrogate_t>(&inGroup, &outGroup);
			swap<surrogate_t>(&inOrder, &outOrder);
		}
	}

	assert(inGroup != NULL);

#if 0
	//output the outGroup for debugging
	for (uint64_t i = 0; i < num_rows_; i++) {
		std::cout << outGroup[i] << std::endl;
	}
#endif

	//destroy
	free(inGroup);
	free(outGroup);
	free(inOrder);
	free(outOrder);
	free(intermediate_group_idx);

	free(hashPtr->Link);	//link is the start memory position of both hash table and collision list
	free(hashPtr);

	if (histogram != NULL) {	//only when # of rounds = 1, histogram is null
		free(histogram);
	}
}

void Composer::multiRoundHashing_non_reordered() {
	std::cout << "start multiRoundHashing_non_reordered.." << std::endl;

	surrogate_t *outGroup = (surrogate_t *)malloc_aligned(num_rows_ * sizeof(surrogate_t));
	surrogate_t *inGroup = NULL;
	uint32_t columnid;
	surrogate_t outGroupNum = 0;
	surrogate_t inGroupNum = 0;

	bool isLastRound = false;	//if the last round, no need to record order and group information
	for (columnid = 0; columnid < num_columns_; ++columnid) {

		isLastRound = (columnid == (num_columns_ - 1));

		//printf("[INFO ] Sort %u-th column...\n", columnid);
		switch(column_values_[columnid]->GetWidth()) {
			case 1:
				hashByGroup_non_reordered<uint8_t>(outGroup, &outGroupNum,
						(uint8_t *)(column_values_[columnid]->GetColumn()), inGroup, &inGroupNum);
				break;
			case 2:
				hashByGroup_non_reordered<uint16_t>(outGroup, &outGroupNum,
						(uint16_t *)(column_values_[columnid]->GetColumn()), inGroup, &inGroupNum);
				break;
			case 4:
				hashByGroup_non_reordered<uint32_t>(outGroup, &outGroupNum,
						(uint32_t *)(column_values_[columnid]->GetColumn()), inGroup, &inGroupNum);
				break;
			case 8:
				hashByGroup_non_reordered<uint64_t>(outGroup, &outGroupNum,
						(uint64_t *)(column_values_[columnid]->GetColumn()), inGroup, &inGroupNum);
				break;
		}

		if (NULL == inGroup) {
			inGroup = (surrogate_t *)malloc_aligned(num_rows_ * sizeof(surrogate_t));
		}

		if (!isLastRound) {
			swap<surrogate_t>(&inGroup, &outGroup);

			/** swap inGroupNum and outGroupNum */
			surrogate_t tmp = outGroupNum;
			outGroupNum = inGroupNum;
			inGroupNum = tmp;
		}
	}

	assert(inGroup != NULL);

#if 0
	//output the outGroup for debugging
	for (uint64_t i = 0; i < num_rows_; i++) {
		std::cout << outGroup[i] << std::endl;
	}
#endif

	//destroy
	free(inGroup);
	free(outGroup);
}

#if 0
void ChainComposer::subsort_pack(BAT_pack_t *values, uint64_t start, uint64_t nitems,
		bool asc_desc, uint32_t bitwidth) {
	if (nitems <= 1)	//trivial case
		return;

	//sort the range [values->elements[start], values->elements[start+nitems])
	assert(start < values->num_elements);
	BAT_pack_t rangeValues;
	rangeValues.num_elements = nitems;
	rangeValues.elements = values->elements + start;


	/*
	 * sort the range, call different functions by its size
	 * TODO: have not differentiate ascending/descending yet, i.e., <asc_desc> is not utilized
	 */
	setting_.bitwidth = bitwidth;	//set the bit width of the target column

	//the whole data can be resided in (half of) L2-cache
	//int64_t *output = (int64_t *)malloc_aligned(nitems * sizeof(int64_t));
	//int64_t *input 	= (int64_t *)rangeValues.elements;
	//printf("enter in sorting small blocks, size:%lu\n", nitems);
	//Mergesort::avxmergesort_rem_aligned(&input, &output, (uint32_t)nitems);
	//std::sort(input, input + nitems);
	//printf("exit sorting small blocks?\n");


	if (nitems < (uint64_t)BLOCKSIZE) {
		int64_t *input 	= (int64_t *)rangeValues.elements;
		std::sort(input, input + nitems);
	} else {
		sorter_->do_sort_pack(&rangeValues, &setting_);
	}



#if 0
		//assign the output to replace the original input
		if (output != (int64_t *)(rangeValues.elements)) {	//the memory of input and output did not swap
			for (uint64_t i = 0; i < nitems; ++i) {
				rangeValues.elements[i] = ((element_t *)output)[i];
			}
			free(output);
		} else {	//the memory of input and output switches
			free(input);
		}
#endif
}
#endif

template <class T>
void Composer::subsort_radixsort(BAT_t<T> *subcolumn, uint64_t start, uint64_t nitems,
		bool asc_desc, uint32_t bitwidth) {

	/* since the interface for radix sort requires the
	 * memory for values and oids to be 64byte aligned, do the memory
	 * re-allocation if necessary
	 */
	T *valPtrs = subcolumn->values + start;
	surrogate_t *oidPtrs = subcolumn->oids + start;
	bool isValAligned = (0 == ((uint64_t)valPtrs & 63ULL));
	bool isOidAligned = (0 == ((uint64_t)oidPtrs & 63ULL));

	if (!isValAligned) {
		//re-allocate aligned memory for values, assigned to <valPtrs>
		valPtrs = (T *) malloc_aligned(nitems * sizeof(T));
		memcpy(valPtrs, subcolumn->values + start, nitems * sizeof(T));
	}

	if (!isOidAligned) {
		//re-allocate aligned memory for oids, assigned to <oidPtrs>
		oidPtrs = (surrogate_t *) malloc_aligned(nitems * sizeof(surrogate_t));
		memcpy(oidPtrs, subcolumn->oids + start, nitems * sizeof(surrogate_t));
	}

	/* call the underlying sorting algorithms after align the memory */
	/* sort the the targetColumn with range [start, start+nitems] */
	setting_.bitwidth = bitwidth;	//set the bit width of the target column
#if 0
	BAT_t<T> targetSubCol;
	targetSubCol.num_elements = nitems;
	targetSubCol.values = valPtrs;
	targetSubCol.oids = oidPtrs;
#endif


#if 0
	{
		//radix_sorter_ should be instantialized here
		assert(radix_sorter != NULL);
		radix_sorter.do_sort<T>(&targetSubCol, &setting_);
	}
#endif

	/* if memory is re-alloccated because of alignment, copy the content back */
	if (!isValAligned) {
		memcpy(subcolumn->values + start, valPtrs, nitems * sizeof(T));
		free(oidPtrs);
	}

	if (!isOidAligned) {
		memcpy(subcolumn->oids + start, oidPtrs, nitems * sizeof(surrogate_t));
		free(oidPtrs);
	}
}

template <class T>
void Composer::subsort_mergesort(BAT_t<T> *subcolumn, uint64_t start, uint64_t nitems,
		bool asc_desc, uint32_t bitwidth) {
	/**
	 * We separate subsort_mergesort from subsort_radix sort because we need to re-allocate
	 * the memory no matter the original memory is aligned or not; this is determined by the interface
	 * of avxmergesort_anytype()...
	 */
	T *origin_val = subcolumn->values + start;
	surrogate_t *origin_oid = subcolumn->oids + start;

	/* allocate aligned memory for [in, out] chunks */
	T * inptr_val = (T *) malloc_aligned(sizeof(T) * nitems);
	surrogate_t * inptr_oid = (surrogate_t *) malloc_aligned(
			sizeof(surrogate_t) * nitems);
	T * outptr_val = (T *) malloc_aligned(sizeof(T) * nitems);
	surrogate_t * outptr_oid = (surrogate_t *) malloc_aligned(
			sizeof(surrogate_t) * nitems);

	/* copy the value content to the input memory, piggyback with the flipping */
	{
		for (uint64_t idx = 0; idx < nitems; ++idx) {
			inptr_val[idx] = origin_val[idx] ^ (1ULL << (sizeof(T) * 8 - 1));//can use SIMD to speed up???
		}
	}
	/* copy the oid to the input memory */
	memcpy(inptr_oid, origin_oid, nitems * sizeof(surrogate_t));

	switch (setting_.intrinType) {
	case IntrinsicsType::AVX:
		Mergesort::avxmergesort_anytype<T>(&inptr_val, &inptr_oid, &outptr_val,
				&outptr_oid, nitems);
		break;
	case IntrinsicsType::SSE:
		break;
	case IntrinsicsType::SCALAR:
		break;
	}

	/* after sorting, copy the value and oid back to the original memory */
	{
		for (uint64_t idx = 0; idx < nitems; ++idx) {
			origin_val[idx] = outptr_val[idx] ^ (1ULL << (sizeof(T) * 8 - 1));//can use SIMD to speed up???
		}
	}

	memcpy(origin_oid, outptr_oid, nitems * sizeof(surrogate_t));

	free(inptr_val);
	free(inptr_oid);
	free(outptr_val);
	free(outptr_oid);
}

template <class T>
void Composer::subsort(BAT_t<T> *subcolumn, uint64_t start, uint64_t nitems,
		bool asc_desc, uint32_t bitwidth) {

	if (nitems <= 1)	//trivial case
		return;

	assert(start < subcolumn->num_elements);

	switch(compose_params_.sortalgo) {
		case SortAlgo::mergesort:
			subsort_mergesort<T>(subcolumn, start, nitems, asc_desc, bitwidth);
			break;
		case SortAlgo::radixsort:
			subsort_radixsort<T>(subcolumn, start, nitems, asc_desc, bitwidth);
			break;
	}
}

#if 0
void ChainComposer::BATproject_pack(BAT_pack_t *inValues, surrogate_t *inOrder) {

	element_t *column_reordered = (element_t *)malloc_aligned(num_rows_ * sizeof(element_t));

	/* rearrange the column values as well as the oid */
	for (surrogate_t oid = 0; oid < num_rows_; ++oid) {
		column_reordered[oid].value = (inValues->elements)[inOrder[oid]].value;
		column_reordered[oid].oid 	= (inValues->elements)[inOrder[oid]].oid;
	}

	//destroy the memory for original inValues
	free(inValues->elements);
	inValues->elements = column_reordered;
}
#endif

template <class T>
void Composer::BATproject(BAT_t<T> *inValues, surrogate_t *inOrder) {

	/** it is freed in SortByGroup()**/
	T *values_reordered = (T *)malloc_aligned(num_rows_ * sizeof(T));
	//std::cout << "start projection" << std::endl;
	for (surrogate_t oid = 0; oid < num_rows_; ++oid) {

		values_reordered[oid] = inValues->values[inOrder[oid]];
	}
	//std::cout << "end projection" << std::endl;

	inValues->values = values_reordered;

}

template <class T>
void Composer::extractGroupInfo(surrogate_t * outGroup, BAT_t<T> *inValues, surrogate_t *inGroup) {
	//surrogate_t *group = *outGroup;
	T prev;
	prev = inValues->values[0];
	surrogate_t curSurrogate = 0;
	uint64_t 	curId = 1;
	outGroup[0] = curSurrogate;
	if (inGroup != NULL) {
		/* need to take both <inValues> and <inGroup> to decide the grouping info */
		while (curId < num_rows_) {
#if 0
			if (prev != inValues->values[curId]) {
				curSurrogate++;
				prev = inValues->values[curId];
			} else if (inGroup[curId] != inGroup[curId-1]) {
				curSurrogate++;
			}

			if ((prev != inValues->values[curId]) || (inGroup[curId] != inGroup[curId-1])) {
				curSurrogate++;
				prev = inValues->values[curId];
			}
#endif

			curSurrogate += ((prev != inValues->values[curId]) || (inGroup[curId] != inGroup[curId-1]));
			prev = inValues->values[curId];

			outGroup[curId] = curSurrogate;
			curId++;
		}
	} else {
		/* just calculate the group info according to <inValues>*/
		while (curId < num_rows_) {
#if 0
			if (prev != inValues->values[curId]) {
				curSurrogate++;
				prev = inValues->values[curId];
			}
#endif

			curSurrogate += (prev != inValues->values[curId]);
			prev = inValues->values[curId];

			outGroup[curId] = curSurrogate;
			curId++;
		}
	}
}

#if 0
void ChainComposer::extractGroupInfo_pack(surrogate_t ** outGroup, BAT_pack_t *inValues, surrogate_t *inGroup) {
	surrogate_t *group = *outGroup;
	value_t prev;
	prev = inValues->elements[0].value;
	surrogate_t curSurrogate = 0;
	uint64_t 	curId = 1;
	group[0] = curSurrogate;
	if (inGroup != NULL) {
		/* need to take both <inValues> and <inGroup> to decide the grouping info */
		while (curId < num_rows_) {
			if (prev != inValues->elements[curId].value) {
				curSurrogate++;
				prev = inValues->elements[curId].value;
			} else if (inGroup[curId] != inGroup[curId-1]) {
				curSurrogate++;
			}
			group[curId] = curSurrogate;
			curId++;
		}
	} else {
		/* just calculate the group info according to <inValues>*/
		while (curId < num_rows_) {
			if (prev != inValues->elements[curId].value) {
				curSurrogate++;
				prev = inValues->elements[curId].value;
			}
			group[curId] = curSurrogate;
			curId++;
		}
	}
}
#endif

#if 0
void ChainComposer::sortByGroup_pack(surrogate_t **outOrder, surrogate_t **outGroup,
		BAT_pack_t *inValues, surrogate_t *inOrder, surrogate_t *inGroup,
		bool asc_desc, uint32_t bitwidth) {

	//printf("enter in sortByGroup_Pack()\n");

	HybridTimer timer_reconstruct;
	timer_reconstruct.Start();

	//BAT_pack_t inValuesReordered;
	if(inOrder != NULL) {
		/** if it is not the 1st round sorting,
		 * rearrange <inValues> according to <inOrder>**/
		BATproject_pack(inValues, inOrder);
	}
	timer_reconstruct.Stop();

	setting_.time_tupleReconstruct += (double)timer_reconstruct.GetNumCycles();

	if (inGroup != NULL) {
#if 0
		//output the inGroup information
		printf("The inGroup info:\n");
		for (uint64_t i = 0; i < num_rows_; ++i) {
			printf("%u\n", inGroup[i]);
		}
#endif

		/** do the sorting within the sub-groups **/
		surrogate_t prev = inGroup[0];
		uint64_t start 	= 0;
		uint64_t end 	= 1;
		for (; end < num_rows_; ++end) {
			if (inGroup[end] != prev) {	//identify the subgroup
				/* do the sorting in the sub range [start, end)*/
				//printf("sorting the sub range:[%lu, %lu)\n", start, end);
				subsort_pack(inValues, start, end - start, asc_desc, bitwidth);
				//printf("exit sub sort\n");

				start = end;
				prev = inGroup[end];
			}
		}
		//sort [start, num_rows_)
		//printf("final sub range:[%lu, %lu)\n", start, end);
		subsort_pack(inValues, start, end - start, asc_desc, bitwidth);
		//printf("exit final sub sort\n");

	} else {
		/** sort the whole column (as the first-round sorting) **/
		setting_.bitwidth = bitwidth;

#if 0
		//debug: output the sorted first column
		printf("the first column before sorting: \n");
		for (uint64_t i = 0; i < num_rows_; ++i) {
			printf("%u\n", inValuesReordered.elements[i].value);
		}
#endif

		HybridTimer timer_firstround;
		timer_firstround.Start();

		sorter_->do_sort_pack(inValues, &setting_);

		timer_firstround.Stop();
		setting_.time_firstRound += timer_firstround.GetNumCycles();

#if 0
		//debug: output the sorted first column
		printf("the first column after sorting: \n");
		for (uint64_t i = 0; i < num_rows_; ++i) {
			printf("%u\n", inValuesReordered.elements[i].value);
		}
#endif
	}

	HybridTimer timer_ordergroup;
	timer_ordergroup.Start();

	/** after sorting, extract the sorted oid (as output) **/
	surrogate_t *order = *outOrder;
	uint64_t i;

	for (i = 0; i < num_rows_; ++i) {
		order[i] = inValues->elements[i].oid;
	}

	/** after sorting, calculate the grouping (tied values) information (as output) **/

	extractGroupInfo_pack(outGroup, inValues, inGroup);


	timer_ordergroup.Stop();
	setting_.time_orderGroupInfo += timer_ordergroup.GetNumCycles();
}
#endif


template <class T>
void Composer::hashByGroup_reordered(surrogate_t *outGroup, surrogate_t *outOrder,
		T *inValues, surrogate_t *inOrder, surrogate_t *inGroup,
		bool isLastRound, Hash_t *hashPtr,
		surrogate_t *intermediate_group_idx, surrogate_t *histogram) {

	assert(NULL != inValues);
	assert(NULL != hashPtr);

	/** Step 1: reconstruct the column by re-ordered oids **/
	HybridTimer timer_reconstruct;
	timer_reconstruct.Start();

	BAT_t<T> targetColumn;
	targetColumn.num_elements = this->num_rows_;
	targetColumn.values = inValues;
	targetColumn.oids = inOrder;

	if(inGroup != NULL) {
		/** if it is not the 1st round hashing,
		 * rearrange <inValues> according to <inOrder>**/
		BATproject<T>(&targetColumn, inOrder);
	}
	timer_reconstruct.Stop();

	setting_.time_tupleReconstruct += (double)timer_reconstruct.GetNumCycles();

	/** Step 2: hashing and histogram-based re-ordering **/
	HybridTimer timer_hashing;
	timer_hashing.Start();

	HybridTimer timer_reordering;

	const size_t block_size = L2_CACHE_SIZE / (2 * sizeof(T));
	surrogate_t* hashed_per_block = (surrogate_t *) malloc_aligned(block_size * sizeof(surrogate_t));

	surrogate_t startGrpIdx = 0;
	surrogate_t grp_num = 0;

	if (inGroup != NULL) {	//other round of hashing
#if 0
		//output the inGroup information
		printf("The inGroup info:\n");
		for (uint64_t i = 0; i < num_rows_; ++i) {
			printf("%u\n", inGroup[i]);
		}
#endif

		/** do the hashing within the sub-groups **/
		surrogate_t prev = inGroup[0];
		uint64_t start 	= 0;
		uint64_t end 	= 1;
		//uint32_t groupsz = 0;
		//uint8_t cmp;
		for (; end < num_rows_; ++end) {
			if (inGroup[end] != prev) {	//identify the subgroup
				//printf("hashing the sub range:[%lu, %lu)\n", start, end);

				grp_num = 0;
				Hashing::subgroup_range<T>(intermediate_group_idx, &grp_num, &targetColumn,
						start, end - start, startGrpIdx, hashPtr, hashed_per_block);

				//std::cout << "group_num: " << grp_num << std::endl;
				assert(grp_num > 0);

				/** reorder the column based on intermediate_group_idx **/
				if (!isLastRound) {	//not last round, do the reordering

					memset((void *)histogram, 0, grp_num * sizeof(surrogate_t));

					Partition::histogram_based_reordering<T>(outGroup, outOrder, &targetColumn,
							intermediate_group_idx, start, end - start, histogram, startGrpIdx, grp_num);
				}

				startGrpIdx += grp_num;

				start = end;
				prev =  inGroup[end];
			}
		}

		//to hash the final sub range (if exists)
		if ((end - start) > 0) {

			grp_num = 0;
			Hashing::subgroup_range<T>(intermediate_group_idx, &grp_num, &targetColumn,
					start, end - start, startGrpIdx, hashPtr, hashed_per_block);


			/** reorder the column based on intermediate_group_idx **/
			if (!isLastRound) {	//not last round, do the reordering

				memset((void *)histogram, 0, grp_num * sizeof(surrogate_t));

				Partition::histogram_based_reordering<T>(outGroup, outOrder, &targetColumn,
						intermediate_group_idx, start, end - start, histogram, startGrpIdx, grp_num);
			}

			startGrpIdx += grp_num;
		}

#if 1
		std::cout << "accum # of groups: " << startGrpIdx << std::endl;
#endif

	} else {	//first round hashing

		HybridTimer timer_firstpass;
		timer_firstpass.Start();

		grp_num = 0;

		/** hash the whole column (as the first-round hashing) **/
		Hashing::subgroup_range<T>(intermediate_group_idx, &grp_num,
				&targetColumn, 0, num_rows_, startGrpIdx, hashPtr, hashed_per_block);

		timer_reordering.Start();

		/** reorder the column based on intermediate_group_idx **/
		if (!isLastRound) {	//not last round, do the reordering

			//initialize the histogram
			assert(histogram == NULL);
			histogram = (surrogate_t *)malloc_aligned(grp_num * sizeof(surrogate_t));
			memset((void *)histogram, 0, grp_num * sizeof(surrogate_t));

#if 0
			std::cout << "# of groups after 1-st round hashing: " << grp_num << std::endl;
#endif

			Partition::histogram_based_reordering<T>(outGroup, outOrder, &targetColumn,
					intermediate_group_idx, 0, num_rows_, histogram, 0, grp_num);
		}

		timer_reordering.Stop();
		setting_.time_hashReordering += timer_reordering.GetNumCycles();

		timer_firstpass.Stop();
		setting_.time_firstpass_sort = timer_firstpass.GetNumCycles();
#if 0
		//debug: output the sorted first column
		printf("column vales after sorting: \n");
		for (uint64_t i = 0; i < num_rows_; ++i) {
			std::cout << targetColumn.values[i] << std::endl;
		}
#endif

#if 1
		std::cout << "accum # of groups: " << grp_num << std::endl;
#endif
	}

	timer_hashing.Stop();
	setting_.time_multipass_sort += timer_hashing.GetNumCycles();

	//important: we need to free the memory in targetColumn.values *if it is not the first column*,
	//because the memory for the *re-ordered* column values is allocated in BATproject() but not released
	if (inGroup != NULL) {	//not the first column
		free(targetColumn.values);
	}

	free(hashed_per_block);
}


template <class T>
void Composer::hashByGroup_partitioned_reordered(surrogate_t *outGroup, surrogate_t *outOrder,
		T *inValues, surrogate_t *inOrder, surrogate_t *inGroup,
		bool isLastRound, Hash_t *hashPtr,
		surrogate_t *intermediate_group_idx, surrogate_t *histogram) {

	assert(NULL != inValues);
	assert(NULL != hashPtr);

	/** Step 1: reconstruct the column by re-ordered oids **/
	HybridTimer timer_reconstruct;
	timer_reconstruct.Start();

	BAT_t<T> targetColumn;
	targetColumn.num_elements = this->num_rows_;
	targetColumn.values = inValues;
	targetColumn.oids = inOrder;

	if(inGroup != NULL) {
		/** if it is not the 1st round hashing,
		 * rearrange <inValues> according to <inOrder>**/
		BATproject<T>(&targetColumn, inOrder);
	}
	timer_reconstruct.Stop();

	setting_.time_tupleReconstruct += (double)timer_reconstruct.GetNumCycles();

	/** Step 2: hashing and histogram-based re-ordering **/
	HybridTimer timer_hashing;
	timer_hashing.Start();

	HybridTimer timer_reordering;

	/* decide how many tuples fit in L2 cache as a partition (including those auxiliary structures) */
	/* x <T>, x surrogate_t as oids, 2x surrogate_t entries for hash/collision list */
	//uint32_t cardinality_per_partition = L2_CACHE_SIZE / (sizeof(T) + 3 * sizeof(surrogate_t));

	const size_t block_size = L2_CACHE_SIZE / (2 * sizeof(T));
	surrogate_t* hashed_per_block = (surrogate_t *) malloc_aligned(block_size * sizeof(surrogate_t));

	surrogate_t startGrpIdx = 0;
	surrogate_t grp_num = 0;

	if (inGroup != NULL) {	//other round of hashing
#if 0
		//output the inGroup information
		printf("The inGroup info:\n");
		for (uint64_t i = 0; i < num_rows_; ++i) {
			printf("%u\n", inGroup[i]);
		}
#endif

		/** do the hashing within the sub-groups **/
		surrogate_t prev = inGroup[0];
		uint64_t start 	= 0;
		uint64_t end 	= 1;
		//uint32_t groupsz = 0;
		//uint8_t cmp;
		for (; end < num_rows_; ++end) {
			if (inGroup[end] != prev) {	//identify the subgroup
				//printf("hashing the sub range:[%lu, %lu)\n", start, end);

				grp_num = 0;
				Hashing::subgroup_range<T>(intermediate_group_idx, &grp_num, &targetColumn,
						start, end - start, startGrpIdx, hashPtr, hashed_per_block);

				//std::cout << "group_num: " << grp_num << std::endl;
				assert(grp_num > 0);

				/** reorder the column based on intermediate_group_idx **/
				if (!isLastRound) {	//not last round, do the reordering

					memset((void *)histogram, 0, grp_num * sizeof(surrogate_t));

					Partition::histogram_based_reordering<T>(outGroup, outOrder, &targetColumn,
							intermediate_group_idx, start, end - start, histogram, startGrpIdx, grp_num);
				}

				startGrpIdx += grp_num;

				start = end;
				prev =  inGroup[end];
			}
		}

		//to hash the final sub range (if exists)
		if ((end - start) > 0) {

			grp_num = 0;
			Hashing::subgroup_range<T>(intermediate_group_idx, &grp_num, &targetColumn,
					start, end - start, startGrpIdx, hashPtr, hashed_per_block);


			/** reorder the column based on intermediate_group_idx **/
			if (!isLastRound) {	//not last round, do the reordering

				memset((void *)histogram, 0, grp_num * sizeof(surrogate_t));

				Partition::histogram_based_reordering<T>(outGroup, outOrder, &targetColumn,
						intermediate_group_idx, start, end - start, histogram, startGrpIdx, grp_num);
			}

			startGrpIdx += grp_num;
		}

#if 1
		std::cout << "accum # of groups: " << startGrpIdx << std::endl;
#endif

	} else {	//first round hashing

		HybridTimer timer_firstpass;
		timer_firstpass.Start();

		/** Calculate how many partitions required as first round hashing **/
		//uint32_t partition_num = num_rows_ / cardinality_per_partition;

		//uint32_t num_radix_bits = ceil(log2(partition_num));


		grp_num = 0;

		/** hash the whole column (as the first-round hashing) **/
		Hashing::subgroup_range<T>(intermediate_group_idx, &grp_num,
				&targetColumn, 0, num_rows_, startGrpIdx, hashPtr, hashed_per_block);

		timer_reordering.Start();

		/** reorder the column based on intermediate_group_idx **/
		if (!isLastRound) {	//not last round, do the reordering

			//initialize the histogram
			assert(histogram == NULL);
			histogram = (surrogate_t *)malloc_aligned(grp_num * sizeof(surrogate_t));
			memset((void *)histogram, 0, grp_num * sizeof(surrogate_t));

#if 0
			std::cout << "# of groups after 1-st round hashing: " << grp_num << std::endl;
#endif

			Partition::histogram_based_reordering<T>(outGroup, outOrder, &targetColumn,
					intermediate_group_idx, 0, num_rows_, histogram, 0, grp_num);
		}

		timer_reordering.Stop();
		setting_.time_hashReordering += timer_reordering.GetNumCycles();

		timer_firstpass.Stop();
		setting_.time_firstpass_sort = timer_firstpass.GetNumCycles();
#if 0
		//debug: output the sorted first column
		printf("column vales after sorting: \n");
		for (uint64_t i = 0; i < num_rows_; ++i) {
			std::cout << targetColumn.values[i] << std::endl;
		}
#endif

#if 1
		std::cout << "accum # of groups: " << grp_num << std::endl;
#endif
	}

	timer_hashing.Stop();
	setting_.time_multipass_sort += timer_hashing.GetNumCycles();

	//important: we need to free the memory in targetColumn.values *if it is not the first column*,
	//because the memory for the *re-ordered* column values is allocated in BATproject() but not released
	if (inGroup != NULL) {	//not the first column
		free(targetColumn.values);
	}

	free(hashed_per_block);
}

template <class T>
void Composer::hashByGroup_non_reordered(surrogate_t *outGroup, surrogate_t *outGroupNum,
		T *inValues, surrogate_t *inGroup, surrogate_t *inGroupNum) {

	assert(NULL != inValues);

	/**Initialize hash table and collision list**/

	/* # of entries for hash table */
	size_t mask = max<size_t>(Hashing::HASHmask(num_rows_), 1 << 16);
	//size_t mask = (1ULL << 30);

	surrogate_t *mem_pool = (surrogate_t *)
			malloc_aligned((mask + num_rows_) * sizeof(surrogate_t));

	Hash_t * hashPtr = (Hash_t *) malloc_aligned(sizeof(Hash_t));

	//hashPtr->width = 8;	/* regard entries of hash table and collision as uint32_t*/
	hashPtr->mask = static_cast<surrogate_t>(mask - 1);
	hashPtr->nil = ENTRY_NONE;
	hashPtr->lim = static_cast<surrogate_t>(num_rows_);	//tuple number; equivalently, size of collision list
	hashPtr->Link = mem_pool;
	hashPtr->Hash = mem_pool + hashPtr->lim;

	//clear the hash table
	memset((void*)hashPtr->Hash, 0xFF, (hashPtr->mask + 1) * sizeof(surrogate_t));

	/* start the grouping operation by hashing*/
	HybridTimer timer_hashing;
	timer_hashing.Start();

	if (inGroup != NULL) {
#if 0
		//output the inGroup information
		printf("The inGroup info:\n");
		for (uint64_t i = 0; i < num_rows_; ++i) {
			printf("%u\n", inGroup[i]);
		}
#endif
		Hashing::subgroup_non_first_round<T>(outGroup, outGroupNum, inValues, inGroup, inGroupNum, hashPtr);

	} else {	//first round hashing

		HybridTimer timer_firstpass;
		timer_firstpass.Start();

		/** sort the whole column (as the first-round sorting) **/
		Hashing::subgroup_first_round<T>(outGroup, outGroupNum, inValues, hashPtr);

		timer_firstpass.Stop();
		setting_.time_firstpass_sort = timer_firstpass.GetNumCycles();
#if 0
		//debug: output the sorted first column
		printf("column vales after sorting: \n");
		for (uint64_t i = 0; i < num_rows_; ++i) {
			std::cout << targetColumn.values[i] << std::endl;
		}
#endif
	}

	timer_hashing.Stop();
	setting_.time_multipass_sort += timer_hashing.GetNumCycles();

	/**free the memory**/
	free(mem_pool);
	free(hashPtr);
}


template <class T>
void Composer::sortByGroup(surrogate_t *outGroup,
		T *inValues, surrogate_t *order, surrogate_t *inGroup,
		bool asc_desc, uint32_t bitwidth, bool isLastRound) {

	//printf("enter in sortByGroup_Pack()\n");
	assert(NULL != inValues);

	HybridTimer timer_reconstruct;
	timer_reconstruct.Start();

	BAT_t<T> targetColumn;
	targetColumn.num_elements = this->num_rows_;
	targetColumn.values = inValues;

	targetColumn.oids = order;

	//BAT_pack_t inValuesReordered;
	if(inGroup != NULL) {
		/** if it is not the 1st round sorting,
		 * rearrange <inValues> according to <inOrder>**/
		BATproject<T>(&targetColumn, order);
	}
	timer_reconstruct.Stop();

	setting_.time_tupleReconstruct += (double)timer_reconstruct.GetNumCycles();

	HybridTimer timer_sorting;
	timer_sorting.Start();

	if (inGroup != NULL) {

		/** do the sorting within the sub-groups **/
		surrogate_t prev = inGroup[0];
		uint64_t start 	= 0;
		uint64_t end 	= 1;
		uint32_t groupsz = 0;
		uint8_t cmp;
		for (; end < num_rows_; ++end) {
			if (inGroup[end] != prev) {	//identify the subgroup
				/* do the sorting in the sub range [start, end)*/
				//printf("sorting the sub range:[%lu, %lu)\n", start, end);

				subsort<T>(&targetColumn, start, end - start, asc_desc, bitwidth);

				assert(std::is_sorted(targetColumn.values+start, targetColumn.values+end));

				groupsz = end - start;
				group_info_.ngroups++;
				group_info_.sum_group_size += groupsz;
				cmp = groupsz < BLOCKSIZE<T>();
				group_info_.ngroups_incache += cmp;
				group_info_.ngroups_outofcache += !cmp;

				group_info_.ngroups_gtOne += (groupsz > 1);

				start = end;
				prev = inGroup[end];
			}
		}
		//sort [start, num_rows_)
		//printf("final sub range:[%lu, %lu)\n", start, end);
		subsort<T>(&targetColumn, start, end - start, asc_desc, bitwidth);

		groupsz = end - start;
		group_info_.ngroups++;
		group_info_.sum_group_size += groupsz;
		cmp = groupsz < BLOCKSIZE<T>();
		group_info_.ngroups_incache += cmp;
		group_info_.ngroups_outofcache += !cmp;

		group_info_.ngroups_gtOne += (groupsz > 1);

		assert(std::is_sorted(targetColumn.values+start, targetColumn.values+end));
		//printf("exit final sub sort\n");

	} else {
		HybridTimer timer_firstpass;
		timer_firstpass.Start();

		for (uint32_t i = 0; i < num_rows_; ++i) {
			assert(targetColumn.oids[i] < num_rows_);
		}

		/** sort the whole column (as the first-round sorting) **/
		subsort<T>(&targetColumn, 0, num_rows_, asc_desc, bitwidth);
		assert(std::is_sorted(targetColumn.values, targetColumn.values+num_rows_));

		for (uint32_t i = 0; i < num_rows_; ++i) {
			assert(targetColumn.oids[i] < num_rows_);
		}

		timer_firstpass.Stop();
		setting_.time_firstpass_sort = timer_firstpass.GetNumCycles();
	}

	timer_sorting.Stop();
	setting_.time_multipass_sort += timer_sorting.GetNumCycles();

	HybridTimer timer_ordergroup;
	timer_ordergroup.Start();

	/** after sorting, calculate the grouping (tied values) information (as output) **/
	if (!isLastRound) {	//optimize: do not extract grouping if it is the last round
		extractGroupInfo<T>(outGroup, &targetColumn, inGroup);
	}

	timer_ordergroup.Stop();
	setting_.time_extractGrop += timer_ordergroup.GetNumCycles();
	setting_.time_orderGroupInfo += timer_ordergroup.GetNumCycles();

	//important: the re-assignment of targetColumn.oids is in BATproject()
	//also, we need to free the memory in targetColumn.values *if it is not the first column*,
	//because the memory for the *re-ordered* column values is allocated in BATproject() but not released
	if (inGroup != NULL) {	//not the first column
		free(targetColumn.values);
	}
}

#if 0
void ChainComposer::SortAllColumns_Pack() {
	//printf("enter in SortAllColumns_Pack()\n");

	/** create columns with BAT_pack_t representation (assume 4B+4B first) from column_values_**/
	HybridTimer timer_packoid;
	timer_packoid.Start();

	BAT_pack_t* columns;
	columns = (BAT_pack_t*) malloc_aligned(num_columns_*sizeof(BAT_pack_t));

	uint32_t columnid;
	for (columnid = 0; columnid < num_columns_; ++columnid) {
		columns[columnid].num_elements = num_rows_;
		columns[columnid].elements = (element_t *)malloc_aligned(num_rows_ * sizeof(element_t));

		/* Packing the oid with column value for each column*/
		for (surrogate_t oid = 0; oid < num_rows_; ++oid) {
			columns[columnid].elements[oid].value 	= column_values_[columnid][oid];
			columns[columnid].elements[oid].oid 	= oid;
		}
	}
	timer_packoid.Stop();
	setting_.time_packingOID += (double)timer_packoid.GetNumCycles();

	/** do the multiple round sorting as the baseline**/
	surrogate_t *outOrder = (surrogate_t *)malloc_aligned(num_rows_ * sizeof(surrogate_t));
	surrogate_t *outGroup = (surrogate_t *)malloc_aligned(num_rows_ * sizeof(surrogate_t));
	surrogate_t *inOrder = NULL;
	surrogate_t *inGroup = NULL;
	for (columnid = 0; columnid < num_columns_; ++columnid) {
		//setting_.bitwidth = compose_params_.bitwidth[columnid];
		sortByGroup_pack(&outOrder, &outGroup,
				&columns[columnid], inOrder, inGroup,
				compose_params_.column_asc_desc[columnid],
				compose_params_.bitwidth[columnid]);

		if (NULL == inOrder) {
			inOrder = (surrogate_t *)malloc_aligned(num_rows_ * sizeof(surrogate_t));
			inGroup = (surrogate_t *)malloc_aligned(num_rows_ * sizeof(surrogate_t));
		}
		memcpy(inOrder, outOrder, num_rows_ * sizeof(surrogate_t));
		memcpy(inGroup, outGroup, num_rows_ * sizeof(surrogate_t));
		//inOrder = outOrder;
		//inGroup = outGroup;
	}

	assert(inOrder != NULL);
	//the final output is in <inOrder/outOrder>
#if 0
	for (uint64_t i = 0; i < 20; ++i) {
		printf("%u\t", inOrder[i]);
	}
#endif

	//destroy
	if (inOrder != NULL) {
		free(inOrder);
		free(inGroup);
	}

	free(outGroup);
	free(outOrder);

	for (columnid = 0; columnid < num_columns_; ++columnid) {
		free(columns[columnid].elements);
	}
	free(columns);

}
#endif



}

