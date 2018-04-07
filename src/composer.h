/*******************************************************************************
 * Copyright (c) 2016
 * The Hong Kong Polytechnic University, Database Group
 *
 * Author: Wenjian Xu (cswxu AT comp DOT polyu.edu.hk)
 *
 * See file LICENSE.md for details.
 *******************************************************************************/

#ifndef COMPOSER_H
#define COMPOSER_H

//#include	"common.h"
#include	"types.h"
//#include	"sorter.h"
#include	"mergesort.h"
#include	"column.h"
#include	"hashing.h"
#include	"partition.h"


namespace multiAttrSort{

class Composer{
public:
    virtual ~Composer() {
    	delete sorter_;
    }

    void SortAllColumns();

    Composer(setting_t setting, uint32_t ncolumns, uint64_t nrows,
		compose_params_t params, Column **column_values):num_columns_(ncolumns),
		num_rows_(nrows), column_values_(column_values), setting_(setting), compose_params_(params) {

    	group_info_.ngroups = 0;
    	group_info_.ngroups_incache = 0;
    	group_info_.ngroups_outofcache = 0;
    	group_info_.sum_group_size = 0;

    	/**** initialize the sorter ****/
    	switch(compose_params_.sortalgo) {
    		case SortAlgo::mergesort:
    			sorter_ = new Mergesort();
    			break;
    		case SortAlgo::radixsort:
    			//sorter_ = new RadixSort();
    			break;
    	}
    }

    Composer(): num_columns_(0), num_rows_(0){
    	column_values_ = NULL;
    	sorter_ = NULL;
    }

    setting_t getSetting() {return setting_;}
    group_info_t getGroupInfo() {return group_info_;}
//protected:
	Mergesort* sorter_;			/* the sorting implementation */
    uint32_t num_columns_;		/* number of columns to be sorted*/
    uint64_t num_rows_;	/* size of the relation */
	Column ** column_values_;
	setting_t setting_;
	compose_params_t compose_params_;
	group_info_t group_info_;

	virtual void SortAllColumns_Pack() = 0;
	virtual void SortAllColumns_Nonpack() = 0;
	virtual void ExhaustiveSearch() = 0;
	virtual void TwoRoundsExhaustive() = 0;
	virtual void RunAnInstance(std::string intance) = 0;
	virtual void CostModelComparison() = 0;

	virtual void RunAnInstance_hash(std::string intance) = 0;
	virtual void TwoRoundsExhaustive_hash() = 0;
	virtual void HashAllColumns() = 0;

	void multiRoundSorting();
	void multiRoundHashing_non_reordered();
	void multiRoundHashing_reordered();
	void multiRoundHashing_partitioned_reordered();

	/**
	 * Corresponds to MonetDB's BATsubsort function
	 *
	 * @param outOrder 	the object-oid after single-round sorting
	 * @param outGroup 	the tied objects are grouped in the same group, e.g, [0,0,1,1,1,2,...]
	 * @param inValues	the column values to be sorted (as input) <BAT type can be BAT_t or BAT_pack_t>
	 * @param inOrder	the object-oid after previous round (as input)
	 * @param inGroup	the group info after previous round (as input)
	 * @param asc_desc	indicate the ascending or descending order
	 */
	template <class T>
	void sortByGroup(surrogate_t *outGroup,
			T *inValues, surrogate_t *order, surrogate_t *inGroup,
			bool asc_desc, uint32_t bitwidth, bool isLastRound);

	/**
	 * Corresponds to MonetDB's BATgroup_internal function
	 *
	 * @param outGroup 	the group mapping vector after current round of hashing (as output)
	 * @param inValues	the column values to be hashed (as input)
	 * @param inGroup	the group mapping vector after previous round of hashing (as input)
	 */
	template <class T>
	void hashByGroup_non_reordered(surrogate_t *outGroup, surrogate_t *outGroupNum,
			T *inValues, surrogate_t *inGroup, surrogate_t *inGroupNum);

	/**
	 * reordered version
	 *
	 * INPUT:
	 * inValues: input column values
	 * inOrder:  used to reconstruct/lookup the reordered column values
	 * inGroup:  we only need to hash within each group ==> to reduce cache miss hopefully
	 * isLastRound: if it is the last round, only need to generate intermediate_group_idx,
	 * 				no need to reorder the tuples (i.e., generate *outGroup* and *outOrder*)
	 *
	 *
	 * OUTPUT:
	 * outGroup: clustered grouping information
	 * outOrder: since this is not in-place reordering, we need this buffer as output
	 *
	 * AUXILIARY:
	 * hashPtr: hash tables
	 * intermediate_group_idx: store the group index (need to be re-ordered)
	 *
	 */
	template <class T>
	void hashByGroup_reordered(surrogate_t *outGroup, surrogate_t *outOrder,
			T *inValues, surrogate_t *inOrder, surrogate_t *inGroup,
			bool isLastRound, Hash_t *hashPtr,
			surrogate_t *intermediate_group_idx, surrogate_t *histogram);

	template <class T>
	void hashByGroup_partitioned_reordered(surrogate_t *outGroup, surrogate_t *outOrder,
			T *inValues, surrogate_t *inOrder, surrogate_t *inGroup,
			bool isLastRound, Hash_t *hashPtr,
			surrogate_t *intermediate_group_idx, surrogate_t *histogram);

	/**
	 * initialize the hash table and the collision list, shared by different rounds of hashing
	 */
	void createHashTable(Hash_t *hashPtr);

	/**
	 * rearrange <inValues> according to <inOrder>
	 *
	 */
	//void BATproject_pack(BAT_pack_t *inValues, surrogate_t *inOrder);
	template <class T>
	void BATproject(BAT_t<T> *inValues, surrogate_t *inOrder);

	/**
	 * extract group information in <sortByGroup> function
	 *
	 * @param outGroup 	as return value
	 * @param inValues 	the sorted column values, the extraction is based on these values
	 */
	//void extractGroupInfo_pack(surrogate_t ** outGroup, BAT_pack_t *inValues, surrogate_t *inGroup);
	template <class T>
	void extractGroupInfo(surrogate_t * outGroup, BAT_t<T> *inValues, surrogate_t *inGroup);

	/**
	 * sort the range [values->elements[start], values->elements[start+nitems])
	 * invoke different functions according to <nitems>
	 *
	 */
	//void subsort_pack(BAT_pack_t *values, uint64_t start, uint64_t nitems,
	//		bool asc_desc, uint32_t bitwidth);
	template <class T>
	void subsort(BAT_t<T> *values, uint64_t start, uint64_t nitems,
			bool asc_desc, uint32_t bitwidth);

	template <class T>
	void subsort_mergesort(BAT_t<T> *values, uint64_t start, uint64_t nitems,
			bool asc_desc, uint32_t bitwidth);

	template <class T>
	void subsort_radixsort(BAT_t<T> *values, uint64_t start, uint64_t nitems,
			bool asc_desc, uint32_t bitwidth);

};




}

#endif  //COMPOSER_H
