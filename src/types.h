/**
 * @file   types.h
 * @author XU Wenjian (zeroxwj@gmail.com)
 * @date   August 11 2015
 *
 * @brief  Provides general type definitions.
 *
 *
 */
#ifndef TYPES_H
#define TYPES_H

#include <cstdint>
#include <cinttypes>
#include "params.h"
#include <cstdlib>
#include "barrier.h"
#include <cassert>
#include "hybrid_timer.h"
#include <string>
#include <functional>

/**
 * Common type definitions used by all sort/composition implementations.
 */

//forward delcaration
//struct element_t;
//struct BAT_pack_t;
enum class IntrinsicsType;
struct dataset_t;
struct compose_params_t;
struct setting_t;
enum class SortAlgo;
enum class PackOIDType;
enum class ComposerType;


/**
 * A unit of shift operation for constructing a new column
 */
struct construct_entry_t {
	uint32_t finalIdx;	//this index is already transferred by bestOrdering
	uint32_t start;
	uint32_t end;
};

/** an column value with physical surrogate id 
 *  it can be 8B or 16B depending on COLUMN_WIDTH_8B
 *
 *  IMPORTANT: the order of these two attribute matters.
 *  Our machine is little-endian, so in order to make <value> as Most Significant Digits (for sorting as uint64_t),
 *  it should be placed LATTER
 */
#if 1
template <class T>
class element_t {
public:
	surrogate_t oid;
	T value;
};
#endif

template <class T>
class mycomparison {
public:
	bool operator()(element_t<T> i, element_t<T> j) {return i.value < j.value;}
};

/** a cache line to store values **/
template <class T>
union cacheline_val_t {
	struct {
		T values[VALUESSPERCACHELINE(surrogate_t)];	//here its size is the # of oids in a cacheline for consistency issue
	} values;
	struct {
		T values[VALUESSPERCACHELINE(surrogate_t)-1];
		uint32_t slot;
	} data;
};

/** a cache line to store oids **/
struct cacheline_oid_t {
	surrogate_t oids[VALUESSPERCACHELINE(surrogate_t)];
};

/** used in histogram_based_reordering **/
typedef union cacheline {
	struct {
		surrogate_t oids[VALUESSPERCACHELINE(surrogate_t)];	//here its size is the # of oids in a cacheline for consistency issue
	} oids;
	struct {
		surrogate_t oids[VALUESSPERCACHELINE(surrogate_t)-1];
		uint32_t slot;
	} data;
} cacheline_oid_enhanced_t;

/**
 * Assemble MonetDB's BAT representation
 */
template <class T>
class BAT_t {
public:
	surrogate_t *oids;
	T		*values;
	uint64_t	num_elements;
};




/**
  * a column representation with # infomation
  * an array of pair<value, oid>
  */
#if 0
struct BAT_pack_t {
	element_t* elements;
	uint64_t num_elements;
};
#endif


/** the intrinsics type for some kernals/functions **/
enum class IntrinsicsType {
	AVX,
	SSE,
	SCALAR
};

/** two versions of baseline **/
enum class Baseline {
	aligned32,		//always use 32 bit to store each column
	adaptive		//fit to the smallest simd bank
	///+++MonetDB's timsort???
};

/** hold the parameters for the dataset **/
struct dataset_t {
	uint64_t nrows;
	uint32_t ncolumns;
	std::string column_bitwidth;
	//uint32_t column_asc_desc;
	uint32_t zipf;
	uint32_t cardinality;
	uint32_t group_num = 2;	//number of groups (with tied values) for each column
	//uint32_t groupsize = L3_CACHE_SIZE;
};

/** the underlying sorting implementation  **/
enum class SortAlgo {
	mergesort,
	radixsort
};

/** indicate how OID are processed during the underlying sorting**/
enum class PackOIDType {
	pack,
	nonpack
};

/** the approaches of composer to compose sorting of multiple columns **/
enum class ComposerType {
	chaining,
	stitching
};

/** hold the parameters for the composition  **/ 
struct compose_params_t {
	SortAlgo sortalgo = SortAlgo::mergesort;
	PackOIDType packtype = PackOIDType::nonpack;
	bool is_oid_encoded = false;
	bool *column_asc_desc = NULL;	/* true means ascending; false means descending*/
	uint32_t *bitwidth = NULL;		/* bitwidth of each column*/
	Baseline baseline = Baseline::adaptive;

	//new
	uint32_t colwidth_sum = 0;		/* the sum of all column width*/
	uint32_t ordered = 1;			/* 0 indicates unordered */
	uint32_t hash_type = 1;			/* 0 indicates not use re-ordering*/
	//uint32_t *cardinality = NULL;	/* number of distinct values for each column => this information is attached to each Column object */

	//for stitching idea
	//uint32_t stitch_nbits = 11;
	//uint8_t stitch_direction = 0;	//0 is left
};

struct group_info_t {
	uint32_t ngroups = 0;	//number of groups
	uint32_t sum_group_size = 0;	//sum of all group size
	uint32_t ngroups_incache = 0;	//number of groups that are cache-resident (smaller than half of L2-cache)
	uint32_t ngroups_outofcache = 0;

	uint32_t ngroups_gtOne = 0;	//number of groups larger than one
};

/** hold the settings for the composer/sorter and the timing **/
struct setting_t {
	/* general */
	uint32_t 			nthreads = 1;
	IntrinsicsType 		intrinType = IntrinsicsType::AVX;
	//int32_t 			my_tid = 0;		//thread logical id
	//pthread_barrier_t*	barrier = NULL;		//for threads synchronization
	uint32_t			bitwidth = 0;

	/* Timing for chain_composer*/
	//double				time_packingOID	= 0.0; /*not used anymore because packOID style is not used*/
	//double				time_assembleBAT = 0.0;
	double	 			time_tupleReconstruct = 0.0;
	double				time_orderGroupInfo = 0.0;
	double				time_multipass_sort	= 0.0;

	double				time_firstpass_sort = 0.0;

	double 				time_mmcpyOrderGroup = 0.0;
	double 				time_extractGrop = 0.0;	//for detailed profiling

	double				time_stitching = 0.0;
	double				time_optimization = 0.0;

	double				time_hashReordering = 0.0;

	/* for merge sort */
	uint32_t		 	partition_fanout = 1;	//same as merge fan-in
#if 0
	element_t* 			columnchunk = NULL;
	element_t*			tmp_partition = NULL;
	element_t*			tmp_sort = NULL;
	uint64_t			chunksize = 0;
	BAT_pack_t**				threadcolchunks = NULL;
	element_t*			sharedmergebuffer = NULL;	//normally size of L3-cache, shared by all threads, used in the multiway merge phase
	BAT_pack_t*				result = NULL;
#endif
	/* for radix sort */

};
//} __attribute__((aligned(CACHE_LINE_SIZE)));

/**parameters for each thread of merge sort**/
template <class T>
class mergesort_params_t {
public:
	uint32_t 			nthreads = 1;
	IntrinsicsType 		intrinType = IntrinsicsType::AVX;
	uint32_t			bitwidth = 0;
	int32_t 			my_tid = 0;			//thread logical id
	pthread_barrier_t*	barrier = NULL;		//for threads synchronization
	uint32_t		 	partition_fanout = 1;	//same as merge fan-in

	T* 					columnchunk = NULL;
	surrogate_t*		oidchunk = NULL;
	T*					tmp_partition = NULL;
	surrogate_t*		tmp_partition_oid = NULL;
	T*					tmp_sort = NULL;
	surrogate_t*		tmp_sort_oid = NULL;
	uint64_t			chunksize = 0;
	BAT_t<T>**			threadcolchunks = NULL;
	T*					sharedmergebuffer = NULL;	//normally size of L3-cache, shared by all threads, used in the multiway merge phase
	BAT_t<T>*			result = NULL;
}__attribute__((aligned(CACHE_LINE_SIZE)));


/** @} */

#endif /* TYPES_H */
