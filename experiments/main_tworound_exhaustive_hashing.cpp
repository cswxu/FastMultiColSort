/**
  * 
  */

#include <stdio.h>
#include <stdint.h>
#include <sys/time.h>           /* gettimeofday */
#include <stdlib.h>
#include <limits.h>             /* LONG_MAX */
#include <math.h>
#include <cpuid.h>              /* for check_avx() */

/* #include <assert.h> */

#include "src/generator.h"
#include "src/params.h"             /* macro parameters */
#include "src/chain_composer.h"
#include "src/stitch_composer.h"
#include "src/composer.h"

#include "src/types.h"
#include "src/column.h"
#include <getopt.h>             /* getopt */

#include <sstream>
#include <cstring>             /* strlen(), memcpy() */
#include <string>
#include <map>
#include <iostream>

#include "src/hybrid_timer.h"
//#include "pcm/system_perf_monitor.h"

using namespace multiAttrSort;

/******************************************************************************
 *                                                                            *
 *                 Command-line handling & Driver Program                     *
 *                                                                            *
 ******************************************************************************/
struct cmdparam_t {
    
	/**************** Merge sort Related*********************/
	/* partitioning fanout, or merge fanin, e.g., 2^rdxbits */
    uint32_t part_fanout = 1;

	
	/******************** General ***************************/
    uint32_t nthreads = 1;
	/* the approach to compose sorting of multiple comlumns */
	ComposerType comptype = ComposerType::chaining;
	/* the underlying sort type */
	SortAlgo sortalgo = SortAlgo::mergesort;
	Baseline baseline = Baseline::adaptive;			//not used yet
	IntrinsicsType intrinsics_type = IntrinsicsType::AVX;
	PackOIDType packtype = PackOIDType::nonpack;	//always use nonpack currently
	int encode_oid = 0;
	/* 121 means asc order for 1st column, desc order for 2nd column, asc order for 3rd column */
	uint32_t column_asc_desc = 12;					//
	
	uint32_t repeat = 1;

	uint8_t stitch_direction = 0;	//0 is left, 1 is right
	uint32_t stitch_nbits = 0;

	uint32_t ordered = 1;	//1 is ordered (for ORDER-BY case), 0 is unordered (for GROUP-BY case)

	uint32_t hashtype = 1; 	//1: use reordering; 0: not use reordering, 2: partitioned reordering

	/******************** Generate Data **********************/
	uint64_t nrows = 64*1024*1024;
	uint32_t ncolumns = 2;
	/* 11021 means 1st column has 11 bit-width, 2nd column has 21 column width*/
	std::string column_bitwidth = "25B21";

	/* NO NEED TO SET THE FOLLOWING THREE WHEN RAED DATA AND DO THE SORTING*/
	/* data skew: determined by the zipf/10.0, i.e., [0, 2] -- 0 indicate uniform distribution */
	uint32_t zipf = 0;
	/* cardinality, i.e., # of distinct values, 
	 * represented as percentage [0,100] of maximum value (determined by bit width) */
	uint32_t cardinality = 100;
	uint32_t group_num = 2;
	
	/* NOTE: if gendata is not set, zipf and cardinality is not required*/
	int gendata = 0;
	std::string filename;
	
	/*TODO: 1) more data set parameters: entropy by Kurt92
	 *      2) zipf and cardinality info are for ALL columns now, diff conf for diff columns?*/
};

std::map<std::string, Baseline> baselineMap = {
	{"aligned32",	Baseline::aligned32},
	{"adaptive", 	Baseline::adaptive}
};

std::map<std::string, SortAlgo> sortalgoMap = {
	{"ms",	SortAlgo::mergesort},
	{"rs", 	SortAlgo::radixsort}
};

std::map<std::string, PackOIDType> packTypeMap = {
	{"pack", 	PackOIDType::pack},
	{"nonpack", 	PackOIDType::nonpack}
};

std::map<std::string, ComposerType> comptypeMap = {
	{"ch",	ComposerType::chaining},
	{"st", 	ComposerType::stitching}
};

std::map<std::string, IntrinsicsType> intrinsicsMap = {
	{"avx", 	IntrinsicsType::AVX},
	{"sse", 	IntrinsicsType::SSE},
	{"scalar", 	IntrinsicsType::SCALAR}
};

//extern char * optarg;
//extern int    optind, opterr, optopt;

/* command line handling functions */
void print_help(char * progname);

void parse_args(int argc, char ** argv, cmdparam_t& cmd_params);

void print_args(cmdparam_t& cmd_params);

void parse_columnwidths(bool *asc_desc_outp, uint32_t *widths_outp,
		uint32_t asc_desc_inp, std::string widths_inp, uint32_t num_columns);

int main(int argc, char *argv[])
{
    //struct timeval start, end;

    /* Command line parameters */
    cmdparam_t cmd_params;

    parse_args(argc, argv, cmd_params);


	dataset_t dataset;
	DataGenerator generator;
	if (cmd_params.gendata)	{ //just generate data and save to files
		/* set the data set parameters */
		dataset.nrows = cmd_params.nrows;
		dataset.ncolumns = cmd_params.ncolumns;
		dataset.column_bitwidth = cmd_params.column_bitwidth;
		dataset.zipf = cmd_params.zipf;
		dataset.cardinality = cmd_params.cardinality;
		dataset.group_num = cmd_params.group_num;

		generator.GenAndSaveAsFile(&dataset);

	} else {

#if 0
		B *column_values;	//to store all column values
		column_values = (value_t **) malloc_aligned(cmd_params.ncolumns * sizeof(value_t *));
		uint32_t i;
		for (i = 0; i < cmd_params.ncolumns; ++i) {
			column_values[i] = (value_t *) malloc_aligned(cmd_params.nrows * sizeof(value_t));
		}
#endif

		//parsing asc_desc and widths
		bool *column_asc_desc = (bool*) malloc_aligned(cmd_params.ncolumns * sizeof(bool));
		uint32_t *column_widths = (uint32_t*) malloc_aligned(cmd_params.ncolumns * sizeof(uint32_t));

		parse_columnwidths(column_asc_desc, column_widths,
				cmd_params.column_asc_desc, cmd_params.column_bitwidth, cmd_params.ncolumns);

		uint32_t i;
		for (i = 0; i < cmd_params.ncolumns; ++i) {
			printf("[INFO ] %u-th column with %u-bit width, to be sorted in %s order\n",
					i, column_widths[i],
					(column_asc_desc[i]) ? std::string("ascending").c_str() : std::string("descending").c_str()
					);
		}

		/** initialize the #ncolumns columns, each with #nrow values **/
		Column **columns;
		columns = (Column **) malloc_aligned(cmd_params.ncolumns * sizeof(Column *));

	    std::cout << "[INFO ] Repeat " << cmd_params.repeat << " times, start running..." << std::endl;
	    double totalSec = 0;
	    double totalCycle = 0;
	    double time_stitching = 0;
	    double time_tupleReconstruct = 0;
	    double time_orderGroupInfo = 0;
	    double time_firstRound = 0;

#if 1
		uint32_t ngroups = 0;	//number of groups
		float avg_group_size = 0.0;	//avg of group size
		uint32_t ngroups_incache = 0;	//number of groups that are cache-resident (smaller than half of L2-cache)
		uint32_t ngroups_outofcache = 0;
#endif

	    for (uint32_t tmp = 0; tmp < cmd_params.repeat; ++tmp) {

			/* allocate memory for each column; read column values from file */
	    	//printf("[INFO ] Start reading data...\n");
			generator.ReadFromFile(columns, cmd_params.filename,
					cmd_params.ncolumns, cmd_params.nrows, column_widths);

#if 0
			std::cout << "first column first 20 elements: " << std::endl;
			for (int tmpidx = 0; tmpidx < 20; ++tmpidx) {
				std::cout << ((uint64_t*)(columns[0]->GetColumn()))[tmpidx] << std::endl;
			}
#endif
			/* assign settings for sort algorithms */
			setting_t setting;
			setting.nthreads = cmd_params.nthreads;
			setting.intrinType = cmd_params.intrinsics_type;
			setting.partition_fanout = cmd_params.part_fanout;

			/* assign params for composition */
			compose_params_t compose_params;
			compose_params.sortalgo = cmd_params.sortalgo;
			compose_params.packtype = cmd_params.packtype;
			compose_params.is_oid_encoded = (cmd_params.encode_oid == 1);
			compose_params.ordered = cmd_params.ordered;
			compose_params.hash_type = cmd_params.hashtype;

			//deep copy
			compose_params.column_asc_desc =
					(bool*) malloc_aligned(cmd_params.ncolumns * sizeof(bool));
			compose_params.bitwidth =
					(uint32_t*) malloc_aligned(cmd_params.ncolumns * sizeof(uint32_t));
			for (i = 0; i < cmd_params.ncolumns; ++i) {
				(compose_params.column_asc_desc)[i] = column_asc_desc[i];
				(compose_params.bitwidth)[i] = column_widths[i];
				compose_params.colwidth_sum += column_widths[i];	//calculate the sum
			}
			std::cout << "[INFO ] Sum of all column width: " << compose_params.colwidth_sum << std::endl;

			//compose_params.stitch_nbits = cmd_params.stitch_nbits;
			//compose_params.stitch_direction = cmd_params.stitch_direction;

			//seed_generator(cmd_params.seed);
			Composer* composer = NULL;
			switch(cmd_params.comptype) {
				case ComposerType::chaining:
					composer = new ChainComposer(setting, cmd_params.ncolumns, cmd_params.nrows,
						compose_params, columns);
					break;
				case ComposerType::stitching:
					composer = new StitchComposer(setting, cmd_params.ncolumns, cmd_params.nrows,
							compose_params, columns);
					break;
			}

			//SystemPerfMonitor pm;
			//pm.Init();
			//printf("[INFO ] Start Sorting...\n");
		    HybridTimer totalTimer;
			totalTimer.Start();
			//pm.Start();

			//composer->SortAllColumns();
			composer->TwoRoundsExhaustive_hash();
			std::cout << "[INFO ] Finish multi-round hashing, clear memory..." << std::endl;

			totalTimer.Stop();

			//pm.Stop();
			totalSec += totalTimer.GetSeconds();
			totalCycle += double(totalTimer.GetNumCycles());
			time_stitching += composer->getSetting().time_stitching;
			time_tupleReconstruct += composer->getSetting().time_tupleReconstruct;
			time_orderGroupInfo += composer->getSetting().time_orderGroupInfo;
			time_firstRound += composer->getSetting().time_multipass_sort;

			ngroups = composer->getGroupInfo().ngroups;
			avg_group_size = (float)(composer->getGroupInfo().sum_group_size) / (float)(ngroups);
			ngroups_incache = composer->getGroupInfo().ngroups_incache;
			ngroups_outofcache = composer->getGroupInfo().ngroups_outofcache;

			//note: columns may be changed, so use the updated value (in composer) here
#if 1
			for (i = 0; i < composer->num_columns_; ++i) {
				switch(composer->column_values_[i]->GetWidth()) {
					case 1:
						//free((uint8_t *)composer->column_values_[i]->GetColumn());
					{
						std::cout << "[Error ] Not allow 8-bank size..." << std::endl;
						exit(EXIT_SUCCESS);
					}
						break;
					case 2:
						free((uint16_t *)composer->column_values_[i]->GetColumn());
						break;
					case 4:
						free((uint32_t *)composer->column_values_[i]->GetColumn());
						break;
					case 8:
						free((uint64_t *)composer->column_values_[i]->GetColumn());
						break;
				}
				delete composer->column_values_[i];
			}
#endif
			free(composer->compose_params_.column_asc_desc);	//note: updated because of the switched memory address
			free(composer->compose_params_.bitwidth);
			delete composer;

#if 0
			std::cout << tmp << "-th repeat: " << std::endl
					<<"Group info <#group, avg_group_sz, #groups_incache, #group_outofcache>:" << std::endl
					<< "<" << ngroups << ", " << avg_group_size << ", " << ngroups_incache<< ", " << ngroups_outofcache << ">"
					<< std::endl;
#endif
		}

#if 1
	    std::cout << "wall time (sec), total CPU cost (cycle/value), stitching, reconstruct, ordergroup, firstround, remainround" << std::endl;
	    std::cout << "do not make filter it as output," << totalSec/cmd_params.repeat << ", "
	                << totalCycle/cmd_params.nrows/cmd_params.repeat << ", "
					<< time_stitching/cmd_params.nrows/cmd_params.repeat << ", "
					<< time_tupleReconstruct/cmd_params.nrows/cmd_params.repeat << ", "
					<< time_orderGroupInfo/cmd_params.nrows/cmd_params.repeat << ", "
					<< time_firstRound/cmd_params.nrows/cmd_params.repeat << ", "
					<< (totalCycle-time_stitching-time_tupleReconstruct-time_orderGroupInfo-time_firstRound)/cmd_params.nrows/cmd_params.repeat
	                //<< std::endl << std::endl
					//<<"Group info <#group, avg_group_sz, #groups_incache, #group_outofcache>:" << std::endl
					<< ", "
					<< ngroups << ", " << avg_group_size << ", " << ngroups_incache<< ", " << ngroups_outofcache
					<< std::endl;
#endif

	
		/** destroy the column_values **/
#if 0
		for (i = 0; i < cmd_params.ncolumns; ++i) {
			free(column_values[i]);
		}
		free(column_values);
#endif
		free(column_asc_desc);
		free(column_widths);

		free(columns);
	}


    return 0;
}

void parse_columnwidths(bool *asc_desc_outp, uint32_t *widths_outp,
		uint32_t asc_desc_inp, std::string widths_inp, uint32_t num_columns) {
	/**** parse the (1) ascending/descending feature, (2) column bitwidth ****/
	uint32_t bw = 0;
	uint32_t desc_asc = 0;
	uint32_t tmpRemained_asc_desc = asc_desc_inp;

	std::string bitwidth_str = widths_inp;

	char bw_str_tmp[50];
	strcpy(bw_str_tmp, bitwidth_str.c_str());
	char *bw_str = strtok(bw_str_tmp, "B");
	int colIdx = 0;
	while (bw_str != NULL) {
		bw = atoi(bw_str);
		widths_outp[colIdx] = bw;

		desc_asc = tmpRemained_asc_desc % 10;
		asc_desc_outp[num_columns-colIdx-1] = (desc_asc==1);

#ifdef DEBUG
		printf("[INFO ] The %d-th column has bit_width=%d\n", colIdx, bw);
#endif

		tmpRemained_asc_desc /= 10;

		bw_str = strtok(NULL, "B");
		colIdx++;
	}
}

/* command line handling functions */
void
print_help(char * progname)
{
    printf("Usage: %s [options]\n", progname);

    printf("\n");

    printf("\
\tData generation options, with default values in []:      		           			\n\
\t\t--gen-data			Just generate data and save as files; do not sort [NO] 		\n\
\t\t-r --nrows=<R>      Number of rows in each column [64M]							\n\
\t\t-c --ncolumns=<C>   Number of columns to be sorted [2]							\n\
\t\t-w --bitwidth=<B>   the bit width of each column; separated by `B' [25B21]		\n\
\t\t-z --zipf=<Z>	    data skew: [0, 20], 0 indicate uniform distribution [0]		\n\
\t\t-C --cardinality=<> # of distinct values, as perc [0,100] of max value [100]	\n\
\t\t-g --ngroups=<G> 	Number of groups in each column [2]							\n\
\t\t-i --inputfile=<F>  read column values from files (--gen-data is not set)		\n\
	   																				\n\
\tGeneral options:																	\n\
\t\t--encode-oid		Encode OID with necessary bits; otherwise use 32|64 [NO]	\n\
\t\t-n --nthreads=<N>   Number of threads to use [1]	                    		\n\
\t\t-a --asc-desc=<A>   Indicate ascending(1)|descending(2) of each column [12]		\n\
\t\t-s --sortalgo=<..>  Implementation of underlying sort algorithms. [ms]			\n\
						   Acceptable sort algorithms are:							\n\
						   ms: Merge sort											\n\
						   rs: Radix sort											\n\
\t\t-p --comptype=<..>  Approaches to compose sorting of multiple columns. [ch]		\n\
	   					   Acceptable compose types are:							\n\
						   ch: chaining-based; baseline as in MonetDB				\n\
						   st: stitching-based; 									\n\
\t\t-P --oid-pack=<..>  Indicate how OID are processed during sorting. [nonpack]	\n\
						   Acceptable pack types are:								\n\
						   pack: 	pack the oid with (stitched) column values		\n\
						   nonpack:	OID stored as separated array	   				\n\
\t\t-b --baseline=<..>	Indicate the versions of baseline	[Disabled...]			\n\
    						Acceptable baseline type are:							\n\
    						aligned32: always stored as 32bit for each column		\n\
    						adaptive:  fit the lowest bank size						\n\
																					\n\
\t\t-x --stitchstyle=<>	Indicate # of bits to stitch left or right [Disabled...]	\n\
    						Acceptable examples are:								\n\
    						l3: stitch to left 3 bits								\n\
    						r4: stitch to right 4 bits								\n\
\t\t-o --ordered=<>		Indicate Group-by case (0) or Order-by case (1) [1]			\n\
																					\n\
\tMerge sort related options:														\n\
\t\t-f --partfanout=<F> Fanout of partition phase in merge sort. [1]				\n\
	   					   (also equivalent to fanin of multiway merge phase)		\n\
																					\n\
\tRadix sort related options:														\n\
                                                                               		\n\
																				    \n\
\tHashAggr related options:															\n\
\t\t-H --hashtype=<H>	Indicate whether multi-round hash apply reorder strategy	\n\
    		            Acceptable examples are: 									\n\
                        0: not reordering 											\n\
    			        1: reordering												\n\
\tBasic user options:                                                         		\n\
\t\t-h --help           Show this message                                   		\n");
}

void print_args(cmdparam_t& cmd_params) {

    printf("[INFO ] #threads=%d\n[INFO ] comptype=%s\n[INFO ] sortalgo=%s\n\[INFO ] #columns=%d\n[INFO ] #rows=%lu\n",
			cmd_params.nthreads,
			(cmd_params.comptype == ComposerType::chaining) ? 
			std::string("chaining").c_str() : std::string("stitching").c_str(),
			(cmd_params.sortalgo == SortAlgo::mergesort) ?
			std::string("mergesort").c_str() : std::string("radixsort").c_str(),
			//(cmd_params.packtype == PackOIDType::pack) ?
			//std::string("pack").c_str() : std::string("nonpack").c_str(),
			cmd_params.ncolumns,
            cmd_params.nrows);
    //std::cout << "[INFO ] Stitch " << cmd_params.stitch_nbits << " bits to the " <<  ((cmd_params.stitch_direction) ?
	//		std::string("right").c_str() : std::string("left").c_str()) << std::endl;
    std::cout << "[INFO ] In " << ((cmd_params.ordered == 1) ?
    		std::string("ORDER-BY").c_str() : std::string("GROUP-BY").c_str()) << " case" << std::endl;
    std::cout << "[INFO ] Hash type: " << ((cmd_params.hashtype == 1) ?
    		std::string("reorder based").c_str() : std::string("not reorder based").c_str()) << std::endl;
}

void parse_args(int argc, char ** argv, cmdparam_t& cmd_params)
{
    int c;
	std::string s;
	
	static int gendata_flag = 0;
	static int encode_oid_flag = 0;

    while(1) {
		static struct option long_options[] =
			{
                /* These options set a flag. Do not require arguments */
                {"gen-data",   no_argument,    &gendata_flag,   1},
                {"encode-oid", no_argument,    &encode_oid_flag, 1},
                {"help",       no_argument,    0, 'h'},
                /* These options require arguments. */
                {"nthreads",	required_argument, 0, 'n'},
                {"nrows",  		required_argument, 0, 'r'},
                {"ncolumns",  	required_argument, 0, 'c'},
                {"bitwidth", 	required_argument, 0, 'w'},
				{"zipf", 		required_argument, 0, 'z'},
				{"cardinality", required_argument, 0, 'C'},
				{"inputfile", 	required_argument, 0, 'i'},
                {"asc-desc",  	required_argument, 0, 'a'},
                {"sortalgo",  	required_argument, 0, 's'},
                {"comptype",    required_argument, 0, 'p'},
                {"partfanout",  required_argument, 0, 'f'},
                {"oid-pack",	required_argument, 0, 'P'},
				{"baseline", 	required_argument, 0, 'b'},
				{"stitchstyle", required_argument, 0, 'x'},
				{"ngroups", 	required_argument, 0, 'g'},
				{"ordered", 	required_argument, 0, 'o'},
				{"hashtype", 	required_argument, 0, 'H'}
            };
        
		int option_index = 0;

        c = getopt_long(argc, argv, "h:n:r:c:w:z:C:f:a:s:p:i:P:b:x:g:o:H",
						long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;
        switch (c)
        {
          case 'h':
              print_help(argv[0]);
              exit(EXIT_SUCCESS);
              break;
          
		  case 'n':
              cmd_params.nthreads = atoi(optarg);
              break;

		  case 'r':
			  cmd_params.nrows = (uint64_t)(atol(optarg));
			  break;

          case 'c':
              cmd_params.ncolumns = atoi(optarg);
              break;

          case 'w':
			  /* overflow may happen depend on ncolumns*/
              cmd_params.column_bitwidth = std::string(optarg);
              break;

          case 'a':
              cmd_params.column_asc_desc = atoi(optarg);
              break;

          case 's':
			  s = std::string(optarg);
			  if (sortalgoMap.find(s) == sortalgoMap.end()) {
			      printf("[ERROR ] Sort Algorithm `%s' does not exist.\n ", s.c_str());
			  	  exit(EXIT_SUCCESS);	  
			  } else {
			   	  cmd_params.sortalgo = sortalgoMap[s];
			  }
              break;

          case 'x':
          {
		  	  s = std::string(optarg);
		  	  //std::cout << s << std::endl;
		  	  cmd_params.stitch_direction = (s[s.length()-1] == 'r');	//0 is left, 1 is right
		  	  cmd_params.stitch_nbits = atoi(s.c_str());
		  	  //std::cout << cmd_params.stitch_direction << "\t" << cmd_params.stitch_nbits << std::endl;

		  	  break;
          }

          case 'p':
			  s = std::string(optarg);
			  if (comptypeMap.find(s) == comptypeMap.end()) {
			      printf("[ERROR ] Composition type `%s' does not exist.\n ", s.c_str());
			  	  exit(EXIT_SUCCESS);	  
			  } else {
			   	  cmd_params.comptype = comptypeMap[s];
			  }
              break;

          case 'b':
			  s = std::string(optarg);
			  if (baselineMap.find(s) == baselineMap.end()) {
			      printf("[ERROR ] Baseline type `%s' does not exist.\n ", s.c_str());
			  	  exit(EXIT_SUCCESS);
			  } else {
			   	  cmd_params.baseline = baselineMap[s];
			  }
              break;

          case 'f':
              cmd_params.part_fanout = atoi(optarg);
              break;

		  case 'P':
			  s = std::string(optarg);
			  if (packTypeMap.find(s) == packTypeMap.end()) {
			      printf("[ERROR ] OID pack type `%s' does not exist.\n ", s.c_str());
			  	  exit(EXIT_SUCCESS);	  
			  } else {
			   	  cmd_params.packtype = packTypeMap[s];
			  }
              break;
			
		  case 'z':
			  cmd_params.zipf = atoi(optarg);
			  break;
		  case 'C':
			  cmd_params.cardinality = atoi(optarg);
			  break;	  
		  case 'i':
			  cmd_params.filename = std::string(optarg);	  
			  break;
		  case 'g':
			  cmd_params.group_num = atoi(optarg);
			  break;
		  case 'o':
			  cmd_params.ordered = atoi(optarg);
			  break;
		  case 'H':
			  cmd_params.hashtype = atoi(optarg);
			  break;
          default:
              break;
        }
    }

    /* make sure part_fanout is greater the num-threads */
    if(cmd_params.part_fanout < cmd_params.nthreads){
        printf("[ERROR] Partitioning fan-out of merge sort must be >= num-threads.\n");
        exit(EXIT_SUCCESS);
    }
		
	cmd_params.gendata 		= gendata_flag;
	cmd_params.encode_oid 	= encode_oid_flag;

	print_args(cmd_params);
}
