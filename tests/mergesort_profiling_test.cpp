#include 	"src/composer.h"
#include	"src/chain_composer.h"
#include    <cstdlib>
#include    "gtest/gtest.h"
#include	<cstdint>
#include	<ctime>
#include	<vector>
#include	<algorithm>
#include    <random>
#include    <functional>
#include    <unistd.h>
#include	<iostream>
#include 	"src/hybrid_timer.h"

namespace multiAttrSort{

class MergeSortPorfiling: public ::testing::Test{
public:
    virtual void SetUp(){
    	std::srand(std::time(0));
        num_ = 16*1024*1024;

        bank64Col_val = (uint64_t *)malloc_aligned(num_ * sizeof(uint64_t));
        bank64Col_oid = (surrogate_t *)malloc_aligned(num_ * sizeof(surrogate_t));


        bank32Col_val = (uint32_t *)malloc_aligned(num_ * sizeof(uint32_t));
        bank32Col_oid = (surrogate_t *)malloc_aligned(num_ * sizeof(surrogate_t));

        bank16Col_val = (uint16_t *)malloc_aligned(num_ * sizeof(uint16_t));
        bank16Col_oid = (surrogate_t *)malloc_aligned(num_ * sizeof(surrogate_t));

        uint64_t i;

        auto dice = std::bind(std::uniform_int_distribution<uint64_t>(
                                std::numeric_limits<uint64_t>::min(),
                                std::numeric_limits<uint64_t>::max()),
                                std::default_random_engine(std::time(0)));
        for (i = 0; i < num_; ++i) {
        	bank64Col_val[i] = dice();
        	bank64Col_oid[i] = (surrogate_t)i;
        }

        auto dice2 = std::bind(std::uniform_int_distribution<uint32_t>(
                                std::numeric_limits<uint32_t>::min(),
                                std::numeric_limits<uint32_t>::max()),
                                std::default_random_engine(std::time(0)));
        for (i = 0; i < num_; ++i) {
        	bank32Col_val[i] = dice2();
        	bank32Col_oid[i]	= (surrogate_t)i;
        }

        auto dice3 = std::bind(std::uniform_int_distribution<uint16_t>(
                                std::numeric_limits<uint16_t>::min(),
                                std::numeric_limits<uint16_t>::max()),
                                std::default_random_engine(std::time(0)));
        for (i = 0; i < num_; ++i) {
        	bank16Col_val[i] = dice3();
        	bank16Col_oid[i] = (surrogate_t)i;
        }
    }

    virtual void TearDown(){
        free(bank64Col_val);
        free(bank64Col_oid);
        free(bank32Col_val);
        free(bank32Col_oid);
        free(bank16Col_val);
        free(bank16Col_oid);
    }

protected:
    struct element_64_t {
    	surrogate_t oid;
    	uint64_t value;
    };
    struct Compare {
    	bool operator()(element_64_t i, element_64_t j) {return i.value < j.value;}
    } mycomparison;
    uint64_t *bank64Col_val;
    surrogate_t *bank64Col_oid;
    uint32_t *bank32Col_val;
    surrogate_t *bank32Col_oid;
    uint16_t * bank16Col_val;
    surrogate_t *bank16Col_oid;

    uint64_t num_;
};


#if 0
TEST_F(MergeSortPorfiling, subsort_mergesort_uint32_profile){

		uint32_t repeat = 1000;
		uint32_t idx;
		uint32_t testnum = 2;

		ChainComposer* composer = new ChainComposer();
		composer->setting_.nthreads = 1;
		composer->setting_.intrinType = IntrinsicsType::AVX;

		BAT_t<uint32_t> subcolumn;
		subcolumn.num_elements = num_;
		subcolumn.values = bank32Col_val;
		subcolumn.oids = bank32Col_oid;

		//allocate for original values
		uint32_t *origin_val;
		origin_val = (uint32_t *)malloc_aligned(num_ * sizeof(uint32_t));
		memcpy(origin_val, bank32Col_val, num_ * sizeof(uint32_t));

		while (testnum < 8192) {
			HybridTimer timer1;
			double sumcycles = 0.0;
			//recover the values


			for (uint32_t rep_idx = 0; rep_idx < repeat; ++rep_idx) {
				for (idx = 0; idx < testnum; ++idx) {
					subcolumn.values[idx] = origin_val[idx];
				}
				timer1.Start();
				composer->subsort_mergesort<uint32_t>(&subcolumn, 0, testnum, true, 32);
				timer1.Stop();
				sumcycles += timer1.GetNumCycles();
			}

			std::cout << testnum << "\t" << sumcycles / repeat / (testnum * log(testnum)) << std::endl;


			if (testnum < 784) {
				testnum++;
			} else {
				testnum += 128;
			}

			//testnum += 1024*1024;
		}

		delete composer;
		free(origin_val);
}
#endif

#if 0
TEST_F(MergeSortPorfiling, subsort_mergesort_uint64_profile){

		uint32_t repeat = 1;
		uint32_t idx;
		uint32_t testnum = 2;

		ChainComposer* composer = new ChainComposer();
		composer->setting_.nthreads = 1;
		composer->setting_.intrinType = IntrinsicsType::AVX;

		BAT_t<uint64_t> subcolumn;
		subcolumn.num_elements = num_;
		subcolumn.values = bank64Col_val;
		subcolumn.oids = bank64Col_oid;

		//allocate for original values
		uint64_t *origin_val;
		origin_val = (uint64_t *)malloc_aligned(num_ * sizeof(uint64_t));
		memcpy(origin_val, bank64Col_val, num_ * sizeof(uint64_t));

		while (testnum < num_) {
			HybridTimer timer1;
			double sumcycles = 0.0;
			//recover the values


			for (uint32_t rep_idx = 0; rep_idx < repeat; ++rep_idx) {
				for (idx = 0; idx < testnum; ++idx) {
					subcolumn.values[idx] = origin_val[idx];
				}
				timer1.Start();
				composer->subsort_mergesort<uint64_t>(&subcolumn, 0, testnum, true, 64);
				timer1.Stop();
				sumcycles += timer1.GetNumCycles();
			}

			std::cout << testnum << "\t" << sumcycles / repeat / (testnum * log(testnum)) << std::endl;


			//if (testnum < 512) {
			//	testnum++;
			//} else {
			//	testnum += 128;
			//}

			testnum += 1024*1024;
		}

		delete composer;
		free(origin_val);
}
#endif

#if 1
TEST_F(MergeSortPorfiling, subsort_mergesort_uint16_profile){

		uint32_t repeat = 1;
		uint32_t idx;
		uint32_t testnum = 2;

		ChainComposer* composer = new ChainComposer();
		composer->setting_.nthreads = 1;
		composer->setting_.intrinType = IntrinsicsType::AVX;

		BAT_t<uint16_t> subcolumn;
		subcolumn.num_elements = num_;
		subcolumn.values = bank16Col_val;
		subcolumn.oids = bank16Col_oid;

		//allocate for original values
		uint16_t *origin_val;
		origin_val = (uint16_t *)malloc_aligned(num_ * sizeof(uint16_t));
		memcpy(origin_val, bank16Col_val, num_ * sizeof(uint16_t));

		while (testnum < num_) {
			HybridTimer timer1;
			double sumcycles = 0.0;
			//recover the values


			for (uint32_t rep_idx = 0; rep_idx < repeat; ++rep_idx) {
				for (idx = 0; idx < testnum; ++idx) {
					subcolumn.values[idx] = origin_val[idx];
				}
				timer1.Start();
				composer->subsort_mergesort<uint16_t>(&subcolumn, 0, testnum, true, 16);
				timer1.Stop();
				sumcycles += timer1.GetNumCycles();
			}

			std::cout << testnum << "\t" << sumcycles / repeat / (testnum * log(testnum)) << std::endl;

			/*
			if (testnum < 512) {
				testnum++;
			} else {
				testnum += 128;
			}
			*/

			testnum += 1024*1024;
		}

		delete composer;
		free(origin_val);
}
#endif

}   //namespace
