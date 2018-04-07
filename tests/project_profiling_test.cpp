#include	<cstdlib>
#include    "gtest/gtest.h"
#include	"src/types.h"
#include	<ctime>
#include	<iostream>
#include	<cstdint>
#include 	<string>
#include	<functional>
#include	<vector>
#include	"src/common.h"
#include 	<algorithm>
#include	"src/column.h"

namespace multiAttrSort{

class ProjectProfileTest: public ::testing::Test{
public:
    virtual void SetUp(){
    	std::srand(std::time(0));
    	num_rows_ = 16*1024*1024;
    	uint64_t i;

    	column64 = (uint64_t *)malloc_aligned(num_rows_ * sizeof(uint64_t));
        auto dice = std::bind(std::uniform_int_distribution<uint64_t>(
                                std::numeric_limits<uint64_t>::min(),
                                std::numeric_limits<uint64_t>::max()),
                                std::default_random_engine(std::time(0)));
        for (i = 0; i < num_rows_; ++i) {
        	column64[i] = dice();
        }


		column32 = (uint32_t *)malloc_aligned(num_rows_ * sizeof(uint32_t));
        auto dice2 = std::bind(std::uniform_int_distribution<uint32_t>(
                                std::numeric_limits<uint32_t>::min(),
                                std::numeric_limits<uint32_t>::max()),
                                std::default_random_engine(std::time(0)));
        for (i = 0; i < num_rows_; ++i) {
        	column32[i] = dice2();
        }

		column16 = (uint16_t *)malloc_aligned(num_rows_ * sizeof(uint16_t));
        auto dice3 = std::bind(std::uniform_int_distribution<uint16_t>(
                                std::numeric_limits<uint16_t>::min(),
                                std::numeric_limits<uint16_t>::max()),
                                std::default_random_engine(std::time(0)));
        for (i = 0; i < num_rows_; ++i) {
        	column16[i] = dice3();
        }

        oids_ = (surrogate_t *)malloc_aligned(num_rows_ * sizeof(surrogate_t));
        auto dice_oid = std::bind(std::uniform_int_distribution<uint32_t>(
                                0,
                                num_rows_-1),
                                std::default_random_engine(std::time(0)));
        for (i = 0; i < num_rows_; ++i) {
        	oids_[i] = dice_oid();	//random oids within range [0, num_rows-1]
        }
    }

    virtual void TearDown(){

    	free(column16);
    	free(column32);
    	free(column64);
    	free(oids_);
    }

protected:
    uint16_t *column16;
	uint32_t *column32;
	uint64_t *column64;
	surrogate_t *oids_;
	uint32_t num_rows_;

	template <class T>
	void project(BAT_t<T> *inValues, surrogate_t *inOrder, uint32_t num);
};

template <class T>
void ProjectProfileTest::project(BAT_t<T> *inValues, surrogate_t *inOrder, uint32_t num) {

	T *values_reordered = (T *)malloc_aligned(num * sizeof(T));
	for (surrogate_t oid = 0; oid < num; ++oid) {
		values_reordered[oid] = inValues->values[inOrder[oid]];
	}

	free(values_reordered);
}

#if 1
TEST_F(ProjectProfileTest, project_int64){
	uint32_t testnum = 200 * BLOCKSIZE<uint64_t>();
	uint32_t repeat = 10;
	BAT_t<uint64_t> *inBAT = new BAT_t<uint64_t>();
	uint32_t i;
	//inValues->num_elements = 0;
	//inValues->values = column64;
	//inValues->oids = oids;
	//uint32_t accu_count = 0;

	while (testnum > 0) {
		HybridTimer timer1;
		double sumcycles = 0.0;

		for (uint32_t rep_idx = 0; rep_idx < repeat; ++rep_idx) {

			uint64_t * inValues = (uint64_t *)malloc_aligned(testnum * sizeof(uint64_t));
			memcpy(inValues, column64, testnum * sizeof(uint64_t));
			surrogate_t *oids = (surrogate_t *)malloc_aligned(testnum * sizeof(surrogate_t));
			for (i = 0; i < testnum; ++i) {
				oids[i] = this->oids_[i] % testnum;
			}
			inBAT->num_elements = 0;
			inBAT->oids = NULL;
			inBAT->values = inValues;

			timer1.Start();
			project<uint64_t>(inBAT, oids, testnum);
			timer1.Stop();
			sumcycles += timer1.GetNumCycles();

			free(inValues);
			free(oids);
		}

		//total cost = sequential access + random access
		uint32_t COST_PER_SEQ_ACCESS_test = 4;
		uint32_t COST_PER_RAND_ACCESS_test = 16;

		uint32_t seq_cache_num = (double)(testnum * (sizeof(surrogate_t) + 2*sizeof(uint64_t))) / CACHE_LINE_SIZE;
		double seq_access_cost = seq_cache_num * COST_PER_SEQ_ACCESS_test;
		uint32_t L2_cache_block_num = testnum * (sizeof(surrogate_t) + 2*sizeof(uint64_t)) / L2_CACHE_SIZE + 1;
		double rand_access_perc = (double)(L2_cache_block_num-1)/(double)L2_cache_block_num;
		double rand_access_cost = testnum*COST_PER_RAND_ACCESS_test*rand_access_perc;
		double totalCost = seq_access_cost + rand_access_cost;

		std::cout << testnum << "\t" << sumcycles / repeat / testnum
				<< "\t" << totalCost / testnum << std::endl;


		if (testnum > BLOCKSIZE<uint64_t>()) {
			testnum -= BLOCKSIZE<uint64_t>();
		} else {
			repeat = 1000;
			testnum -= 256;
		}
	}

	free(inBAT);
}
#endif

#if 0
TEST_F(ProjectProfileTest, project_int32){
	uint32_t testnum = 200 * BLOCKSIZE<uint32_t>();
	uint32_t repeat = 10;
	BAT_t<uint32_t> *inBAT = new BAT_t<uint32_t>();
	uint32_t i;
	//inValues->num_elements = 0;
	//inValues->values = column32;
	//inValues->oids = oids;
	//uint32_t accu_count = 0;

	while (testnum > 0) {
		HybridTimer timer1;
		double sumcycles = 0.0;

		for (uint32_t rep_idx = 0; rep_idx < repeat; ++rep_idx) {

			uint32_t * inValues = (uint32_t *)malloc_aligned(testnum * sizeof(uint32_t));
			memcpy(inValues, column32, testnum * sizeof(uint32_t));
			surrogate_t *oids = (surrogate_t *)malloc_aligned(testnum * sizeof(surrogate_t));
			for (i = 0; i < testnum; ++i) {
				oids[i] = this->oids_[i] % testnum;
			}
			inBAT->num_elements = 0;
			inBAT->oids = NULL;
			inBAT->values = inValues;

			timer1.Start();
			project<uint32_t>(inBAT, oids, testnum);
			timer1.Stop();
			sumcycles += timer1.GetNumCycles();

			free(inValues);
			free(oids);
		}

		std::cout << testnum << "\t" << sumcycles / repeat / testnum << std::endl;


		if (testnum > BLOCKSIZE<uint32_t>()) {
			testnum -= BLOCKSIZE<uint32_t>();
		} else {
			repeat = 1000;
			testnum -= 256;
		}
	}

	free(inBAT);
}
#endif

#if 0
TEST_F(ProjectProfileTest, project_int16){
	uint32_t testnum = 200 * BLOCKSIZE<uint16_t>();
	uint32_t repeat = 10;
	BAT_t<uint16_t> *inBAT = new BAT_t<uint16_t>();
	uint32_t i;
	//inValues->num_elements = 0;
	//inValues->values = column16;
	//inValues->oids = oids;
	//uint32_t accu_count = 0;

	while (testnum > 0) {
		HybridTimer timer1;
		double sumcycles = 0.0;

		for (uint32_t rep_idx = 0; rep_idx < repeat; ++rep_idx) {

			uint16_t * inValues = (uint16_t *)malloc_aligned(testnum * sizeof(uint16_t));
			memcpy(inValues, column16, testnum * sizeof(uint16_t));
			surrogate_t *oids = (surrogate_t *)malloc_aligned(testnum * sizeof(surrogate_t));
			for (i = 0; i < testnum; ++i) {
				oids[i] = this->oids_[i] % testnum;
			}
			inBAT->num_elements = 0;
			inBAT->oids = NULL;
			inBAT->values = inValues;

			timer1.Start();
			project<uint16_t>(inBAT, oids, testnum);
			timer1.Stop();
			sumcycles += timer1.GetNumCycles();

			free(inValues);
			free(oids);
		}

		std::cout << testnum << "\t" << sumcycles / repeat / testnum << std::endl;


		if (testnum > BLOCKSIZE<uint16_t>()) {
			testnum -= BLOCKSIZE<uint16_t>();
		} else {
			repeat = 1000;
			testnum -= 256;
		}
	}

	free(inBAT);
}
#endif

}   //namespace
