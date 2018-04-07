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

class SubgroupTest: public ::testing::Test{
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


#if 1
TEST_F(SubgroupTest, subgroup_test){

		ChainComposer* composer = new ChainComposer();
		composer->setting_.nthreads = 1;
		composer->setting_.intrinType = IntrinsicsType::AVX;

		const size_t tuple_num = 9;
		uint32_t *col1_values = (uint32_t *) malloc_aligned(tuple_num * sizeof(uint32_t));
		uint32_t col1_origin_values[tuple_num] = {100, 200, 200, 200, 100, 100, 300, 100, 100};
		memcpy(col1_values, col1_origin_values, tuple_num * sizeof(uint32_t));

		uint32_t *col2_values = (uint32_t *) malloc_aligned(tuple_num * sizeof(uint32_t));
		uint32_t col2_origin_values[tuple_num] = {1000, 2000, 2000, 1000, 2000, 1000, 3000, 4000, 4000};
		memcpy(col2_values, col2_origin_values, tuple_num * sizeof(uint32_t));

		surrogate_t *outGroup = (surrogate_t *)malloc_aligned(tuple_num * sizeof(surrogate_t));
		surrogate_t *inGroup = NULL;


		composer->num_rows_ = tuple_num;

		composer->hashByGroup<uint32_t>(outGroup, col1_values, inGroup);

		inGroup = (surrogate_t *)malloc_aligned(tuple_num * sizeof(surrogate_t));

		swap<surrogate_t>(&inGroup, &outGroup);

		composer->hashByGroup<uint32_t>(outGroup, col2_values, inGroup);

		//output group info
		size_t i = 0;
		for (i = 0; i < tuple_num; ++i) {
			std::cout << outGroup[i] << std::endl;
		}

		std::cout << "stitching" << std::endl;

		/** stitch together**/
		surrogate_t *outGroup_stitched = (surrogate_t *)malloc_aligned(tuple_num * sizeof(surrogate_t));
		surrogate_t *inGroup_stitched = NULL;

		uint32_t *col_values_stitched = (uint32_t *) malloc_aligned(tuple_num * sizeof(uint32_t));

		for (i = 0; i < tuple_num; ++i) {
			col_values_stitched[i] = col1_values[i] * 8192 + col2_values[i];
		}

		composer->hashByGroup<uint32_t>(outGroup_stitched, col_values_stitched, inGroup_stitched);

		//output group info
		for (i = 0; i < tuple_num; ++i) {
			std::cout << outGroup_stitched[i] << std::endl;
			EXPECT_EQ(outGroup_stitched[i], outGroup[i]);
		}

		free(outGroup);
		free(outGroup_stitched);
		free(inGroup);
		free(col1_values);
		free(col2_values);
		free(col_values_stitched);

		delete composer;
}
#endif

}   //namespace
