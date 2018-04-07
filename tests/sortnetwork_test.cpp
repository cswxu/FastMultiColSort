#include 	"src/mergesort.h"
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
#include	"src/hybrid_timer.h"

namespace multiAttrSort{

class SortnetworkTest: public ::testing::Test{
public:
    virtual void SetUp(){
    	std::srand(std::time(0));
        //small_num_ = 0.25 * BLOCKSIZE;
    	small_num_ = 16;
        large_num_ = 64;
        num_bank16_ = 256;


        smallCol_val = (int64_t *)malloc_aligned(small_num_ * sizeof(int64_t));
        smallCol_oid = (surrogate_t *)malloc_aligned(small_num_ * sizeof(surrogate_t));


        largeCol_val = (int32_t *)malloc_aligned(large_num_ * sizeof(int32_t));
        largeCol_oid = (surrogate_t *)malloc_aligned(large_num_ * sizeof(surrogate_t));

        bank16Col_val = (int16_t *)malloc_aligned(num_bank16_ * sizeof(int16_t));
        bank16Col_oid = (surrogate_t *)malloc_aligned(num_bank16_ * sizeof(surrogate_t));

        uint64_t i;

        auto dice = std::bind(std::uniform_int_distribution<int64_t>(
                                std::numeric_limits<int64_t>::min(),
                                std::numeric_limits<int64_t>::max()),
                                std::default_random_engine(std::time(0)));
        uint64_t mask = (1ULL << 32) - 1;
        for (i = 0; i < small_num_; ++i) {
        	smallCol_val[i] = dice() & mask;
        	//smallCol_val[i] = dice();
        	smallCol_oid[i]	= (surrogate_t)i;
        }

        auto dice2 = std::bind(std::uniform_int_distribution<int32_t>(
                                std::numeric_limits<int32_t>::min(),
                                std::numeric_limits<int32_t>::max()),
                                std::default_random_engine(std::time(0)));
        for (i = 0; i < large_num_; ++i) {
        	largeCol_val[i] = dice2();
        	largeCol_oid[i]	= (surrogate_t)i;
        }

        auto dice3 = std::bind(std::uniform_int_distribution<int16_t>(
                                std::numeric_limits<int16_t>::min(),
                                std::numeric_limits<int16_t>::max()),
                                std::default_random_engine(std::time(0)));
        for (i = 0; i < num_bank16_; ++i) {
        	bank16Col_val[i] = dice3();
        	bank16Col_oid[i] = (surrogate_t)i;
        }
    }

    virtual void TearDown(){
        free(smallCol_val);
        free(smallCol_oid);
        free(largeCol_val);
        free(largeCol_oid);
        free(bank16Col_val);
        free(bank16Col_oid);
        //free(largeCol);
    }

protected:
    struct element_64_t {
    	surrogate_t oid;
    	int64_t value;
    };
    struct Compare {
    	bool operator()(element_64_t i, element_64_t j) {return i.value < j.value;}
    } mycomparison;
    int64_t *smallCol_val;
    surrogate_t *smallCol_oid;
    int32_t *largeCol_val;
    surrogate_t *largeCol_oid;
    int16_t * bank16Col_val;
    surrogate_t *bank16Col_oid;
    //element_t *reverseCol;
    //element_t *largeCol;
    uint32_t small_num_;
    uint64_t large_num_;
    uint32_t num_bank16_;
};

TEST_F(SortnetworkTest, inregister_sort_int16_aligned){
	uint32_t i;
	//original values and oids
#if 0
	for (i = 0; i < large_num_; i++) {
		std::cout << largeCol_val[i] << "\t" << largeCol_oid[i] << std::endl;
	}
#endif
	//copy the original values for the checking of the oids

	int16_t *temp_val;
	temp_val = (int16_t *)malloc_aligned(num_bank16_ * sizeof(int16_t));
	memcpy(temp_val, bank16Col_val, num_bank16_ * sizeof(int16_t));

	//do the sorting using sorting network
	HybridTimer timer1;
	timer1.Start();

	Mergesort::inregister_sort_int16_aligned(bank16Col_val, bank16Col_val,
	        		(int32_t*)bank16Col_oid, (int32_t*)bank16Col_oid);
	timer1.Stop();
	std::cout << "sorting 16*16 elements: " << timer1.GetNumCycles()/num_bank16_ << "cycles/tuple" << std::endl;

	//check value equivalence
	//uint32_t j;
	for (i = 0; i < num_bank16_; i += 16) {
#if 0
		std::cout << "run " << i << std::endl;
		for (j = 0; j < 16; j++) {
			std::cout << bank16Col_val[i+j] << "\t" << bank16Col_oid[i+j] << std::endl;
		}
#endif
		EXPECT_TRUE(std::is_sorted(bank16Col_val+i, bank16Col_val+i+16));
	}

	//check oid equivalence
	for (i = 0; i < num_bank16_; ++i) {
		EXPECT_EQ(bank16Col_val[i], temp_val[bank16Col_oid[i]]);
		//std::cout << largeCol_val[i] << "\t" << temp_val[largeCol_oid[i]] << "\t" << largeCol_oid[i] << std::endl;
		//std::cout << "oid: " << i << std::endl;
	}

	free(temp_val);
}

TEST_F(SortnetworkTest, inregister_sort_int32_aligned){
	uint32_t i;
	//original values and oids
#if 0
	for (i = 0; i < large_num_; i++) {
		std::cout << largeCol_val[i] << "\t" << largeCol_oid[i] << std::endl;
	}
#endif
	//copy the original values for the checking of the oids
	int32_t *temp_val;
	temp_val = (int32_t *)malloc_aligned(large_num_ * sizeof(int32_t));
	memcpy(temp_val, largeCol_val, large_num_ * sizeof(int32_t));

	//do the sorting using sorting network
	HybridTimer timer1;
	timer1.Start();
	Mergesort::inregister_sort_int32_aligned(largeCol_val, largeCol_val,
	        		(int32_t*)largeCol_oid, (int32_t*)largeCol_oid);
	timer1.Stop();
	std::cout << "sorting 8*8 elements: " << timer1.GetNumCycles()/large_num_ << "cycles/tuple" << std::endl;

	//check value equivalence
	//uint32_t j;
	for (i = 0; i < large_num_; i += 8) {
#if 0
		std::cout << "run " << i << std::endl;
		for (j = 0; j < 8; j++) {
			std::cout << largeCol_val[i+j] << "\t" << largeCol_oid[i+j] << std::endl;
		}
#endif
		EXPECT_TRUE(std::is_sorted(largeCol_val+i, largeCol_val+i+8));
	}

	//check oid equivalence
	for (i = 0; i < large_num_; ++i) {
		EXPECT_EQ(largeCol_val[i], temp_val[largeCol_oid[i]]);
		//std::cout << largeCol_val[i] << "\t" << temp_val[largeCol_oid[i]] << "\t" << largeCol_oid[i] << std::endl;
	}


	free(temp_val);
}

TEST_F(SortnetworkTest, inregister_sort_int64_aligned){
	uint32_t i;
	//original values and oids
#if 0
	for (i = 0; i < small_num_; i++) {
		std::cout << smallCol_val[i] << "\t" << smallCol_oid[i] << std::endl;
	}
#endif
	//copy the original values for the checking of the oids
	int64_t *temp_val;
	temp_val = (int64_t *)malloc_aligned(small_num_ * sizeof(int64_t));
	memcpy(temp_val, smallCol_val, small_num_ * sizeof(int64_t));

	//do the sorting using sorting network
	HybridTimer timer1;
	timer1.Start();
	Mergesort::inregister_sort_int64_aligned(smallCol_val, smallCol_val,
	        		(int32_t*)smallCol_oid, (int32_t*)smallCol_oid);
	timer1.Stop();
	std::cout << "sorting 4*4 elements: " << timer1.GetNumCycles()/small_num_ << "cycles/tuple" << std::endl;

	//check value equivalence
	//uint32_t j;
	for (i = 0; i < small_num_; i += 4) {
#if 0
		std::cout << "run " << i << std::endl;
		for (j = 0; j < 4; j++) {
			std::cout << smallCol_val[i+j] << "\t" << smallCol_oid[i+j] << std::endl;
		}
#endif
		EXPECT_TRUE(std::is_sorted(smallCol_val+i, smallCol_val+i+4));
	}

	//check oid equivalence
	for (i = 0; i < small_num_; ++i) {
		EXPECT_EQ(smallCol_val[i], temp_val[smallCol_oid[i]]);
		//std::cout << "transposed oid:" << smallCol_oid[i] << std::endl;
	}


	free(temp_val);
}

/**
 * TODO: CANNOT pass this test, when <small_num_> is small (e.g., hundreds), passed;
 * When <small_num_> is close to BLOCKSIZE, almost cannot be passed
 */
#if 0
TEST_F(SortnetworkTest, avxmergesort_rem_aligned){
	//do the sorting with STL
	std::vector<element_t> comparedCol;
	uint64_t i;
	element_t tmpElement;
	for (i = 0; i < small_num_; ++i) {
		tmpElement.value 	= smallCol[i].value;
		tmpElement.oid		= smallCol[i].oid;
		comparedCol.push_back(tmpElement);
	}

#if 1
	std::cout << "List before sorting:" << std::endl;
	for (i = 0; i < 10; ++i) {
		std::cout << comparedCol[i].value << "[" << comparedCol[i].oid << "]" << "\t";
	}
	std::cout << std::endl;
#endif

	std::sort(comparedCol.begin(), comparedCol.end(), mycomparison);

#if 1
	std::cout << "List after STL sorting:" << std::endl;
	for (i = 0; i < 10; ++i) {
		std::cout << comparedCol[i].value << "[" << comparedCol[i].oid << "]" << "\t";
	}
	std::cout << std::endl;
#endif


	//do the sorting with *avxmergesort_rem_aligned*
	int64_t *output = (int64_t *)malloc_aligned(small_num_ * sizeof(int64_t));
	int64_t *input 	= (int64_t *)smallCol;
	Mergesort::avxmergesort_rem_aligned(&input, &output, small_num_);

	//test whether the result is correct
	element_t *outputElement = (element_t *)output;
	for (i = 0; i < small_num_; ++i) {
		EXPECT_EQ(outputElement[i].value, comparedCol[i].value);
	}
#if 0
	for (i = 0; i < small_num_; ++i) {
		EXPECT_EQ(outputElement[i].oid, comparedCol[i].oid);
	}
#endif

#if 1
	std::cout << "List after avxmergesort_rem_aligned:" << std::endl;
	for (i = 0; i < 10; ++i) {
		std::cout << outputElement[i].value << "[" << outputElement[i].oid << "]" << "\t";
	}
	std::cout << std::endl;
#endif

	free(output);//something wrong here; may not clear about the pointer switch in avxmergesort_rem_aligned()
}
#endif

}   //namespace
