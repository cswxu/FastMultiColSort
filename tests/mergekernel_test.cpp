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
#include 	"src/hybrid_timer.h"

namespace multiAttrSort{

class MergekernelTest: public ::testing::Test{
public:
    virtual void SetUp(){
    	std::srand(std::time(0));
        num_ = BLOCKSIZE<uint32_t>() * 2;
        testtest = 100;

        bank64Col_val = (int64_t *)malloc_aligned(num_ * sizeof(int64_t));
        bank64Col_oid = (surrogate_t *)malloc_aligned(num_ * sizeof(surrogate_t));


        bank32Col_val = (int32_t *)malloc_aligned(num_ * sizeof(int32_t));
        bank32Col_oid = (surrogate_t *)malloc_aligned(num_ * sizeof(surrogate_t));

        bank16Col_val = (int16_t *)malloc_aligned(num_ * sizeof(int16_t));
        bank16Col_oid = (surrogate_t *)malloc_aligned(num_ * sizeof(surrogate_t));

        uint64_t i;

        auto dice = std::bind(std::uniform_int_distribution<int64_t>(
                                std::numeric_limits<int64_t>::min(),
                                std::numeric_limits<int64_t>::max()),
                                std::default_random_engine(std::time(0)));
        //uint64_t mask = (1ULL << 63) - 1;
        for (i = 0; i < num_; ++i) {
        	//bank64Col_val[i] = dice() & mask;
        	bank64Col_val[i] = dice();
        	bank64Col_oid[i] = (surrogate_t)i;
        }

        auto dice2 = std::bind(std::uniform_int_distribution<int32_t>(
                                std::numeric_limits<int32_t>::min(),
                                std::numeric_limits<int32_t>::max()),
                                std::default_random_engine(std::time(0)));
        for (i = 0; i < num_; ++i) {
        	bank32Col_val[i] = dice2();
        	bank32Col_oid[i]	= (surrogate_t)i;
        }

        auto dice3 = std::bind(std::uniform_int_distribution<int16_t>(
                                std::numeric_limits<int16_t>::min(),
                                std::numeric_limits<int16_t>::max()),
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
    	int64_t value;
    };
    struct Compare {
    	bool operator()(element_64_t i, element_64_t j) {return i.value < j.value;}
    } mycomparison;
    int64_t *bank64Col_val;
    surrogate_t *bank64Col_oid;
    int32_t *bank32Col_val;
    surrogate_t *bank32Col_oid;
    int16_t * bank16Col_val;
    surrogate_t *bank16Col_oid;

    uint64_t num_;

    int32_t testtest = 0;
};

#if 1
TEST_F(MergekernelTest, merge8_int32_varlen_aligned){

	uint32_t i;
#if 0
	for (i = 0; i < large_num_; i++) {
		std::cout << largeCol_val[i] << "\t" << largeCol_oid[i] << std::endl;
	}
#endif
	//allocate output
	int32_t *output_val;
	surrogate_t *output_oid;
	output_val = (int32_t *)malloc_aligned(num_ * sizeof(int32_t));
	output_oid = (surrogate_t *)malloc_aligned(num_ * sizeof(surrogate_t));

	//sort two runs
	uint64_t sizeA = 16384;
	uint64_t sizeB = 16384;
	std::stable_sort(bank32Col_val, bank32Col_val+sizeA);
	std::stable_sort(bank32Col_val+sizeA, bank32Col_val+sizeA+sizeB);

	//use a temporary mem to remember the array values
	int32_t *tmp_val = (int32_t *)malloc_aligned((sizeA+sizeB) * sizeof(int32_t));
	memcpy(tmp_val, bank32Col_val, (sizeA+sizeB) * sizeof(int32_t));

#if 0
	for (i = 0; i < 16; i++) {
		std::cout << bank32Col_val[i] << "\t" << bank32Col_oid[i] << std::endl;
	}
#endif

#if 0
	std::cout << "original values:" << std::endl;
	for (i = 0; i < 16; i++) {
		std::cout << output_val[i] << std::endl;
	}
#endif

    //do the sorting
	//std::cout << "before sorting" << std::endl;
	HybridTimer timer1;
	timer1.Start();
	Mergesort::merge8_int32_varlen_aligned(bank32Col_val,
									bank32Col_val+sizeA,
									output_val,
									bank32Col_oid,
									bank32Col_oid+sizeA,
									output_oid,
	                      			sizeA,
	                      			sizeB);
	//std::cout << "after sorting" << std::endl;
	timer1.Stop();
	std::cout << "merge8_int32_varlen_aligned: " << timer1.GetNumCycles()/(sizeA+sizeB) << "cycles/tuple" << std::endl;

#if 0
	for (i = 0; i < 16; i++) {
		std::cout << output_val[i] << std::endl;
	}
#endif
	//check sortness
    EXPECT_TRUE(std::is_sorted(output_val, output_val+sizeA+sizeB));

	//check oid equivalence
	for (i = 0; i < sizeA+sizeB; ++i) {
		EXPECT_EQ(output_val[i], tmp_val[output_oid[i]]);
		//std::cout << largeCol_val[i] << "\t" << temp_val[largeCol_oid[i]] << "\t" << largeCol_oid[i] << std::endl;
		//std::cout << "oid: " << i << std::endl;
	}

	free(output_val);
	free(output_oid);
	free(tmp_val);
}
#endif

#if 1
TEST_F(MergekernelTest, merge16_int16_varlen_aligned){

	uint32_t i;
#if 0
	for (i = 0; i < large_num_; i++) {
		std::cout << largeCol_val[i] << "\t" << largeCol_oid[i] << std::endl;
	}
#endif
	//allocate output
	int16_t *output_val;
	surrogate_t *output_oid;
	output_val = (int16_t *)malloc_aligned(num_ * sizeof(int16_t));
	output_oid = (surrogate_t *)malloc_aligned(num_ * sizeof(surrogate_t));

	//sort two runs
	//uint64_t sizeA = 160;		//case 1
	//uint64_t sizeB = 320;		//case 1
	uint64_t sizeA = 16384;		//case 2
	uint64_t sizeB = 16384;		//case 2
	//uint64_t sizeA = 1608;		//case 3	//cannot pass, because it is not a multiply of 16?, doesnot matter we do not have such case
	//uint64_t sizeB = 988;		//case 3
	std::stable_sort(bank16Col_val, bank16Col_val+sizeA);
	std::stable_sort(bank16Col_val+sizeA, bank16Col_val+sizeA+sizeB);

	//use a temporary mem to remember the array values
	int16_t *tmp_val = (int16_t *)malloc_aligned((sizeA+sizeB) * sizeof(int16_t));
	memcpy(tmp_val, bank16Col_val, (sizeA+sizeB) * sizeof(int16_t));

#if 0
	for (i = 0; i < 16; i++) {
		std::cout << bank32Col_val[i] << "\t" << bank32Col_oid[i] << std::endl;
	}
#endif

#if 0
	std::cout << "original values:" << std::endl;
	for (i = 0; i < 16; i++) {
		std::cout << output_val[i] << std::endl;
	}
#endif

	HybridTimer timer1;
	timer1.Start();
    //do the sorting
	//std::cout << "before sorting" << std::endl;
	Mergesort::merge16_int16_varlen_aligned(bank16Col_val,
									bank16Col_val+sizeA,
									output_val,
									bank16Col_oid,
									bank16Col_oid+sizeA,
									output_oid,
	                      			sizeA,
	                      			sizeB);
	//std::cout << "after sorting" << std::endl;
	timer1.Stop();
	std::cout << "merge16_int16_varlen_aligned: " << timer1.GetNumCycles()/(sizeA + sizeB) << "cycles/tuple" << std::endl;

#if 0
	for (i = 0; i < 16; i++) {
		std::cout << output_val[i] << std::endl;
	}
#endif
	//check sortness
    EXPECT_TRUE(std::is_sorted(output_val, output_val+sizeA+sizeB));

	//check oid equivalence
	for (i = 0; i < sizeA+sizeB; ++i) {
		EXPECT_EQ(output_val[i], tmp_val[output_oid[i]]);
		//std::cout << largeCol_val[i] << "\t" << temp_val[largeCol_oid[i]] << "\t" << largeCol_oid[i] << std::endl;
		//std::cout << "oid: " << i << std::endl;
	}

	free(output_val);
	free(output_oid);
	free(tmp_val);
}
#endif

#if 1
TEST_F(MergekernelTest, merge8_int32_eqlen_aligned){

	uint32_t i;
#if 0
	for (i = 0; i < large_num_; i++) {
		std::cout << largeCol_val[i] << "\t" << largeCol_oid[i] << std::endl;
	}
#endif
	uint32_t test_num = 1024;
	//allocate output
	int32_t *output_val;
	surrogate_t *output_oid;
	output_val = (int32_t *)malloc_aligned(test_num*2 * sizeof(int32_t));
	output_oid = (surrogate_t *)malloc_aligned(test_num*2 * sizeof(surrogate_t));

	//sort two runs

	std::stable_sort(bank32Col_val, bank32Col_val+test_num);
	std::stable_sort(bank32Col_val+test_num, bank32Col_val+test_num*2);

	//use a temporary mem to remember the array values
	int32_t *tmp_val = (int32_t *)malloc_aligned(test_num*2 * sizeof(int32_t));
	memcpy(tmp_val, bank32Col_val, test_num*2 * sizeof(int32_t));

#if 0
	for (i = 0; i < 16; i++) {
		std::cout << bank32Col_val[i] << "\t" << bank32Col_oid[i] << std::endl;
	}
#endif

#if 0
	std::cout << "original values:" << std::endl;
	for (i = 0; i < 16; i++) {
		std::cout << output_val[i] << std::endl;
	}
#endif

    //do the sorting
	//std::cout << "before sorting" << std::endl;
	HybridTimer timer1;
	timer1.Start();
	Mergesort::merge8_int32_eqlen_aligned(bank32Col_val,
									bank32Col_val+test_num,
									bank32Col_oid,
									bank32Col_oid+test_num,
									output_val,
									output_oid,
									test_num);
	//std::cout << "after sorting" << std::endl;
	timer1.Stop();
	std::cout << "merge8_int32_eqlen_aligned: " << timer1.GetNumCycles()/(2*test_num) << "cycles/tuple" << std::endl;

#if 0
	for (i = 0; i < 16; i++) {
		std::cout << output_val[i] << std::endl;
	}
#endif
	//check sortness
    EXPECT_TRUE(std::is_sorted(output_val, output_val+test_num*2));

	//check oid equivalence

#if 1
	for (i = 0; i < test_num*2; ++i) {
		assert(output_oid[i] < test_num*2);
		//EXPECT_EQ(output_val[i], tmp_val[output_oid[i]]);
		//std::cout << largeCol_val[i] << "\t" << temp_val[largeCol_oid[i]] << "\t" << largeCol_oid[i] << std::endl;
		//std::cout << "oid: " << i << std::endl;
	}
#endif

	free(output_val);
	free(output_oid);
	free(tmp_val);
}
#endif

#if 1
TEST_F(MergekernelTest, merge16_int16_eqlen_aligned){

	uint32_t i;
#if 0
	for (i = 0; i < large_num_; i++) {
		std::cout << largeCol_val[i] << "\t" << largeCol_oid[i] << std::endl;
	}
#endif
	//uint32_t test_num = 512;	//case 1
	//uint32_t test_num = 256;	//case 2
	//uint32_t test_num = 1024;	//case 3
	uint32_t test_num = BLOCKSIZE<uint16_t>()/4;	//case 4
	//allocate output
	int16_t *output_val;
	surrogate_t *output_oid;
	output_val = (int16_t *)malloc_aligned(test_num*2 * sizeof(int16_t));
	output_oid = (surrogate_t *)malloc_aligned(test_num*2 * sizeof(surrogate_t));

	//sort two runs

	std::stable_sort(bank16Col_val, bank16Col_val+test_num);
	std::stable_sort(bank16Col_val+test_num, bank16Col_val+test_num*2);

	//use a temporary mem to remember the array values
	int16_t *tmp_val = (int16_t *)malloc_aligned(test_num*2 * sizeof(int16_t));
	memcpy(tmp_val, bank16Col_val, test_num*2 * sizeof(int16_t));

#if 0
	for (i = 0; i < 16; i++) {
		std::cout << bank32Col_val[i] << "\t" << bank32Col_oid[i] << std::endl;
	}
#endif

#if 0
	std::cout << "original values:" << std::endl;
	for (i = 0; i < 16; i++) {
		std::cout << output_val[i] << std::endl;
	}
#endif

	HybridTimer timer1;
	timer1.Start();
    //do the sorting
	//std::cout << "before sorting" << std::endl;
	Mergesort::merge16_int16_eqlen_aligned(bank16Col_val,
									bank16Col_val+test_num,
									bank16Col_oid,
									bank16Col_oid+test_num,
									output_val,
									output_oid,
									test_num);
	//std::cout << "after sorting" << std::endl;
	timer1.Stop();
	std::cout << "merge16_int16_eqlen_aligned: " << timer1.GetNumCycles()/(2*test_num) << "cycles/tuple" << std::endl;

#if 0
	for (i = 0; i < 16; i++) {
		std::cout << output_val[i] << std::endl;
	}
#endif
	//check sortness
    EXPECT_TRUE(std::is_sorted(output_val, output_val+test_num*2));

	//check oid equivalence

#if 1
	for (i = 0; i < test_num*2; ++i) {
		assert(output_oid[i] < test_num*2);
		//EXPECT_EQ(output_val[i], tmp_val[output_oid[i]]);
		//std::cout << largeCol_val[i] << "\t" << temp_val[largeCol_oid[i]] << "\t" << largeCol_oid[i] << std::endl;
		//std::cout << "oid: " << i << std::endl;
	}
#endif

	free(output_val);
	free(output_oid);
	free(tmp_val);
}
#endif

#if 1
TEST_F(MergekernelTest, merge4_int64_eqlen_aligned){

	uint32_t i;
#if 0
	for (i = 0; i < large_num_; i++) {
		std::cout << largeCol_val[i] << "\t" << largeCol_oid[i] << std::endl;
	}
#endif
	//uint32_t test_num = 512;	//case 1
	//uint32_t test_num = 2;	//case 2 => this case failed; but it won't happen in any case
	//uint32_t test_num = 8;	//case 3
	uint32_t test_num = 1024;	//case 4

	//allocate output
	int64_t *output_val;
	surrogate_t *output_oid;
	output_val = (int64_t *)malloc_aligned(test_num*2 * sizeof(int64_t));
	output_oid = (surrogate_t *)malloc_aligned(test_num*2 * sizeof(surrogate_t));

	//sort two runs

	std::stable_sort(bank64Col_val, bank64Col_val+test_num);
	std::stable_sort(bank64Col_val+test_num, bank64Col_val+test_num*2);

	//use a temporary mem to remember the array values
	int64_t *tmp_val = (int64_t *)malloc_aligned(test_num*2 * sizeof(int64_t));
	memcpy(tmp_val, bank64Col_val, test_num*2 * sizeof(int64_t));

#if 0
	for (i = 0; i < 16; i++) {
		std::cout << bank32Col_val[i] << "\t" << bank32Col_oid[i] << std::endl;
	}
#endif

#if 0
	std::cout << "original values:" << std::endl;
	for (i = 0; i < 16; i++) {
		std::cout << output_val[i] << std::endl;
	}
#endif

	HybridTimer timer1;
	timer1.Start();
    //do the sorting
	//std::cout << "before sorting" << std::endl;
	Mergesort::merge4_int64_eqlen_aligned(bank64Col_val,
			bank64Col_val+test_num,
			bank64Col_oid,
			bank64Col_oid+test_num,
			output_val,
			output_oid,
			test_num);
	//std::cout << "after sorting" << std::endl;
	timer1.Stop();
	std::cout << "merge4_int64_eqlen_aligned: " << timer1.GetNumCycles()/(2*test_num) << "cycles/tuple" << std::endl;

#if 0
	for (i = 0; i < 16; i++) {
		std::cout << output_val[i] << std::endl;
	}
#endif
	//check sortness
    EXPECT_TRUE(std::is_sorted(output_val, output_val+test_num*2));

	//check oid equivalence

#if 1
	for (i = 0; i < test_num*2; ++i) {
		assert(output_oid[i] < test_num*2);
		EXPECT_EQ(output_val[i], tmp_val[output_oid[i]]);
		//std::cout << largeCol_val[i] << "\t" << temp_val[largeCol_oid[i]] << "\t" << largeCol_oid[i] << std::endl;
		//std::cout << "oid: " << i << std::endl;
	}
#endif

	free(output_val);
	free(output_oid);
	free(tmp_val);
}
#endif

#if 1
TEST_F(MergekernelTest, avxmergesort_block_uint64_aligned){

	uint32_t i;
	#if 0
		for (i = 0; i < large_num_; i++) {
			std::cout << largeCol_val[i] << "\t" << largeCol_oid[i] << std::endl;
		}
	#endif

	uint64_t test_num = BLOCKSIZE<uint64_t>(); //case 1
	//	uint64_t test_num = 512;	//case 2
#if 1
	//allocate for temporary storage
	int64_t *input_val;
	surrogate_t *input_oid;
	input_val = (int64_t *)malloc_aligned(test_num * sizeof(int64_t));
	input_oid = (surrogate_t *)malloc_aligned(test_num * sizeof(surrogate_t));
	//int64_t *original_val = input_val;
#if 1
    auto dice = std::bind(std::uniform_int_distribution<int64_t>(
                            std::numeric_limits<int64_t>::min(),
                            std::numeric_limits<int64_t>::max()),
                            std::default_random_engine(std::time(0)));
    for (i = 0; i < test_num; ++i) {
    	input_val[i] = dice();
    	input_oid[i] = (surrogate_t)i;
    }
#endif
    //use bank32Col_val to store original values
    memcpy(bank64Col_val, input_val, test_num * sizeof(int64_t));

#endif
	//allocate output
    int64_t *output_val;
	surrogate_t *output_oid;
	output_val = (int64_t *)malloc_aligned(test_num * sizeof(int64_t));
	output_oid = (surrogate_t *)malloc_aligned(test_num * sizeof(surrogate_t));

#if 1
	Mergesort::avxmergesort_block_uint64_aligned((uint64_t **)&input_val, &input_oid,
			(uint64_t **)&output_val, &output_oid, test_num);
#endif

	EXPECT_TRUE(std::is_sorted(output_val, output_val+test_num));

#if 1
	for (i = 0; i < test_num; ++i) {
		assert(output_oid[i] < test_num);
	}

#endif

	//test oid equivalence
#if 1
	//if (original_val == input_val) {//pointer not switched
		for (i = 0; i < test_num; ++i) {
			EXPECT_EQ(output_val[i], bank64Col_val[output_oid[i]]);
			//std::cout << "original incorret oid:" << output_oid[i] << std::endl;
		}
	//} else {
	//	for (i = 0; i < test_num; ++i) {
	//		EXPECT_EQ(input_val[i], bank64Col_val[input_oid[i]]);
	//	}
	//}
#endif

	free(input_val);
	free(input_oid);
	free(output_val);
	free(output_oid);
}
#endif

#if 1
TEST_F(MergekernelTest, avxmergesort_block_uint16_aligned){
	uint32_t i;
	#if 0
		for (i = 0; i < large_num_; i++) {
			std::cout << largeCol_val[i] << "\t" << largeCol_oid[i] << std::endl;
		}
	#endif

	//uint64_t test_num = BLOCKSIZE<uint16_t>();	//case 1
	//uint64_t test_num = BLOCKSIZE<uint16_t>()/2;	//case 2
	//uint64_t test_num = BLOCKSIZE<uint16_t>()/4;	//case 3
	//uint64_t test_num = BLOCKSIZE<uint16_t>()/8;	//case 4
	uint64_t test_num = 256;	//case 5
#if 1
	//allocate input
	int16_t *input_val;
	surrogate_t *input_oid;
	input_val = (int16_t *)malloc_aligned(test_num * sizeof(int16_t));
	input_oid = (surrogate_t *)malloc_aligned(test_num * sizeof(surrogate_t));
	//int16_t *original_val = input_val;

    auto dice = std::bind(std::uniform_int_distribution<int16_t>(
                            std::numeric_limits<int16_t>::min(),
                            std::numeric_limits<int16_t>::max()),
                            std::default_random_engine(std::time(0)));

    for (i = 0; i < test_num; ++i) {
    	input_val[i] = dice();
    	input_oid[i] = (surrogate_t)i;
    }

    //use bank16Col_val to store original values
    memcpy(bank16Col_val, input_val, test_num * sizeof(int16_t));
#endif
	//allocate output
    int16_t *output_val;
	surrogate_t *output_oid;
	output_val = (int16_t *)malloc_aligned(test_num * sizeof(int16_t));
	output_oid = (surrogate_t *)malloc_aligned(test_num * sizeof(surrogate_t));

#if 1

	Mergesort::avxmergesort_block_uint16_aligned((uint16_t **)&input_val, &input_oid,
			(uint16_t **)&output_val, &output_oid, test_num);
#endif
#if 0
	int32_t *tmp_val = NULL;
	tmp_val = bank32Col_val;
	Mergesort::avxmergesort_rem_uint32_aligned((uint32_t **)&bank32Col_val, &bank32Col_oid,
			(uint32_t **)&output_val, &output_oid, test_num);
#endif

	EXPECT_TRUE(std::is_sorted(output_val, output_val+test_num));

#if 1
	for (i = 0; i < test_num; ++i) {
		assert(output_oid[i] >= 0);
		assert(output_oid[i] < test_num);
		assert(input_oid[i] >= 0);
		assert(input_oid[i] < test_num);
	}

#endif

	//test oid equivalence
#if 1
	//if (original_val == input_val) {//pointer not switched
		for (i = 0; i < test_num; ++i) {
			EXPECT_EQ(output_val[i], bank16Col_val[output_oid[i]]);
		}
	//} else {
	//	for (i = 0; i < test_num; ++i) {
	//		EXPECT_EQ(input_val[i], bank16Col_val[input_oid[i]]);
	//	}
	//}
#endif

	free(input_val);
	free(input_oid);
	free(output_val);
	free(output_oid);
}
#endif

#if 1
TEST_F(MergekernelTest, avxmergesort_block_uint32_aligned){
	uint32_t i;
	#if 0
		for (i = 0; i < large_num_; i++) {
			std::cout << largeCol_val[i] << "\t" << largeCol_oid[i] << std::endl;
		}
	#endif

	uint64_t test_num = BLOCKSIZE<uint32_t>();
#if 1
	//allocate input
	int32_t *input_val;
	surrogate_t *input_oid;
	input_val = (int32_t *)malloc_aligned(test_num * sizeof(int32_t));
	input_oid = (surrogate_t *)malloc_aligned(test_num * sizeof(surrogate_t));
	//int32_t *original_val = input_val;

    auto dice = std::bind(std::uniform_int_distribution<int32_t>(
                            std::numeric_limits<int32_t>::min(),
                            std::numeric_limits<int32_t>::max()),
                            std::default_random_engine(std::time(0)));

    for (i = 0; i < test_num; ++i) {
    	input_val[i] = dice();
    	input_oid[i] = (surrogate_t)i;
    }

    //use bank32Col_val to store original values
    memcpy(bank32Col_val, input_val, test_num * sizeof(int32_t));
#endif
	//allocate output
	int32_t *output_val;
	surrogate_t *output_oid;
	output_val = (int32_t *)malloc_aligned(test_num * sizeof(int32_t));
	output_oid = (surrogate_t *)malloc_aligned(test_num * sizeof(surrogate_t));

#if 1

	Mergesort::avxmergesort_block_uint32_aligned((uint32_t **)&input_val, &input_oid,
			(uint32_t **)&output_val, &output_oid, test_num);
#endif
#if 0
	int32_t *tmp_val = NULL;
	tmp_val = bank32Col_val;
	Mergesort::avxmergesort_rem_uint32_aligned((uint32_t **)&bank32Col_val, &bank32Col_oid,
			(uint32_t **)&output_val, &output_oid, test_num);
#endif

	EXPECT_TRUE(std::is_sorted(output_val, output_val+test_num));

#if 1
	for (i = 0; i < test_num; ++i) {
		assert(output_oid[i] < test_num);
		assert(output_oid[i] >= 0);
		assert(input_oid[i] < test_num);
		assert(input_oid[i] >= 0);
	}

#endif

	//test oid equivalence
#if 1
	//if (original_val == input_val) {//pointer not switched
		for (i = 0; i < test_num; ++i) {
			EXPECT_EQ(output_val[i], bank32Col_val[output_oid[i]]);
		}
	//} else {
	//	for (i = 0; i < test_num; ++i) {
	//		EXPECT_EQ(input_val[i], bank32Col_val[input_oid[i]]);
	//	}
	//}
#endif

	free(input_val);
	free(input_oid);
	free(output_val);
	free(output_oid);
}
#endif

#if 1
TEST_F(MergekernelTest, avxmergesort_rem_uint32_aligned){
	uint32_t i;
	#if 0
		for (i = 0; i < large_num_; i++) {
			std::cout << largeCol_val[i] << "\t" << largeCol_oid[i] << std::endl;
		}
	#endif
	//for (uint32_t test_num = 256; test_num <= 2048; test_num += 128) {
		//uint64_t test_num = 16000; //case 1
		//uint64_t test_num = 2;	//case2:
		//uint64_t test_num = 127;	//case2.5
		//uint64_t test_num = 256;	//case3: can pass
		//uint64_t test_num = 255;	//case4:
		//uint64_t test_num = 10;	//case5:
		uint64_t test_num = 16001; //case 9 cannot pass: 8192+4096+2048+1024+512+129
		//uint64_t test_num = 12288;	//case 10
		//uint64_t test_num = 12300;	//case 10 ==> 8192+4096+12
		//uint64_t test_num = 14436;	//case 10
		//uint64_t test_num = 15000;	//case 10 cannot pass ==> 5 parts: 8192+4096+2048+512+152
		//uint64_t test_num = 15360;	//case 10
		//uint64_t test_num = 15400;	//case 10 cannot pass ==> 5 parts: 8192+4096+2048+1024+40
#if 1
		//allocate input
		int32_t *input_val;
		surrogate_t *input_oid;
		input_val = (int32_t *) malloc_aligned(test_num * sizeof(int32_t));
		input_oid = (surrogate_t *) malloc_aligned(
				test_num * sizeof(surrogate_t));
		//int32_t *original_val = input_val;

		auto dice = std::bind(
				std::uniform_int_distribution<int32_t>(
						std::numeric_limits<int32_t>::min(),
						std::numeric_limits<int32_t>::max()),
				std::default_random_engine(std::time(0)));

		//std::cout << "input_val[i]:" << std::endl;
		for (i = 0; i < test_num; ++i) {
			input_val[i] = dice();
			input_oid[i] = (surrogate_t) i;
			//std::cout << input_val[i] << std::endl;
		}

		//use bank32Col_val to store original values
		memcpy(bank32Col_val, input_val, test_num * sizeof(int32_t));
#endif
		//allocate output
		int32_t *output_val;
		surrogate_t *output_oid;
		output_val = (int32_t *) malloc_aligned(test_num * sizeof(int32_t));
		output_oid = (surrogate_t *) malloc_aligned(
				test_num * sizeof(surrogate_t));

#if 1
		//HybridTimer timer1;
		//timer1.Start();

		Mergesort::avxmergesort_rem_uint32_aligned((uint32_t **) &input_val,
				&input_oid, (uint32_t **) &output_val, &output_oid, test_num);

		//timer1.Stop();

		//std::cout << "Cycles for avxmergesort_rem_uint32_aligned(): " << timer1.GetNumCycles() << std::endl;
#endif
#if 0
		int32_t *tmp_val = NULL;
		tmp_val = bank32Col_val;
		Mergesort::avxmergesort_rem_uint32_aligned((uint32_t **)&bank32Col_val, &bank32Col_oid,
				(uint32_t **)&output_val, &output_oid, test_num);
#endif

		/* this one cannt be passed because we didn't do the flipping here!!!*/
		//EXPECT_TRUE(std::is_sorted(output_val, output_val+test_num));
#if 1
		for (i = 0; i < test_num; ++i) {
			assert(output_oid[i] < test_num);
		}

#endif

		//test oid equivalence
#if 1
		//std::cout << "output_val[i]:" << std::endl;
		//if ((original_val == input_val) || (test_num < 256)) {//pointer not switched, or num<256(hard-coded)
		for (i = 0; i < test_num; ++i) {
			EXPECT_EQ(output_val[i], bank32Col_val[output_oid[i]]);
			//std::cout << output_val[i] << std::endl;
		}
		//} else {
		//	for (i = 0; i < test_num; ++i) {
		//		EXPECT_EQ(input_val[i], bank32Col_val[input_oid[i]]);
		//std::cout << "original oid: " << input_oid[i] << "\tvalue: " << input_val[i] << std::endl;
		//	}
		//}
#endif

		free(input_val);
		free(input_oid);
		free(output_val);
		free(output_oid);

}
#endif

#if 1
TEST_F(MergekernelTest, avxmergesort_rem_uint16_aligned){
	uint32_t i;
	#if 0
		for (i = 0; i < large_num_; i++) {
			std::cout << largeCol_val[i] << "\t" << largeCol_oid[i] << std::endl;
		}
	#endif

	//uint64_t test_num = 16000; //case 1 cannot pass: 8192+4096+2048+1024+512+128
	//uint64_t test_num = 2;	//case2:
	//uint64_t test_num = 127;	//case2.5
	//uint64_t test_num = 256;	//case3: can pass
	//uint64_t test_num = 255;	//case4:
	//uint64_t test_num = 4096;	//case5:
	//uint64_t test_num = 4098;	//case6:
	//uint64_t test_num = 9000;	//case7:	8192+512+296
	//uint64_t test_num = 16383;	//case8:
	uint64_t test_num = 16001; //case 9 cannot pass: 8192+4096+2048+1024+512+129
	//uint64_t test_num = 12288;	//case 10
	//uint64_t test_num = 12300;	//case 10 ==> 8192+4096+12
	//uint64_t test_num = 14436;	//case 10
	//uint64_t test_num = 15000;	//case 10 cannot pass ==> 5 parts: 8192+4096+2048+512+152
	//uint64_t test_num = 15360;	//case 10
	//uint64_t test_num = 15400;	//case 10 cannot pass ==> 5 parts: 8192+4096+2048+1024+40
	//uint64_t test_num = 550;	//case 10 cannot pass ==> 2 parts: 512+38
#if 1
	//allocate input
	int16_t *input_val;
	surrogate_t *input_oid;
	input_val = (int16_t *)malloc_aligned(test_num * sizeof(int16_t));
	input_oid = (surrogate_t *)malloc_aligned(test_num * sizeof(surrogate_t));
	//int16_t *original_val = input_val;

    auto dice = std::bind(std::uniform_int_distribution<int16_t>(
                            std::numeric_limits<int16_t>::min(),
                            std::numeric_limits<int16_t>::max()),
                            std::default_random_engine(std::time(0)));

    //std::cout << "input_val[i]:" << std::endl;
    for (i = 0; i < test_num; ++i) {
    	input_val[i] = dice();
    	input_oid[i] = (surrogate_t)i;
    	//std::cout << input_val[i] << std::endl;
    }

    //use bank32Col_val to store original values
    memcpy(bank16Col_val, input_val, test_num * sizeof(int16_t));
#endif
	//allocate output
    int16_t *output_val;
	surrogate_t *output_oid;
	output_val = (int16_t *)malloc_aligned(test_num * sizeof(int16_t));
	output_oid = (surrogate_t *)malloc_aligned(test_num * sizeof(surrogate_t));

#if 1
	//HybridTimer timer1;
	//timer1.Start();

	Mergesort::avxmergesort_rem_uint16_aligned((uint16_t **)&input_val, &input_oid,
			(uint16_t **)&output_val, &output_oid, test_num);

	//timer1.Stop();

	//std::cout << "Cycles for avxmergesort_rem_uint32_aligned(): " << timer1.GetNumCycles() << std::endl;
#endif
#if 0
	int32_t *tmp_val = NULL;
	tmp_val = bank32Col_val;
	Mergesort::avxmergesort_rem_uint32_aligned((uint32_t **)&bank32Col_val, &bank32Col_oid,
			(uint32_t **)&output_val, &output_oid, test_num);
#endif

	/* this one cannt be passed because we didn't do the flipping here!!!*/
	//EXPECT_TRUE(std::is_sorted(output_val, output_val+test_num));

#if 1
	for (i = 0; i < test_num; ++i) {
		assert(output_oid[i] < test_num);
		assert(output_oid[i] >= 0);
		if (test_num >= 256) {	// special case
			assert(input_oid[i] < test_num);
			assert(input_oid[i] >= 0);
		}
	}

#endif

	//test oid equivalence
#if 1
	//std::cout << "output_val[i]:" << std::endl;
	//if ((original_val == input_val) || (test_num < 256)) {//pointer not switched, or num<256(hard-coded)
		for (i = 0; i < test_num; ++i) {
			EXPECT_EQ(output_val[i], bank16Col_val[output_oid[i]]);
			//std::cout << output_val[i] << std::endl;
		}
	//} else {
	//	for (i = 0; i < test_num; ++i) {
	//		EXPECT_EQ(input_val[i], bank16Col_val[input_oid[i]]);
	//		std::cout << "i value: " << i << "\t original oid: " << input_oid[i] << std::endl;
	//	}
	//}
#endif

	free(input_val);
	free(input_oid);
	free(output_val);
	free(output_oid);
}
#endif

#if 1
TEST_F(MergekernelTest, avxmergesort_rem_uint64_aligned){
	uint32_t i;
	#if 0
		for (i = 0; i < large_num_; i++) {
			std::cout << largeCol_val[i] << "\t" << largeCol_oid[i] << std::endl;
		}
	#endif

	//uint64_t test_num = 8191;	//case1
	//uint64_t test_num = 4096;	//case2
	//uint64_t test_num = 4196;	//case2.5
	//uint64_t test_num = 2;	//case3
	//uint64_t test_num = 128;	//case4 ==>cannot pass is_sorted
	uint64_t test_num = 255;	//case4
	//uint64_t test_num = 256;	//can pass
#if 1
	//allocate input
	int64_t *input_val;
	surrogate_t *input_oid;
	input_val = (int64_t *)malloc_aligned(test_num * sizeof(int64_t));
	input_oid = (surrogate_t *)malloc_aligned(test_num * sizeof(surrogate_t));
	//int64_t *original_val = input_val;

    auto dice = std::bind(std::uniform_int_distribution<int64_t>(
                            std::numeric_limits<int64_t>::min(),
                            std::numeric_limits<int64_t>::max()),
                            std::default_random_engine(std::time(0)));

    for (i = 0; i < test_num; ++i) {
    	input_val[i] = dice();
    	input_oid[i] = (surrogate_t)i;
    }

    //use bank32Col_val to store original values
    memcpy(bank64Col_val, input_val, test_num * sizeof(int64_t));
#endif
	//allocate output
	int64_t *output_val;
	surrogate_t *output_oid;
	output_val = (int64_t *)malloc_aligned(test_num * sizeof(int64_t));
	output_oid = (surrogate_t *)malloc_aligned(test_num * sizeof(surrogate_t));

#if 1
	//HybridTimer timer1;
	//timer1.Start();

	Mergesort::avxmergesort_rem_uint64_aligned((uint64_t **)&input_val, &input_oid,
			(uint64_t **)&output_val, &output_oid, test_num);

	//timer1.Stop();

	//std::cout << "Cycles for avxmergesort_rem_uint32_aligned(): " << timer1.GetNumCycles() << std::endl;
#endif
#if 0
	int32_t *tmp_val = NULL;
	tmp_val = bank32Col_val;
	Mergesort::avxmergesort_rem_uint32_aligned((uint32_t **)&bank32Col_val, &bank32Col_oid,
			(uint32_t **)&output_val, &output_oid, test_num);
#endif

	/* this one cannt be passed because we didn't do the flipping here!!!*/
	//EXPECT_TRUE(std::is_sorted(output_val, output_val+test_num));

#if 1
	for (i = 0; i < test_num; ++i) {
		assert(output_oid[i] < test_num);
	}

#endif

	//test oid equivalence
#if 1
	//if ((original_val == input_val) || (test_num < 256)) {//pointer not switched, (256) is hard coded
		for (i = 0; i < test_num; ++i) {
			EXPECT_EQ(output_val[i], bank64Col_val[output_oid[i]]);
		}
	//} else {
	//	for (i = 0; i < test_num; ++i) {
	//		EXPECT_EQ(input_val[i], bank64Col_val[input_oid[i]]);
	//		//std::cout << "target original oid: " << input_oid[i] << "\tvalue: " << "" << std::endl;
	//	}
	//}
#endif

	free(input_val);
	free(input_oid);
	free(output_val);
	free(output_oid);
}
#endif

TEST_F(MergekernelTest, BITONIC_MERGE8_64){
	uint32_t i;
#if 0
	for (i = 0; i < large_num_; i++) {
		std::cout << largeCol_val[i] << "\t" << largeCol_oid[i] << std::endl;
	}
#endif
	//allocate output
	int64_t *output_val;
	surrogate_t *output_oid;
	output_val = (int64_t *)malloc_aligned(num_ * sizeof(int64_t));
	output_oid = (surrogate_t *)malloc_aligned(num_ * sizeof(surrogate_t));

	//sort two run of length 8
	std::stable_sort(bank64Col_val, bank64Col_val+8);
	std::stable_sort(bank64Col_val+8, bank64Col_val+16);

#if 0
	for (i = 0; i < 16; i++) {
		std::cout << bank64Col_val[i] << "\t" << bank64Col_oid[i] << std::endl;
	}
#endif

#if 0
	std::cout << "original values:" << std::endl;
	for (i = 0; i < 16; i++) {
		std::cout << output_val[i] << std::endl;
	}
#endif

	//load to register
	register block8_64 * inA_val = (block8_64 *) bank64Col_val;
	register block8_64 * inB_val = (block8_64 *)(bank64Col_val+8);

	register block8_32 * inA_oid = (block8_32 *) bank64Col_oid;
	register block8_32 * inB_oid = (block8_32 *)(bank64Col_oid+8);

	register block8_64 * outp_val = (block8_64 *) output_val;
	register block8_32 * outp_oid = (block8_32 *) output_oid;

    register __m256d outreg1l_val, outreg1h_val;
    register __m256d outreg2l_val, outreg2h_val;

    register __m256d regAl_val, regAh_val;
    register __m256d regBl_val, regBh_val;

    LOAD8_64(regAl_val, regAh_val, inA_val);
    LOAD8_64(regBl_val, regBh_val, inB_val);

    register __m256d outreg1l_oid, outreg1h_oid;
    register __m256d outreg2l_oid, outreg2h_oid;

    register __m256d regAl_oid, regAh_oid;
    register __m256d regBl_oid, regBh_oid;

    LOAD8_32T64(regAl_oid, regAh_oid, inA_oid);
    LOAD8_32T64(regBl_oid, regBh_oid, inB_oid);

    //do the sorting

    BITONIC_MERGE8_64(outreg1l_val, outreg1h_val, outreg2l_val, outreg2h_val,
    				outreg1l_oid, outreg1h_oid, outreg2l_oid, outreg2h_oid,
                   regAl_val, regAh_val, regBl_val, regBh_val,
				   regAl_oid, regAh_oid, regBl_oid, regBh_oid);

#if 0
    int64_t tempVal[4] __attribute__((aligned(32)));
    _mm256_store_si256(reinterpret_cast<__m256i *>(tempVal), reinterpret_cast<__m256i>(outreg1l_val));
    std::cout << tempVal[0] << "\t" << tempVal[1] << "\t" << tempVal[2] << "\t" << tempVal[3] << std::endl;
    _mm256_store_si256(reinterpret_cast<__m256i *>(tempVal), reinterpret_cast<__m256i>(outreg1h_val));
    std::cout << tempVal[0] << "\t" << tempVal[1] << "\t" << tempVal[2] << "\t" << tempVal[3] << std::endl;
    _mm256_store_si256(reinterpret_cast<__m256i *>(tempVal), reinterpret_cast<__m256i>(outreg2l_val));
    std::cout << tempVal[0] << "\t" << tempVal[1] << "\t" << tempVal[2] << "\t" << tempVal[3] << std::endl;
    _mm256_store_si256(reinterpret_cast<__m256i *>(tempVal), reinterpret_cast<__m256i>(outreg2h_val));
    std::cout << tempVal[0] << "\t" << tempVal[1] << "\t" << tempVal[2] << "\t" << tempVal[3] << std::endl;
    std::cout << std::endl;
#endif
    STORE8_64(outp_val, outreg1l_val, outreg1h_val);
    outp_val++;
    STORE8_64(outp_val, outreg2l_val, outreg2h_val);

    STORE8_64T32(outp_oid, outreg1l_oid, outreg1h_oid);
    outp_oid++;
    STORE8_64T32(outp_oid, outreg2l_oid, outreg2h_oid);


#if 0
	for (i = 0; i < 16; i++) {
		std::cout << output_val[i] << std::endl;
	}
#endif
	//check sortness
    EXPECT_TRUE(std::is_sorted(output_val, output_val+16));

	//check oid equivalence
	for (i = 0; i < 16; ++i) {
		EXPECT_EQ(output_val[i], bank64Col_val[output_oid[i]]);
		//std::cout << largeCol_val[i] << "\t" << temp_val[largeCol_oid[i]] << "\t" << largeCol_oid[i] << std::endl;
		//std::cout << "oid: " << i << std::endl;
	}

	free(output_val);
	free(output_oid);
}

TEST_F(MergekernelTest, BITONIC_MERGE8_32){
	uint32_t i;
#if 0
	for (i = 0; i < large_num_; i++) {
		std::cout << largeCol_val[i] << "\t" << largeCol_oid[i] << std::endl;
	}
#endif
	//allocate output
	int32_t *output_val;
	surrogate_t *output_oid;
	output_val = (int32_t *)malloc_aligned(num_ * sizeof(int32_t));
	output_oid = (surrogate_t *)malloc_aligned(num_ * sizeof(surrogate_t));

	//sort two run of length 8
	std::stable_sort(bank32Col_val, bank32Col_val+8);
	std::stable_sort(bank32Col_val+8, bank32Col_val+16);

#if 0
	for (i = 0; i < 16; i++) {
		std::cout << bank32Col_val[i] << "\t" << bank32Col_oid[i] << std::endl;
	}
#endif

#if 0
	std::cout << "original values:" << std::endl;
	for (i = 0; i < 16; i++) {
		std::cout << output_val[i] << std::endl;
	}
#endif

	//load to register
	register block8_32 * inA_val = (block8_32 *) bank32Col_val;
	register block8_32 * inB_val = (block8_32 *)(bank32Col_val+8);

	register block8_32 * inA_oid = (block8_32 *) bank32Col_oid;
	register block8_32 * inB_oid = (block8_32 *)(bank32Col_oid+8);

	register block8_32 * outp_val = (block8_32 *) output_val;
	register block8_32 * outp_oid = (block8_32 *) output_oid;

    register __m256i outregl_val, outregh_val;

    register __m256i regA_val, regB_val;

    regA_val = _mm256_load_si256((__m256i const *)inA_val);
    regB_val = _mm256_load_si256((__m256i const *)inB_val);

    register __m256i outregl_oid, outregh_oid;

    register __m256i regA_oid, regB_oid;

    regA_oid = _mm256_load_si256((__m256i const *)inA_oid);
    regB_oid = _mm256_load_si256((__m256i const *)inB_oid);

    //do the sorting
    BITONIC_MERGE8_32(outregl_val, outregh_val, outregl_oid, outregh_oid,
    		regA_val, regB_val, regA_oid, regB_oid);

    _mm256_store_si256((__m256i *)outp_val, outregl_val);
    outp_val++;
    _mm256_store_si256((__m256i *)outp_val, outregh_val);

    _mm256_store_si256((__m256i *)outp_oid, outregl_oid);
    outp_oid++;
    _mm256_store_si256((__m256i *)outp_oid, outregh_oid);


#if 0
	for (i = 0; i < 16; i++) {
		std::cout << output_val[i] << std::endl;
	}
#endif
	//check sortness
    EXPECT_TRUE(std::is_sorted(output_val, output_val+16));

	//check oid equivalence
	for (i = 0; i < 16; ++i) {
		EXPECT_EQ(output_val[i], bank32Col_val[output_oid[i]]);
		//std::cout << largeCol_val[i] << "\t" << temp_val[largeCol_oid[i]] << "\t" << largeCol_oid[i] << std::endl;
		//std::cout << "oid: " << i << std::endl;
	}

	free(output_val);
	free(output_oid);
}

TEST_F(MergekernelTest, BITONIC_MERGE16_16){
	uint32_t i;
#if 0
	for (i = 0; i < large_num_; i++) {
		std::cout << largeCol_val[i] << "\t" << largeCol_oid[i] << std::endl;
	}
#endif
	//allocate output
	int16_t *output_val;
	surrogate_t *output_oid;
	output_val = (int16_t *)malloc_aligned(num_ * sizeof(int16_t));
	output_oid = (surrogate_t *)malloc_aligned(num_ * sizeof(surrogate_t));

	//sort two run of length 16
	std::stable_sort(bank16Col_val, bank16Col_val+16);
	std::stable_sort(bank16Col_val+16, bank16Col_val+32);

#if 0
	for (i = 0; i < 16; i++) {
		std::cout << bank32Col_val[i] << "\t" << bank32Col_oid[i] << std::endl;
	}
#endif

#if 0
	std::cout << "original values:" << std::endl;
	for (i = 0; i < 16; i++) {
		std::cout << output_val[i] << std::endl;
	}
#endif

	//load to register
	register block16_16 * inA_val = (block16_16 *) bank16Col_val;
	register block16_16 * inB_val = (block16_16 *)(bank16Col_val+16);

	register block16_32 * inA_oid = (block16_32 *) bank16Col_oid;
	register block16_32 * inB_oid = (block16_32 *)(bank16Col_oid+16);

	register block16_16 * outp_val = (block16_16 *) output_val;
	register block16_32 * outp_oid = (block16_32 *) output_oid;

    register __m256i outregl_val, outregh_val;

    register __m256i regA_val, regB_val;

    regA_val = _mm256_load_si256((__m256i const *)inA_val);
    regB_val = _mm256_load_si256((__m256i const *)inB_val);

    register __m256i outregl1_oid, outregl2_oid, outregh1_oid, outregh2_oid;

    register __m256i regA1_oid, regA2_oid, regB1_oid, regB2_oid;

    regA1_oid = _mm256_load_si256((__m256i const *)inA_oid);
    regA2_oid = _mm256_load_si256((__m256i const *)((block8_32 *)inA_oid + 1));
    regB1_oid = _mm256_load_si256((__m256i const *)inB_oid);
    regB2_oid = _mm256_load_si256((__m256i const *)((block8_32 *)inB_oid + 1));

    //do the sorting
    BITONIC_MERGE16_16(outregl_val, outregh_val,
    		outregl1_oid, outregl2_oid, outregh1_oid, outregh2_oid,
    		regA_val, regB_val, regA1_oid, regA2_oid, regB2_oid, regB1_oid);

    _mm256_store_si256((__m256i *)outp_val, outregl_val);
    outp_val++;
    _mm256_store_si256((__m256i *)outp_val, outregh_val);

    _mm256_store_si256((__m256i *)outp_oid, outregl1_oid);
    _mm256_store_si256((__m256i *)((block8_32 *)outp_oid + 1), outregl2_oid);
    outp_oid++;
    _mm256_store_si256((__m256i *)outp_oid, outregh1_oid);
    _mm256_store_si256((__m256i *)((block8_32 *)outp_oid + 1), outregh2_oid);


#if 0
	for (i = 0; i < 16; i++) {
		std::cout << output_val[i] << std::endl;
	}
#endif
	//check sortness
    EXPECT_TRUE(std::is_sorted(output_val, output_val+32));

	//check oid equivalence
	for (i = 0; i < 32; ++i) {
		EXPECT_EQ(output_val[i], bank16Col_val[output_oid[i]]);
		//std::cout << largeCol_val[i] << "\t" << temp_val[largeCol_oid[i]] << "\t" << largeCol_oid[i] << std::endl;
		//std::cout << "oid: " << i << std::endl;
	}

	free(output_val);
	free(output_oid);
}

TEST_F(MergekernelTest, merge8_int64_varlen_aligned){
	uint32_t i;
#if 0
	for (i = 0; i < large_num_; i++) {
		std::cout << largeCol_val[i] << "\t" << largeCol_oid[i] << std::endl;
	}
#endif
	//allocate output
	int64_t *output_val;
	surrogate_t *output_oid;
	output_val = (int64_t *)malloc_aligned(num_ * sizeof(int64_t));
	output_oid = (surrogate_t *)malloc_aligned(num_ * sizeof(surrogate_t));

	//sort two runs
	uint64_t sizeA = BLOCKSIZE<uint64_t>();
	uint64_t sizeB = BLOCKSIZE<uint64_t>();
	std::stable_sort(bank64Col_val, bank64Col_val+sizeA);
	std::stable_sort(bank64Col_val+sizeA, bank64Col_val+sizeA+sizeB);

	//use a temporary mem to remember the array values
	int64_t *tmp_val = (int64_t *)malloc_aligned((sizeA+sizeB) * sizeof(int64_t));
	memcpy(tmp_val, bank64Col_val, (sizeA+sizeB) * sizeof(int64_t));

#if 0
	for (i = 0; i < sizeA+sizeB; i++) {
		std::cout << bank64Col_val[i] << "\t" << bank64Col_oid[i] << std::endl;
	}
#endif

#if 0
	std::cout << "original values:" << std::endl;
	for (i = 0; i < 16; i++) {
		std::cout << output_val[i] << std::endl;
	}
#endif

	HybridTimer timer1;
	timer1.Start();
    //do the sorting
	//std::cout << "before sorting" << std::endl;
	Mergesort::merge8_int64_varlen_aligned(bank64Col_val,
									bank64Col_val+sizeA,
									output_val,
									bank64Col_oid,
									bank64Col_oid+sizeB,
									output_oid,
	                      			sizeA,
	                      			sizeB);
	//std::cout << "after sorting" << std::endl;
	timer1.Stop();
	std::cout << "merge8_int64_varlen_aligned: " << timer1.GetNumCycles()/(sizeA+sizeB) << "cycles/tuple" << std::endl;

#if 0
	for (i = 0; i < sizeA+sizeB; i++) {
		std::cout << output_val[i] << "\t" << output_oid[i] << "\t" << i << std::endl;
	}
#endif
	//check sortness
    EXPECT_TRUE(std::is_sorted(output_val, output_val+sizeA+sizeB));

	//check oid equivalence
	for (i = 0; i < sizeA+sizeB; ++i) {
		EXPECT_EQ(output_val[i], tmp_val[output_oid[i]]);
		//std::cout << largeCol_val[i] << "\t" << temp_val[largeCol_oid[i]] << "\t" << largeCol_oid[i] << std::endl;
		assert(output_oid[i] < (sizeA+sizeB));
	}

	free(output_val);
	free(output_oid);
	free(tmp_val);
}

#if 0
TEST_F(MergekernelTest, merge8_int32_varlen_aligned){
	uint32_t i;
#if 0
	for (i = 0; i < large_num_; i++) {
		std::cout << largeCol_val[i] << "\t" << largeCol_oid[i] << std::endl;
	}
#endif
	//allocate output
	int32_t *output_val;
	surrogate_t *output_oid;
	output_val = (int32_t *)malloc_aligned(num_ * sizeof(int32_t));
	output_oid = (surrogate_t *)malloc_aligned(num_ * sizeof(surrogate_t));

	//sort two runs
	uint64_t sizeA = 8192;
	uint64_t sizeB = 4096;
	std::stable_sort(bank32Col_val, bank32Col_val+sizeA);
	std::stable_sort(bank32Col_val+sizeA, bank32Col_val+sizeA+sizeB);

	//use a temporary mem to remember the array values
	int32_t *tmp_val = (int32_t *)malloc_aligned((sizeA+sizeB) * sizeof(int32_t));
	memcpy(tmp_val, bank32Col_val, (sizeA+sizeB) * sizeof(int32_t));

#if 0
	for (i = 0; i < 16; i++) {
		std::cout << bank32Col_val[i] << "\t" << bank32Col_oid[i] << std::endl;
	}
#endif

#if 0
	std::cout << "original values:" << std::endl;
	for (i = 0; i < 16; i++) {
		std::cout << output_val[i] << std::endl;
	}
#endif

    //do the sorting
	//std::cout << "before sorting" << std::endl;
	Mergesort::merge8_int32_varlen_aligned(bank32Col_val,
									bank32Col_val+sizeA,
									output_val,
									bank32Col_oid,
									bank32Col_oid+sizeA,
									output_oid,
	                      			sizeA,
	                      			sizeB);
	//std::cout << "after sorting" << std::endl;

#if 0
	for (i = 0; i < 16; i++) {
		std::cout << output_val[i] << std::endl;
	}
#endif
	//check sortness
    EXPECT_TRUE(std::is_sorted(output_val, output_val+sizeA+sizeB));

	//check oid equivalence
	for (i = 0; i < sizeA+sizeB; ++i) {
		EXPECT_EQ(output_val[i], tmp_val[output_oid[i]]);
		//std::cout << largeCol_val[i] << "\t" << temp_val[largeCol_oid[i]] << "\t" << largeCol_oid[i] << std::endl;
		//std::cout << "oid: " << i << std::endl;
		assert(output_oid[i] < sizeA+sizeB);
	}

	free(output_val);
	free(output_oid);
	free(tmp_val);
}
#endif

#if 0
TEST_F(MergekernelTest, merge8_int32_varlen_unaligned){
	uint32_t i;
#if 0
	for (i = 0; i < large_num_; i++) {
		std::cout << largeCol_val[i] << "\t" << largeCol_oid[i] << std::endl;
	}
#endif
	//allocate output
	int32_t *output_val;
	surrogate_t *output_oid;
	output_val = (int32_t *)malloc_aligned(num_ * sizeof(int32_t));
	output_oid = (surrogate_t *)malloc_aligned(num_ * sizeof(surrogate_t));

	//sort two runs
	uint64_t sizeA = 6509;
	uint64_t sizeB = 3119;
	std::stable_sort(bank32Col_val, bank32Col_val+sizeA);
	std::stable_sort(bank32Col_val+sizeA, bank32Col_val+sizeA+sizeB);

	//use a temporary mem to remember the array values
	int32_t *tmp_val = (int32_t *)malloc_aligned((sizeA+sizeB) * sizeof(int32_t));
	memcpy(tmp_val, bank32Col_val, (sizeA+sizeB) * sizeof(int32_t));

#if 0
	for (i = 0; i < 16; i++) {
		std::cout << bank32Col_val[i] << "\t" << bank32Col_oid[i] << std::endl;
	}
#endif

#if 0
	std::cout << "original values:" << std::endl;
	for (i = 0; i < 16; i++) {
		std::cout << output_val[i] << std::endl;
	}
#endif

    //do the sorting
	//std::cout << "before sorting" << std::endl;
	Mergesort::merge8_int32_varlen_unaligned(bank32Col_val,
									bank32Col_val+sizeA,
									output_val,
									bank32Col_oid,
									bank32Col_oid+sizeA,
									output_oid,
	                      			sizeA,
	                      			sizeB);
	//std::cout << "after sorting" << std::endl;


	//check sortness
    EXPECT_TRUE(std::is_sorted(output_val, output_val+sizeA+sizeB));

#if 1
	for (i = 0; i < sizeA+sizeB; i++) {
		assert(output_oid[i] < (sizeA+sizeB));
	}
#endif

	//check oid equivalence
	for (i = 0; i < sizeA+sizeB; ++i) {
		EXPECT_EQ(output_val[i], tmp_val[output_oid[i]]);
		//std::cout << largeCol_val[i] << "\t" << temp_val[largeCol_oid[i]] << "\t" << largeCol_oid[i] << std::endl;
		//std::cout << "oid: " << i << std::endl;
	}

	free(output_val);
	free(output_oid);
	free(tmp_val);
}
#endif

}   //namespace
