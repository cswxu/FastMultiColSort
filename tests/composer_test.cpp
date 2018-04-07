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

class ComposerTest: public ::testing::Test{
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
TEST_F(ComposerTest, subsort_mergesort_uint16){

	uint64_t i;
	#if 0
		for (i = 0; i < large_num_; i++) {
			std::cout << largeCol_val[i] << "\t" << largeCol_oid[i] << std::endl;
		}
	#endif

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

	/*test case 1: the whole first half, use whole 16 bits*/
	uint64_t test_num_case1 = 0.5 * num_;

#if 0
	//make MSB=0
	uint64_t mask = (1ULL << 63) - 1;
	for (i = 0; i < test_num_case1; ++i) {
		bank64Col_val[i] &= mask;
		origin_val[i] &= mask;
	}
#endif

	composer->subsort_mergesort<uint16_t>(&subcolumn, 0, test_num_case1, true, 16);

	EXPECT_TRUE(std::is_sorted(subcolumn.values, subcolumn.values+test_num_case1));

#if 1
	for (i = 0; i < test_num_case1; ++i) {
		assert(subcolumn.oids[i] < test_num_case1);
	}
#endif

#if 0
	int64_t * sortedVal = subcolumn.values;
	surrogate_t * sortedOid = subcolumn.oids;
	for (i = 0; i < test_num_case1; ++i) {
		EXPECT_EQ(sortedVal[i], origin_val[sortedOid[i]]);
		//std::cout << "sortedVal[i]: " << sortedVal[i] <<
		//		"\tsortedOid[i]: " << sortedOid[i] <<
		//		"\torigin_val[sortedOid[i]]: " << origin_val[sortedOid[i]] << std::endl;
	}
#endif


#if 1
	/*test case 2: the next unaligned quarter, with MSB=0*/
	uint64_t startpoint = 33 + test_num_case1;
	uint64_t test_num_case2 = 0.25 * num_;

	//make MSB=0
	uint64_t mask = (1ULL << 15) - 1;
	for (i = startpoint; i < startpoint + test_num_case2; ++i) {
		bank16Col_val[i] &= mask;
		origin_val[i] &= mask;
	}

	//std::cout << "start subsort [" << startpoint << ", " << startpoint + test_num_case2 << "]" << std::endl;
	composer->subsort_mergesort<uint16_t>(&subcolumn, startpoint, test_num_case2, true, 15);
	//std::cout << "end subsort..." << std::endl;

	EXPECT_TRUE(std::is_sorted(subcolumn.values+startpoint, subcolumn.values+startpoint+test_num_case2));

#if 1
	for (i = startpoint; i < startpoint + test_num_case2; ++i) {
		assert(subcolumn.oids[i] < (startpoint + test_num_case2));
	}
#endif

	//test oid equivalence
#if 1
	uint16_t * sortedVal = subcolumn.values;
	surrogate_t * sortedOid = subcolumn.oids;
	for (i = startpoint; i < startpoint + test_num_case2; ++i) {
		EXPECT_EQ(sortedVal[i], origin_val[sortedOid[i]]);
		//std::cout << "sortedVal[i]: " << sortedVal[i] <<
		//		"\tsortedOid[i]: " << sortedOid[i] <<
		//		"\torigin_val[sortedOid[i]]: " << origin_val[sortedOid[i]] << std::endl;
	}
#endif
#endif

	/*test case 3: small size = 2*/
	uint64_t startpoint3 = startpoint + test_num_case2;
	uint64_t test_num_case3 = 2;

	//make MSB=0
#if 1
	for (i = startpoint3; i < startpoint3 + test_num_case3; ++i) {
		bank16Col_val[i] &= mask;
		origin_val[i] &= mask;
	}
#endif

	//std::cout << "start subsort [" << startpoint << ", " << startpoint + test_num_case2 << "]" << std::endl;
	composer->subsort_mergesort<uint16_t>(&subcolumn, startpoint3, test_num_case3, true, 15);
	//std::cout << "end subsort..." << std::endl;

	EXPECT_TRUE(std::is_sorted(subcolumn.values+startpoint3, subcolumn.values+startpoint3+test_num_case3));

	for (i = startpoint3; i < startpoint3 + test_num_case3; ++i) {
		EXPECT_EQ(sortedVal[i], origin_val[sortedOid[i]]);
	}



	free(origin_val);
	free(composer);
}
#endif

#if 0
TEST_F(ComposerTest, subsort_mergesort_uint32_profile){

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
TEST_F(ComposerTest, subsort_mergesort_uint64_profile){

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
TEST_F(ComposerTest, subsort_mergesort_uint16_profile){

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

#if 0
TEST_F(ComposerTest, subsort_mergesort_uint32){

	uint64_t i;
	#if 0
		for (i = 0; i < large_num_; i++) {
			std::cout << largeCol_val[i] << "\t" << largeCol_oid[i] << std::endl;
		}
	#endif

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

	/*test case 1: the whole first half, use whole 32 bits*/
	uint64_t test_num_case1 = 0.5 * num_;

#if 0
	//make MSB=0
	uint64_t mask = (1ULL << 63) - 1;
	for (i = 0; i < test_num_case1; ++i) {
		bank64Col_val[i] &= mask;
		origin_val[i] &= mask;
	}
#endif

	composer->subsort_mergesort<uint32_t>(&subcolumn, 0, test_num_case1, true, 32);

	EXPECT_TRUE(std::is_sorted(subcolumn.values, subcolumn.values+test_num_case1));

#if 1
	for (i = 0; i < test_num_case1; ++i) {
		assert(subcolumn.oids[i] < test_num_case1);
	}
#endif

#if 0
	int64_t * sortedVal = subcolumn.values;
	surrogate_t * sortedOid = subcolumn.oids;
	for (i = 0; i < test_num_case1; ++i) {
		EXPECT_EQ(sortedVal[i], origin_val[sortedOid[i]]);
		//std::cout << "sortedVal[i]: " << sortedVal[i] <<
		//		"\tsortedOid[i]: " << sortedOid[i] <<
		//		"\torigin_val[sortedOid[i]]: " << origin_val[sortedOid[i]] << std::endl;
	}
#endif


#if 1
	/*test case 2: the next unaligned quarter, with MSB=0*/
	uint64_t startpoint = 33 + test_num_case1;
	uint64_t test_num_case2 = 0.25 * num_;

	//make MSB=0
	uint64_t mask = (1ULL << 31) - 1;
	for (i = startpoint; i < startpoint + test_num_case2; ++i) {
		bank32Col_val[i] &= mask;
		origin_val[i] &= mask;
	}

	//std::cout << "start subsort [" << startpoint << ", " << startpoint + test_num_case2 << "]" << std::endl;
	composer->subsort_mergesort<uint32_t>(&subcolumn, startpoint, test_num_case2, true, 31);
	//std::cout << "end subsort..." << std::endl;

	EXPECT_TRUE(std::is_sorted(subcolumn.values+startpoint, subcolumn.values+startpoint+test_num_case2));

#if 1
	for (i = startpoint; i < startpoint + test_num_case2; ++i) {
		assert(subcolumn.oids[i] < (startpoint + test_num_case2));
	}
#endif

	//test oid equivalence
#if 1
	uint32_t * sortedVal = subcolumn.values;
	surrogate_t * sortedOid = subcolumn.oids;
	for (i = startpoint; i < startpoint + test_num_case2; ++i) {
		EXPECT_EQ(sortedVal[i], origin_val[sortedOid[i]]);
		//std::cout << "sortedVal[i]: " << sortedVal[i] <<
		//		"\tsortedOid[i]: " << sortedOid[i] <<
		//		"\torigin_val[sortedOid[i]]: " << origin_val[sortedOid[i]] << std::endl;
	}
#endif
#endif

	/*test case 3: small size = 2*/
	uint64_t startpoint3 = startpoint + test_num_case2;
	uint64_t test_num_case3 = 2;

	//make MSB=0
#if 1
	for (i = startpoint3; i < startpoint3 + test_num_case3; ++i) {
		bank32Col_val[i] &= mask;
		origin_val[i] &= mask;
	}
#endif

	//std::cout << "start subsort [" << startpoint << ", " << startpoint + test_num_case2 << "]" << std::endl;
	composer->subsort_mergesort<uint32_t>(&subcolumn, startpoint3, test_num_case3, true, 31);
	//std::cout << "end subsort..." << std::endl;

	EXPECT_TRUE(std::is_sorted(subcolumn.values+startpoint3, subcolumn.values+startpoint3+test_num_case3));

	for (i = startpoint3; i < startpoint3 + test_num_case3; ++i) {
		EXPECT_EQ(sortedVal[i], origin_val[sortedOid[i]]);
	}



	free(origin_val);
	free(composer);
}
#endif

#if 0
TEST_F(ComposerTest, subsort_mergesort_uint64){

	uint64_t i;
	#if 0
		for (i = 0; i < large_num_; i++) {
			std::cout << largeCol_val[i] << "\t" << largeCol_oid[i] << std::endl;
		}
	#endif

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

	/*test case 1: the whole first half, use whole 64 bits*/
	uint64_t test_num_case1 = 0.5 * num_;

#if 0
	//make MSB=0
	uint64_t mask = (1ULL << 63) - 1;
	for (i = 0; i < test_num_case1; ++i) {
		bank64Col_val[i] &= mask;
		origin_val[i] &= mask;
	}
#endif

	composer->subsort_mergesort<uint64_t>(&subcolumn, 0, test_num_case1, true, 64);

	EXPECT_TRUE(std::is_sorted(subcolumn.values, subcolumn.values+test_num_case1));

#if 0
	//output the first 500 and last 500 elements
	std::cout << "first 500 elements:" << std::endl;
	for (i = 0; i < 500; i++) {
		std::cout << subcolumn.values[i] << std::endl;
	}
	std::cout << "last 500 elements:" << std::endl;
	for (i = test_num_case1 - 500; i < test_num_case1; i++) {
		std::cout << subcolumn.values[i] << std::endl;
	}
#endif
#if 1
	for (i = 0; i < test_num_case1; ++i) {
		assert(subcolumn.oids[i] < test_num_case1);
	}
#endif

#if 0
	int64_t * sortedVal = subcolumn.values;
	surrogate_t * sortedOid = subcolumn.oids;
	for (i = 0; i < test_num_case1; ++i) {
		EXPECT_EQ(sortedVal[i], origin_val[sortedOid[i]]);
		//std::cout << "sortedVal[i]: " << sortedVal[i] <<
		//		"\tsortedOid[i]: " << sortedOid[i] <<
		//		"\torigin_val[sortedOid[i]]: " << origin_val[sortedOid[i]] << std::endl;
	}
#endif


#if 1
	/*test case 2: the next unaligned quarter, with MSB=0*/
	uint64_t startpoint = 33 + test_num_case1;
	uint64_t test_num_case2 = 0.25 * num_;

	//make MSB=0
	uint64_t mask = (1ULL << 63) - 1;
	for (i = startpoint; i < startpoint + test_num_case2; ++i) {
		bank64Col_val[i] &= mask;
		origin_val[i] &= mask;
	}

	//std::cout << "start subsort [" << startpoint << ", " << startpoint + test_num_case2 << "]" << std::endl;
	composer->subsort_mergesort<uint64_t>(&subcolumn, startpoint, test_num_case2, true, 63);
	//std::cout << "end subsort..." << std::endl;

	EXPECT_TRUE(std::is_sorted(subcolumn.values+startpoint, subcolumn.values+startpoint+test_num_case2));

#if 1
	for (i = startpoint; i < startpoint + test_num_case2; ++i) {
		assert(subcolumn.oids[i] < (startpoint + test_num_case2));
	}
#endif

	//test oid equivalence
#if 1
	uint64_t * sortedVal = subcolumn.values;
	surrogate_t * sortedOid = subcolumn.oids;
	for (i = startpoint; i < startpoint + 10; ++i) {
		EXPECT_EQ(sortedVal[i], origin_val[sortedOid[i]]);
		//std::cout << "sortedVal[i]: " << sortedVal[i] <<
		//		"\tsortedOid[i]: " << sortedOid[i] <<
		//		"\torigin_val[sortedOid[i]]: " << origin_val[sortedOid[i]] << std::endl;
	}
#endif
#endif

	/*test case 3: small size = 2*/
	uint64_t startpoint3 = startpoint + test_num_case2;
	uint64_t test_num_case3 = 2;

	//make MSB=0
#if 1
	for (i = startpoint3; i < startpoint3 + test_num_case3; ++i) {
		bank64Col_val[i] &= mask;
		origin_val[i] &= mask;
	}
#endif

	//std::cout << "start subsort [" << startpoint << ", " << startpoint + test_num_case2 << "]" << std::endl;
	composer->subsort_mergesort<uint64_t>(&subcolumn, startpoint3, test_num_case3, true, 63);
	//std::cout << "end subsort..." << std::endl;

	EXPECT_TRUE(std::is_sorted(subcolumn.values+startpoint3, subcolumn.values+startpoint3+test_num_case3));

	for (i = startpoint3; i < startpoint3 + test_num_case3; ++i) {
		EXPECT_EQ(sortedVal[i], origin_val[sortedOid[i]]);
	}



	free(origin_val);
	free(composer);
}
#endif

}   //namespace
