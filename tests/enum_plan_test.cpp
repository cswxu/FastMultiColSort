#include	<cstdlib>
#include    "gtest/gtest.h"
#include	"src/types.h"
#include	<ctime>
#include	<iostream>
#include	<cstdint>
#include 	<string>
#include	<vector>
#include	"src/common.h"
#include 	<algorithm>
#include	"src/column.h"

namespace multiAttrSort{

class EnumPlanTest: public ::testing::Test{
public:
    virtual void SetUp(){
     }

    virtual void TearDown(){
    }

protected:
    int count;
    int partitions[64];

    /** Given a round number, say 3. Return ALL massaging plans with that round number
     *
     * @param n  			remained value
     * @param curRound		filling the integer for which round
     * @param round_num   	how many rounds do we need to fill up
     * @param maxValue		the largest eligible value allowed for one round (i.e., 64 in our scenario)
     * @param split_enum	plans to be stored
     */
	void split_exhaustive(int n, int curRound, const int round_num,
			const int maxValue, std::vector<std::vector<uint32_t> > &split_enum);

	/**
	 * Given a round number, say 3. Return ALL possible banksize categories
	 *
	 * for the value in cate_enum, 0 stands for 16-banksize, 1 stands for 32-banksize, 2 stands for 64-banksize
	 */
	void enum_banksize_category(int curRound, const int round_num, std::vector<std::vector<uint32_t> > &cate_enum);
};

void EnumPlanTest::enum_banksize_category(int curRound, const int round_num, std::vector<std::vector<uint32_t> > &cate_enum) {

	int i;
	if (curRound == round_num) {	//surpass the last round
		std::vector<uint32_t> banksize_cate;
		for (i = 0; i < round_num; ++i) {
			banksize_cate.push_back(partitions[i]);
		}
		cate_enum.push_back(banksize_cate);
	} else {
		for (i = 0; i < 3; ++i) {
			partitions[curRound] = i;

			enum_banksize_category(curRound+1, round_num, cate_enum);
		}
	}
}

TEST_F(EnumPlanTest, enum_banksize_category){
	uint32_t test_time = 10;
	uint32_t round_num;
	std::vector<std::vector<uint32_t> > cate_enum;
	while (--test_time > 0) {
		std::cin >> round_num;
		std::vector<std::vector<uint32_t> >().swap(cate_enum);
		assert(cate_enum.empty());

		enum_banksize_category(0, round_num, cate_enum);

		//for display
		for (uint32_t i = 0; i < cate_enum.size(); ++i) {
			std::vector<uint32_t> &eachEnum = cate_enum[i];

			for (uint32_t j = 0; j < eachEnum.size(); ++j) {
				std::cout << eachEnum[j] << "\t";
			}
			std::cout << std::endl;

		}
		std::cout << cate_enum.size() << std::endl;
	}
}

void EnumPlanTest::split_exhaustive(int n, int curRound, const int round_num,
		const int maxValue, std::vector<std::vector<uint32_t> > &split_enum) {

	int i;
	if (curRound == (round_num-1)) {	//filling the last round
		assert(n > 0);
		if (n > maxValue){
			return;			//not eligible
		}
		partitions[curRound] = n;	//assign for the last round

		std::vector < uint32_t > aSplit;
		for (i = 0; i <= curRound; i++) {
			aSplit.push_back(partitions[i]);
		}
		split_enum.push_back(aSplit);

	} else {
		for ((n > maxValue) ? (i = maxValue) : (i = n); i > 0; i--) {
			partitions[curRound] = i;	//写入数组
			if ((n-i) > 0) {//still have values to be filled
				split_exhaustive(n - i, curRound + 1, round_num, maxValue, split_enum);
			}
		}
	}
}

#if 0
TEST_F(EnumPlanTest, split_exhaustive){
	uint32_t test_time = 10;
	uint32_t round_num;
	uint32_t col_sum;
	std::vector<std::vector<uint32_t> > split_enum;
	while (--test_time > 0) {
		std::cin >> round_num >> col_sum;
		std::vector<std::vector<uint32_t> >().swap(split_enum);
		assert(split_enum.empty());

		split_exhaustive(col_sum, 0, round_num, 64, split_enum);

		//for display
		for (uint32_t i = 0; i < split_enum.size(); ++i) {
			std::vector<uint32_t> &eachSplit = split_enum[i];

			for (uint32_t j = 0; j < eachSplit.size(); ++j) {
				std::cout << eachSplit[j] << "\t";
			}
			std::cout << std::endl;
		}
	}
}
#endif

}   //namespace
