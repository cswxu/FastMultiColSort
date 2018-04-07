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

class OptimizerTest: public ::testing::Test{
public:
    virtual void SetUp(){
    	std::srand(std::time(0));
    	count = 0;
    }

    virtual void TearDown(){
    }

protected:
    int count;
    int partitions[64];

    /** for the integer division algorithm **/
    void split(int n, int curLevel, const int maxLevel,
    		const int maxValue, const int threshold_loc, const int sum);
    void display(int curLevel);

    /** for the all combinations algorithm **/
    /**
     * k-th number to swap, totally (m+1) numbers
     */
    void allRange(uint32_t * arr, uint32_t k, uint32_t m, std::vector<uint32_t *> &comb);
    void swap(uint32_t * a, uint32_t * b);

    void stitch(uint32_t * oldCol, size_t oldSize, uint32_t * newCol, size_t newSize);

    void assemble(Column** column_values_, Column** column_values_new,
    		uint32_t new_col_num, uint64_t num_rows_, std::vector < construct_entry_t > *constructRecord);

    void getHeuristics(std::vector < std::vector<uint32_t> > & plans, uint32_t colwidth_sum);
};

void OptimizerTest::getHeuristics(std::vector < std::vector<uint32_t> > & plans, uint32_t colwidth_sum) {
	std::vector<uint32_t> eachPlan;

	//64 as a boundary
	uint32_t temp_sum = colwidth_sum;
	uint32_t i;
	for (i = 0; i < colwidth_sum/64; ++i) {
		eachPlan.push_back(64);
		temp_sum -= 64;
	}
	if (temp_sum > 0) {
		eachPlan.push_back(temp_sum);
	}
	if (colwidth_sum > 64) {
		plans.push_back(eachPlan);
	}

	//32 as a boundary,
	temp_sum = colwidth_sum;
	eachPlan.clear();
	for (i = 0; i < colwidth_sum/32; ++i) {
		eachPlan.push_back(32);
		temp_sum -= 32;
	}
	if (temp_sum > 0) {
		eachPlan.push_back(temp_sum);
	}
	if (colwidth_sum > 32) {
		plans.push_back(eachPlan);
	}

	//16 as a boundary not a choice

	//stitching together
	if (colwidth_sum < 64) {
		plans.push_back(std::vector<uint32_t>(1, colwidth_sum));
	}
}

TEST_F(OptimizerTest, getHeuristicsTest){
	std::vector < std::vector<uint32_t> > plans;
	uint32_t i, j;

	getHeuristics(plans, 13);


	for (i = 0; i < plans.size(); ++i) {
		for (j = 0; j < plans[i].size(); ++j) {
			std::cout << plans[i][j] << "\t";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl << std::endl;

	plans.clear();
	getHeuristics(plans, 25);
	for (i = 0; i < plans.size(); ++i) {
		for (j = 0; j < plans[i].size(); ++j) {
			std::cout << plans[i][j] << "\t";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl << std::endl;


	plans.clear();
	getHeuristics(plans, 32);
	for (i = 0; i < plans.size(); ++i) {
		for (j = 0; j < plans[i].size(); ++j) {
			std::cout << plans[i][j] << "\t";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;


	plans.clear();
	getHeuristics(plans,39);
	for (i = 0; i < plans.size(); ++i) {
		for (j = 0; j < plans[i].size(); ++j) {
			std::cout << plans[i][j] << "\t";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

	plans.clear();
	getHeuristics(plans,64);
	for (i = 0; i < plans.size(); ++i) {
		for (j = 0; j < plans[i].size(); ++j) {
			std::cout << plans[i][j] << "\t";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl << std::endl;


	plans.clear();
	getHeuristics(plans,70);
	for (i = 0; i < plans.size(); ++i) {
		for (j = 0; j < plans[i].size(); ++j) {
			std::cout << plans[i][j] << "\t";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl << std::endl;


	plans.clear();
	getHeuristics(plans,64);
	for (i = 0; i < plans.size(); ++i) {
		for (j = 0; j < plans[i].size(); ++j) {
			std::cout << plans[i][j] << "\t";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl << std::endl;
}

TEST_F(OptimizerTest, assembleTest1){
	uint64_t num_rows_ = 16*1048576;
	auto dice = std::bind(std::uniform_int_distribution<int32_t>(
									std::numeric_limits<int32_t>::min(),
									std::numeric_limits<int32_t>::max()),
									std::default_random_engine(std::time(0)));


	//construct two columns(17/32+17/32)
	uint32_t colWidthA = 17;
	uint32_t colWidthB = 17;
	Column* column_values_[2];
	uint32_t *values = (uint32_t *)malloc_aligned(num_rows_ * sizeof(uint32_t));
	uint32_t mask = (1 << colWidthA) - 1;
	for (uint32_t idx = 0; idx < num_rows_; ++idx) {
		values[idx] = dice() & mask;

	}
	SpecialiedColumn<uint32_t> *col = new SpecialiedColumn<uint32_t>();
	col->SetColumn(values);
	col->SetWidth(sizeof(uint32_t));
	column_values_[0] = col;

	values = (uint32_t *)malloc_aligned(num_rows_ * sizeof(uint32_t));
	mask = (1 << colWidthB) - 1;
	for (uint32_t idx = 0; idx < num_rows_; ++idx) {
		values[idx] = dice() & mask;

	}
	col = new SpecialiedColumn<uint32_t>();
	col->SetColumn(values);
	col->SetWidth(sizeof(uint32_t));
	column_values_[1] = col;
	{
		/** test case 1, stitch together to 34/64**/
		Column* column_values_new_case1[1];

		uint64_t *values_case1 = (uint64_t *)malloc_aligned(num_rows_ * sizeof(uint64_t));
		SpecialiedColumn<uint64_t> *col_case1 = new SpecialiedColumn<uint64_t>();
		col_case1->SetColumn(values_case1);
		col_case1->SetWidth(sizeof(uint64_t));
		column_values_new_case1[0] = col_case1;

		std::vector < construct_entry_t > constructRecord_case1[1];
		construct_entry_t entry1;
		entry1.finalIdx = 0;
		entry1.start = 1;
		entry1.end = 17;
		constructRecord_case1[0].push_back(entry1);

		construct_entry_t entry2;
		entry2.finalIdx = 1;
		entry2.start = 1;
		entry2.end = 17;
		constructRecord_case1[0].push_back(entry2);

		HybridTimer timer;
		timer.Start();

		assemble(column_values_, column_values_new_case1, 1, num_rows_, constructRecord_case1);

		timer.Stop();
		std::cout << timer.GetNumCycles()/num_rows_ << "cycles/tuple" << std::endl;

		free((uint64_t *)column_values_new_case1[0]->GetColumn());
		free(column_values_new_case1[0]);
	}


	free((uint32_t *)column_values_[0]->GetColumn());
	free((uint32_t *)column_values_[1]->GetColumn());
	free(column_values_[0]);
	free(column_values_[1]);
}

TEST_F(OptimizerTest, assembleTest2){
	uint64_t num_rows_ = 16*1048576;
	auto dice = std::bind(std::uniform_int_distribution<int32_t>(
									std::numeric_limits<int32_t>::min(),
									std::numeric_limits<int32_t>::max()),
									std::default_random_engine(std::time(0)));


	//construct two columns(17/32+17/32)
	uint32_t colWidthA = 17;
	uint32_t colWidthB = 17;
	Column* column_values_[2];
	uint32_t *values = (uint32_t *)malloc_aligned(num_rows_ * sizeof(uint32_t));
	uint32_t mask = (1 << colWidthA) - 1;
	for (uint32_t idx = 0; idx < num_rows_; ++idx) {
		values[idx] = dice() & mask;

	}
	SpecialiedColumn<uint32_t> *col = new SpecialiedColumn<uint32_t>();
	col->SetColumn(values);
	col->SetWidth(sizeof(uint32_t));
	column_values_[0] = col;

	values = (uint32_t *)malloc_aligned(num_rows_ * sizeof(uint32_t));
	mask = (1 << colWidthB) - 1;
	for (uint32_t idx = 0; idx < num_rows_; ++idx) {
		values[idx] = dice() & mask;

	}
	col = new SpecialiedColumn<uint32_t>();
	col->SetColumn(values);
	col->SetWidth(sizeof(uint32_t));
	column_values_[1] = col;
	{
		/** test case 2, stitch together to 32/32+2/16**/

		HybridTimer timer;
		timer.Start();


		Column* column_values_new_case2[2];

		uint32_t *values_case2 = (uint32_t *)malloc_aligned(num_rows_ * sizeof(uint32_t));
		SpecialiedColumn<uint32_t> *col_case2 = new SpecialiedColumn<uint32_t>();
		col_case2->SetColumn(values_case2);
		col_case2->SetWidth(sizeof(uint32_t));
		column_values_new_case2[0] = col_case2;

		uint16_t *values = (uint16_t *)malloc_aligned(num_rows_ * sizeof(uint16_t));
		SpecialiedColumn<uint16_t> *col = new SpecialiedColumn<uint16_t>();
		col->SetColumn(values);
		col->SetWidth(sizeof(uint16_t));
		column_values_new_case2[1] = col;

		std::vector < construct_entry_t > constructRecord_case2[2];
		construct_entry_t entry1;
		entry1.finalIdx = 0;
		entry1.start = 1;
		entry1.end = 17;
		constructRecord_case2[0].push_back(entry1);

		construct_entry_t entry2;
		entry2.finalIdx = 1;
		entry2.start = 1;
		entry2.end = 15;
		constructRecord_case2[0].push_back(entry2);

		construct_entry_t entry3;
		entry3.finalIdx = 1;
		entry3.start = 16;
		entry3.end = 17;
		constructRecord_case2[1].push_back(entry3);



		assemble(column_values_, column_values_new_case2, 2, num_rows_, constructRecord_case2);

		timer.Stop();
		std::cout << timer.GetNumCycles()/num_rows_ << "cycles/tuple" << std::endl;

		free((uint32_t *)column_values_new_case2[0]->GetColumn());
		free((uint16_t *)column_values_new_case2[1]->GetColumn());
		free(column_values_new_case2[0]);
		free(column_values_new_case2[1]);
	}

	free((uint32_t *)column_values_[0]->GetColumn());
	free((uint32_t *)column_values_[1]->GetColumn());
	free(column_values_[0]);
	free(column_values_[1]);
}

void OptimizerTest::assemble(Column** column_values_, Column** column_values_new,
		uint32_t new_col_num, uint64_t num_rows_, std::vector < construct_entry_t > *constructRecord) {
	uint32_t colIdx, entryIdx;
	uint32_t finalIdx, start, end;	//start is already declared
	//uint64_t assemble_val = 0;
	uint32_t portion_width = 0;
	//uint32_t kPrefetchDistance = 0;
#if 0	//version 1
	std::vector < construct_entry_t > records;
	for (colIdx = 0; colIdx < new_col_num; ++colIdx) {
		for (rowIdx = 0; rowIdx < num_rows_; ++rowIdx) {
			records = constructRecord[colIdx];
			assert(!records.empty());
			assemble_val = 0;
			for (entryIdx = 0; entryIdx < records.size(); ++entryIdx) {
				finalIdx = records[entryIdx].finalIdx;
				start = records[entryIdx].start;
				end = records[entryIdx].end;
				portion_width = end - start + 1;

				assemble_val <<= portion_width;	//shift enough bits for next portion

				//retrieve the value for column[finalIdx]'s <start, end> value
				assemble_val |= (column_values_[finalIdx]->GetValPortion(rowIdx, start, end));

			}
			//assign assemble_val to the target column
			column_values_new[colIdx]->SetValueAt(assemble_val, rowIdx);

			//cardinality extraction (not extract for the LAST column)
			//if (colIdx != (new_col_num - 1)) {
			//	cardInfo[colIdx].insert(assemble_val);
			//}
		}
	}
#endif

#if 0	//version 2: switch for loop
	std::vector < construct_entry_t > records;
	for (colIdx = 0; colIdx < new_col_num; ++colIdx) {
		records = constructRecord[colIdx];
		assert(!records.empty());
		for (entryIdx = 0; entryIdx < records.size(); ++entryIdx) {
			finalIdx = records[entryIdx].finalIdx;
			start = records[entryIdx].start;
			end = records[entryIdx].end;
			portion_width = end - start + 1;

			for (rowIdx = 0; rowIdx < num_rows_; ++rowIdx) {

				//shift enough bits for next portion
				column_values_new[colIdx]->ShiftValueLeft(portion_width * (0 != entryIdx), rowIdx);


				//retrieve the value for column[finalIdx]'s <start, end> value
				column_values_new[colIdx]->LogicalOr((column_values_[finalIdx]->GetValPortion(rowIdx, start, end)), rowIdx);

			}
		}
	}
#endif

#if 0	//version 3: partial SIMD
	std::vector < construct_entry_t > records;
	for (colIdx = 0; colIdx < new_col_num; ++colIdx) {
		records = constructRecord[colIdx];
		assert(!records.empty());
		for (entryIdx = 0; entryIdx < records.size(); ++entryIdx) {
			finalIdx = records[entryIdx].finalIdx;
			start = records[entryIdx].start;
			end = records[entryIdx].end;
			portion_width = end - start + 1;

			if (0 != entryIdx) {
				column_values_new[colIdx]->ShiftLeftBatch(num_rows_, portion_width);
			}

			for (rowIdx = 0; rowIdx < num_rows_; ++rowIdx) {

				//retrieve the value for column[finalIdx]'s <start, end> value
				column_values_new[colIdx]->LogicalOr((column_values_[finalIdx]->GetValPortion(rowIdx, start, end)), rowIdx);
			}
		}
	}
#endif

#if 1	//version 3: full SIMD
	std::vector < construct_entry_t > records;
	for (colIdx = 0; colIdx < new_col_num; ++colIdx) {
		records = constructRecord[colIdx];
		assert(!records.empty());
		for (entryIdx = 0; entryIdx < records.size(); ++entryIdx) {
			finalIdx = records[entryIdx].finalIdx;
			start = records[entryIdx].start;
			end = records[entryIdx].end;
			portion_width = end - start + 1;

			if (0 != entryIdx) {
				column_values_new[colIdx]->ShiftLeftBatch(num_rows_, portion_width);
			}
			switch (column_values_[finalIdx]->GetWidth()) {
			case 1:
				break;
			case 2:
				column_values_new[colIdx]->AssembleFromColumn(
						reinterpret_cast<SpecialiedColumn<uint16_t> *>(column_values_[finalIdx]), start, end, num_rows_);
				break;
			case 4:
				column_values_new[colIdx]->AssembleFromColumn(
						reinterpret_cast<SpecialiedColumn<uint32_t> *>(column_values_[finalIdx]), start, end, num_rows_);
				break;
			case 8:
				break;

			}

		}
	}
#endif

}


#if 0
TEST_F(MiscTest, stitchTest){

	uint32_t oldCol1[3] = {14, 16, 20};
	uint32_t newCol1[1] = {50};
	stitch(oldCol1, 3, newCol1, 1);

	std::cout << std::endl << std::endl;
	uint32_t oldCol2[1] = {63};
	uint32_t newCol2[3] = {21, 21, 21};
	stitch(oldCol2, 1, newCol2, 3);

	std::cout << std::endl << std::endl;
	uint32_t oldCol3[2] = {63, 20};
	uint32_t newCol3[4] = {31, 21, 21, 10};
	stitch(oldCol3, 2, newCol3, 4);

	std::cout << std::endl << std::endl;
	uint32_t oldCol4[4] = {16, 20, 30, 1};
	uint32_t newCol4[2] = {30, 37};
	stitch(oldCol4, 4, newCol4, 2);

	std::cout << std::endl << std::endl;
	uint32_t oldCol5[3] = {14, 16, 20};
	uint32_t newCol5[3] = {14, 16, 20};
	stitch(oldCol5, 3, newCol5, 3);

	std::cout << std::endl << std::endl;
	uint32_t oldCol6[3] = {16, 20, 30};
	uint32_t newCol6[3] = {16, 32, 18};
	stitch(oldCol6, 3, newCol6, 3);
}
#endif

void OptimizerTest::stitch(uint32_t * oldCol, size_t oldSize, uint32_t * newCol, size_t newSize) {
	std::vector<construct_entry_t> constructRecord[newSize];
	uint32_t curOldColIdx = 0;
	uint32_t curNewColIdx = 0;
	uint32_t start = 1; 	//[start, end] is from the MSB to LSB, as the location pointer for current old column
	uint32_t minVal = 0;
	uint32_t curOldColWidth = oldCol[curOldColIdx];
	uint32_t curNewColWidth = newCol[curNewColIdx];
	while (curNewColIdx < newSize) {
		//curOldColWidth = bw[curOldColIdx];
		//curNewColWidth = bestGlobalSplit[curNewColIdx];
		minVal = (curOldColWidth < curNewColWidth) ? curOldColWidth : curNewColWidth;

		curOldColWidth -= minVal;
		curNewColWidth -= minVal;

		//insert a new entry
		construct_entry_t entry;
		entry.finalIdx = curOldColIdx;
		entry.start = start;
		entry.end = start + minVal - 1;
		constructRecord[curNewColIdx].push_back(entry);
		start += minVal;

		if (0 == curOldColWidth) {	//an old column is consumed up
			curOldColIdx++;
			if (curOldColIdx < oldSize) {
				curOldColWidth = oldCol[curOldColIdx];
			}
			start = 1;
			//std::cout << "to proceed to next old column: " << curOldColIdx << std::endl;
		}

		if (0 == curNewColWidth) {	//a new column is consumed up
			curNewColIdx++;
			if (curNewColIdx < newSize) {
				curNewColWidth = newCol[curNewColIdx];
			}
		}
	}
	assert(curOldColIdx == oldSize);
	assert(curNewColIdx == newSize);

	//print the assemble result
	//std::cout << "to print result: " << std::endl;
	for (uint32_t i = 0; i < newSize; ++i) {
		//print constructRecord[i]
		std::cout << "For the new column " << i << std::endl;
		for (uint32_t j = 0; j < constructRecord[i].size(); ++j) {
			std::cout << constructRecord[i][j].finalIdx << "\t[" <<
					constructRecord[i][j].start << ", " << constructRecord[i][j].end << "]\n";
		}
	}
}


#if 0
TEST_F(MiscTest, combinations){
	uint32_t n;
	std::cout << "To calculate all combinations of #: ";
	std::cin >> n;
	uint32_t* arr = (uint32_t*) malloc_aligned(n * sizeof(uint32_t));
	uint32_t i;
	for (i = 0; i < n; ++i) {
		arr[i] = i;
	}

	std::vector<uint32_t *> combinations;
	allRange(arr, 0, n-1, combinations);	//

	free(arr);
	for (i = 0; i < combinations.size(); ++i) {
		free(combinations[i]);
	}
}
#endif

void OptimizerTest::swap(uint32_t * a, uint32_t * b) {
	uint32_t temp = *a;
	*a = *b;
	*b = temp;
}

void OptimizerTest::allRange(uint32_t * arr, uint32_t k, uint32_t m, std::vector<uint32_t *> &comb) {
	uint32_t i;
	if (k == m) {
		//deep copy of elements
		uint32_t* newArr = (uint32_t *)malloc_aligned((m+1) * sizeof(uint32_t));
		for (i = 0; i <= m; ++i) {
			newArr[i] = arr[i];
			std::cout << arr[i] << "\t";
		}
		std::cout << std::endl;
		comb.push_back(newArr);
	} else {
		for (i = k; i <= m; ++i) {	//第i个数分别与它后面的数字交换就能得到新的排列
			swap(arr + k, arr + i);
			allRange(arr, k + 1, m, comb);
			swap(arr + k, arr + i);
		}
	}
}

void OptimizerTest::split(int n, int curLevel, const int maxLevel,
		const int maxValue, const int threshold_loc, const int sum)//curLevel是数组中的位置，亦是递归层数
{
	//int128_t a;
    int i;
    if(curLevel > maxLevel) {
        return;
    } else if (0 == n){
    	display(curLevel);//分解完成，输出结果
    	count++;
    } else {
        for((n > maxValue) ? (i = maxValue) : (i = n); i > 0; i--)
            //if((0 == curLevel) || (i<=x[k-1]))
            {
                partitions[curLevel] = i;//写入数组

                if ((n - i) > 0) {	//i is not the last round width
                	if (((sum-n) <= threshold_loc) && ((sum-n+i) >= threshold_loc)) { //if i is overlapped with threshold_loc, do not require the alignment
                		split(n - i, curLevel + 1, maxLevel, maxValue, threshold_loc, sum);
                	} else if ((0 == (i % 16)) && (i != 48)) {
                		split(n - i, curLevel + 1, maxLevel, maxValue, threshold_loc, sum);
                	}
                } else {
                	split(n - i, curLevel + 1, maxLevel, maxValue, threshold_loc, sum);
                }
            }
    }
}

void OptimizerTest::display(int curLevel)
{
    int i;
    for(i = 0; i< curLevel; i++) {
    	std::cout << partitions[i] << " ";
    }
    std::cout << std::endl;
}

#if 0
TEST_F(MiscTest, integerDivide){
	int sum = 0;
	//int nparts = 0;	//divide <sum> into <nparts> with each one >=1 and <=64
	int maxRounds = 0;
	int threshold_loc = 0;
	std::cout << "sum\t maxrounds \t threshold_loc:" << std::endl;
	std::cin >> sum >> maxRounds >> threshold_loc;

	//for (nparts = 2; nparts <= maxRounds; ++nparts) {
		split(sum, 0, maxRounds, 64, threshold_loc, sum);
	//}

	std::cout << "Total division #" << count << std::endl;
}
#endif

#if 0
TEST_F(MiscTest, constAsArraySize){
	int fanout = 0;
	fanout = std::rand() % 256;

	const int arraysz = fanout;
	int array[arraysz];
	array[0] = 0;
	std::cout << fanout << std::endl << array[0] << std::endl;
}


TEST_F(MiscTest, unsignedIntConversion){

	uint64_t wider = 16000ULL;
	uint8_t shorter = 100;

	std::cout << "uint64_t to uint32_t: " << (uint32_t)(wider) << std::endl;
	std::cout << "uint64_t to uint16_t: " << (uint16_t)(wider) << std::endl;
	std::cout << "uint8_t to uint16_t: " << (uint16_t)(shorter) << std::endl;
	std::cout << "uint8_t to uint32_t: " << (uint32_t)(shorter) << std::endl;
	std::cout << "uint8_t to uint64_t: " << (uint64_t)(shorter) << std::endl;
}

TEST_F(MiscTest, movebitsleft){
	uint64_t val1 = 0x123442311234;
	uint64_t val2 = 0x5678876556789;
	//uint32_t bitw1 = 48;
	uint32_t bitw2 = 52;

	uint32_t nbits = 4;

	uint64_t val1_moved = (val1 << nbits) | (val2 >> (bitw2 - nbits));
	uint64_t val2_moved = ((1ULL << (bitw2-nbits)) - 1) & val2;

	uint64_t val1_moved_r = (val1 >> nbits);
	uint64_t val2_moved_r = ((((1ULL << nbits) - 1) & val1) << (bitw2)) | val2;

	std::cout << std::hex << val1_moved << "\t" << val2_moved << std::endl;
	std::cout << std::hex << val1_moved_r << "\t" << val2_moved_r << std::endl;


	std::string s = std::string("12l");
	std::cout << std::dec << (s[s.length()-1] == 'r') << "\t"	//0 is left, 1 is right
			<< (int32_t)(atoi(s.c_str()));
}
#endif
#if 0
TEST_F(MiscTest, memalloc_free_cost){
	uint64_t alloc_cost_free[30];
	uint64_t alloc_free[30];

	uint32_t byte_num = 0;
	uint32_t idx = 0;
	HybridTimer timer1;
	for (byte_num = 2; byte_num <= 1024*1024*1024; byte_num *= 2) {
		assert(idx < 30);
		timer1.Start();

		uint32_t *values = (uint32_t *)malloc_aligned(byte_num);
		timer1.Stop();
		alloc_cost[idx] = timer1.GetNumCycles();

		timer1.Start();
		free(values);
		timer1.Stop();
		alloc_free[idx] = timer1.GetNumCycles();

		idx++;
	}

	//output:
	for (idx = 0; idx < 30; ++idx) {
		std::cout <<
	}
}
#endif

#if 0
TEST_F(MiscTest, profiling){
	uint32_t num_rows = 64*1024*1024;

	HybridTimer timer1;
	timer1.Start();

	uint32_t *values = (uint32_t *)malloc_aligned(num_rows * sizeof(uint32_t));
	timer1.Stop();
	//std::cout << "memory allocation: " << ((double)timer1.GetNumCycles())/(num_rows * sizeof(uint32_t)) << "cycles" << std::endl;
	std::cout << "memory allocation: " << timer1.GetNumCycles() << "cycles" << std::endl;

	timer1.Start();
	for (uint32_t i = 0; i < num_rows; ++i) {
		values[i] = i;
	}
	timer1.Stop();
	std::cout << "sequential access: " << timer1.GetNumCycles()/(num_rows * sizeof(uint32_t)/CACHE_LINE_SIZE) << "cycles" << std::endl;

    auto dice = std::bind(std::uniform_int_distribution<int32_t>(
                            0,
							num_rows-1),
                            std::default_random_engine(std::time(0)));

	timer1.Start();
	for (uint32_t i = 0; i < num_rows; ++i) {
		values[dice()]++;
	}
	timer1.Stop();
	std::cout << "random access: " << timer1.GetNumCycles()/num_rows << "cycles" << std::endl;

    auto dice2 = std::bind(std::uniform_int_distribution<int32_t>(
                            std::numeric_limits<int32_t>::min(),
                            std::numeric_limits<int32_t>::max()),
                            std::default_random_engine(std::time(0)));

    timer1.Start();
	std::vector<element_t<uint32_t> > comparedCol;
	element_t<uint32_t> tmpElement;
	uint32_t idx;
	for (idx = 0; idx < 256; ++idx) {
		tmpElement.value 	= dice2();
		tmpElement.oid		= idx;
		comparedCol.push_back(tmpElement);
	}

	std::sort(comparedCol.begin(), comparedCol.end(), mycomparison<uint32_t>());
	timer1.Stop();
	std::cout << "sorting 256 elements: " << timer1.GetNumCycles() << "cycles" << std::endl;


	timer1.Start();

	free(values);
	timer1.Stop();
	//std::cout << "memory free: " << ((double)timer1.GetNumCycles())/(num_rows * sizeof(uint32_t)) << "cycles/Byte" << std::endl;
	std::cout << "memory free: " << timer1.GetNumCycles() << "cycles" << std::endl;


}
#endif

#if 0
TEST_F(MiscTest, STLsort) {

	HybridTimer timer1;
	{
		auto dice1 = std::bind(
				std::uniform_int_distribution<int32_t>(
						std::numeric_limits<int32_t>::min(),
						std::numeric_limits<int32_t>::max()),
				std::default_random_engine(std::time(0)));
		uint32_t values[256];
		uint32_t oids[256];
		for (uint32_t i = 0; i < 256; ++i) {
			values[i] = dice1();
			oids[i] = i;
		}

		for (uint32_t size = 256; size > 1; size--) {
			int repeat = 1000;
			timer1.Start();
			for (int tmp = 0; tmp < repeat; ++tmp) {

				std::vector<element_t<int32_t> > comparedCol;
				element_t<int32_t> tmpElement;
				uint32_t idx;

				for (idx = 0; idx < size; ++idx) {
					tmpElement.value = values[idx];
					tmpElement.oid = oids[idx];
					comparedCol.push_back(tmpElement);
				}

				std::sort(comparedCol.begin(), comparedCol.end(),
						mycomparison<int32_t>());
			}
			timer1.Stop();
			std::cout << "sorting " << size << " int32_t elements: "
					<< ((double) timer1.GetNumCycles()) / repeat
							/ (size * log(size)) << "cycles/tuple" << std::endl;
		}
	}
#if 0
	{
	auto dice2 = std::bind(
			std::uniform_int_distribution<int64_t>(
					std::numeric_limits<int64_t>::min(),
					std::numeric_limits<int64_t>::max()),
			std::default_random_engine(std::time(0)));

	timer1.Start();
	std::vector<element_t<int64_t> > comparedCol2;
	element_t<int64_t> tmpElement2;
	uint32_t idx;
	for (idx = 0; idx < 256; ++idx) {
		tmpElement2.value = dice2();
		tmpElement2.oid = idx;
		comparedCol2.push_back(tmpElement2);
	}

	std::sort(comparedCol2.begin(), comparedCol2.end(), mycomparison<int64_t>());
	timer1.Stop();
	std::cout << "sorting 256 int64_t elements: " << timer1.GetNumCycles() << "cycles"
			<< std::endl;
	}

	auto dice3 = std::bind(
			std::uniform_int_distribution<int16_t>(
					std::numeric_limits<int16_t>::min(),
					std::numeric_limits<int16_t>::max()),
			std::default_random_engine(std::time(0)));

	timer1.Start();
	std::vector<element_t<int16_t> > comparedCol3;
	element_t<int16_t> tmpElement3;
	uint32_t idx;
	for (idx = 0; idx < 256; ++idx) {
		tmpElement3.value = dice3();
		tmpElement3.oid = idx;
		comparedCol3.push_back(tmpElement3);
	}

	std::sort(comparedCol3.begin(), comparedCol3.end(), mycomparison<int16_t>());
	timer1.Stop();
	std::cout << "sorting 256 int16_t elements: " << timer1.GetNumCycles() << "cycles"
			<< std::endl;
#endif
}
#endif

#if 0
TEST_F(MiscTest, endianess){
	element_t A;
	A.value = 100;
	A.oid = 11;

	element_t *ptr_A = &A;
	uint64_t *ptr_int64 = (uint64_t *)ptr_A;
	uint8_t *ptr_char = (uint8_t *)ptr_int64;
	std::cout << std::hex << *ptr_int64 << std::endl;
	std::cout << "[0]: " << (unsigned int)ptr_char[0] << std::endl
		<< "[1]: " << (unsigned int)ptr_char[1] << std::endl
		<< "[2]: " << (unsigned int)ptr_char[2] << std::endl
		<< "[3]: " << (unsigned int)ptr_char[3] << std::endl
		<< "[4]: " << (unsigned int)ptr_char[4] << std::endl
		<< "[5]: " << (unsigned int)ptr_char[5] << std::endl
		<< "[6]: " << (unsigned int)ptr_char[6] << std::endl
		<< "[7]: " << (unsigned int)ptr_char[7] << std::endl;
}
#endif

}   //namespace
