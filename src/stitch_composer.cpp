/*******************************************************************************
 * Copyright (c) 2016
 * The Hong Kong Polytechnic University, Database Group
 *
 * Author: Wenjian Xu (cswxu AT comp DOT polyu.edu.hk)
 *
 * See file LICENSE.md for details.
 *******************************************************************************/

#include "stitch_composer.h"
#include "chain_composer.h"
#include "plan_enum_rrs.h"
#include <algorithm>
#include <unistd.h>

namespace multiAttrSort {

void StitchComposer::swap_val(uint32_t * a, uint32_t * b) {
	uint32_t temp = *a;
	*a = *b;
	*b = temp;
}

void StitchComposer::permutations(uint32_t * arr, uint32_t k, uint32_t m,
		std::vector<uint32_t *> &comb) {
	uint32_t i;
	if (k == m) {
		//deep copy of elements
		uint32_t* newArr = (uint32_t *) malloc_aligned(
				(m + 1) * sizeof(uint32_t));
		for (i = 0; i <= m; ++i) {
			newArr[i] = arr[i];
			//std::cout << arr[i] << "\t";
		}
		//std::cout << std::endl;
		comb.push_back(newArr);
	} else {
		for (i = k; i <= m; ++i) {	//第i个数分别与它后面的数字交换就能得到新的排列
			swap_val(arr + k, arr + i);
			permutations(arr, k + 1, m, comb);
			swap_val(arr + k, arr + i);
		}
	}
}

#if 0
uint32_t StitchComposer::estimateThresholdLoc(uint32_t *targetComb) {

	uint64_t threshold_ngroups_lb = num_rows_ / THRESHOLD_AVG_GROUP_SIZE_UB;
	uint64_t threshold_ngroups_ub = num_rows_ / THRESHOLD_AVG_GROUP_SIZE_LB;

	uint32_t virtualColIdx = 0;
	uint32_t realColIdx = targetComb[virtualColIdx];
	uint64_t cur_accu_ngroups = 1;
	assert(realColIdx >= 0);
	assert(realColIdx < num_columns_);
	uint64_t next_accu_ngroups = column_values_[realColIdx]->GetCardinality();
	uint32_t cur_accu_columwidth_sum = 0;	//may be better renamed as accu_factor!!

	/** find the threshold in which **/
	while (next_accu_ngroups < threshold_ngroups_lb) {

		cur_accu_ngroups = next_accu_ngroups;
		cur_accu_columwidth_sum += compose_params_.bitwidth[realColIdx];

		virtualColIdx++;
		if (virtualColIdx == num_columns_) {
			break;
		}
		realColIdx = targetComb[virtualColIdx];
		assert(realColIdx >= 0);
		assert(realColIdx < num_columns_);
		next_accu_ngroups *= column_values_[realColIdx]->GetCardinality();
	}

	if (virtualColIdx == num_columns_) {	//multiplication of all column cardinality cannot hit the threshold
		return compose_params_.colwidth_sum;
	} else if ((next_accu_ngroups >= threshold_ngroups_lb) && (next_accu_ngroups <= threshold_ngroups_ub)) {	//the threhshold is just in the <realColIdx>-th column's division

		assert((cur_accu_columwidth_sum + compose_params_.bitwidth[realColIdx] - 1) < compose_params_.colwidth_sum);
		return cur_accu_columwidth_sum + compose_params_.bitwidth[realColIdx] - 1;
	}
	else {
		//the concrete threshold is within <realColIdx>-th column
		uint32_t nbits_shifted_right = 0;
		switch(column_values_[realColIdx]->GetWidth()) {
		case 1:
			nbits_shifted_right = estimateDetailedLoc<uint8_t>(
					cur_accu_ngroups, threshold_ngroups_lb, threshold_ngroups_ub,
					(uint8_t *)(column_values_[realColIdx]->GetColumn()), compose_params_.bitwidth[realColIdx]);
			break;
		case 2:
			nbits_shifted_right = estimateDetailedLoc<uint16_t>(
					cur_accu_ngroups, threshold_ngroups_lb, threshold_ngroups_ub,
					(uint16_t *)(column_values_[realColIdx]->GetColumn()), compose_params_.bitwidth[realColIdx]);
			break;
		case 4:
			nbits_shifted_right = estimateDetailedLoc<uint32_t>(
					cur_accu_ngroups, threshold_ngroups_lb, threshold_ngroups_ub,
					(uint32_t *)(column_values_[realColIdx]->GetColumn()), compose_params_.bitwidth[realColIdx]);
			break;
		case 8:
			nbits_shifted_right = estimateDetailedLoc<uint64_t>(
					cur_accu_ngroups, threshold_ngroups_lb, threshold_ngroups_ub,
					(uint64_t *)(column_values_[realColIdx]->GetColumn()), compose_params_.bitwidth[realColIdx]);
			break;
		}
		assert(nbits_shifted_right > 0);
		assert(nbits_shifted_right < compose_params_.bitwidth[realColIdx]);

		return cur_accu_columwidth_sum + (compose_params_.bitwidth[realColIdx]-nbits_shifted_right) - 1;
	}
}
#endif

#if 0
template <class T>
uint32_t StitchComposer::estimateDetailedLoc(uint64_t cur_accu_ngroups, uint64_t threshold_lb, uint64_t threshold_ub,
					T* values, uint32_t width) {
	assert(threshold_lb < threshold_ub);
	uint32_t nbits_right_shifted = 1;
	width--;
	//copy the content with right shifted 1 bit
	T* copiedValues = (T *)malloc_aligned(num_rows_ * sizeof(T));
	std::set<T> container;
	uint64_t i;
	for (i = 0; i < num_rows_; ++i) {
		copiedValues[i] = values[i] >> 1;
		container.insert(copiedValues[i]);
	}

	uint32_t updated_ngroups = container.size();

	while ((cur_accu_ngroups * updated_ngroups > threshold_ub)
			|| (cur_accu_ngroups * updated_ngroups < threshold_lb)) {

		std::set<T>().swap(container);
		assert(container.empty());
		for (i = 0; i < num_rows_; ++i) {
			copiedValues[i] >>= 1;
			container.insert(copiedValues[i]);
		}
		updated_ngroups = container.size();

		nbits_right_shifted++;
		width--;
		assert(width > 0);
	}

	free(copiedValues);

	return nbits_right_shifted;
}
#endif

void StitchComposer::split_exhaustive_all_round(int n, int curLevel, const int maxLevel,
		const int maxValue, const int sum,
		std::vector<std::vector<uint32_t> > &split_enum) {

	int i;
	if (curLevel > maxLevel) {
		return;
	} else if (0 == n) {
		//display(curLevel);//分解完成，输出结果
		std::vector < uint32_t > aSplit;
		for (i = 0; i < curLevel; i++) {
			aSplit.push_back(partitions[i]);
		}
		split_enum.push_back(aSplit);
		//count++;
	} else {
		for ((n > maxValue) ? (i = maxValue) : (i = n); i > 0; i--)
		//if((0 == curLevel) || (i<=x[k-1]))
				{
			partitions[curLevel] = i;	//写入数组

			split_exhaustive_all_round(n - i, curLevel + 1, maxLevel, maxValue, sum, split_enum);

#if 0
			if ((n - i) > 0) {	//i is not the last round width
				if (((sum - n) <= threshold_loc)
						&& ((sum - n + i) >= threshold_loc)) { //if i is overlapped with threshold_loc, do not require the alignment
					split_exhaustive_all_round(n - i, curLevel + 1, maxLevel, maxValue,
							threshold_loc, sum, split_enum);
				} else if ((0 == (i % 16)) && (i != 48)) {
					split_exhaustive_all_round(n - i, curLevel + 1, maxLevel, maxValue,
							threshold_loc, sum, split_enum);
				}
			} else {
				split_exhaustive_all_round(n - i, curLevel + 1, maxLevel, maxValue, threshold_loc,
						sum, split_enum);
			}
#endif

#if 0
			if ((n - i) > 0) {	//i is not the last round width
				if ((sum - n + i) >= threshold_loc) { //revised: if i is after threshold_loc, do not require the alignment
					split_exhaustive_all_round(n - i, curLevel + 1, maxLevel, maxValue,
							threshold_loc, sum, split_enum);
				} else if ((0 == (i % 16)) && (i != 48)) {
					split_exhaustive_all_round(n - i, curLevel + 1, maxLevel, maxValue,
							threshold_loc, sum, split_enum);
				}
			} else {
				split_exhaustive_all_round(n - i, curLevel + 1, maxLevel, maxValue, threshold_loc,
						sum, split_enum);
			}
#endif

		}
	}
}

void StitchComposer::split_exhaustive_certain_round(int n, int curRound, const int round_num,
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
				split_exhaustive_certain_round(n - i, curRound + 1, round_num, maxValue, split_enum);
			}
		}
	}
}

#if 0
void StitchComposer::split_applyRule1(int n, int curLevel, const int maxLevel,
		const int maxValue, std::vector<std::vector<uint32_t> > &split_enum, bool isBroken) {

	int i;
	if (curLevel > maxLevel) {
		return;
	} else if (0 == n) {
		//display(curLevel);//分解完成，输出结果
		std::vector < uint32_t > aSplit;
		for (i = 0; i < curLevel; i++) {
			aSplit.push_back(partitions[i]);
		}
		split_enum.push_back(aSplit);
		//count++;
	} else {
		for ((n > maxValue) ? (i = maxValue) : (i = n); i > 0; i--)
		{
			bool memoIsBroken = isBroken;
			partitions[curLevel] = i;	//写入数组

			//split(n - i, curLevel + 1, maxLevel, maxValue, threshold_loc, sum, split_enum);
#if 1
			if ((n - i) > 0) {	//i is not the last round width
				if ((0 == (i % 16)) && (i != 48)) {
					split_applyRule1(n - i, curLevel + 1, maxLevel, maxValue, split_enum, isBroken);
				} else if (!isBroken) {
					isBroken = true;
					split_applyRule1(n - i, curLevel + 1, maxLevel, maxValue, split_enum, isBroken);
				}
			} else {
				split_applyRule1(n - i, curLevel + 1, maxLevel, maxValue, split_enum, isBroken);
			}
			isBroken = memoIsBroken;
#endif

#if 0
			if ((n - i) > 0) {	//i is not the last round width
				if ((sum - n + i) >= threshold_loc) { //revised: if i is after threshold_loc, do not require the alignment
					split_exhaustive_all_round(n - i, curLevel + 1, maxLevel, maxValue,
							threshold_loc, sum, split_enum);
				} else if ((0 == (i % 16)) && (i != 48)) {
					split_exhaustive_all_round(n - i, curLevel + 1, maxLevel, maxValue,
							threshold_loc, sum, split_enum);
				}
			} else {
				split_exhaustive_all_round(n - i, curLevel + 1, maxLevel, maxValue, threshold_loc,
						sum, split_enum);
			}
#endif

		}
	}
}
#endif

void StitchComposer::getCardinalityFromStatistics(uint32_t *columnOrder, std::vector<uint32_t> &curSplit,
		double* cardinality, uint32_t assembleStatistics[][3], uint32_t shiftStatistics[]) {

	/** stitching statistics, for costing stitching/shifting overhead
	 *
	 * for matrix assembleStatistics, x axis is ToColumns, y axis is FromColumns
	 */

	/** original bitwidth **/
	uint32_t *bw = compose_params_.bitwidth;

	const uint32_t new_col_num = curSplit.size();

	/**
	 * Step a) Find the shifting operations for each new column
	 */
	std::vector < construct_entry_t > constructRecord[new_col_num];
	uint32_t curOldColIdx = 0;
	uint32_t curNewColIdx = 0;
	uint32_t start = 1; //[start, end] is from the MSB to LSB, as the location pointer for current old column
	uint32_t minVal = 0;
	uint32_t curOldColWidth = bw[columnOrder[curOldColIdx]];	//important!!!
	uint32_t curNewColWidth = curSplit[curNewColIdx];
	while (curNewColIdx < new_col_num) {
		//curOldColWidth = bw[curOldColIdx];
		//curNewColWidth = bestGlobalSplit[curNewColIdx];
		minVal =
				(curOldColWidth < curNewColWidth) ?
						curOldColWidth : curNewColWidth;

		curOldColWidth -= minVal;
		curNewColWidth -= minVal;

		//insert a new entry
		construct_entry_t entry;
		entry.finalIdx = columnOrder[curOldColIdx];	//the index is transferred by columnOrder
		entry.start = start;
		entry.end = start + minVal - 1;
		constructRecord[curNewColIdx].push_back(entry);
		start += minVal;

		/*********************update stitching/shifting statistics************************/
		if (curSplit[curNewColIdx] <= 16) {
			shiftStatistics[0]++;

			if (bw[columnOrder[curOldColIdx]] <= 16) {
				assembleStatistics[0][0]++;
			} else if (bw[columnOrder[curOldColIdx]] <= 32) {
				assembleStatistics[1][0]++;
			} else {
				assembleStatistics[2][0]++;
			}

		} else if (curSplit[curNewColIdx] <= 32) {
			shiftStatistics[1]++;

			if (bw[columnOrder[curOldColIdx]] <= 16) {
				assembleStatistics[0][1]++;
			} else if (bw[columnOrder[curOldColIdx]] <= 32) {
				assembleStatistics[1][1]++;
			} else {
				assembleStatistics[2][1]++;
			}
		} else {	//	<=64
			shiftStatistics[2]++;

			if (bw[columnOrder[curOldColIdx]] <= 16) {
				assembleStatistics[0][2]++;
			} else if (bw[columnOrder[curOldColIdx]] <= 32) {
				assembleStatistics[1][2]++;
			} else {
				assembleStatistics[2][2]++;
			}
		}
		/*********************************************************************************/

		if (0 == curOldColWidth) {	//an old column is consumed up
			curOldColIdx++;
			if (curOldColIdx < num_columns_) {
				curOldColWidth = bw[columnOrder[curOldColIdx]];
			}
			start = 1;
		}

		if (0 == curNewColWidth) {	//a new column is consumed up
			curNewColIdx++;
			if (curNewColIdx < new_col_num) {
				curNewColWidth = curSplit[curNewColIdx];
			}
		}
	}
	assert(curOldColIdx == num_columns_);
	assert(curNewColIdx == new_col_num);

	/**
	 * Step b) Estimate the cardinality info from cardinality of base columns
	 *
	 * use this info: uint32_t *bw = compose_params_.bitwidth;
	 * NOTE: Below commented operation is done by constructor of StitchComposer
	 */
#if 0
	const uint32_t origin_col_num = num_columns_;
	uint64_t base_col_card[origin_col_num];
	double card_per_bit[origin_col_num];

	/**calculate the bit info carried by each bit in each base column**/
	for (uint32_t colIdx = 0; colIdx < origin_col_num; ++colIdx) {
		base_col_card[colIdx] = column_values_[columnOrder[colIdx]]->GetCardinality();
		card_per_bit[colIdx] = pow((double)(base_col_card[colIdx]), 1.0/(double)(bw[columnOrder[colIdx]]));
	}
#endif

	/** use this:
	 * std::vector < construct_entry_t > constructRecord[new_col_num];
	 * to *assemble* the cardinality info
	 */
	std::vector < construct_entry_t > records;
	uint32_t finalIdx;
	uint32_t end;
	uint32_t portion_width;
	//std::set<uint32_t> columns_crossed;
	//uint32_t crossed_num = 0;

	for (uint32_t new_colIdx = 0; new_colIdx < (new_col_num -1); ++new_colIdx) {
		records = constructRecord[new_colIdx];
		assert(!records.empty());

		//columns_crossed.clear();
		for (uint32_t entryIdx = 0; entryIdx < records.size(); ++entryIdx) {
			finalIdx = records[entryIdx].finalIdx;
			start = records[entryIdx].start;
			end = records[entryIdx].end;
			portion_width = end - start + 1;

			cardinality[new_colIdx] *= pow(card_per_bit[finalIdx], (double)portion_width);

			//columns_crossed.insert(finalIdx);
#if 0
			std::cout << "assembled column "
					<< new_colIdx << " cardinality: " << cardinality[new_colIdx] << std::endl;
#endif
		}
#if 0
		crossed_num = columns_crossed.size();
		if (crossed_num > 1) {

			cardinality[new_colIdx] *= pow(shrink_factor, (double)(crossed_num-1));
		}
#endif
	}

#if 0
	/**
	 * Step b) Re-construct each column according to ===>std::vector < construct_entry_t > constructRecord[new_col_num];
	 */

	/** Step b1) allocate the memories for all new columns: **/
	Column ** column_values_new = (Column **) malloc_aligned(new_col_num * sizeof(Column *));
	uint32_t idx;
	uint32_t nbits = 0;
	//std::cout << "#column: " << new_col_num << std::endl;
	for (idx = 0; idx < new_col_num; ++idx) {
		nbits = curSplit[idx];
		assert(nbits > 0 && nbits <= 64);

		if (nbits <= 16) { /*IMPORTANT: use uint16_t to deal with values less than 8 bits*/
			uint16_t *values = (uint16_t *)malloc_aligned(num_rows_ * sizeof(uint16_t));
			memset(values, 0, num_rows_ * sizeof(uint16_t));
			SpecialiedColumn<uint16_t> *col = new SpecialiedColumn<uint16_t>();
			col->SetColumn(values);
			col->SetWidth(sizeof(uint16_t));
			column_values_new[idx] = col;
		} else if (nbits > 16 && nbits <= 32) {
			uint32_t *values = (uint32_t *)malloc_aligned(num_rows_ * sizeof(uint32_t));
			memset(values, 0, num_rows_ * sizeof(uint32_t));
			SpecialiedColumn<uint32_t> *col = new SpecialiedColumn<uint32_t>();
			col->SetColumn(values);
			col->SetWidth(sizeof(uint32_t));
			column_values_new[idx] = col;
		} else {	//width>32 && width<=64
			uint64_t *values = (uint64_t *)malloc_aligned(num_rows_ * sizeof(uint64_t));
			memset(values, 0, num_rows_ * sizeof(uint64_t));
			SpecialiedColumn<uint64_t> *col = new SpecialiedColumn<uint64_t>();
			col->SetColumn(values);
			col->SetWidth(sizeof(uint64_t));
			column_values_new[idx] = col;
		}
	}

	/** Step b2) re-construct the values **/
	/** Note about the column index issue, i.e., whether it is after transferred by columnOrder **/
	uint32_t colIdx, entryIdx;
	uint32_t finalIdx, end;	//start is already declared
	uint32_t portion_width = 0;

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
				std::cout << "[Error ] Not Allow 8b column width" << std::endl;
				exit(EXIT_SUCCESS);
				break;
			case 2:
				column_values_new[colIdx]->AssembleFromColumn(
						reinterpret_cast<SpecialiedColumn<uint16_t> *>(column_values_[finalIdx]), start, end, num_rows_, bw[finalIdx]);
				break;
			case 4:
				column_values_new[colIdx]->AssembleFromColumn(
						reinterpret_cast<SpecialiedColumn<uint32_t> *>(column_values_[finalIdx]), start, end, num_rows_, bw[finalIdx]);
				break;
			case 8:
				column_values_new[colIdx]->AssembleFromColumn(
						reinterpret_cast<SpecialiedColumn<uint64_t> *>(column_values_[finalIdx]), start, end, num_rows_, bw[finalIdx]);
				break;
			}

		}
	}

	//cardinality extraction
	std::set<uint64_t> cardInfo[new_col_num - 1];
	for (colIdx = 0; colIdx < new_col_num - 1; ++colIdx) {
		for (uint32_t rowIdx = 0; rowIdx < num_rows_; ++rowIdx) {
				cardInfo[colIdx].insert((uint64_t)column_values_new[colIdx]->GetValueAt(rowIdx));
			}
	}

	for (idx = 0; idx < (new_col_num -1); ++idx) {
		cardinality[idx] = cardInfo[idx].size();
	}

	/**
	 * free the memory for current get statistics
	 */
	uint32_t i;
	for (i = 0; i < new_col_num; ++i) {
		switch(column_values_new[i]->GetWidth()) {
			case 1:
				free((uint8_t *)column_values_new[i]->GetColumn());
				break;
			case 2:
				free((uint16_t *)column_values_new[i]->GetColumn());
				break;
			case 4:
				free((uint32_t *)column_values_new[i]->GetColumn());
				break;
			case 8:
				free((uint64_t *)column_values_new[i]->GetColumn());
				break;
		}
		delete column_values_new[i];
	}
#endif

}

double StitchComposer::run_single_plan(uint32_t *columnOrder, std::vector<uint32_t> &curSplit, double *breakdown) {

	/**
	 * construct new column values and feed them to a chain composer
	 * NOTE: use columnOrder and curSplit, copy most code from stitching() function
	 */
	setting_t setting_new = this->setting_;

	/** original bitwidth **/
	uint32_t *bw = compose_params_.bitwidth;

	const uint32_t new_col_num = curSplit.size();

	/**
	 * Step a) Find the shifting operations for each new column
	 */
	std::vector < construct_entry_t > constructRecord[new_col_num];
	uint32_t curOldColIdx = 0;
	uint32_t curNewColIdx = 0;
	uint32_t start = 1; //[start, end] is from the MSB to LSB, as the location pointer for current old column
	uint32_t minVal = 0;
	uint32_t curOldColWidth = bw[columnOrder[curOldColIdx]];	//important!!!
	uint32_t curNewColWidth = curSplit[curNewColIdx];
	while (curNewColIdx < new_col_num) {
		minVal =
				(curOldColWidth < curNewColWidth) ?
						curOldColWidth : curNewColWidth;

		curOldColWidth -= minVal;
		curNewColWidth -= minVal;

		//insert a new entry
		construct_entry_t entry;
		entry.finalIdx = columnOrder[curOldColIdx];	//the index is transferred by bestOrdering
		entry.start = start;
		entry.end = start + minVal - 1;
		constructRecord[curNewColIdx].push_back(entry);
		start += minVal;

		if (0 == curOldColWidth) {	//an old column is consumed up
			curOldColIdx++;
			if (curOldColIdx < num_columns_) {
				curOldColWidth = bw[columnOrder[curOldColIdx]];
			}
			start = 1;
		}

		if (0 == curNewColWidth) {	//a new column is consumed up
			curNewColIdx++;
			if (curNewColIdx < new_col_num) {
				curNewColWidth = curSplit[curNewColIdx];
			}
		}
	}
	assert(curOldColIdx == num_columns_);
	assert(curNewColIdx == new_col_num);
	/**
	 * Step b) Re-construct each column according to ===>std::vector < construct_entry_t > constructRecord[new_col_num];
	 */

	/** Step b1) allocate the memories for all new columns: **/
	Column ** column_values_new = (Column **) malloc_aligned(new_col_num * sizeof(Column *));
	uint32_t idx;
	uint32_t nbits = 0;
	for (idx = 0; idx < new_col_num; ++idx) {
		nbits = curSplit[idx];
		assert(nbits > 0 && nbits <= 64);

		//if (nbits <= 8) {
		//
		if (nbits <= 16) { /*IMPORTANT: use uint16_t to deal with values less than 8 bits*/
			uint16_t *values = (uint16_t *)malloc_aligned(num_rows_ * sizeof(uint16_t));
			memset(values, 0, num_rows_ * sizeof(uint16_t));
			SpecialiedColumn<uint16_t> *col = new SpecialiedColumn<uint16_t>();
			col->SetColumn(values);
			col->SetWidth(sizeof(uint16_t));
			column_values_new[idx] = col;
		} else if (nbits > 16 && nbits <= 32) {
			uint32_t *values = (uint32_t *)malloc_aligned(num_rows_ * sizeof(uint32_t));
			memset(values, 0, num_rows_ * sizeof(uint32_t));
			SpecialiedColumn<uint32_t> *col = new SpecialiedColumn<uint32_t>();
			col->SetColumn(values);
			col->SetWidth(sizeof(uint32_t));
			column_values_new[idx] = col;
		} else {	//width>32 && width<=64
			uint64_t *values = (uint64_t *)malloc_aligned(num_rows_ * sizeof(uint64_t));
			memset(values, 0, num_rows_ * sizeof(uint64_t));
			SpecialiedColumn<uint64_t> *col = new SpecialiedColumn<uint64_t>();
			col->SetColumn(values);
			col->SetWidth(sizeof(uint64_t));
			column_values_new[idx] = col;
		}
	}

	HybridTimer timer_stitching;
	timer_stitching.Start();

	/** Step b2) re-construct the values **/
	/** Note about the column index issue, i.e., whether it is after transferred by columnOrder **/
	uint32_t colIdx, entryIdx;
	uint32_t finalIdx, end;	//start is already declared
	uint32_t portion_width = 0;

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
				std::cout << "[Error ] Not Allow 8b column width" << std::endl;
				exit(EXIT_SUCCESS);
				break;
			case 2:
				column_values_new[colIdx]->AssembleFromColumn(
						reinterpret_cast<SpecialiedColumn<uint16_t> *>(column_values_[finalIdx]), start, end, num_rows_, bw[finalIdx]);
				break;
			case 4:
				column_values_new[colIdx]->AssembleFromColumn(
						reinterpret_cast<SpecialiedColumn<uint32_t> *>(column_values_[finalIdx]), start, end, num_rows_, bw[finalIdx]);
				break;
			case 8:
				column_values_new[colIdx]->AssembleFromColumn(
						reinterpret_cast<SpecialiedColumn<uint64_t> *>(column_values_[finalIdx]), start, end, num_rows_, bw[finalIdx]);
				break;
			}

		}
	}

	timer_stitching.Stop();
	setting_new.time_stitching = (double) timer_stitching.GetNumCycles();
	breakdown[0] = setting_new.time_stitching;

	/**
	 * Step c) feed new column values and other parameters to a chain composer (different from original stitching function)
	 **/
	setting_new.nthreads = setting_.nthreads;
	setting_new.intrinType = setting_.intrinType;
	setting_new.partition_fanout = setting_.partition_fanout;

	compose_params_t compose_params_new;
	compose_params_new.sortalgo = compose_params_.sortalgo;
	compose_params_new.packtype = compose_params_.packtype;
	compose_params_new.is_oid_encoded = compose_params_.is_oid_encoded;
	compose_params_new.colwidth_sum = compose_params_.colwidth_sum;
	compose_params_new.ordered = compose_params_.ordered;

	compose_params_new.column_asc_desc = (bool*) malloc_aligned(new_col_num * sizeof(bool));
	compose_params_new.bitwidth = (uint32_t*) malloc_aligned(new_col_num * sizeof(uint32_t));

	uint32_t i;
	for (i = 0; i < new_col_num; ++i) {
		compose_params_new.column_asc_desc[i] = true;
		compose_params_new.bitwidth[i] = curSplit[i];
	}

	/** note: use the chain composer **/
	Composer* composer_new = new ChainComposer(setting_new, new_col_num, num_rows_,
			compose_params_new, column_values_new);

    HybridTimer totalTimer;
	totalTimer.Start();

	composer_new->SortAllColumns();	//NOTE: this is ChainComposer's SortAllColumns, i.e., equivalent to multiRoundSorting

	totalTimer.Stop();

	setting_t updated = composer_new->getSetting();
	breakdown[1] = updated.time_tupleReconstruct;
	breakdown[2] = updated.time_orderGroupInfo;
	breakdown[3] = updated.time_multipass_sort;
	//breakdown[4] = updated.time_mmcpyOrderGroup;//this cost is obselete
	breakdown[4] = updated.time_firstpass_sort;
	breakdown[5] = updated.time_extractGrop;
	assert(breakdown[0] > 0);
	assert(breakdown[1] > 0);
	assert(breakdown[2] > 0);
	assert(breakdown[3] > 0);
	assert(breakdown[4] > 0);
	assert(breakdown[5] > 0);

#if 1
	//output group information:
	std::cout << "[INFO ] 2nd round group num: " << composer_new->getGroupInfo().ngroups
			<< ", group num larger than 2: " << composer_new->getGroupInfo().ngroups_gtOne
			<< ", avg group size: "
			<< (double)composer_new->getGroupInfo().sum_group_size/(double)composer_new->getGroupInfo().ngroups
			<< std::endl;
#endif

	/**
	 * Step d) free the memory for current sumilation
	 */
	for (i = 0; i < composer_new->num_columns_; ++i) {
		switch(composer_new->column_values_[i]->GetWidth()) {
			case 1:
				free((uint8_t *)composer_new->column_values_[i]->GetColumn());
				break;
			case 2:
				free((uint16_t *)composer_new->column_values_[i]->GetColumn());
				break;
			case 4:
				free((uint32_t *)composer_new->column_values_[i]->GetColumn());
				break;
			case 8:
				free((uint64_t *)composer_new->column_values_[i]->GetColumn());
				break;
		}
		delete composer_new->column_values_[i];
	}
	free(composer_new->compose_params_.column_asc_desc);
	free(composer_new->compose_params_.bitwidth);

	free(composer_new->column_values_);
	free(composer_new);

	//return totalTimer.GetNumCycles() + breakdown[0];
	return totalTimer.GetNumCycles() + breakdown[0];
}

double StitchComposer::run_single_plan_hash(uint32_t *columnOrder, std::vector<uint32_t> &curSplit, double *breakdown) {

	HybridTimer timer_stitching;
	timer_stitching.Start();

	/**
	 * construct new column values and feed them to a chain composer (HASHING)
	 */
	//setting_t setting_new = this->setting_;
	setting_t setting_new;

	/** original bitwidth **/
	uint32_t *bw = compose_params_.bitwidth;

	const uint32_t new_col_num = curSplit.size();

	/**
	 * Step a) Find the shifting operations for each new column
	 */
	std::vector < construct_entry_t > constructRecord[new_col_num];
	uint32_t curOldColIdx = 0;
	uint32_t curNewColIdx = 0;
	uint32_t start = 1; //[start, end] is from the MSB to LSB, as the location pointer for current old column
	uint32_t minVal = 0;
	uint32_t curOldColWidth = bw[columnOrder[curOldColIdx]];	//important!!!
	uint32_t curNewColWidth = curSplit[curNewColIdx];
	while (curNewColIdx < new_col_num) {
		//curOldColWidth = bw[curOldColIdx];
		//curNewColWidth = bestGlobalSplit[curNewColIdx];
		minVal =
				(curOldColWidth < curNewColWidth) ?
						curOldColWidth : curNewColWidth;

		curOldColWidth -= minVal;
		curNewColWidth -= minVal;

		//insert a new entry
		construct_entry_t entry;
		entry.finalIdx = columnOrder[curOldColIdx];	//the index is transferred by bestOrdering
		entry.start = start;
		entry.end = start + minVal - 1;
		constructRecord[curNewColIdx].push_back(entry);
		start += minVal;

		if (0 == curOldColWidth) {	//an old column is consumed up
			curOldColIdx++;
			if (curOldColIdx < num_columns_) {
				curOldColWidth = bw[columnOrder[curOldColIdx]];
			}
			start = 1;
		}

		if (0 == curNewColWidth) {	//a new column is consumed up
			curNewColIdx++;
			if (curNewColIdx < new_col_num) {
				curNewColWidth = curSplit[curNewColIdx];
			}
		}
	}
	assert(curOldColIdx == num_columns_);
	assert(curNewColIdx == new_col_num);

	/**
	 * Step b) Re-construct each column according to ===>std::vector < construct_entry_t > constructRecord[new_col_num];
	 */

	/** Step b1) allocate the memories for all new columns: **/
	Column ** column_values_new = (Column **) malloc_aligned(new_col_num * sizeof(Column *));
	uint32_t idx;
	uint32_t nbits = 0;
	//std::cout << "#column: " << new_col_num << std::endl;
	for (idx = 0; idx < new_col_num; ++idx) {
		nbits = curSplit[idx];
		assert(nbits > 0 && nbits <= 64);

		if (nbits <= 8) {	/*IMPORTANT: now we can have type uint8_t*/
			uint8_t *values = (uint8_t *)malloc_aligned(num_rows_ * sizeof(uint8_t));
			memset(values, 0, num_rows_ * sizeof(uint8_t));
			SpecialiedColumn<uint8_t> *col = new SpecialiedColumn<uint8_t>();
			col->SetColumn(values);
			col->SetWidth(sizeof(uint8_t));
			column_values_new[idx] = col;
		}
		else if (nbits > 8 && nbits <= 16) {
			uint16_t *values = (uint16_t *)malloc_aligned(num_rows_ * sizeof(uint16_t));
			memset(values, 0, num_rows_ * sizeof(uint16_t));
			SpecialiedColumn<uint16_t> *col = new SpecialiedColumn<uint16_t>();
			col->SetColumn(values);
			col->SetWidth(sizeof(uint16_t));
			column_values_new[idx] = col;
		} else if (nbits > 16 && nbits <= 32) {
			uint32_t *values = (uint32_t *)malloc_aligned(num_rows_ * sizeof(uint32_t));
			memset(values, 0, num_rows_ * sizeof(uint32_t));
			SpecialiedColumn<uint32_t> *col = new SpecialiedColumn<uint32_t>();
			col->SetColumn(values);
			col->SetWidth(sizeof(uint32_t));
			column_values_new[idx] = col;
		} else {	//width>32 && width<=64
			uint64_t *values = (uint64_t *)malloc_aligned(num_rows_ * sizeof(uint64_t));
			memset(values, 0, num_rows_ * sizeof(uint64_t));
			SpecialiedColumn<uint64_t> *col = new SpecialiedColumn<uint64_t>();
			col->SetColumn(values);
			col->SetWidth(sizeof(uint64_t));
			column_values_new[idx] = col;
		}
		//std::cout << "after allocating memory" << std::endl;
	}
	//printf("arrive here3\n");
	/** Step b2) re-construct the values **/
	/** Note about the column index issue, i.e., whether it is after transferred by columnOrder **/
	uint32_t colIdx, entryIdx;
	uint32_t finalIdx, end;	//start is already declared
	//uint64_t assemble_val = 0;
	uint32_t portion_width = 0;
	//uint32_t kPrefetchDistance = 0;

	std::vector < construct_entry_t > records;
	for (colIdx = 0; colIdx < new_col_num; ++colIdx) {
		records = constructRecord[colIdx];
		assert(!records.empty());
		for (entryIdx = 0; entryIdx < records.size(); ++entryIdx) {
			finalIdx = records[entryIdx].finalIdx;
			start = records[entryIdx].start;
			end = records[entryIdx].end;
			portion_width = end - start + 1;
#if 0
			std::cout << "finalIdx, start, end: " << finalIdx << ", " << start << ", " << end << std::endl;
#endif

			if (0 != entryIdx) {
				column_values_new[colIdx]->ShiftLeftBatch(num_rows_, portion_width);
			}
#if 0
			std::cout << "first value: " << column_values_[finalIdx]->GetValueAt(0);
			std::cout << "second value: " << column_values_[finalIdx]->GetValueAt(1);
#endif
			switch (column_values_[finalIdx]->GetWidth()) {
			case 1:
				//std::cout << "[Error ] Not Allow 8b column width" << std::endl;
				//exit(EXIT_SUCCESS);
				column_values_new[colIdx]->AssembleFromColumn(
						reinterpret_cast<SpecialiedColumn<uint8_t> *>(column_values_[finalIdx]), start, end, num_rows_, bw[finalIdx]);
				break;
			case 2:
				column_values_new[colIdx]->AssembleFromColumn(
						reinterpret_cast<SpecialiedColumn<uint16_t> *>(column_values_[finalIdx]), start, end, num_rows_, bw[finalIdx]);
				break;
			case 4:
				column_values_new[colIdx]->AssembleFromColumn(
						reinterpret_cast<SpecialiedColumn<uint32_t> *>(column_values_[finalIdx]), start, end, num_rows_, bw[finalIdx]);
				break;
			case 8:
				column_values_new[colIdx]->AssembleFromColumn(
						reinterpret_cast<SpecialiedColumn<uint64_t> *>(column_values_[finalIdx]), start, end, num_rows_, bw[finalIdx]);
				break;
			}

		}
	}

#if 0	//debug: output the cardinality of first column:
	std::set<uint64_t> cardInfo;
	for (uint32_t rowIdx = 0; rowIdx < 100; ++rowIdx) {
			cardInfo.insert((uint64_t)column_values_new[0]->GetValueAt(rowIdx));
			std::cout << std::hex << column_values_new[0]->GetValueAt(rowIdx) << std::endl;
	}
	std::cout << "cardinality of 1st column: " << cardInfo.size() << std::endl;

#endif
	//printf("arrive here2\n");
	timer_stitching.Stop();
	//setting_new.time_stitching = (double) timer_stitching.GetNumCycles();
	setting_new.time_stitching = (double) timer_stitching.GetNumCycles();	//TODO: TO REVISE HERE
	//setting_new.time_stitching = (double) timer_stitching.GetSeconds();	//TODO: TO REVISE HERE
	breakdown[0] = setting_new.time_stitching;

	/**
	 * Step c) feed new column values and other parameters to a chain composer (different from original stitching function)
	 **/
	setting_new.nthreads = setting_.nthreads;
	setting_new.intrinType = setting_.intrinType;
	setting_new.partition_fanout = setting_.partition_fanout;

	compose_params_t compose_params_new;
	compose_params_new.sortalgo = compose_params_.sortalgo;
	compose_params_new.packtype = compose_params_.packtype;
	compose_params_new.is_oid_encoded = compose_params_.is_oid_encoded;
	compose_params_new.colwidth_sum = compose_params_.colwidth_sum;
	compose_params_new.ordered = compose_params_.ordered;
	compose_params_new.hash_type = compose_params_.hash_type;

	compose_params_new.column_asc_desc = (bool*) malloc_aligned(new_col_num * sizeof(bool));
	compose_params_new.bitwidth = (uint32_t*) malloc_aligned(new_col_num * sizeof(uint32_t));

	uint32_t i;
	for (i = 0; i < new_col_num; ++i) {
		compose_params_new.column_asc_desc[i] = true;
		compose_params_new.bitwidth[i] = curSplit[i];
	}

	/** note: use the chain composer **/
	Composer* composer_new = new ChainComposer(setting_new, new_col_num, num_rows_,
			compose_params_new, column_values_new);

    HybridTimer totalTimer;
	totalTimer.Start();

	//printf("arrive here\n");
	composer_new->HashAllColumns();	//NOTE: this is ChainComposer's HashAllColumns, i.e., equivalent to multiRoundHashing

	totalTimer.Stop();

	setting_t updated = composer_new->getSetting();
	breakdown[1] = updated.time_tupleReconstruct;
	breakdown[2] = updated.time_hashReordering;
	breakdown[3] = updated.time_multipass_sort;
	//breakdown[4] = updated.time_mmcpyOrderGroup;//this cost is obselete
	breakdown[4] = updated.time_firstpass_sort;
	breakdown[5] = updated.time_extractGrop;
	assert(breakdown[0] > 0);
	assert(breakdown[1] > 0);
	//assert(breakdown[2] > 0);
	assert(breakdown[3] > 0);
	assert(breakdown[4] > 0);
	//assert(breakdown[5] > 0);

#if 0
	std::cout << "breakdown[1]: " << updated.time_tupleReconstruct << std::endl;
	std::cout << "breakdown[2]: " << updated.time_orderGroupInfo << std::endl;
	std::cout << "breakdown[3]: " << updated.time_multipass_sort << std::endl;
	std::cout << "breakdown[4]: " << updated.time_firstpass_sort << std::endl;
	std::cout << "breakdown[5]: " << updated.time_extractGrop << std::endl;
#endif

#if 0
	//output group information:
	std::cout << "[INFO ] 2nd round group num: " << composer_new->getGroupInfo().ngroups
			<< ", group num larger than 2: " << composer_new->getGroupInfo().ngroups_gtOne
			<< ", avg group size: "
			<< (double)composer_new->getGroupInfo().sum_group_size/(double)composer_new->getGroupInfo().ngroups
			<< std::endl;
#endif

	/**
	 * Step d) free the memory for current sumilation
	 */
	for (i = 0; i < composer_new->num_columns_; ++i) {
		switch(composer_new->column_values_[i]->GetWidth()) {
			case 1:
				free((uint8_t *)composer_new->column_values_[i]->GetColumn());
				break;
			case 2:
				free((uint16_t *)composer_new->column_values_[i]->GetColumn());
				break;
			case 4:
				free((uint32_t *)composer_new->column_values_[i]->GetColumn());
				break;
			case 8:
				free((uint64_t *)composer_new->column_values_[i]->GetColumn());
				break;
		}
		delete composer_new->column_values_[i];
	}
	free(composer_new->compose_params_.column_asc_desc);
	free(composer_new->compose_params_.bitwidth);

	free(composer_new->column_values_);
	free(composer_new);

	return totalTimer.GetNumCycles() + breakdown[0];
	//return totalTimer.GetSeconds() + breakdown[0];
}

#if 0
double StitchComposer::do_simulation_fake(uint32_t *columnOrder,
		std::vector<uint32_t> &curSplit, uint64_t* cardinality, double * breakdown, uint32_t &entry_num,
		uint32_t assembleStatistics[][3], uint32_t shiftStatistics[]) {

	/**
	 * construct new column values and feed them to a chain composer
	 * NOTE: use columnOrder and curSplit, copy most code from stitching() function
	 */
	setting_t setting_new = this->setting_;

	HybridTimer timer_stitching;
	timer_stitching.Start();

	/** original bitwidth **/
	uint32_t *bw = compose_params_.bitwidth;

	const uint32_t new_col_num = curSplit.size();

	HybridTimer timer_part1;
	timer_part1.Start();
	/**
	 * Step a) Find the shifting operations for each new column
	 */
	std::vector < construct_entry_t > constructRecord[new_col_num];
	uint32_t curOldColIdx = 0;
	uint32_t curNewColIdx = 0;
	uint32_t start = 1; //[start, end] is from the MSB to LSB, as the location pointer for current old column
	uint32_t minVal = 0;
	uint32_t curOldColWidth = bw[columnOrder[curOldColIdx]];	//important!!!
	uint32_t curNewColWidth = curSplit[curNewColIdx];
	while (curNewColIdx < new_col_num) {
		//curOldColWidth = bw[curOldColIdx];
		//curNewColWidth = bestGlobalSplit[curNewColIdx];
		minVal =
				(curOldColWidth < curNewColWidth) ?
						curOldColWidth : curNewColWidth;

		curOldColWidth -= minVal;
		curNewColWidth -= minVal;

		//insert a new entry
		construct_entry_t entry;
		entry.finalIdx = columnOrder[curOldColIdx];	//the index is transferred by bestOrdering
		entry.start = start;
		entry.end = start + minVal - 1;
		constructRecord[curNewColIdx].push_back(entry);
		entry_num++;
		start += minVal;

		/*********************update stitching/shifting statistics************************/
		if (curSplit[curNewColIdx] <= 16) {
			shiftStatistics[0]++;

			if (bw[columnOrder[curOldColIdx]] <= 16) {
				assembleStatistics[0][0]++;
			} else if (bw[columnOrder[curOldColIdx]] <= 32) {
				assembleStatistics[1][0]++;
			} else {
				assembleStatistics[2][0]++;
			}

		} else if (curSplit[curNewColIdx] <= 32) {
			shiftStatistics[1]++;

			if (bw[columnOrder[curOldColIdx]] <= 16) {
				assembleStatistics[0][1]++;
			} else if (bw[columnOrder[curOldColIdx]] <= 32) {
				assembleStatistics[1][1]++;
			} else {
				assembleStatistics[2][1]++;
			}
		} else {	//	<=64
			shiftStatistics[2]++;

			if (bw[columnOrder[curOldColIdx]] <= 16) {
				assembleStatistics[0][2]++;
			} else if (bw[columnOrder[curOldColIdx]] <= 32) {
				assembleStatistics[1][2]++;
			} else {
				assembleStatistics[2][2]++;
			}
		}
		/*********************************************************************************/

		if (0 == curOldColWidth) {	//an old column is consumed up
			curOldColIdx++;
			if (curOldColIdx < num_columns_) {
				curOldColWidth = bw[columnOrder[curOldColIdx]];
			}
			start = 1;
		}

		if (0 == curNewColWidth) {	//a new column is consumed up
			curNewColIdx++;
			if (curNewColIdx < new_col_num) {
				curNewColWidth = curSplit[curNewColIdx];
			}
		}
	}
	assert(curOldColIdx == num_columns_);
	assert(curNewColIdx == new_col_num);

	timer_part1.Stop();
	//std::cout << "Step1 cost: " << ((double)timer_part1.GetNumCycles())/num_rows_ << std::endl;

	HybridTimer timer_part2;
	timer_part2.Start();
	/**
	 * Step b) Re-construct each column according to ===>std::vector < construct_entry_t > constructRecord[new_col_num];
	 */

	/** Step b1) allocate the memories for all new columns: **/
	Column ** column_values_new = (Column **) malloc_aligned(new_col_num * sizeof(Column *));
	uint32_t idx;
	uint32_t nbits = 0;
	//std::cout << "#column: " << new_col_num << std::endl;
	for (idx = 0; idx < new_col_num; ++idx) {
		nbits = curSplit[idx];
		assert(nbits > 0 && nbits <= 64);

		//if (nbits <= 8) {
		//
		if (nbits <= 16) { /*IMPORTANT: use uint16_t to deal with values less than 8 bits*/
			uint16_t *values = (uint16_t *)malloc_aligned(num_rows_ * sizeof(uint16_t));
			memset(values, 0, num_rows_ * sizeof(uint16_t));
			SpecialiedColumn<uint16_t> *col = new SpecialiedColumn<uint16_t>();
			col->SetColumn(values);
			col->SetWidth(sizeof(uint16_t));
			column_values_new[idx] = col;
		} else if (nbits > 16 && nbits <= 32) {
			uint32_t *values = (uint32_t *)malloc_aligned(num_rows_ * sizeof(uint32_t));
			memset(values, 0, num_rows_ * sizeof(uint32_t));
			SpecialiedColumn<uint32_t> *col = new SpecialiedColumn<uint32_t>();
			col->SetColumn(values);
			col->SetWidth(sizeof(uint32_t));
			column_values_new[idx] = col;
		} else {	//width>32 && width<=64
			uint64_t *values = (uint64_t *)malloc_aligned(num_rows_ * sizeof(uint64_t));
			memset(values, 0, num_rows_ * sizeof(uint64_t));
			SpecialiedColumn<uint64_t> *col = new SpecialiedColumn<uint64_t>();
			col->SetColumn(values);
			col->SetWidth(sizeof(uint64_t));
			column_values_new[idx] = col;
		}
	}

	timer_part2.Stop();
	//std::cout << "Step2 cost: " << ((double)timer_part2.GetNumCycles())/num_rows_ << std::endl;

	HybridTimer timer_part3;
	timer_part3.Start();

	/** Step b2) re-construct the values **/
	/** Note about the column index issue, i.e., whether it is after transferred by columnOrder **/
	uint32_t colIdx, entryIdx;
	uint32_t finalIdx, end;	//start is already declared
	//uint64_t assemble_val = 0;
	uint32_t portion_width = 0;
	//uint32_t kPrefetchDistance = 0;

	std::vector < construct_entry_t > records;
	for (colIdx = 0; colIdx < new_col_num; ++colIdx) {
		records = constructRecord[colIdx];
		assert(!records.empty());
		for (entryIdx = 0; entryIdx < records.size(); ++entryIdx) {
			finalIdx = records[entryIdx].finalIdx;
			start = records[entryIdx].start;
			end = records[entryIdx].end;
			portion_width = end - start + 1;

			//std::cout << "finalIdx, start, end: " << finalIdx << ", " << start << ", " << end << std::endl;

			if (0 != entryIdx) {
				column_values_new[colIdx]->ShiftLeftBatch(num_rows_, portion_width);
			}
			switch (column_values_[finalIdx]->GetWidth()) {
			case 1:
				std::cout << "[Error ] Not Allow 8b column width" << std::endl;
				exit(EXIT_SUCCESS);
				break;
			case 2:
				column_values_new[colIdx]->AssembleFromColumn(
						reinterpret_cast<SpecialiedColumn<uint16_t> *>(column_values_[finalIdx]), start, end, num_rows_, bw[finalIdx]);
				break;
			case 4:
				column_values_new[colIdx]->AssembleFromColumn(
						reinterpret_cast<SpecialiedColumn<uint32_t> *>(column_values_[finalIdx]), start, end, num_rows_, bw[finalIdx]);
				break;
			case 8:
				column_values_new[colIdx]->AssembleFromColumn(
						reinterpret_cast<SpecialiedColumn<uint64_t> *>(column_values_[finalIdx]), start, end, num_rows_, bw[finalIdx]);
				break;
			}

		}
	}

	timer_part3.Stop();
	//std::cout << "Step3 cost: " << ((double)timer_part3.GetNumCycles())/num_rows_ << std::endl;

	timer_stitching.Stop();
	//setting_new.time_stitching = (double) timer_stitching.GetNumCycles();
	setting_new.time_stitching = (double) timer_part3.GetNumCycles();	//TODO: TO REVISE HERE
	breakdown[0] = setting_new.time_stitching;

	//cardinality extraction
	std::set<uint64_t> cardInfo[new_col_num - 1];
	for (colIdx = 0; colIdx < new_col_num - 1; ++colIdx) {
		for (uint32_t rowIdx = 0; rowIdx < num_rows_; ++rowIdx) {
				cardInfo[colIdx].insert((uint64_t)column_values_new[colIdx]->GetValueAt(rowIdx));
			}
	}

	for (idx = 0; idx < (new_col_num -1); ++idx) {
		cardinality[idx] = cardInfo[idx].size();
	}

	/**
	 * Step c) feed new column values and other parameters to a chain composer (different from original stitching function)
	 **/
	setting_new.nthreads = setting_.nthreads;
	setting_new.intrinType = setting_.intrinType;
	setting_new.partition_fanout = setting_.partition_fanout;

	compose_params_t compose_params_new;
	compose_params_new.sortalgo = compose_params_.sortalgo;
	compose_params_new.packtype = compose_params_.packtype;
	compose_params_new.is_oid_encoded = compose_params_.is_oid_encoded;
	compose_params_new.colwidth_sum = compose_params_.colwidth_sum;
	compose_params_new.ordered = compose_params_.ordered;

	compose_params_new.column_asc_desc = (bool*) malloc_aligned(new_col_num * sizeof(bool));
	compose_params_new.bitwidth = (uint32_t*) malloc_aligned(new_col_num * sizeof(uint32_t));

	uint32_t i;
	for (i = 0; i < new_col_num; ++i) {
		compose_params_new.column_asc_desc[i] = true;
		compose_params_new.bitwidth[i] = curSplit[i];
	}

	/** note: use the chain composer **/
	Composer* composer_new = new ChainComposer(setting_new, new_col_num, num_rows_,
			compose_params_new, column_values_new);

    HybridTimer totalTimer;
	totalTimer.Start();

	composer_new->SortAllColumns();

	totalTimer.Stop();

	setting_t updated = composer_new->getSetting();
	breakdown[1] = updated.time_tupleReconstruct;
	breakdown[2] = updated.time_orderGroupInfo;
	breakdown[3] = updated.time_multipass_sort;
	breakdown[4] = updated.time_mmcpyOrderGroup;
	breakdown[5] = updated.time_extractGrop;
	assert(breakdown[0] > 1);
	assert(breakdown[1] > 1);
	assert(breakdown[2] > 1);
	assert(breakdown[3] > 1);
	assert(breakdown[4] > 1);
	assert(breakdown[5] > 1);

	/**
	 * Step d) free the memory for current sumilation
	 */
	for (i = 0; i < composer_new->num_columns_; ++i) {
		switch(composer_new->column_values_[i]->GetWidth()) {
			case 1:
				free((uint8_t *)composer_new->column_values_[i]->GetColumn());
				break;
			case 2:
				free((uint16_t *)composer_new->column_values_[i]->GetColumn());
				break;
			case 4:
				free((uint32_t *)composer_new->column_values_[i]->GetColumn());
				break;
			case 8:
				free((uint64_t *)composer_new->column_values_[i]->GetColumn());
				break;
		}
		delete composer_new->column_values_[i];
	}
	free(composer_new->compose_params_.column_asc_desc);
	free(composer_new->compose_params_.bitwidth);

	free(composer_new->column_values_);
	free(composer_new);

	return totalTimer.GetNumCycles() + setting_new.time_stitching;
}
#endif

#if 0
void StitchComposer::stitching() {
	/**
	 * do the stitching based on <bestOrdering> and <bestGlobalSplit>
	 * IMPORTANT: the [bestGlobalSplit] is under the context the [bestOrdering] as the ordering of columns
	 */
	HybridTimer timer_stitching;
	timer_stitching.Start();

	/** original bitwidth **/
	uint32_t *bw = compose_params_.bitwidth;

	const uint32_t new_col_num = bestGlobalSplit.size();

	/**
	 * Step a) Find the shifting operations for each new column
	 */
	std::vector < construct_entry_t > constructRecord[new_col_num];
	uint32_t curOldColIdx = 0;
	uint32_t curNewColIdx = 0;
	uint32_t start = 1; //[start, end] is from the MSB to LSB, as the location pointer for current old column
	uint32_t minVal = 0;
	uint32_t curOldColWidth = bw[bestOrdering[curOldColIdx]];	//important!!!
	uint32_t curNewColWidth = bestGlobalSplit[curNewColIdx];
	while (curNewColIdx < new_col_num) {
		//curOldColWidth = bw[curOldColIdx];
		//curNewColWidth = bestGlobalSplit[curNewColIdx];
		minVal =
				(curOldColWidth < curNewColWidth) ?
						curOldColWidth : curNewColWidth;

		curOldColWidth -= minVal;
		curNewColWidth -= minVal;

		//insert a new entry
		construct_entry_t entry;
		entry.finalIdx = bestOrdering[curOldColIdx];	//the index is transferred by bestOrdering
		entry.start = start;
		entry.end = start + minVal - 1;
		constructRecord[curNewColIdx].push_back(entry);
		start += minVal;

		if (0 == curOldColWidth) {	//an old column is consumed up
			curOldColIdx++;
			if (curOldColIdx < num_columns_) {
				curOldColWidth = bw[bestOrdering[curOldColIdx]];
			}
			start = 1;
		}

		if (0 == curNewColWidth) {	//a new column is consumed up
			curNewColIdx++;
			if (curNewColIdx < new_col_num) {
				curNewColWidth = bestGlobalSplit[curNewColIdx];
			}
		}
	}
	assert(curOldColIdx == num_columns_);
	assert(curNewColIdx == new_col_num);

	/**
	 * Step b) Re-construct each column according to ===>std::vector < construct_entry_t > constructRecord[new_col_num];
	 */

	/** Step b1) allocate the memories for all new columns: **/
	Column ** column_values_new = (Column **) malloc_aligned(new_col_num * sizeof(Column *));
	uint32_t idx;
	uint32_t nbits = 0;
	for (idx = 0; idx < new_col_num; ++idx) {
		nbits = bestGlobalSplit[idx];
		assert(nbits > 0 && nbits <= 64);

		//if (nbits <= 8) {
		//
		if (nbits <= 16) { /*IMPORTANT: use uint16_t to deal with values less than 8 bits*/
			uint16_t *values = (uint16_t *)malloc_aligned(num_rows_ * sizeof(uint16_t));
			memset(values, 0, num_rows_ * sizeof(uint16_t));
			SpecialiedColumn<uint16_t> *col = new SpecialiedColumn<uint16_t>();
			col->SetColumn(values);
			col->SetWidth(sizeof(uint16_t));
			column_values_new[idx] = col;
		} else if (nbits > 16 && nbits <= 32) {
			uint32_t *values = (uint32_t *)malloc_aligned(num_rows_ * sizeof(uint32_t));
			memset(values, 0, num_rows_ * sizeof(uint32_t));
			SpecialiedColumn<uint32_t> *col = new SpecialiedColumn<uint32_t>();
			col->SetColumn(values);
			col->SetWidth(sizeof(uint32_t));
			column_values_new[idx] = col;
		} else {	//width>32 && width<=64
			uint64_t *values = (uint64_t *)malloc_aligned(num_rows_ * sizeof(uint64_t));
			memset(values, 0, num_rows_ * sizeof(uint64_t));
			SpecialiedColumn<uint64_t> *col = new SpecialiedColumn<uint64_t>();
			col->SetColumn(values);
			col->SetWidth(sizeof(uint64_t));
			column_values_new[idx] = col;
		}
	}

	/** Step b2) re-construct the values **/
	/** Note about the column index issue, i.e., whether it is after transferred by bestOrdering **/
	uint32_t colIdx, entryIdx;
	uint32_t finalIdx, end;	//start is already declared
	uint32_t portion_width = 0;

	std::vector < construct_entry_t > records;
	for (colIdx = 0; colIdx < new_col_num; ++colIdx) {
		records = constructRecord[colIdx];
		assert(!records.empty());
		for (entryIdx = 0; entryIdx < records.size(); ++entryIdx) {
			finalIdx = records[entryIdx].finalIdx;
			start = records[entryIdx].start;
			end = records[entryIdx].end;
			portion_width = end - start + 1;

			std::cout << "finalIdx, start, end: " << finalIdx << ", " << start << ", " << end << std::endl;

			if (0 != entryIdx) {
				column_values_new[colIdx]->ShiftLeftBatch(num_rows_, portion_width);
			}
			switch (column_values_[finalIdx]->GetWidth()) {
			case 1:
				std::cout << "[Error ] Not Allow 8b column width" << std::endl;
				exit(EXIT_SUCCESS);
				break;
			case 2:
				column_values_new[colIdx]->AssembleFromColumn(
						reinterpret_cast<SpecialiedColumn<uint16_t> *>(column_values_[finalIdx]), start, end, num_rows_, bw[finalIdx]);
				break;
			case 4:
				column_values_new[colIdx]->AssembleFromColumn(
						reinterpret_cast<SpecialiedColumn<uint32_t> *>(column_values_[finalIdx]), start, end, num_rows_, bw[finalIdx]);
				break;
			case 8:
				column_values_new[colIdx]->AssembleFromColumn(
						reinterpret_cast<SpecialiedColumn<uint64_t> *>(column_values_[finalIdx]), start, end, num_rows_, bw[finalIdx]);
				break;

			}

		}
	}

	/**
	 * Step c) delete original column values, and use the new one; also, updated all information of class <Composer>
	 */
	uint32_t i;
	for (i = 0; i < num_columns_; ++i) {
		switch(column_values_[i]->GetWidth()) {
			case 1:
				free((uint8_t *)column_values_[i]->GetColumn());
				break;
			case 2:
				free((uint16_t *)column_values_[i]->GetColumn());
				break;
			case 4:
				free((uint32_t *)column_values_[i]->GetColumn());
				break;
			case 8:
				free((uint64_t *)column_values_[i]->GetColumn());
				break;
		}
		delete column_values_[i];
	}
	free(column_values_);

	column_values_ = column_values_new;

	num_columns_ = new_col_num;

	free(compose_params_.column_asc_desc);
	free(compose_params_.bitwidth);

	compose_params_.column_asc_desc = (bool*) malloc_aligned(num_columns_ * sizeof(bool));
	compose_params_.bitwidth = (uint32_t*) malloc_aligned(num_columns_ * sizeof(uint32_t));

	for (i = 0; i < num_columns_; ++i) {
		compose_params_.column_asc_desc[i] = true;
		compose_params_.bitwidth[i] = bestGlobalSplit[i];
	}


#if 0
	uint32_t nbits = compose_params_.stitch_nbits;
	uint32_t *bw = compose_params_.bitwidth;
	uint64_t idx;
	/**
	 * note: first, deal with simple case with TWO columns
	 * case 1: 51+51, won't change the # of bytes
	 * case 2: 51+35, may change the bank size when shifting 3 bits left or 19 bit right
	 */
	uint64_t * col1 = (uint64_t *) column_values_[0]->GetColumn();
	uint64_t * col2 = (uint64_t *) column_values_[1]->GetColumn();
	if (0 == compose_params_.stitch_direction) { //stitch to left

		for (idx = 0; idx < num_rows_; ++idx) {
			col1[idx] = (col1[idx] << nbits) | (col2[idx] >> (bw[1] - nbits));
			col2[idx] = ((1ULL << (bw[1]-nbits)) - 1) & col2[idx];
		}

		if ((bw[1]-nbits) <= 32) { //need to change the bank size
			//allocate a new column type
			uint32_t *newColVal = (uint32_t *)malloc_aligned(num_rows_ * sizeof(uint32_t));
			for (idx = 0; idx < num_rows_; ++idx) {
				newColVal[idx] = (uint32_t) col2[idx];
			}

			//free original values
			free(col2);
			free(column_values_[1]);

			SpecialiedColumn<uint32_t> *newcol = new SpecialiedColumn<uint32_t>();
			newcol->SetColumn(newColVal);
			newcol->SetWidth(sizeof(uint32_t));
			column_values_[1] = newcol;
		}

		//change compose_params_t.bitwidth accordingly
		bw[0] += nbits;
		bw[1] -= nbits;

	} else {			//stitch to the right
		for (idx = 0; idx < num_rows_; ++idx) {
			col1[idx] = col1[idx] >> nbits;
			col2[idx] = ((((1ULL << nbits) - 1) & col1[idx]) << (bw[1])) | col2[idx];
		}

		if ((bw[0]-nbits) <= 32) { //need to change the bank size for the first column
			//allocate a new column type
			uint32_t *newColVal = (uint32_t *)malloc_aligned(num_rows_ * sizeof(uint32_t));
			for (idx = 0; idx < num_rows_; ++idx) {
				newColVal[idx] = (uint32_t) col1[idx];
			}

			//free original values
			free(col1);
			free(column_values_[0]);

			SpecialiedColumn<uint32_t> *newcol = new SpecialiedColumn<uint32_t>();
			newcol->SetColumn(newColVal);
			newcol->SetWidth(sizeof(uint32_t));
			column_values_[0] = newcol;
		}

		//change compose_params_t.bitwidth accordingly
		bw[0] -= nbits;
		bw[1] += nbits;
	}
#endif

#if 0
	/**
	 *
	 * case 1: 17+17, can be a 64-bank, can be 8-bank+32-bank (or reverse)
	 *
	 * Note: since bank-16 and bank-8 is not available, we still use bank-32 if the bitwidth is smaller than 16
	 */
	uint32_t * col1 = (uint32_t *) column_values_[0]->GetColumn();
	uint32_t * col2 = (uint32_t *) column_values_[1]->GetColumn();
	if (0 == compose_params_.stitch_direction) { //stitch to left

		if ((bw[0]+nbits) > 32) { //need to increase the bank size to 64 for the first column
			//allocate a new column type
			uint64_t *newColVal = (uint64_t *)malloc_aligned(num_rows_ * sizeof(uint64_t));
			for (idx = 0; idx < num_rows_; ++idx) {
				newColVal[idx] = (uint64_t) col1[idx];
			}

			//free original values
			free(col1);
			free(column_values_[0]);

			SpecialiedColumn<uint64_t> *newcol = new SpecialiedColumn<uint64_t>();
			newcol->SetColumn(newColVal);
			newcol->SetWidth(sizeof(uint64_t));
			column_values_[0] = newcol;

			for (idx = 0; idx < num_rows_; ++idx) {
				newColVal[idx] = (newColVal[idx] << nbits) | ((uint64_t)(col2[idx]) >> (bw[1] - nbits));
				col2[idx] = ((1ULL << (bw[1]-nbits)) - 1) & col2[idx];
			}

		} else {
			for (idx = 0; idx < num_rows_; ++idx) {
				col1[idx] = (col1[idx] << nbits) | (col2[idx] >> (bw[1] - nbits));
				col2[idx] = ((1ULL << (bw[1]-nbits)) - 1) & col2[idx];
			}
		}

		//change compose_params_t.bitwidth accordingly
		bw[0] += nbits;
		bw[1] -= nbits;

		//if the second column is unnecessary anymore
		assert((bw[0] >= 0) && (bw[1] >= 0));
		if (0 == bw[1]) {
			num_columns_--;
		}

	} else {			//stitch to the right

		if ((bw[1]+nbits) > 32) { //need to increase the bank size to 64 for the second column
			//allocate a new column type
			uint64_t *newColVal = (uint64_t *)malloc_aligned(num_rows_ * sizeof(uint64_t));
			for (idx = 0; idx < num_rows_; ++idx) {
				newColVal[idx] = (uint64_t) col2[idx];
			}

			//free original values
			free(col2);
			free(column_values_[1]);

			SpecialiedColumn<uint64_t> *newcol = new SpecialiedColumn<uint64_t>();
			newcol->SetColumn(newColVal);
			newcol->SetWidth(sizeof(uint64_t));
			column_values_[1] = newcol;

			for (idx = 0; idx < num_rows_; ++idx) {
				col1[idx] = col1[idx] >> nbits;
				newColVal[idx] = ((((1ULL << nbits) - 1) & (uint64_t)(col1[idx])) << (bw[1])) | newColVal[idx];
			}

		} else {
			for (idx = 0; idx < num_rows_; ++idx) {
				col1[idx] = col1[idx] >> nbits;
				col2[idx] = ((((1ULL << nbits) - 1) & col1[idx]) << (bw[1])) | col2[idx];
			}
		}

		//change compose_params_t.bitwidth accordingly
		bw[0] -= nbits;
		bw[1] += nbits;

		assert((bw[0] > 0) && (bw[1] > 0));//do not make a column blank by shifting right
	}
#endif

	timer_stitching.Stop();
	setting_.time_stitching += (double) timer_stitching.GetNumCycles();
}
#endif

void StitchComposer::enum_banksize_category(int curRound, const int round_num,
		std::vector<std::vector<uint32_t> > &cate_enum) {

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

bool StitchComposer::determineBitRange(const uint32_t banksize, int & bit_range_min,
		int & bit_range_max, const uint32_t bitOffset) {

	uint32_t sum_colWidth = compose_params_.colwidth_sum;
	assert(sum_colWidth > 0);

	switch(banksize) {
	case 0:	//16 banksize
		bit_range_min = bitOffset + 1;
		bit_range_max = bitOffset + 16;
		break;
	case 1: //32 banksize
		bit_range_min = bitOffset + 17;
		bit_range_max = bitOffset + 32;
		break;
	case 2:	//64 banksize
		bit_range_min = bitOffset + 33;
		bit_range_max = bitOffset + 64;
		break;
	}

	if (bit_range_min > (int)sum_colWidth) {	//no bit available for current banksize
		return false;
	} else {
		bit_range_max = min<int>(bit_range_max, (int)sum_colWidth);
		return true;
	}
}

void StitchComposer::determineBowlGroupsz(double & groupsz_bowl, int & bit_position_bowl,
		const uint32_t banksz_type, const double avg_groupsz, uint32_t * columnOrder){

	/**Step 1: estimate groupsz_bowl**/
	switch(banksz_type) {
	case 0:
		groupsz_bowl = calculate_bowl_int16(avg_groupsz);
		break;
	case 1:
		groupsz_bowl = calculate_bowl_int32(avg_groupsz);
		break;
	case 2:
		groupsz_bowl = calculate_bowl_int64(avg_groupsz);
		break;
	}

	/**Step 2: estimate bit_position_bowl according to groupsz_bowl**/
	uint32_t *bw = compose_params_.bitwidth;
#if 0
	const uint32_t origin_col_num = num_columns_;
	uint64_t base_col_card[origin_col_num];
	double card_per_bit[origin_col_num];

	/** calculate the bit info carried by each bit in each base column **/
	uint32_t colIdx;
	for (colIdx = 0; colIdx < origin_col_num; ++colIdx) {
		base_col_card[colIdx] = column_values_[columnOrder[colIdx]]->GetCardinality();
		card_per_bit[colIdx] = pow((double)(base_col_card[colIdx]), 1.0/(double)(bw[columnOrder[colIdx]]));
	}
#endif

	uint32_t colIdx = 0;
	double accuCardinality = 1.0;
	uint32_t bitOffset = 0;
	while (colIdx < num_columns_) {
		accuCardinality *= base_col_card[columnOrder[colIdx]];
		if (accuCardinality > groupsz_bowl) {
			break;
		}

		bitOffset += bw[columnOrder[colIdx]];
		colIdx++;
	}

	if (colIdx == num_columns_) {
		bit_position_bowl = bitOffset;
	} else {
		bit_position_bowl = bitOffset + bw[columnOrder[colIdx]] -
				ceil(log(accuCardinality/groupsz_bowl) / log(card_per_bit[columnOrder[colIdx]]));
	}

}

void StitchComposer::determineGroupszRange(const int bit_range_min, const int bit_range_max,
			double & groupsz_min, double & groupsz_max, uint32_t * columnOrder) {
#if 0
	/** redundant code: get the cardinality info carried per bit **/
	uint32_t *bw = compose_params_.bitwidth;

	const uint32_t origin_col_num = num_columns_;
	uint64_t base_col_card[origin_col_num];
	double card_per_bit[origin_col_num];

	/** calculate the bit info carried by each bit in each base column **/
	for (uint32_t colIdx = 0; colIdx < origin_col_num; ++colIdx) {
		base_col_card[colIdx] = column_values_[columnOrder[colIdx]]->GetCardinality();
		card_per_bit[colIdx] = pow((double)(base_col_card[colIdx]), 1.0/(double)(bw[columnOrder[colIdx]]));
	}
#endif
	/** determine groupsz_min and groupsz_max **/
	double group_num_min = 1.0;	//i.e., cardinality_min/max
	double group_num_max = 1.0;
	group_num_min = estimateGroupNum(bit_range_min, columnOrder);
	group_num_max = estimateGroupNum(bit_range_max, columnOrder);

	groupsz_min = (double)num_rows_ / group_num_max;
	groupsz_max = (double)num_rows_ / group_num_min;

}

double StitchComposer::estimateGroupNum(const uint32_t bit_position, uint32_t * columnOrder){
	uint32_t *bw = compose_params_.bitwidth;

	uint32_t colIdx = 0;
	uint32_t bits_avail = bit_position;
	uint32_t bits_in_each_basecol = 0;
	double group_num = 1.0;

	while (colIdx < num_columns_) {
		bits_in_each_basecol = min<uint32_t>(bits_avail, bw[columnOrder[colIdx]]);
		group_num *= pow(card_per_bit[columnOrder[colIdx]], bits_in_each_basecol);

		bits_avail -= bits_in_each_basecol;
		if (0 == bits_avail) {	//this is modified from (0 == bits_in_each_basecol)
			break;
		}

		//group_num *= shrink_factor;

		colIdx++;
	}
	return group_num;
}


#if 0
void StitchComposer::find_greedy_plans(uint32_t round_num,
		std::vector<std::vector<uint32_t> > &validPlans, uint32_t * columnOrder) {

	uint32_t sum_colWidth = compose_params_.colwidth_sum;
	assert(round_num > 2);	//no need for greedy under round_num=1 and round_num=2

	/** Step a) get all banksize categories under this round_num**/

	std::vector<std::vector<uint32_t> > cate_enum;

	enum_banksize_category(0, round_num, cate_enum);

	/** Step b) enumerate all banksize categories, greedily select one under each category**/
	std::vector<uint32_t> banksize_cate;
	std::vector<uint32_t> plan_cur_cate;
	for (uint32_t i = 0; i < cate_enum.size(); ++i) {
		banksize_cate = cate_enum[i];
		assert(banksize_cate.size() == round_num);

		plan_cur_cate.clear();

		/** enumerate through each round, and determine one plan greedily
		 * NOTE: no need to determine the last round, so *round_num-1* here **/
		uint32_t bitOffset = 0;			//how many bits has been determined for previous round
		//double accumulatedCardi = 1.0;	//what is the cardinality up to previous round
		int bit_range_min = -1;
		int bit_range_max = -1;
		bool isValid = true;
		double groupsz_min = -1.0;
		double groupsz_max = -1.0;

		double groupsz_bowl = -1.0;
		int bit_position_bowl = -1;

		int bit_position_opt = -1;

		double avg_groupsz = 0.0;

		for (uint32_t roundIdx = 0; roundIdx < round_num-1; ++roundIdx) {

			bit_range_min = -1;
			bit_range_max = -1;

			/* Step b1) validate if the banksize is valid; if true, return the bit range */
			isValid = determineBitRange(banksize_cate[roundIdx], bit_range_min, bit_range_max, bitOffset);
			if (!isValid) {
				break;
			}

			assert(bit_range_max > 0);
			assert(bit_range_min > 0);
			assert(bit_range_max >= bit_range_min);

			/** Step b2) determine the group size range (according to bit_range_min/max)
			 *  for the NEXT round sorting  **/
			groupsz_min = -1.0;
			groupsz_max = -1.0;

			determineGroupszRange(bit_range_min, bit_range_max, groupsz_min, groupsz_max, columnOrder);

			assert(groupsz_min > 0);
			assert(groupsz_max > 0);
			assert(groupsz_max >= groupsz_min);

			/** Step b3) Determine groupsz_bowl and bit_position_bowl by differentiate the cost equation for sorting **/
			groupsz_bowl = -1.0;
			bit_position_bowl = -1;

			avg_groupsz = (groupsz_min + groupsz_max)/2.0;

			determineBowlGroupsz(groupsz_bowl, bit_position_bowl,
					banksize_cate[roundIdx+1], avg_groupsz, columnOrder);

			assert(groupsz_bowl > 0);
			assert(bit_position_bowl > 0);
			assert(bit_position_bowl <= (int)sum_colWidth);

			/** Step b4) Get the best decision among (groupsz_min, groupsz_max, groupsz_opt)
			 * return corresponding bit_position_opt**/
			bit_position_opt = -1;
			//NOTE1: *bit_range_min* corresponds to *groupsz_max*!!
			//NOTE2: use banksz for NEXT round to determine the best option
			bit_position_opt = determineOptBitPosition(bit_range_min, bit_range_max, bit_position_bowl,
					groupsz_max, groupsz_min, groupsz_bowl,
					banksize_cate[roundIdx+1]);
			assert(bit_position_opt > 0);
			assert(bit_position_opt > (int)bitOffset);

			/** Step b5) Assign the decision, update bitOffset**/
			plan_cur_cate.push_back(bit_position_opt - bitOffset);
			bitOffset = bit_position_opt;
		}

		/** Step b6) Determine the last round bit num **/
		if (isValid) {
			assert(bitOffset <= sum_colWidth);

			uint32_t last_round_bit_num = sum_colWidth - bitOffset;
			switch(banksize_cate[round_num - 1]) {
			case 0:	//bank 16
				isValid = (last_round_bit_num > 0) && (last_round_bit_num <= 16);
				break;
			case 1:	//bank 32
				isValid = (last_round_bit_num > 16) && (last_round_bit_num <= 32);
				break;
			case 2:	//bank 64
				isValid = (last_round_bit_num > 32) && (last_round_bit_num <= 64);
				break;
			}

			if (isValid) {
				plan_cur_cate.push_back(last_round_bit_num);
				assert(plan_cur_cate.size() == round_num);
				validPlans.push_back(plan_cur_cate);
			}
		}
	}
}
#endif

#if 1
void StitchComposer::find_greedy_plans(uint32_t round_num,
		std::vector<std::vector<uint32_t> > &validPlans, uint32_t * columnOrder) {

	uint32_t sum_colWidth = compose_params_.colwidth_sum;
	assert(round_num > 2);	//no need for greedy under round_num=1 and round_num=2

	/** Step a) get all banksize categories under this round_num**/

	std::vector<std::vector<uint32_t> > cate_enum;

	enum_banksize_category(0, round_num, cate_enum);

	/** Step b) enumerate all banksize categories, greedily select one under each category**/
	std::vector<uint32_t> banksize_cate;
	std::vector<uint32_t> plan_cur_cate;
	for (uint32_t i = 0; i < cate_enum.size(); ++i) {
		banksize_cate = cate_enum[i];
		assert(banksize_cate.size() == round_num);

		plan_cur_cate.clear();

		/** enumerate through each round, and determine one plan greedily
		 * NOTE: no need to determine the last round, so *round_num-1* here **/
		uint32_t bitOffset = 0;			//how many bits has been determined for previous round
		//double accumulatedCardi = 1.0;	//what is the cardinality up to previous round
		int bit_range_min = -1;
		int bit_range_max = -1;
		bool isValid = true;
		//double groupsz_min = -1.0;
		//double groupsz_max = -1.0;

		//double groupsz_bowl = -1.0;
		//int bit_position_bowl = -1;

#if 0	//debug
		std::cout << "Banksize category: [";
		for (uint32_t i = 0; i < banksize_cate.size(); ++i) {
			std::cout << banksize_cate[i] << ", ";
		}
		std::cout << "]\n";

#endif

		int bit_position_opt = -1;

		//double avg_groupsz = 0.0;

		for (uint32_t roundIdx = 0; roundIdx < round_num-1; ++roundIdx) {

			bit_range_min = -1;
			bit_range_max = -1;

			/* Step b1) validate if the banksize is valid; if true, return the bit range */
			isValid = determineBitRange(banksize_cate[roundIdx], bit_range_min, bit_range_max, bitOffset);
			if (!isValid) {
				break;
			}

			assert(bit_range_max > 0);
			assert(bit_range_min > 0);
			assert(bit_range_max >= bit_range_min);

			/** Step b2) determine the group size range (according to bit_range_min/max)
			 *  for the NEXT round sorting  **/
			//groupsz_min = -1.0;
			//groupsz_max = -1.0;

			//determineGroupszRange(bit_range_min, bit_range_max, groupsz_min, groupsz_max, columnOrder);

			//assert(groupsz_min > 0);
			//assert(groupsz_max > 0);
			//assert(groupsz_max >= groupsz_min);

			/** Step b3) Determine groupsz_bowl and bit_position_bowl by differentiate the cost equation for sorting **/
			//groupsz_bowl = -1.0;
			//bit_position_bowl = -1;

			//avg_groupsz = (groupsz_min + groupsz_max)/2.0;

			//determineBowlGroupsz(groupsz_bowl, bit_position_bowl,
			//		banksize_cate[roundIdx+1], avg_groupsz, columnOrder);

			//assert(groupsz_bowl > 0);
			//assert(bit_position_bowl > 0);
			//assert(bit_position_bowl <= (int)sum_colWidth);

			/** Step b4) Get the best decision among (groupsz_min, groupsz_max, groupsz_opt)
			 * return corresponding bit_position_opt**/
			bit_position_opt = -1;
			//NOTE1: *bit_range_min* corresponds to *groupsz_max*!!
			//NOTE2: use banksz for NEXT round to determine the best option
			//bit_position_opt = determineOptBitPosition(bit_range_min, bit_range_max, bit_position_bowl,
			//		groupsz_max, groupsz_min, groupsz_bowl,
			//		banksize_cate[roundIdx+1]);
			bit_position_opt = determineOptBitPosition_v2(bit_range_min,
					bit_range_max, banksize_cate[roundIdx+1], columnOrder);

			assert(bit_position_opt > 0);
			assert(bit_position_opt > (int)bitOffset);

			/** Step b5) Assign the decision, update bitOffset**/
			plan_cur_cate.push_back(bit_position_opt - bitOffset);
			bitOffset = bit_position_opt;

#if 0
			std::cout << "round " << roundIdx << ": " << bitOffset << std::endl;
#endif
		}

		/** Step b6) Determine the last round bit num **/
		if (isValid) {
			assert(bitOffset <= sum_colWidth);

			uint32_t last_round_bit_num = sum_colWidth - bitOffset;
			switch(banksize_cate[round_num - 1]) {
			case 0:	//bank 16
				isValid = (last_round_bit_num > 0) && (last_round_bit_num <= 16);
				break;
			case 1:	//bank 32
				isValid = (last_round_bit_num > 16) && (last_round_bit_num <= 32);
				break;
			case 2:	//bank 64
				isValid = (last_round_bit_num > 32) && (last_round_bit_num <= 64);
				break;
			}

			if (isValid) {
				plan_cur_cate.push_back(last_round_bit_num);
				assert(plan_cur_cate.size() == round_num);
				validPlans.push_back(plan_cur_cate);
			}
		}
	}
}
#endif


int StitchComposer::determineOptBitPosition(const int bit_range_min, const int bit_range_max, const int bit_position_bowl,
					const double groupsz_max, const double groupsz_min, const double groupsz_bowl,
					uint32_t banksz_type) {

	double cost_groupsz_min = 0.0;
	double cost_groupsz_max = 0.0;
	double cost_groupsz_bowl = 0.0;

	cost_groupsz_min = calculate_single_column_cost(groupsz_min, banksz_type);
	cost_groupsz_max = calculate_single_column_cost(groupsz_max, banksz_type);
	if ((bit_position_bowl >= bit_range_min) && (bit_position_bowl <= bit_range_max)) {//see if bowl \in [min, max]
		cost_groupsz_bowl = calculate_single_column_cost(groupsz_bowl, banksz_type);
	} else {
		cost_groupsz_bowl = std::numeric_limits<double>::max();
	}

	/** output the best option **/
	if ((cost_groupsz_min <= cost_groupsz_max) && (cost_groupsz_min <= cost_groupsz_bowl)) {
		return bit_range_max;
	} else if ((cost_groupsz_max <= cost_groupsz_min) && (cost_groupsz_max <= cost_groupsz_bowl)) {
		return bit_range_min;
	} else {
		return bit_position_bowl;
	}
}

int StitchComposer::determineOptBitPosition_v2(const int bit_range_min, const int bit_range_max,
		uint32_t banksz_type, uint32_t * columnOrder) {

	double bestCost = std::numeric_limits<double>::max();
	int optBitPos = -1;

	int curBitPos = bit_range_min;
	double curGroupsz = 0.0;
	double curCost = 0.0;

	while (curBitPos <= bit_range_max) {
		curGroupsz = (double)num_rows_ / (estimateGroupNum(curBitPos, columnOrder));
		curCost = calculate_single_column_cost(curGroupsz, banksz_type);
#if 0	//debug
		std::cout << "curBitPos: " << curBitPos
				<< "; curGroupsz: " << curGroupsz << "; curCost: " << curCost << std::endl;
#endif

		if (curCost <= bestCost) {
			bestCost = curCost;
			optBitPos = curBitPos;
		}

		curBitPos++;
	}

	assert(optBitPos > 0);

	return optBitPos;
}

double StitchComposer::calculate_single_column_cost(const double groupsz, const uint32_t banksz_type) {

	assert(groupsz > 0);
	assert((banksz_type >= 0) && (banksz_type < 3));
	double group_num = ((double)num_rows_) / groupsz;
	double factor = 1.0;
	double alpha = 1.0;

	//assert(avg_group_sz >= 1.0);
#if 0		//the markov inequalities is unable to
	if (avg_group_sz < 2) {
		//percentage of groups with size larger than TWO, use Markov inequality: P(x < a) <= E[x]/a
		factor = avg_group_sz/2;

		//revise the avg group size
		avg_group_sz = ((double)num_rows_ - ((double)cardinality[colIdx])*(1-factor))
				/ (((double)cardinality[colIdx])*factor);

	}
#endif
	if ((groupsz - 1.0) < 0.000001) {
		factor = 0;
	} else {
		factor = 1.0/(1.0+pow(1.0/(groupsz-1.0), alpha));
	}

	switch(banksz_type) {
		case 0:	//banksize 16
			return group_num * factor * estimate_mergesort_uint16(groupsz);
			break;
		case 1: //banksize 32
			return group_num * factor * estimate_mergesort_uint32(groupsz);
			break;
		case 2: //banksize 64
			return group_num * factor * estimate_mergesort_uint64(groupsz);
			break;
		default:
			std::cout << "[ERROR ] No such banksz_type: " << banksz_type << std::endl;
			exit(EXIT_SUCCESS);
	}
}

void StitchComposer::estimate_plan_collection_greedy(uint32_t * columnOrder,
		std::vector<std::vector<uint32_t> > &greedy_plans,
		int & bestLocalSplitIdx, double & optimalLocCost){

#if 1
	optimalLocCost = std::numeric_limits<double>::max();

	/**
	 * Pay much attention to columnOrder!!!
	 */
	//std::cout << "[INFO ] Total number of plans: " << split_enum.size() << std::endl;
	std::vector<uint32_t> curSplit;

	/**aggregate Costing**/
	double estimatedTotalCost = 0.0;

	/** a) stitching/shifting overhead **/
	double stitch_cost = 0.0;

	/** b) tuple reconstruction overhead (i.e., BATproject function)**/
	double tuple_reconstruct_cost = 0.0;
	uint32_t byte_num_sum_remove_1st = 0;

	/** c) Group information extraction **/
	double order_group_cost = 0.0;

	/** d) multi-pass sorting**/
	double multi_pass_sort_cost = 0.0;
	double avg_group_sz = 0.0;

	uint32_t idx;
	uint32_t colIdx;
	//uint32_t i;

	for (uint32_t enum_idx = 0; enum_idx < greedy_plans.size(); ++enum_idx) {

		curSplit = greedy_plans[enum_idx];

		const uint32_t new_column_num = curSplit.size();

		//byte_num_sum = 0;
		byte_num_sum_remove_1st = 0;
		uint32_t bytewidths[new_column_num];
		for (idx = 0; idx < new_column_num; ++idx) {
			assert(curSplit[idx] > 0);
			assert(curSplit[idx] <= 64);
			if (curSplit[idx] <= 16) {

				byte_num_sum_remove_1st += 2;
				bytewidths[idx] = 2;
			} else if (curSplit[idx] <= 32) {

				byte_num_sum_remove_1st += 4;
				bytewidths[idx] = 4;
			} else {

				byte_num_sum_remove_1st += 8;
				bytewidths[idx] = 8;
			}

			if (idx == 0) {
				byte_num_sum_remove_1st = 0;
			}
		}

		//no need info for last column, but allocate it, becuase the factorization needs it
		double cardinality[new_column_num];
		for (idx = 0; idx < new_column_num; ++idx) {
			cardinality[idx] = 1;
		}

		//also need assemble statistics for costing stitching overhead
		uint32_t assembleStatistics[3][3] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

		uint32_t shiftStatistics[3] = {0, 0, 0};

		getCardinalityFromStatistics(columnOrder, curSplit, cardinality, assembleStatistics, shiftStatistics);

		assert(1 == cardinality[new_column_num - 1]);//make sure the last element did not change

		/** d) Get grouping (cardinality) information from simulation and do the costing for multi-pass sorting**/
		/* factorize the cardinality, i.e., x0 = a0, x1 = a0*a1, x2 = a0*a1*a2, ...*/
		double factor_accu = cardinality[0];
		double memo_element = 0;
		//double shrink = pow(shrink_factor, (double)(new_column_num-1));
		cardinality[0] = 1;	//the first column is ONE BIG group
		for (idx = 1; idx < new_column_num; ++idx) {
			memo_element = cardinality[idx];
			cardinality[idx] = factor_accu;
			factor_accu *= memo_element;
		}
		//assert(cardinality[new_column_num-1] == factor_accu);

		multi_pass_sort_cost = 0.0;
		avg_group_sz = 0.0;
		//double factor = 1.0;
		//double alpha = 1.0;
		for (colIdx = 0; colIdx < new_column_num; ++colIdx) {
			//assert(cardinality[colIdx] > 0);
			avg_group_sz = ((double)num_rows_) / ((double)cardinality[colIdx]);
			//std::cout << "Estimated group number for Round " << colIdx << " : " << cardinality[colIdx] << std::endl;
#if 1
			switch(bytewidths[colIdx]) {
				case 2:
					multi_pass_sort_cost += calculate_single_column_cost(avg_group_sz, 0);
 					break;
				case 4:
					multi_pass_sort_cost += calculate_single_column_cost(avg_group_sz, 1);
					break;
				case 8:
					multi_pass_sort_cost += calculate_single_column_cost(avg_group_sz, 2);
					break;
			}
#endif
		}

		/** a) stitching/shifting overhead **/
		stitch_cost = estimate_stitch(assembleStatistics, shiftStatistics, num_rows_);

		/** b) tuple reconstruction overhead (i.e., BATproject function)**/
		/** b1. random access. b2. sequential access over new created values **/
		//tuple_reconstruct_cost = ((new_column_num - 1) * num_rows_ +
		//		byte_num_sum_remove_1st * num_rows_ / CACHE_LINE_SIZE) * COST_PER_SEQ_ACCESS;
		tuple_reconstruct_cost = estimate_tuple_reconstruct(bytewidths, new_column_num, num_rows_);
		//tuple_reconstruct_cost = 0;

		/** c) group information extraction**/
		order_group_cost = 0;	//temporary

		/** aggregation for the total estimated cost**/
		estimatedTotalCost = stitch_cost + tuple_reconstruct_cost + order_group_cost + multi_pass_sort_cost;
#if 0
		uint32_t i;
		std::cout << "[INFO ] The current plan is [";
		for (i = 0; i < curSplit.size(); ++i) {
			std::cout << curSplit[i] << ", ";
		}
		std::cout << "], with cost: " << estimatedTotalCost << std::endl;

#endif

		assert(estimatedTotalCost > 0);

		if (estimatedTotalCost < optimalLocCost) {
			optimalLocCost = estimatedTotalCost;
			bestLocalSplitIdx = enum_idx;
		}
	}

#endif
}

void StitchComposer::estimate_plan_collection(uint32_t * columnOrder,
		std::vector<std::vector<uint32_t> > &split_enum, uint64_t startIndex,
		std::set<std::vector<uint32_t> > &greedyPlans,
		std::set<std::vector<uint32_t> > &randomPlans) {

#if 1
	/**
	 * Pay much attention to columnOrder!!!
	 */
	//std::cout << "[INFO ] Total number of plans: " << split_enum.size() << std::endl;
	std::vector<uint32_t> curSplit;

	/**aggregate Costing**/
	double estimatedTotalCost = 0.0;

	/** a) stitching/shifting overhead **/
	double stitch_cost = 0.0;

	/** b) tuple reconstruction overhead (i.e., BATproject function)**/
	double tuple_reconstruct_cost = 0.0;
	uint32_t byte_num_sum_remove_1st = 0;

	/** c) Group information extraction **/
	double order_group_cost = 0.0;

	/** d) multi-pass sorting**/
	double multi_pass_sort_cost = 0.0;
	double avg_group_sz = 0.0;

	uint32_t idx;
	uint32_t colIdx;
	//uint32_t i;

	for (uint32_t enum_idx = 0; enum_idx < split_enum.size(); ++enum_idx) {

		curSplit = split_enum[enum_idx];

		const uint32_t new_column_num = curSplit.size();

		//byte_num_sum = 0;
		byte_num_sum_remove_1st = 0;
		uint32_t bytewidths[new_column_num];
		for (idx = 0; idx < new_column_num; ++idx) {
			assert(curSplit[idx] > 0);
			assert(curSplit[idx] <= 64);
			if (curSplit[idx] <= 16) {

				byte_num_sum_remove_1st += 2;
				bytewidths[idx] = 2;
			} else if (curSplit[idx] <= 32) {

				byte_num_sum_remove_1st += 4;
				bytewidths[idx] = 4;
			} else {

				byte_num_sum_remove_1st += 8;
				bytewidths[idx] = 8;
			}

			if (idx == 0) {
				byte_num_sum_remove_1st = 0;
			}
		}

		/********************************************************************************
		 *
		 * In order to get precise grouping information, really do the stitching/shifting
		 *
		 ********************************************************************************/
#if 0
		uint32_t i;
		std::cout << "[INFO ] The current plan is [";
		for (i = 0; i < curSplit.size(); ++i) {
			std::cout << curSplit[i] << ", ";
		}
		std::cout << "] under the ordering <";
		for (i = 0; i < num_columns_; ++i) {
			std::cout << columnOrder[i] << ", ";
		}
		std::cout << ">\n";
#endif

		//no need info for last column, but allocate it, becuase the factorization needs it
		double cardinality[new_column_num];
		for (idx = 0; idx < new_column_num; ++idx) {
			cardinality[idx] = 1;
		}

		//also need assemble statistics for costing stitching overhead
		uint32_t assembleStatistics[3][3] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

		uint32_t shiftStatistics[3] = {0, 0, 0};

		getCardinalityFromStatistics(columnOrder, curSplit, cardinality, assembleStatistics, shiftStatistics);

		assert(1 == cardinality[new_column_num - 1]);//make sure the last element did not change

		/** d) Get grouping (cardinality) information from simulation and do the costing for multi-pass sorting**/
		/* factorize the cardinality, i.e., x0 = a0, x1 = a0*a1, x2 = a0*a1*a2, ...*/
		double factor_accu = cardinality[0];
		double memo_element = 0;
		cardinality[0] = 1;	//the first column is ONE BIG group
		for (idx = 1; idx < new_column_num; ++idx) {
			memo_element = cardinality[idx];
			cardinality[idx] = factor_accu;
			factor_accu *= memo_element;
		}
		assert(cardinality[new_column_num-1] == factor_accu);

		multi_pass_sort_cost = 0.0;
		avg_group_sz = 0.0;
		//double factor = 1.0;
		//double alpha = 1.0;
		for (colIdx = 0; colIdx < new_column_num; ++colIdx) {
			//assert(cardinality[colIdx] > 0);
			avg_group_sz = ((double)num_rows_) / ((double)cardinality[colIdx]);
#if 1
			switch(bytewidths[colIdx]) {
				case 2:
					multi_pass_sort_cost += calculate_single_column_cost(avg_group_sz, 0);
					break;
				case 4:
					multi_pass_sort_cost += calculate_single_column_cost(avg_group_sz, 1);
					break;
				case 8:
					multi_pass_sort_cost += calculate_single_column_cost(avg_group_sz, 2);
					break;
			}
#endif
#if 0
			factor = 1.0;

			//assert(avg_group_sz >= 1.0);
#if 0		//the markov inequalities is unable to
			if (avg_group_sz < 2) {
				//percentage of groups with size larger than TWO, use Markov inequality: P(x < a) <= E[x]/a
				factor = avg_group_sz/2;

				//revise the avg group size
				avg_group_sz = ((double)num_rows_ - ((double)cardinality[colIdx])*(1-factor))
						/ (((double)cardinality[colIdx])*factor);

			}
#endif
			if ((avg_group_sz - 1.0) < 0.000001) {
				factor = 0;
			} else {
				factor = 1.0/(1.0+pow(1.0/(avg_group_sz-1.0), alpha));
			}

			switch(bytewidths[colIdx]) {
				case 2:
					multi_pass_sort_cost += cardinality[colIdx] * factor * estimate_mergesort_uint16(avg_group_sz);
					break;
				case 4:
					multi_pass_sort_cost += cardinality[colIdx] * factor * estimate_mergesort_uint32(avg_group_sz);
					break;
				case 8:
					multi_pass_sort_cost += cardinality[colIdx] * factor * estimate_mergesort_uint64(avg_group_sz);
					break;
			}
#endif
		}

		/** a) stitching/shifting overhead **/
		stitch_cost = estimate_stitch(assembleStatistics, shiftStatistics, num_rows_);

		/** b) tuple reconstruction overhead (i.e., BATproject function)**/
		/** b1. random access. b2. sequential access over new created values **/
		//tuple_reconstruct_cost = ((new_column_num - 1) * num_rows_ +
		//		byte_num_sum_remove_1st * num_rows_ / CACHE_LINE_SIZE) * COST_PER_SEQ_ACCESS;
		tuple_reconstruct_cost = estimate_tuple_reconstruct(bytewidths, new_column_num, num_rows_);
		//tuple_reconstruct_cost = 0;

		/** c) group information extraction**/
		order_group_cost = 0;	//temporary

		/** aggregation for the total estimated cost**/
		estimatedTotalCost = stitch_cost + tuple_reconstruct_cost + order_group_cost + multi_pass_sort_cost;

		assert(estimatedTotalCost > 0);

		/*this line of information is used for plotting*/
		std::cout << (enum_idx+startIndex) << ", "
				<< estimatedTotalCost/((double)num_rows_) << ", "
				<< stitch_cost/((double)num_rows_) << ", "
				<< tuple_reconstruct_cost/((double)num_rows_) << ", "
				<< order_group_cost/((double)num_rows_) << ", "
				<< multi_pass_sort_cost/((double)num_rows_)
				<< std::endl;

		if ((new_column_num <= 2) || (greedyPlans.find(curSplit) != greedyPlans.end())) {
			std::cout << "GREEDY," << (enum_idx+startIndex) << ", "
					<< estimatedTotalCost/((double)num_rows_) << ", "
					<< stitch_cost/((double)num_rows_) << ", "
					<< tuple_reconstruct_cost/((double)num_rows_) << ", "
					<< order_group_cost/((double)num_rows_) << ", "
					<< multi_pass_sort_cost/((double)num_rows_)
					<< std::endl;
		}

		if (randomPlans.find(curSplit) != randomPlans.end()) {
			std::cout << "RANDOM," << (enum_idx+startIndex) << ", "
					<< estimatedTotalCost/((double)num_rows_) << ", "
					<< stitch_cost/((double)num_rows_) << ", "
					<< tuple_reconstruct_cost/((double)num_rows_) << ", "
					<< order_group_cost/((double)num_rows_) << ", "
					<< multi_pass_sort_cost/((double)num_rows_)
					<< std::endl;
		}

#if 0
		std::cout << enum_idx << ", " << estimatedTotalCost/((double)num_rows_)
				<< ", " << actualTotalCost/((double)num_rows_) <<  ", "
				<< ((isBaseline) ? 0.0 : (stitch_cost)/((double)num_rows_)) <<  ", "
				<< tuple_reconstruct_cost/((double)num_rows_) << ", "
				<< order_group_cost/((double)num_rows_) << ", "
				<< multi_pass_sort_cost/((double)num_rows_) << ", "
				<< ((isBaseline) ? 0.0 : (cost_breakdown[0]/((double)num_rows_)))<< ", "
				<< cost_breakdown[1]/((double)num_rows_) << ", "
				<< cost_breakdown[2]/((double)num_rows_) << ", "
				<< cost_breakdown[3]/((double)num_rows_) << ", [ "
				<< cost_breakdown[4]/((double)num_rows_) << ", "
				<< cost_breakdown[5]/((double)num_rows_) << ", ]"
				<< std::endl;
#endif
	}

#endif
}

#if 0
void StitchComposer::cost_estimation_plans(uint32_t * columnOrder, std::vector<std::vector<uint32_t> > &split_enum,
		int &optimalIdx, double &optimalCost, double& extractCardinalityNumCycles) {

	/**
	 * Pay much attention to columnOrder!!!
	 */
	//std::cout << "[INFO ] Total number of plans: " << split_enum.size() << std::endl;

	std::vector<uint32_t> curSplit;

	/**aggregate Costing**/
	double estimatedTotalCost = 0.0;

	/** a) stitching/shifting overhead **/
	double stitch_cost = 0.0;

	/** b) tuple reconstruction overhead (i.e., BATproject function)**/
	double tuple_reconstruct_cost = 0.0;
	uint32_t byte_num_sum_remove_1st = 0;

	/** c) Group information extraction **/
	double order_group_cost = 0.0;

	/** d) multi-pass sorting**/
	double multi_pass_sort_cost = 0.0;
	double avg_group_sz = 0.0;

	uint32_t idx;
	uint32_t colIdx;
	//uint32_t i;

	for (uint32_t enum_idx = 0; enum_idx < split_enum.size(); ++enum_idx) {

		curSplit = split_enum[enum_idx];

		const uint32_t new_column_num = curSplit.size();

		//byte_num_sum = 0;
		byte_num_sum_remove_1st = 0;
		//uint32_t * bytewidths = (uint32_t *)malloc_aligned(new_column_num * sizeof(uint32_t));
		uint32_t bytewidths[new_column_num];
		for (idx = 0; idx < new_column_num; ++idx) {
			assert(curSplit[idx] > 0);
			assert(curSplit[idx] <= 64);
			if (curSplit[idx] <= 16) {

				byte_num_sum_remove_1st += 2;
				bytewidths[idx] = 2;
			} else if (curSplit[idx] <= 32) {

				byte_num_sum_remove_1st += 4;
				bytewidths[idx] = 4;
			} else {

				byte_num_sum_remove_1st += 8;
				bytewidths[idx] = 8;
			}

			if (idx == 0) {
				byte_num_sum_remove_1st = 0;
			}
		}

		/********************************************************************************
		 *
		 * In order to get precise grouping information, really do the stitching/shifting
		 *
		 ********************************************************************************/
#if 0
		std::cout << "[INFO ] The current plan is [";
		for (i = 0; i < curSplit.size(); ++i) {
			std::cout << curSplit[i] << ", ";
		}
		std::cout << "] under the ordering <";
		for (i = 0; i < num_columns_; ++i) {
			std::cout << columnOrder[i] << ", ";
		}
		std::cout << ">\n";
#endif

		//no need info for last column, but allocate it, becuase the factorization needs it
		uint64_t cardinality[new_column_num];
		cardinality[new_column_num-1] = 1;

		//also need assemble statistics for costing stitching overhead
		uint32_t assembleStatistics[3][3] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

		uint32_t shiftStatistics[3] = {0, 0, 0};

		extractCardinalityNumCycles +=
				getCardinalityFromStatistics(columnOrder, curSplit, cardinality, assembleStatistics, shiftStatistics);

		assert(1 == cardinality[new_column_num - 1]);//make sure the last element did not change

		/** d) Get grouping (cardinality) information from simulation and do the costing for multi-pass sorting**/
		/* factorize the cardinality, i.e., x0 = a0, x1 = a0*a1, x2 = a0*a1*a2, ...*/
		uint64_t factor_accu = cardinality[0];
		uint64_t memo_element = 0;
		cardinality[0] = 1;	//the first column is ONE BIG group
		for (idx = 1; idx < new_column_num; ++idx) {
			memo_element = cardinality[idx];
			cardinality[idx] = factor_accu;
			factor_accu *= memo_element;
		}
		assert(cardinality[new_column_num-1] == factor_accu);

		multi_pass_sort_cost = 0.0;
		avg_group_sz = 0.0;
		double factor = 1.0;
		for (colIdx = 0; colIdx < new_column_num; ++colIdx) {
			assert(cardinality[colIdx] > 0);
			avg_group_sz = ((double)num_rows_) / ((double)cardinality[colIdx]);

			//assert(avg_group_sz >= 1.0);

			if (avg_group_sz < 2) {
				factor = avg_group_sz - 1.0;	//percentage of groups with size TWO
			}
			switch(bytewidths[colIdx]) {
				case 2:
					multi_pass_sort_cost += cardinality[colIdx] * factor * estimate_mergesort_uint16(avg_group_sz);
					break;
				case 4:
					multi_pass_sort_cost += cardinality[colIdx] * factor * estimate_mergesort_uint32(avg_group_sz);
					break;
				case 8:
					multi_pass_sort_cost += cardinality[colIdx] * factor * estimate_mergesort_uint64(avg_group_sz);
					break;
			}
		}

		/** a) stitching/shifting overhead **/
		stitch_cost = estimate_stitch(assembleStatistics, shiftStatistics, num_rows_);

		/** b) tuple reconstruction overhead (i.e., BATproject function)**/
		/** b1. random access. b2. sequential access over new created values **/
		//tuple_reconstruct_cost = ((new_column_num - 1) * num_rows_ +
		//		byte_num_sum_remove_1st * num_rows_ / CACHE_LINE_SIZE) * COST_PER_SEQ_ACCESS;
		tuple_reconstruct_cost = estimate_tuple_reconstruct(bytewidths, new_column_num, num_rows_);

		/** c) group information extraction**/
		order_group_cost = 0;	//temporary

		/** aggregation for the total estimated cost**/
		estimatedTotalCost = stitch_cost + tuple_reconstruct_cost + order_group_cost + multi_pass_sort_cost;

		assert(estimatedTotalCost > 0);

		/*this line of information is used for plotting*/
#if 0
		std::cout << enum_idx << ", " << estimatedTotalCost/((double)num_rows_)
				<< ", " << actualTotalCost/((double)num_rows_) <<  ", "
				<< ((isBaseline) ? 0.0 : (stitch_cost)/((double)num_rows_)) <<  ", "
				<< tuple_reconstruct_cost/((double)num_rows_) << ", "
				<< order_group_cost/((double)num_rows_) << ", "
				<< multi_pass_sort_cost/((double)num_rows_) << ", "
				<< ((isBaseline) ? 0.0 : (cost_breakdown[0]/((double)num_rows_)))<< ", "
				<< cost_breakdown[1]/((double)num_rows_) << ", "
				<< cost_breakdown[2]/((double)num_rows_) << ", "
				<< cost_breakdown[3]/((double)num_rows_) << ", [ "
				<< cost_breakdown[4]/((double)num_rows_) << ", "
				<< cost_breakdown[5]/((double)num_rows_) << ", ]"
				<< std::endl;
#endif

		/**update the optimalLocCost**/
		if (estimatedTotalCost < optimalCost) {
			optimalCost = estimatedTotalCost;
			optimalIdx = enum_idx;
		}
	}

	assert(optimalIdx >= 0);
	assert(optimalCost < std::numeric_limits<double>::max());


}
#endif

#if 0
double StitchComposer::estimate_tuple_reconstruct(const uint32_t *bytewidths, const uint32_t column_num, const uint64_t rows_num) {

	uint32_t COST_PER_SEQ_ACCESS_test = 4;
	uint32_t COST_PER_RAND_ACCESS_test = 16;

	double totalCost = 0.0;

	for (uint32_t colIdx = 1; colIdx < column_num; ++colIdx) {

		uint32_t seq_cache_num = (double)(rows_num * (sizeof(surrogate_t) + 2*bytewidths[colIdx])) / CACHE_LINE_SIZE;
		double seq_access_cost = seq_cache_num * COST_PER_SEQ_ACCESS_test;
		uint32_t L2_cache_block_num = rows_num * (sizeof(surrogate_t) + 2*bytewidths[colIdx]) / L2_CACHE_SIZE + 1;
		double rand_access_perc = (double)(L2_cache_block_num-1)/(double)L2_cache_block_num;
		double rand_access_cost = rows_num*COST_PER_RAND_ACCESS_test*rand_access_perc;
		totalCost += (seq_access_cost + rand_access_cost);
	}

	return totalCost;
}
#endif

double StitchComposer::estimate_tuple_reconstruct(const uint32_t *bytewidths, const uint32_t column_num, const uint64_t rows_num) {

	const double C_MEM_INT16 = 2.0;
	const double C_MEM_INT32 = 34.0;
	const double C_MEM_INT64 = 35.0;

	const double C_CACHE_INT16 = 2.0;
	const double C_CACHE_INT32 = 12.0;
	const double C_CACHE_INT64 = 11.0;

	double totalCost = 0.0;
	double hit_ratio = 0.0;

	for (uint32_t colIdx = 1; colIdx < column_num; ++colIdx) {

		hit_ratio = (double)L3_CACHE_SIZE / (double)(rows_num * bytewidths[colIdx]);
		if (hit_ratio > 1) {
			hit_ratio = 1.0;
		}

		switch(bytewidths[colIdx]) {
		case 2:
			totalCost += ((double)rows_num * (C_CACHE_INT16 * hit_ratio + C_MEM_INT16 * (1.0-hit_ratio)));
			break;
		case 4:
			totalCost += ((double)rows_num * (C_CACHE_INT32 * hit_ratio + C_MEM_INT32 * (1.0-hit_ratio)));
			break;
		case 8:
			totalCost += ((double)rows_num * (C_CACHE_INT64 * hit_ratio + C_MEM_INT64 * (1.0-hit_ratio)));
			break;
		}

	}

	return totalCost;
}

double StitchComposer::estimate_stitch(uint32_t assembleStatistics[][3], uint32_t shiftStatistics[], uint32_t nitems) {
	double leftShiftCost = 0.0;
	double assembleCost = 0.0;
	assert(nitems > 1);

	if (nitems <= 256)	//small case
	{
		leftShiftCost = (shiftStatistics[0] * SMALL_SHIFT_LEFT_COST_PER_ITEM[nitems-1][0] +
				shiftStatistics[1] * SMALL_SHIFT_LEFT_COST_PER_ITEM[nitems-1][1] +
				shiftStatistics[2] * SMALL_SHIFT_LEFT_COST_PER_ITEM[nitems-1][2]) * nitems;


		assembleCost = (assembleStatistics[0][0] * SMALL_ASSEMBLE_FROM_INT16_COST_PER_ITEM[nitems-1][0] +
				assembleStatistics[0][1] * SMALL_ASSEMBLE_FROM_INT16_COST_PER_ITEM[nitems-1][1] +
				assembleStatistics[0][2] * SMALL_ASSEMBLE_FROM_INT16_COST_PER_ITEM[nitems-1][2] +
				assembleStatistics[1][0] * SMALL_ASSEMBLE_FROM_INT32_COST_PER_ITEM[nitems-1][0] +
				assembleStatistics[1][1] * SMALL_ASSEMBLE_FROM_INT32_COST_PER_ITEM[nitems-1][1] +
				assembleStatistics[1][2] * SMALL_ASSEMBLE_FROM_INT32_COST_PER_ITEM[nitems-1][2] +
				assembleStatistics[2][0] * SMALL_ASSEMBLE_FROM_INT64_COST_PER_ITEM[nitems-1][0] +
				assembleStatistics[2][1] * SMALL_ASSEMBLE_FROM_INT64_COST_PER_ITEM[nitems-1][1] +
				assembleStatistics[2][2] * SMALL_ASSEMBLE_FROM_INT64_COST_PER_ITEM[nitems-1][2]) * nitems;

	} else {	//average case
		leftShiftCost = (shiftStatistics[0] * AVG_SHIFT_LEFT_COST_PER_ITEM[0] +
						shiftStatistics[1] * AVG_SHIFT_LEFT_COST_PER_ITEM[1] +
						shiftStatistics[2] * AVG_SHIFT_LEFT_COST_PER_ITEM[2]) * nitems;

		assembleCost = (assembleStatistics[0][0] * AVG_ASSEMBLE_FROM_INT16_COST_PER_ITEM[0] +
				assembleStatistics[0][1] * AVG_ASSEMBLE_FROM_INT16_COST_PER_ITEM[1] +
				assembleStatistics[0][2] * AVG_ASSEMBLE_FROM_INT16_COST_PER_ITEM[2] +
				assembleStatistics[1][0] * AVG_ASSEMBLE_FROM_INT32_COST_PER_ITEM[0] +
				assembleStatistics[1][1] * AVG_ASSEMBLE_FROM_INT32_COST_PER_ITEM[1] +
				assembleStatistics[1][2] * AVG_ASSEMBLE_FROM_INT32_COST_PER_ITEM[2] +
				assembleStatistics[2][0] * AVG_ASSEMBLE_FROM_INT64_COST_PER_ITEM[0] +
				assembleStatistics[2][1] * AVG_ASSEMBLE_FROM_INT64_COST_PER_ITEM[1] +
				assembleStatistics[2][2] * AVG_ASSEMBLE_FROM_INT64_COST_PER_ITEM[2]) * nitems;
	}

	return leftShiftCost + assembleCost;
}

void StitchComposer::run_plan_collection(uint32_t * columnOrder,
		std::vector<std::vector<uint32_t> > &split_enum, uint64_t startIndex) {

	std::vector<uint32_t> curSplit;

	uint32_t i;
	for (uint32_t enum_idx = 0; enum_idx < split_enum.size(); ++enum_idx) {

		curSplit = split_enum[enum_idx];
#if 1
		std::cout << "[INFO ] The current plan is [";
		for (i = 0; i < curSplit.size(); ++i) {
			std::cout << curSplit[i] << ", ";
		}
		std::cout << "] under the ordering <";
		for (i = 0; i < num_columns_; ++i) {
			std::cout << columnOrder[i] << ", ";
		}
		std::cout << ">\n";
#endif
		double cost_breakdown[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

		double totalCost = run_single_plan(columnOrder, curSplit, cost_breakdown);

		/*this line of information is used for plotting*/
		std::cout << (enum_idx+startIndex) << ", "
				<< totalCost/((double)num_rows_)<< ", "
				<< cost_breakdown[0]/((double)num_rows_)<< ", "
				<< cost_breakdown[1]/((double)num_rows_) << ", "
				<< cost_breakdown[2]/((double)num_rows_) << ", "
				<< cost_breakdown[3]/((double)num_rows_) << ", [ "
				<< cost_breakdown[4]/((double)num_rows_) << ", "
				<< cost_breakdown[5]/((double)num_rows_) << ", ]"
				<< std::endl;
	}
}

#if 0
void StitchComposer::cost_estimation(uint32_t * columnOrder,
		std::vector<std::vector<uint32_t> > &split_enum, int &bestLocalSplitIdx,
		double &optimalLocCost, int &bestLocalSplitIdx_ac, double &optimalLocCost_ac,
		double & baseline_cost, double & bestSplit_es2ac_each_enum) {

	/**
	 * Pay much attention to columnOrder!!!
	 */
	std::cout << "[INFO ] Total number of plans: " << split_enum.size() << std::endl;

	optimalLocCost = std::numeric_limits<double>::max();
	std::vector<uint32_t> curSplit;
	uint32_t new_column_num = 0;

	//auxiliary
	optimalLocCost_ac = std::numeric_limits<double>::max();

	/**aggregate Costing**/
	double actualTotalCost = 0.0;
	double estimatedTotalCost = 0.0;

	/** a) stitching/shifting overhead **/
	//double reallocate_cost = 0.0;
	//double assemble_cost = 0.0;
	//double free_cost = 0.0;
	//uint32_t byte_num_sum = 0;
	//uint32_t origin_byte_num_sum = 0;
	double stitch_cost = 0.0;

	/** b) tuple reconstruction overhead (i.e., BATproject function)**/
	double tuple_reconstruct_cost = 0.0;
	uint32_t byte_num_sum_remove_1st = 0;

	/** c) Group information extraction **/
	double order_group_cost = 0.0;

	/** d) multi-pass sorting**/
	double multi_pass_sort_cost = 0.0;
	double avg_group_sz = 0.0;

	uint32_t idx;
#if 0
	for (idx = 0; idx < num_columns_; ++idx) {
		if (compose_params_.bitwidth[idx] <= 8) {
			origin_byte_num_sum += 1;
		} else if (compose_params_.bitwidth[idx] <= 16) {
			origin_byte_num_sum += 2;
		} else if (compose_params_.bitwidth[idx] <= 32) {
			origin_byte_num_sum += 4;
		} else {
			origin_byte_num_sum += 8;
		}
	}
	assert(origin_byte_num_sum > 0);
#endif
	uint32_t colIdx;
	uint32_t i;
	bool isBaseline = false;
	baseline_cost = 0.0;

	for (uint32_t enum_idx = 0; enum_idx < split_enum.size(); ++enum_idx) {

		//enum_idx = 2;

		isBaseline = (enum_idx == (split_enum.size()-1));
		curSplit = split_enum[enum_idx];

		new_column_num = curSplit.size();

		//byte_num_sum = 0;
		byte_num_sum_remove_1st = 0;
		uint32_t * bytewidths = (uint32_t *)malloc_aligned(new_column_num * sizeof(uint32_t));
		for (idx = 0; idx < new_column_num; ++idx) {
			assert(curSplit[idx] > 0);
			assert(curSplit[idx] <= 64);
			if (curSplit[idx] <= 16) {
				//byte_num_sum += 2;
				byte_num_sum_remove_1st += 2;
				bytewidths[idx] = 2;
			} else if (curSplit[idx] <= 32) {
				//byte_num_sum += 4;
				byte_num_sum_remove_1st += 4;
				bytewidths[idx] = 4;
			} else {
				//byte_num_sum += 8;
				byte_num_sum_remove_1st += 8;
				bytewidths[idx] = 8;
			}

			if (idx == 0) {
				byte_num_sum_remove_1st = 0;
			}
		}
		//assert(byte_num_sum_remove_1st < byte_num_sum);

		/********************************************************************************
		 *
		 * In order to get precise grouping information, really do the stitching/shifting
		 *
		 ********************************************************************************/
#if 1
		std::cout << "[INFO ] The current plan is [";
		for (i = 0; i < curSplit.size(); ++i) {
			std::cout << curSplit[i] << ", ";
		}
		std::cout << "] under the ordering <";
		for (i = 0; i < num_columns_; ++i) {
			std::cout << columnOrder[i] << ", ";
		}
		std::cout << ">\n";
#endif

		//no need info for last column, but allocate it, becuase the factorization needs it
		uint64_t* cardinality = (uint64_t *)malloc_aligned(new_column_num * sizeof(uint64_t));
		cardinality[new_column_num-1] = 1;

		uint32_t assembleStatistics[3][3] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

		uint32_t shiftStatistics[3] = {0, 0, 0};

		double cost_breakdown[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		uint32_t entry_sum = 0;
		actualTotalCost = do_simulation_fake(columnOrder, curSplit, cardinality,
				cost_breakdown, entry_sum, assembleStatistics, shiftStatistics);
		assert(entry_sum > 0);
		assert(1 == cardinality[new_column_num - 1]);//make sure the last element did not change

		/**XU: if it is the baseline (the last plan), the stitching overhead is ZERO)**/
		if (isBaseline) {
			actualTotalCost -= cost_breakdown[0];
			baseline_cost = actualTotalCost;
		}

		/** d) Get grouping (cardinality) information from simulation and do the costing for multi-pass sorting**/
		/* factorize the cardinality, i.e., x0 = a0, x1 = a0*a1, x2 = a0*a1*a2, ...*/
		uint64_t factor_accu = cardinality[0];
		uint64_t memo_element = 0;
		cardinality[0] = 1;	//the first column is ONE BIG group
		for (idx = 1; idx < new_column_num; ++idx) {
			memo_element = cardinality[idx];
			cardinality[idx] = factor_accu;
			factor_accu *= memo_element;
		}
		assert(cardinality[new_column_num-1] == factor_accu);

		multi_pass_sort_cost = 0.0;
		avg_group_sz = 0.0;
		double factor = 1.0;
		for (colIdx = 0; colIdx < new_column_num; ++colIdx) {
			assert(cardinality[colIdx] > 0);
			avg_group_sz = ((double)num_rows_) / ((double)cardinality[colIdx]);

			if (avg_group_sz < 2) {
				factor = avg_group_sz - 1.0;	//percentage of groups with size TWO
			}
			switch(bytewidths[colIdx]) {
				case 2:
					multi_pass_sort_cost += cardinality[colIdx] * factor * estimate_mergesort_uint16(avg_group_sz);
					break;
				case 4:
					multi_pass_sort_cost += cardinality[colIdx] * factor * estimate_mergesort_uint32(avg_group_sz);
					break;
				case 8:
					multi_pass_sort_cost += cardinality[colIdx] * factor * estimate_mergesort_uint64(avg_group_sz);
					break;
			}
		}

		/** a) stitching/shifting overhead **/
#if 0
		reallocate_cost = byte_num_sum * num_rows_ * COST_ALLOC_PER_BYTE;
		//assemble_cost = ((byte_num_sum + origin_byte_num_sum) * num_rows_ / CACHE_LINE_SIZE) * COST_PER_SEQ_ACCESS;
		assemble_cost = num_rows_ * (entry_sum + new_column_num) * COST_PER_RAND_ACCESS;
		free_cost = origin_byte_num_sum * num_rows_ * COST_FREE_PER_BYTE;
#endif
		stitch_cost = estimate_stitch(assembleStatistics, shiftStatistics, num_rows_);;

		/** b) tuple reconstruction overhead (i.e., BATproject function)**/
		/** b1. random access. b2. sequential access over new created values **/
		//tuple_reconstruct_cost = ((new_column_num - 1) * num_rows_ +
		//		byte_num_sum_remove_1st * num_rows_ / CACHE_LINE_SIZE) * COST_PER_SEQ_ACCESS;
		tuple_reconstruct_cost = 0;

		/** c) group information extraction**/
		order_group_cost = 0;	//temporary

		/** aggregation for the total estimated cost**/
		if (isBaseline) {
			estimatedTotalCost = tuple_reconstruct_cost + order_group_cost + multi_pass_sort_cost;
		} else {
			estimatedTotalCost = stitch_cost + tuple_reconstruct_cost + order_group_cost + multi_pass_sort_cost;
		}

		assert(actualTotalCost > 0);
		assert(estimatedTotalCost > 0);

		/*this line of information is used for plotting*/
		std::cout << enum_idx << ", " << estimatedTotalCost/((double)num_rows_)
				<< ", " << actualTotalCost/((double)num_rows_) <<  ", "
				<< ((isBaseline) ? 0.0 : (stitch_cost)/((double)num_rows_)) <<  ", "
				<< tuple_reconstruct_cost/((double)num_rows_) << ", "
				<< order_group_cost/((double)num_rows_) << ", "
				<< multi_pass_sort_cost/((double)num_rows_) << ", "
				<< ((isBaseline) ? 0.0 : (cost_breakdown[0]/((double)num_rows_)))<< ", "
				<< cost_breakdown[1]/((double)num_rows_) << ", "
				<< cost_breakdown[2]/((double)num_rows_) << ", "
				<< cost_breakdown[3]/((double)num_rows_) << ", [ "
				<< cost_breakdown[4]/((double)num_rows_) << ", "
				<< cost_breakdown[5]/((double)num_rows_) << ", ]"
				<< std::endl;

		/**update the optimalLocCost**/
		if (estimatedTotalCost < optimalLocCost) {
			optimalLocCost = estimatedTotalCost;
			bestLocalSplitIdx = enum_idx;

			bestSplit_es2ac_each_enum = actualTotalCost;
		}

		/**auxialiary: upate the actual cost**/
		if (actualTotalCost < optimalLocCost_ac) {
			optimalLocCost_ac = actualTotalCost;
			bestLocalSplitIdx_ac = enum_idx;
		}

		free(cardinality);
		free(bytewidths);
		//usleep(1000*2);
	}

	assert(bestLocalSplitIdx >= 0);
	assert(optimalLocCost < std::numeric_limits<double>::max());

	assert(bestLocalSplitIdx_ac >= 0);
	assert(optimalLocCost_ac < std::numeric_limits<double>::max());
}
#endif

#if 0
void StitchComposer::cost_estimation_profiling(uint32_t * columnOrder,
		std::vector<std::vector<uint32_t> > &split_enum, int &bestLocalSplitIdx,
		double &optimalLocCost, int &bestLocalSplitIdx_ac, double &optimalLocCost_ac,
		double & baseline_cost, double & bestSplit_es2ac_each_enum) {

	/**
	 * Pay much attention to columnOrder!!!
	 */
	//std::cout << "[INFO ] Total number of plans: " << split_enum.size() << std::endl;

	optimalLocCost = std::numeric_limits<double>::max();
	std::vector<uint32_t> curSplit;
	uint32_t new_column_num = 0;

	//auxiliary
	optimalLocCost_ac = std::numeric_limits<double>::max();

	/**aggregate Costing**/
	double actualTotalCost = 0.0;
	double estimatedTotalCost = 0.0;

	/** a) stitching/shifting overhead **/
	//double reallocate_cost = 0.0;
	//double assemble_cost = 0.0;
	//double free_cost = 0.0;
	//uint32_t byte_num_sum = 0;
	//uint32_t origin_byte_num_sum = 0;
	double stitch_cost = 0.0;

	/** b) tuple reconstruction overhead (i.e., BATproject function)**/
	double tuple_reconstruct_cost = 0.0;
	uint32_t byte_num_sum_remove_1st = 0;

	/** c) Group information extraction **/
	double order_group_cost = 0.0;

	/** d) multi-pass sorting**/
	double multi_pass_sort_cost = 0.0;
	double avg_group_sz = 0.0;

	uint32_t idx;
#if 0
	for (idx = 0; idx < num_columns_; ++idx) {
		if (compose_params_.bitwidth[idx] <= 8) {
			origin_byte_num_sum += 1;
		} else if (compose_params_.bitwidth[idx] <= 16) {
			origin_byte_num_sum += 2;
		} else if (compose_params_.bitwidth[idx] <= 32) {
			origin_byte_num_sum += 4;
		} else {
			origin_byte_num_sum += 8;
		}
	}
	assert(origin_byte_num_sum > 0);
#endif
	uint32_t colIdx;
	uint32_t i;
	bool isBaseline = false;
	baseline_cost = 0.0;

	for (uint32_t enum_idx = 0; enum_idx < split_enum.size(); ++enum_idx) {

		//enum_idx = 2;

		isBaseline = (enum_idx == (split_enum.size()-1));
		curSplit = split_enum[enum_idx];

		new_column_num = curSplit.size();

		//byte_num_sum = 0;
		byte_num_sum_remove_1st = 0;
		uint32_t * bytewidths = (uint32_t *)malloc_aligned(new_column_num * sizeof(uint32_t));
		for (idx = 0; idx < new_column_num; ++idx) {
			assert(curSplit[idx] > 0);
			assert(curSplit[idx] <= 64);
			if (curSplit[idx] <= 16) {
				//byte_num_sum += 2;
				byte_num_sum_remove_1st += 2;
				bytewidths[idx] = 2;
			} else if (curSplit[idx] <= 32) {
				//byte_num_sum += 4;
				byte_num_sum_remove_1st += 4;
				bytewidths[idx] = 4;
			} else {
				//byte_num_sum += 8;
				byte_num_sum_remove_1st += 8;
				bytewidths[idx] = 8;
			}

			if (idx == 0) {
				byte_num_sum_remove_1st = 0;
			}
		}
		//assert(byte_num_sum_remove_1st < byte_num_sum);

		/********************************************************************************
		 *
		 * In order to get precise grouping information, really do the stitching/shifting
		 *
		 ********************************************************************************/
#if 0
		std::cout << "[INFO ] The current plan is [";
		for (i = 0; i < curSplit.size(); ++i) {
			std::cout << curSplit[i] << ", ";
		}
		std::cout << "] under the ordering <";
		for (i = 0; i < num_columns_; ++i) {
			std::cout << columnOrder[i] << ", ";
		}
		std::cout << ">\n";
#endif

		//no need info for last column, but allocate it, becuase the factorization needs it
		uint64_t* cardinality = (uint64_t *)malloc_aligned(new_column_num * sizeof(uint64_t));
		cardinality[new_column_num-1] = 1;

		double cost_breakdown[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		//uint32_t entry_sum = 0;
		/***************************************************************************************/
		//actualTotalCost = do_simulation_fake(columnOrder, curSplit, cardinality, cost_breakdown, entry_sum);
		//============>
		//entry_sum = 1;
		for (i = 0; i < new_column_num-1; ++i) {
				cardinality[i] = column_values_[min<uint32_t>(i, num_columns_ - 1)]->GetCardinality();
		}
		/***************************************************************************************/
		//assert(entry_sum > 0);
		assert(1 == cardinality[new_column_num - 1]);//make sure the last element did not change

		/**XU: if it is the baseline (the last plan), the stitching overhead is ZERO)**/
		if (isBaseline) {
			actualTotalCost -= cost_breakdown[0];
			baseline_cost = actualTotalCost;
		}

		/** d) Get grouping (cardinality) information from simulation and do the costing for multi-pass sorting**/
		/* factorize the cardinality, i.e., x0 = a0, x1 = a0*a1, x2 = a0*a1*a2, ...*/
		uint64_t factor_accu = cardinality[0];
		uint64_t memo_element = 0;
		cardinality[0] = 1;	//the first column is ONE BIG group
		for (idx = 1; idx < new_column_num; ++idx) {
			memo_element = cardinality[idx];
			cardinality[idx] = factor_accu;
			factor_accu *= memo_element;
		}
		assert(cardinality[new_column_num-1] == factor_accu);

		multi_pass_sort_cost = 0.0;
		avg_group_sz = 0.0;
		double factor = 1.0;
		for (colIdx = 0; colIdx < new_column_num; ++colIdx) {
			assert(cardinality[colIdx] > 0);
			avg_group_sz = ((double)num_rows_) / ((double)cardinality[colIdx]);

			if (avg_group_sz < 2) {
				factor = 1.0 - 1.0/avg_group_sz;
			}
			switch(bytewidths[colIdx]) {
				case 2:
					multi_pass_sort_cost += cardinality[colIdx] * factor * estimate_mergesort_uint16(avg_group_sz);
					break;
				case 4:
					multi_pass_sort_cost += cardinality[colIdx] * factor * estimate_mergesort_uint32(avg_group_sz);
					break;
				case 8:
					multi_pass_sort_cost += cardinality[colIdx] * factor * estimate_mergesort_uint64(avg_group_sz);
					break;
			}
		}

		/** a) stitching/shifting overhead **/
#if 0
		reallocate_cost = byte_num_sum * num_rows_ * COST_ALLOC_PER_BYTE;
		//assemble_cost = ((byte_num_sum + origin_byte_num_sum) * num_rows_ / CACHE_LINE_SIZE) * COST_PER_SEQ_ACCESS;
		assemble_cost = num_rows_ * (entry_sum + new_column_num) * COST_PER_RAND_ACCESS;
		free_cost = origin_byte_num_sum * num_rows_ * COST_FREE_PER_BYTE;
#endif
		stitch_cost = num_rows_ * STITCH_COST_PER_ITEM;

		/** b) tuple reconstruction overhead (i.e., BATproject function)**/
		/** b1. random access. b2. sequential access over new created values **/
		tuple_reconstruct_cost = ((new_column_num - 1) * num_rows_ +
				byte_num_sum_remove_1st * num_rows_ / CACHE_LINE_SIZE) * COST_PER_SEQ_ACCESS;

		/** c) group information extraction**/
		order_group_cost = tuple_reconstruct_cost / 2;	//temporary


		/** aggregation for the total estimated cost**/
		if (isBaseline) {
			estimatedTotalCost = tuple_reconstruct_cost + order_group_cost + multi_pass_sort_cost;
		} else {
			estimatedTotalCost = stitch_cost + tuple_reconstruct_cost + order_group_cost + multi_pass_sort_cost;
		}

		assert(actualTotalCost > 0);
		assert(estimatedTotalCost > 0);

		/*this line of information is used for plotting*/
#if 0
		std::cout << enum_idx << ", " << estimatedTotalCost/((double)num_rows_)
				<< ", " << actualTotalCost/((double)num_rows_) <<  ", "
				<< ((isBaseline) ? 0.0 : (stitch_cost)/((double)num_rows_)) <<  ", "
				<< tuple_reconstruct_cost/((double)num_rows_) << ", "
				<< order_group_cost/((double)num_rows_) << ", "
				<< multi_pass_sort_cost/((double)num_rows_) << ", "
				<< ((isBaseline) ? 0.0 : (cost_breakdown[0]/((double)num_rows_)))<< ", "
				<< cost_breakdown[1]/((double)num_rows_) << ", "
				<< cost_breakdown[2]/((double)num_rows_) << ", "
				<< cost_breakdown[3]/((double)num_rows_) << ", [ "
				<< cost_breakdown[4]/((double)num_rows_) << ", "
				<< cost_breakdown[5]/((double)num_rows_) << ", ]"
				<< std::endl;
#endif

		/**update the optimalLocCost**/
		if (estimatedTotalCost < optimalLocCost) {
			optimalLocCost = estimatedTotalCost;
			bestLocalSplitIdx = enum_idx;

			bestSplit_es2ac_each_enum = actualTotalCost;
		}

		/**auxialiary: upate the actual cost**/
		if (actualTotalCost < optimalLocCost_ac) {
			optimalLocCost_ac = actualTotalCost;
			bestLocalSplitIdx_ac = enum_idx;
		}

		free(cardinality);
		free(bytewidths);
		//usleep(1000*2);
	}

	assert(bestLocalSplitIdx >= 0);
	assert(optimalLocCost < std::numeric_limits<double>::max());

	assert(bestLocalSplitIdx_ac >= 0);
	assert(optimalLocCost_ac < std::numeric_limits<double>::max());
}
#endif


#if 0
double StitchComposer::estimate_mergesort_uint64(double nitems) {

	if (nitems > BLOCKSIZE<uint64_t>()) {	//out of cache data volume
		return C1_incache_int64 +
				C2_incache_int64 * nitems +
				C3_incache_int64 * nitems * ceil(log(nitems)/log(2));
	} else {
		return 0;
	}

}

double StitchComposer::estimate_mergesort_uint32(double nitems) {
	if (nitems > BLOCKSIZE<uint32_t>()) {	//out of cache data volume
		return C1_incache_int32 +
				C2_incache_int32 * nitems +
				C3_incache_int32 * nitems * ceil(log(nitems)/log(2));
	} else {
		return 0;
	}
}


double StitchComposer::estimate_mergesort_uint16(double nitems) {
	if (nitems > BLOCKSIZE<uint16_t>()) {	//out of cache data volume
		return C1_incache_int16 +
				C2_incache_int16 * nitems +
				C3_incache_int16 * nitems * ceil(log(nitems)/log(2));
	} else {
		return 0;
	}
}
#endif


#if 1
double StitchComposer::estimate_mergesort_uint64(double nitems) {
#if 0
	if (nitems < 1.001) {
		return 0;
#endif
	if (nitems <= 2) {
		return STL_SORT_COST_INT64_PER_NLOGN[1] * 2 * log(2);
	} else if (nitems <= 512){
		uint64_t roundedItemNum = (uint64_t)(round(nitems));
		return STL_SORT_COST_INT64_PER_NLOGN[roundedItemNum - 1] * roundedItemNum * log(roundedItemNum);
	} else {
		uint64_t roundedItemNum = (uint64_t)(round(nitems));
		return AVG_SORT_COST_INT64_PER_NLOGN * roundedItemNum * log(roundedItemNum);
	}
}

double StitchComposer::estimate_mergesort_uint32(double nitems) {
#if 0
	if (nitems < 1.001) {
		return 0;
#endif
	if (nitems <= 2) {
		return STL_SORT_COST_INT32_PER_NLOGN[1] * 2 * log(2);
	} else if (nitems <= 256){
		uint64_t roundedItemNum = (uint64_t)(round(nitems));
		return STL_SORT_COST_INT32_PER_NLOGN[roundedItemNum - 1] * roundedItemNum * log(roundedItemNum);
	} else {
		uint64_t roundedItemNum = (uint64_t)(round(nitems));
		return AVG_SORT_COST_INT32_PER_NLOGN * roundedItemNum * log(roundedItemNum);
	}
}


double StitchComposer::estimate_mergesort_uint16(double nitems) {
#if 0
	if (nitems < 1.001) {
		return 0;
#endif
	if (nitems <= 2) {
		return STL_SORT_COST_INT16_PER_NLOGN[1] * 2 * log(2);
	} else if (nitems <= 256){
		uint64_t roundedItemNum = (uint64_t)(round(nitems));
		return STL_SORT_COST_INT16_PER_NLOGN[roundedItemNum - 1] * roundedItemNum * log(roundedItemNum);
	} else {
		uint64_t roundedItemNum = (uint64_t)(round(nitems));
		return AVG_SORT_COST_INT16_PER_NLOGN * roundedItemNum * log(roundedItemNum);
	}
}
#endif

#if 0
double StitchComposer::estimate_mergesort_uint64(double nitems) {

	constexpr double ln2 = 0.693147;

	//fix cost for function calls and memory allocation
	constexpr double fix_cost = SORT_FIX_COST;

	constexpr double C1 = log(32.0/L2_CACHE_SIZE)/ln2;

	double n2waymerge_rounds = C1 + log(nitems)/ln2;

	constexpr double C2 = 48.0 / CACHE_LINE_SIZE  * COST_PER_SEQ_ACCESS;

	double mm_cpy_cost = C2 * nitems;

	constexpr double C3 = 24.0 / CACHE_LINE_SIZE * COST_PER_SEQ_ACCESS;

	double scan1_cost = C3 * nitems;

	constexpr double C4 = C3;

	double scan2_cost = n2waymerge_rounds * C4 * nitems;

	constexpr double C5 = SORT_NETWORK_COST_PER_ITEM_64;

	double inregister_sort_cost = C5 * nitems;

	constexpr double C6 = MERGE_VALEN_COST_PER_ITEM_64;

	double outofcache_merge_cost = n2waymerge_rounds * C6 * nitems;

	constexpr double C7 = log(BLOCKSIZE<uint64_t>()/4.0)/ln2  * MERGE_EQLEN_COST_PER_ITEM_64;

	double incache_merge_cost_for_outofcachesize = C7 * nitems;

	constexpr double C8 = MERGE_EQLEN_COST_PER_ITEM_64;

	double incache_merge_cost_for_incachesize = C8 * nitems * log(nitems)/ln2;

	if (nitems > BLOCKSIZE<uint64_t>()) {
		return fix_cost + mm_cpy_cost + scan1_cost + scan2_cost + inregister_sort_cost +
				outofcache_merge_cost + incache_merge_cost_for_outofcachesize;
	} else {
		return fix_cost + mm_cpy_cost + inregister_sort_cost + incache_merge_cost_for_incachesize;
	}

#if 0
	//double total_cost = 0.0;
	double value_width = 8.0;
	double oid_width = 4.0;
	double ln2 = 0.693147;

	if (nitems <= 256) {
		return SORT_FIX_COST/nitems;
	} else {
		double n2waymerge_rounds = log(ceil(4.0*nitems*max<double>(value_width, oid_width)/L2_CACHE_SIZE))/ln2;

		/**a) memory copy before and after sorting**/
		double mm_cpy_cost = nitems * (value_width + oid_width) / CACHE_LINE_SIZE  * 4.0 * COST_PER_SEQ_ACCESS;

		/**b) scan sorting blocks and 2-way merge blocks**/
		double sort_scan_cost = (2.0 + 2.0 * n2waymerge_rounds) *
				nitems * (value_width + oid_width) / CACHE_LINE_SIZE * COST_PER_SEQ_ACCESS;

		double inregister_sort_cost = nitems * SORT_NETWORK_COST_PER_ITEM_64;

		if (nitems > BLOCKSIZE<uint64_t>()) {	//depends on the data volume is larger than one block or not
			/****/
			return SORT_FIX_COST + mm_cpy_cost + sort_scan_cost + inregister_sort_cost
					+ log(BLOCKSIZE<uint64_t>()/4.0)/ln2 * nitems * MERGE_EQLEN_COST_PER_ITEM_64
					+ n2waymerge_rounds * nitems * MERGE_VALEN_COST_PER_ITEM_64;
		} else {	//smaller than a block, just perform in-block sorting
			/****/
			return SORT_FIX_COST + mm_cpy_cost + inregister_sort_cost + log(nitems)/ln2 * nitems * MERGE_EQLEN_COST_PER_ITEM_64;
		}
	}
#endif
}

double StitchComposer::estimate_mergesort_uint32(double nitems) {

	constexpr double ln2 = 0.693147;

	//fix cost for function calls and memory allocation
	constexpr double fix_cost = SORT_FIX_COST;

	constexpr double C1 = log(16.0/L2_CACHE_SIZE)/ln2;

	double n2waymerge_rounds = C1 + log(nitems)/ln2;

	constexpr double C2 = 32.0 / CACHE_LINE_SIZE  * COST_PER_SEQ_ACCESS;

	double mm_cpy_cost = C2 * nitems;

	constexpr double C3 = 16.0 / CACHE_LINE_SIZE * COST_PER_SEQ_ACCESS;

	double scan1_cost = C3 * nitems;

	constexpr double C4 = C3;

	double scan2_cost = n2waymerge_rounds * C4 * nitems;

	constexpr double C5 = SORT_NETWORK_COST_PER_ITEM_32;

	double inregister_sort_cost = C5 * nitems;

	constexpr double C6 = MERGE_VALEN_COST_PER_ITEM_32;

	double outofcache_merge_cost = n2waymerge_rounds * C6 * nitems;

	constexpr double C7 = log(BLOCKSIZE<uint32_t>()/8.0)/ln2  * MERGE_EQLEN_COST_PER_ITEM_32;

	double incache_merge_cost_for_outofcachesize = C7 * nitems;

	constexpr double C8 = MERGE_EQLEN_COST_PER_ITEM_32;

	double incache_merge_cost_for_incachesize = C8 * nitems * log(nitems)/ln2;

	if (nitems > BLOCKSIZE<uint32_t>()) {
		return fix_cost + mm_cpy_cost + scan1_cost + scan2_cost + inregister_sort_cost +
				outofcache_merge_cost + incache_merge_cost_for_outofcachesize;
	} else {
		return fix_cost + mm_cpy_cost + inregister_sort_cost + incache_merge_cost_for_incachesize;
	}


#if 0
	//double total_cost = 0.0;
	double value_width = 4.0;
	double oid_width = 4.0;
	double ln2 = 0.693147;

	if (nitems <= 256) {
		return SORT_FIX_COST/nitems;
	} else {
		double n2waymerge_rounds = log(ceil(4.0*nitems*max<double>(value_width, oid_width)/L2_CACHE_SIZE))/ln2;

		/**a) memory copy before and after sorting**/
		double mm_cpy_cost = nitems * (value_width + oid_width) / CACHE_LINE_SIZE  * 4.0 * COST_PER_SEQ_ACCESS;

		/**b) scan sorting blocks and 2-way merge blocks**/
		double sort_scan_cost = (2.0 + 2.0 * n2waymerge_rounds) *
				nitems * (value_width + oid_width) / CACHE_LINE_SIZE * COST_PER_SEQ_ACCESS;

		double inregister_sort_cost = nitems * SORT_NETWORK_COST_PER_ITEM_32;

		if (nitems > BLOCKSIZE<uint32_t>()) {	//depends on the data volumn is larger than one block or not
			/****/
			return SORT_FIX_COST + mm_cpy_cost + sort_scan_cost + inregister_sort_cost
					+ log(BLOCKSIZE<uint32_t>()/8.0)/ln2 * nitems * MERGE_EQLEN_COST_PER_ITEM_32
					+ n2waymerge_rounds * nitems * MERGE_VALEN_COST_PER_ITEM_32;
		} else {	//smaller than a block, just perform in-block sorting
			/****/
			return SORT_FIX_COST + mm_cpy_cost + inregister_sort_cost + log(nitems)/ln2 * nitems * MERGE_EQLEN_COST_PER_ITEM_32;
		}
	}
#endif
}

double StitchComposer::estimate_mergesort_uint16(double nitems) {

	constexpr double ln2 = 0.693147;

	//fix cost for function calls and memory allocation
	constexpr double fix_cost = SORT_FIX_COST;

	constexpr double C1 = log(16.0/L2_CACHE_SIZE)/ln2;

	double n2waymerge_rounds = C1 + log(nitems)/ln2;

	constexpr double C2 = 24.0 / CACHE_LINE_SIZE  * COST_PER_SEQ_ACCESS;

	double mm_cpy_cost = C2 * nitems;

	constexpr double C3 = 12.0 / CACHE_LINE_SIZE * COST_PER_SEQ_ACCESS;

	double scan1_cost = C3 * nitems;

	constexpr double C4 = C3;

	double scan2_cost = n2waymerge_rounds * C4 * nitems;

	constexpr double C5 = SORT_NETWORK_COST_PER_ITEM_16;

	double inregister_sort_cost = C5 * nitems;

	constexpr double C6 = MERGE_VALEN_COST_PER_ITEM_16;

	double outofcache_merge_cost = n2waymerge_rounds * C6 * nitems;

	constexpr double C7 = log(BLOCKSIZE<uint16_t>()/16.0)/ln2  * MERGE_EQLEN_COST_PER_ITEM_16;

	double incache_merge_cost_for_outofcachesize = C7 * nitems;

	constexpr double C8 = MERGE_EQLEN_COST_PER_ITEM_16;

	double incache_merge_cost_for_incachesize = C8 * nitems * log(nitems)/ln2;

	if (nitems > BLOCKSIZE<uint16_t>()) {
		return fix_cost + mm_cpy_cost + scan1_cost + scan2_cost + inregister_sort_cost +
				outofcache_merge_cost + incache_merge_cost_for_outofcachesize;
	} else {
		return fix_cost + mm_cpy_cost + inregister_sort_cost + incache_merge_cost_for_incachesize;
	}


#if 0
	//double total_cost = 0.0;
	double value_width = 2.0;
	double oid_width = 4.0;
	double ln2 = 0.693147;

	if (nitems <= 256) {
		return SORT_FIX_COST/nitems;
	} else {
		double n2waymerge_rounds = log(ceil(4.0*nitems*max<double>(value_width, oid_width)/L2_CACHE_SIZE))/ln2;

		/**a) memory copy before and after sorting**/
		double mm_cpy_cost = nitems * (value_width + oid_width) / CACHE_LINE_SIZE  * 4.0 * COST_PER_SEQ_ACCESS;

		/**b) scan sorting blocks and 2-way merge blocks**/
		double sort_scan_cost = (2.0 + 2.0 * n2waymerge_rounds) *
				nitems * (value_width + oid_width) / CACHE_LINE_SIZE * COST_PER_SEQ_ACCESS;

		double inregister_sort_cost = nitems * SORT_NETWORK_COST_PER_ITEM_16;

		if (nitems > BLOCKSIZE<uint16_t>()) {	//depends on the data volumn is larger than one block or not
			/****/
			return SORT_FIX_COST + mm_cpy_cost + sort_scan_cost + inregister_sort_cost
					+ log(BLOCKSIZE<uint16_t>()/16.0)/ln2 * nitems * MERGE_EQLEN_COST_PER_ITEM_16
					+ n2waymerge_rounds * nitems * MERGE_VALEN_COST_PER_ITEM_16;
		} else {	//smaller than a block, just perform in-block sorting
			/****/
			return SORT_FIX_COST + mm_cpy_cost + inregister_sort_cost + log(nitems)/ln2 * nitems * MERGE_EQLEN_COST_PER_ITEM_16;
		}
	}
#endif
}

#endif


double StitchComposer::calculate_bowl_int16(const double avg_groupsz) {

	double ln2 = 0.693147;

	constexpr double D = SORT_FIX_COST;

	constexpr double C4 = 12.0 / CACHE_LINE_SIZE * COST_PER_SEQ_ACCESS;

	constexpr double C6 = MERGE_VALEN_COST_PER_ITEM_16;

	constexpr double C8 = MERGE_EQLEN_COST_PER_ITEM_16;

	if (avg_groupsz > BLOCKSIZE<uint16_t>()) {	//out of cache
		return ln2 * D / (C4 + C6);
	} else {
		return ln2 * D / C8;
	}
}

double StitchComposer::calculate_bowl_int32(const double avg_groupsz) {

	double ln2 = 0.693147;

	constexpr double D = SORT_FIX_COST;

	constexpr double C4 = 16.0 / CACHE_LINE_SIZE * COST_PER_SEQ_ACCESS;

	constexpr double C6 = MERGE_VALEN_COST_PER_ITEM_32;

	constexpr double C8 = MERGE_EQLEN_COST_PER_ITEM_32;

	if (avg_groupsz > BLOCKSIZE<uint32_t>()) {	//out of cache
		return ln2 * D / (C4 + C6);
	} else {
		return ln2 * D / C8;
	}
}

double StitchComposer::calculate_bowl_int64(const double avg_groupsz) {

	double ln2 = 0.693147;

	constexpr double D = SORT_FIX_COST;

	constexpr double C4 = 24.0 / CACHE_LINE_SIZE * COST_PER_SEQ_ACCESS;

	constexpr double C6 = MERGE_VALEN_COST_PER_ITEM_64;

	constexpr double C8 = MERGE_EQLEN_COST_PER_ITEM_64;

	if (avg_groupsz > BLOCKSIZE<uint64_t>()) {	//out of cache
		return ln2 * D / (C4 + C6);
	} else {
		return ln2 * D / C8;
	}
}

#if 0
void StitchComposer::getHeuristics(std::vector < std::vector<uint32_t> > & plans, uint32_t colwidth_sum) {
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
#endif

void StitchComposer::optimizeAndRun() {

#if 0
	uint32_t arr[num_columns_];
	uint32_t i;
	for (i = 0; i < num_columns_; ++i) {
		arr[i] = i;
	}

	/** Apply aggressive heuristics to quickly get some plans **/
	std::vector < std::vector<uint32_t> > heuristicPlans;
	getHeuristics(heuristicPlans, compose_params_.colwidth_sum);
	assert(heuristicPlans.size() > 0);

	/** get best plan over these heuristic plans **/
	int heuOptimalIdx = -1;
	double heuOptimalCost = std::numeric_limits<double>::max();
	cost_estimation_plans(arr, heuristicPlans, heuOptimalIdx,
			heuOptimalCost, extractCardinalityNumCycles);
	assert(heuOptimalIdx >= 0);
	assert(heuOptimalIdx < heuristicPlans.size());

	/** test whether it is time to stop optimization **/
	if (timer_optimization.GetAccuCycles() > (THRESHOLD_PERC * heuOptimalCost)) {

	}
#endif
	/** find all factorial of elements (0, 1, ..., n-1), i.e., n! **/
	uint32_t* arr = (uint32_t *) malloc_aligned(num_columns_ * sizeof(uint32_t));
	uint32_t i;
	for (i = 0; i < num_columns_; ++i) {
		arr[i] = i;
	}

	//uint32_t factorial = factor(num_columns_);
	//const uint32_t permu_num = factorial;
	//uint32_t * permutations = (uint32_t *) malloc_aligned(permu_num * num_columns_ * sizeof(uint32_t));

	std::vector<uint32_t *> combinations;
	if (0 == compose_params_.ordered) {	//unordered (GROUP-BY case), need to shuffle the column
		permutations(arr, 0, num_columns_ - 1, combinations);
	} else {
		combinations.push_back(arr);
	}

	assert(combinations.size() > 0);
	std::vector<std::vector<uint32_t> > greedy_plans;
	//std::set<std::vector<uint32_t> > temp_plans;
	uint32_t round_num = 0;

	uint32_t max_round_num = ceil(
			((double) (2 * compose_params_.colwidth_sum - 2)) / ((double) MIN_BANK_SIZE));

	//const uint32_t original_col_num = num_columns_;

	int best_column_order_index = -1;
	std::vector<uint32_t> best_massaging_plan;
	double best_cost = std::numeric_limits<double>::max();

	HybridTimer timer_optimization;
	timer_optimization.Start();

	double accumulatedCycles = 0;
	double checkpointCycles = 0;

	//under each combination, find all massage plans with incremental #rounds
	for (uint32_t order_idx = 0; order_idx < combinations.size(); ++order_idx) {

		uint32_t * targetComb = combinations[order_idx];

		//find all greedy plans
		for (round_num = 1; round_num <= max_round_num; ++round_num) {

			std::vector<std::vector<uint32_t> >().swap(greedy_plans);	//remove the content;
			assert(greedy_plans.empty());

			/** exhaustive search for round_num <= 2 **/
			if (round_num <= 2) {
				split_exhaustive_certain_round(compose_params_.colwidth_sum, 0,
						round_num, 64, greedy_plans);
				//std::cout << "exhaustive plans: " << greedy_plans.size() << std::endl;
			} else { /**greedy search for round_num > 2**/
				//std::set<std::vector<uint32_t> >().swap(temp_plans);
				//assert(temp_plans.empty());
				find_greedy_plans(round_num, greedy_plans, targetComb);
				//greedy_plans.insert(greedy_plans.end(), temp_plans.begin(), temp_plans.end());
			}

			/**cost these plans**/
			if (!greedy_plans.empty()) {
				int bestLocalSplitIdx = -1;		//best plan index for current column order
				double optimalLocCost = -1.0;	//best plan cost for current column order

				/** find the best greedy plans for current column ordering **/
				estimate_plan_collection_greedy(targetComb, greedy_plans, bestLocalSplitIdx, optimalLocCost);

				assert(bestLocalSplitIdx >= 0);
				assert(optimalLocCost > 0);

				/** update the global optimal plan**/
				if (optimalLocCost < best_cost) {
					best_column_order_index = order_idx;
					best_massaging_plan = greedy_plans[bestLocalSplitIdx];
					best_cost = optimalLocCost;
				}

				//early stop: to prevent too long optimization time.
				accumulatedCycles = (timer_optimization.GetAccuCycles() - checkpointCycles);

				if (accumulatedCycles > (THRESHOLD_PERC * best_cost)) {

					checkpointCycles = accumulatedCycles;
					std::cout << "[INFO ] Early Stop at Round # " << round_num << std::endl;
					break;
				}
			}
		}
#if 0
		assert(!greedy_plans.empty());

		int bestLocalSplitIdx = -1;		//best plan index for current column order
		double optimalLocCost = -1.0;	//best plan cost for current column order

		/** find the best greedy plans for current column ordering **/
		estimate_plan_collection_greedy(targetComb, greedy_plans, bestLocalSplitIdx, optimalLocCost);

		assert(bestLocalSplitIdx >= 0);
		assert(optimalLocCost > 0);

		/** update the global optimal plan**/
		if (optimalLocCost < best_cost) {
			best_column_order_index = order_idx;
			best_massaging_plan = greedy_plans[bestLocalSplitIdx];
			best_cost = optimalLocCost;
		}
#endif
	}
	assert(best_column_order_index >= 0);
	assert(!best_massaging_plan.empty());

	timer_optimization.Stop();
	setting_.time_optimization = (double) timer_optimization.GetNumCycles();

	/**Run the plan suggested by the optimizer**/
	double breakdown[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	double total_time = run_single_plan(combinations[best_column_order_index],
			best_massaging_plan, breakdown);

	/** Output the time breakdown**/
	std::cout << "[INFO ] The optimized plan is [";
	for (i = 0; i < best_massaging_plan.size(); ++i) {
		std::cout << best_massaging_plan[i] << ", ";
	}
	std::cout << "] under the ordering <";
	for (i = 0; i < num_columns_; ++i) {
		std::cout << combinations[best_column_order_index][i] << ", ";
	}
	std::cout << ">\n";

	std::cout << "[INFO ] Breakdown <optimization time, other time, [massaging time, execution time]> in ms:" << std::endl;
	//std::cout << setting_.time_optimization << ", " << total_time << ", ["
	//		<< breakdown[0] << ", " << (total_time - breakdown[0]) << "]" << std::endl;
	//std::cout << setting_.time_optimization/3400000.0 << ", " << total_time/3400000.0 << ", ["
	//		<< breakdown[0]/3400000.0 << ", " << (total_time - breakdown[0])/3400000.0 << "]" << std::endl;
	std::cout << setting_.time_optimization/3400000.0/1000.0 << ", " << total_time << ", ["
			<< breakdown[0] << ", " << (total_time - breakdown[0]) << "]" << std::endl;

	//Destroy
	if (0 == compose_params_.ordered) {
		free(arr);
	}
	for (i = 0; i < combinations.size(); ++i) {
		free(combinations[i]);
	}
}

void StitchComposer::SortAllColumns_Pack() {

}

void StitchComposer::HashAllColumns() {

}
void StitchComposer::SortAllColumns_Nonpack() {

	/**running the baseline**/
	HybridTimer timer_baseline;
	timer_baseline.Start();

	multiRoundSorting();

	timer_baseline.Stop();

	/**output the result for baseline**/
	std::cout << "[INFO ] The baseline plan is [";
	for (uint32_t i = 0; i < num_columns_; ++i) {
		std::cout << compose_params_.bitwidth[i] << ", ";
	}
	std::cout << "]\n";

	std::cout << "[INFO ] Baseline output <optimization time, other time, [massaging time, execution time]> in ms:" << std::endl;
	//std::cout << 0.0 << ", " << timer_baseline.GetNumCycles() << ", [0.0, "
	//		<< timer_baseline.GetNumCycles() << "] " << std::endl;
	std::cout << 0.0 << ", " << timer_baseline.GetNumCycles()/3400000.0 << ", [0.0, "
			<< timer_baseline.GetNumCycles()/3400000.0 << "] " << std::endl;
	//std::cout << 0.0 << ", " << timer_baseline.GetSeconds() << ", [0.0, "
	//		<< timer_baseline.GetSeconds() << "] " << std::endl;

	/**running the optimization**/
	optimizeAndRun();
}

void StitchComposer::RunAnInstance(std::string intance) {

	uint32_t nthreads = setting_.nthreads;
	int firstColWidth = (int)compose_params_.bitwidth[0];
	int secondColWidth = (int)compose_params_.bitwidth[1];

	uint32_t columnOrder[2];
	columnOrder[0] = 0;
	columnOrder[1] = 1;
	std::vector<uint32_t> curSplit;

	//for warmup
	int repeat = 5;
	for (int i = 0; i < repeat; ++i) {
		double temp_warmup[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		curSplit.push_back(firstColWidth+secondColWidth);
		run_single_plan(columnOrder, curSplit, temp_warmup);
		curSplit.clear();
	}

	//run baseline
	curSplit.clear();
	double breakdown_baseline[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

	curSplit.push_back(firstColWidth);
	curSplit.push_back(secondColWidth);

	run_single_plan(columnOrder, curSplit, breakdown_baseline);

	double speedup = 1.0;
	if (nthreads > 1) {
		speedup = speedup * (double)nthreads * 0.8;
	}

	std::cout << "0, 0, " << breakdown_baseline[4]/num_rows_/speedup << ", "
			<< (breakdown_baseline[3]-breakdown_baseline[4])/num_rows_/speedup << ", "
			<< breakdown_baseline[1]/num_rows_ << ", "
			<< breakdown_baseline[2]/num_rows_
			<< std::endl;


	//run a specific plan
	curSplit.clear();
	double breakdown_plan[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

	//parse the plan
	uint32_t bw = 0;
	char bw_str_tmp[50];
	strcpy(bw_str_tmp, intance.c_str());
	char *bw_str = strtok(bw_str_tmp, "B");
	while (bw_str != NULL) {
		bw = atoi(bw_str);
		curSplit.push_back(bw);

		bw_str = strtok(NULL, "B");
	}

	assert(!curSplit.empty());

	run_single_plan(columnOrder, curSplit, breakdown_plan);

	std::cout << "1, " << breakdown_plan[0]/num_rows_ << ", "
			<< breakdown_plan[4]/num_rows_/speedup << ", "
			<< (breakdown_plan[3]-breakdown_plan[4])/num_rows_/speedup << ", "
			<< breakdown_plan[1]/num_rows_ << ", "
			<< breakdown_plan[2]/num_rows_
			<< std::endl;
}

void StitchComposer::RunAnInstance_hash(std::string intance) {

	//uint32_t nthreads = setting_.nthreads;
	int firstColWidth = (int)compose_params_.bitwidth[0];
	int secondColWidth = (int)compose_params_.bitwidth[1];

	uint32_t columnOrder[2];
	columnOrder[0] = 0;
	columnOrder[1] = 1;
	std::vector<uint32_t> curSplit;

	//for warmup
#if 0
	int repeat = 5;
	for (int i = 0; i < repeat; ++i) {
		double temp_warmup[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		curSplit.push_back(firstColWidth+secondColWidth);
		run_single_plan_hash(columnOrder, curSplit, temp_warmup);
		curSplit.clear();
	}
#endif

	//run baseline
	curSplit.clear();
	double breakdown_baseline[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

	curSplit.push_back(firstColWidth);
	curSplit.push_back(secondColWidth);

	double total_baseline = run_single_plan_hash(columnOrder, curSplit, breakdown_baseline);

#if 0
	double speedup = 1.0;
	if (nthreads > 1) {
		speedup = speedup * (double)nthreads * 0.8;
	}
#endif

	std::cout << "Summary of baseline <total, massage, first round hash, other round hash>"
			<< std::endl;
	std::cout << (total_baseline - breakdown_baseline[0])/num_rows_ << ", "
			<< "0, " << breakdown_baseline[4]/num_rows_ << ", "
			<< (breakdown_baseline[3]-breakdown_baseline[4])/num_rows_
			<< std::endl;


	//run a specific plan
	curSplit.clear();
	double breakdown_plan[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

	//parse the plan
	uint32_t bw = 0;
	char bw_str_tmp[50];
	strcpy(bw_str_tmp, intance.c_str());
	char *bw_str = strtok(bw_str_tmp, "B");
	while (bw_str != NULL) {
		bw = atoi(bw_str);
		curSplit.push_back(bw);

		bw_str = strtok(NULL, "B");
	}

	assert(!curSplit.empty());

	double total_specific = run_single_plan_hash(columnOrder, curSplit, breakdown_plan);

	std::cout << "Summary of specific instance <total, massage, first round hash, other round hash>"
			<< std::endl;
	std::cout << total_specific/num_rows_ << ", "
			<< breakdown_plan[0]/num_rows_ << ", "
			<< breakdown_plan[4]/num_rows_ << ", "
			<< (breakdown_plan[3]-breakdown_plan[4])/num_rows_
			<< std::endl;
}

void StitchComposer::TwoRoundsExhaustive() {

	//get all plans for baseline with two rounds
	assert(num_columns_ == 2);

	int firstColWidth = (int)compose_params_.bitwidth[0];
	int secondColWidth = (int)compose_params_.bitwidth[1];

	uint32_t columnOrder[2];
	columnOrder[0] = 0;
	columnOrder[1] = 1;
	std::vector<uint32_t> curSplit;

	int new1stColWidth;
	int new2ndColWidth;
	bool isBaseline = false;

	//warm up the cache:
#if 1
	int repeat = 5;
	for (int i = 0; i < repeat; ++i) {
		double temp_warmup[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		curSplit.push_back(firstColWidth+secondColWidth);
		run_single_plan(columnOrder, curSplit, temp_warmup);
		curSplit.clear();
	}
#endif

	int shiftBits = -secondColWidth;	//original column: 10|17, range of shiftBits: [-17, 10]
	for (; shiftBits <= firstColWidth; shiftBits++) {

		curSplit.clear();

		new1stColWidth = firstColWidth - shiftBits;
		new2ndColWidth = secondColWidth + shiftBits;
		assert(new1stColWidth >= 0);
		assert(new2ndColWidth >= 0);

		isBaseline = (shiftBits == 0);

		if (new1stColWidth > 0) {
			curSplit.push_back(new1stColWidth);
		}
		if (new2ndColWidth > 0) {
			curSplit.push_back(new2ndColWidth);
		}
		assert((curSplit.size() > 0) && (curSplit.size() < 3));

		double cost_breakdown[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

		/**
		 * cost_breakdown[3] is all rounds of sorting cost; cost_breakdown[4] is first round sorting cost
		 * other cost can be calculated as (totalcost - cost_breakdown[3])
		 */
		double totalcost = run_single_plan(columnOrder, curSplit, cost_breakdown);

		std::cout << shiftBits << ", " << cost_breakdown[4]/num_rows_ << ", "
				<< (cost_breakdown[3]-cost_breakdown[4])/num_rows_ << ", "
				<< ((isBaseline) ? (totalcost - cost_breakdown[3] - cost_breakdown[0])/num_rows_ : (totalcost - cost_breakdown[3])/num_rows_)
				<< std::endl;
	}
}

void StitchComposer::TwoRoundsExhaustive_hash() {

	//get all plans for baseline with two rounds
	assert(num_columns_ == 2);

	int firstColWidth = (int)compose_params_.bitwidth[0];
	int secondColWidth = (int)compose_params_.bitwidth[1];

	uint32_t columnOrder[2];
	columnOrder[0] = 0;
	columnOrder[1] = 1;
	std::vector<uint32_t> curSplit;

	int new1stColWidth;
	int new2ndColWidth;
	bool isBaseline = false;

	//warm up the cache:
	std::cout << "[INFO ] Warming up the cache..." << std::endl;
#if 1
	int repeat = 5;
	for (int i = 0; i < repeat; ++i) {
		double temp_warmup[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		curSplit.push_back(firstColWidth+secondColWidth);
		run_single_plan_hash(columnOrder, curSplit, temp_warmup);
		curSplit.clear();
	}
#endif

	std::cout << "[INFO ] Find all combinations for two rounds." << std::endl;
#if 0
	int shiftBits = firstColWidth;	//original column: 10|17, range of shiftBits: [-17, 10] ==> [10, -17]
	for (; shiftBits >= -secondColWidth; shiftBits--) {
#endif
	int shiftBits = -secondColWidth;	//original column: 10|17, range of shiftBits: [-17, 10]
	for (; shiftBits <= firstColWidth; shiftBits++) {
		curSplit.clear();

		new1stColWidth = firstColWidth - shiftBits;
		new2ndColWidth = secondColWidth + shiftBits;
		assert(new1stColWidth >= 0);
		assert(new2ndColWidth >= 0);

		isBaseline = (shiftBits == 0);

		if (new1stColWidth > 0) {
			//std::cout << new1stColWidth << ", ";
			curSplit.push_back(new1stColWidth);
		}
		if (new2ndColWidth > 0) {
			//std::cout << new2ndColWidth << std::endl;
			curSplit.push_back(new2ndColWidth);
		}
		assert((curSplit.size() > 0) && (curSplit.size() < 3));

		double cost_breakdown[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

		/**
		 * cost_breakdown[3] is all rounds of sorting cost; cost_breakdown[4] is first round sorting cost
		 * other cost can be calculated as (totalcost - cost_breakdown[3])
		 */
		run_single_plan_hash(columnOrder, curSplit, cost_breakdown);

		//std::cout << shiftBits << ", " << cost_breakdown[4]/num_rows_ << ", "
		//		<< (cost_breakdown[3]-cost_breakdown[4])/num_rows_ << ", "
		//		<< ((isBaseline) ? (totalcost - cost_breakdown[3] - cost_breakdown[0])/num_rows_ : (totalcost - cost_breakdown[3])/num_rows_)
		//		<< std::endl;

		//[first pass hashing, remaining pass hashing, massage, reordering]
		std::cout << shiftBits << ", " << cost_breakdown[4]/num_rows_ << ", "
				<< (cost_breakdown[3]-cost_breakdown[4])/num_rows_ << ", "
				<< ((isBaseline) ? 0.0 : cost_breakdown[0]/num_rows_) << ", "
				<< cost_breakdown[2]/num_rows_
				<< std::endl;
	}
}

void StitchComposer::CostModelComparison() {

#if 1
	/** find all factorial of elements (0, 1, ..., n-1), i.e., n! **/
	uint32_t* arr = (uint32_t*) malloc_aligned(num_columns_ * sizeof(uint32_t));
	uint32_t i;
	for (i = 0; i < num_columns_; ++i) {
		arr[i] = i;
	}

	std::vector<uint32_t *> combinations;
	if (0 == compose_params_.ordered) {	//unordered (GROUP-BY case), need to shuffle the column
		permutations(arr, 0, num_columns_ - 1, combinations);
	} else {
		combinations.push_back(arr);
	}

	assert(combinations.size() > 0);
	std::vector < std::vector<uint32_t> > split_enum;
	uint32_t round_num = 0;

	//calculate method ONE for max_round_num
#if 0
	uint32_t max_round_num = floor(
			((double) (compose_params_.colwidth_sum)) / ((double) MIN_BANK_SIZE + 1))*2+1;
#endif
	//calculate method TWO for max_round_num
#if 0
	uint32_t max_round_num = ceil(
					((double) (compose_params_.colwidth_sum)) / ((double) MIN_BANK_SIZE));
#endif
	//calculate method THREE for max_round_num
	uint32_t max_round_num = ceil(
			((double) (2 * compose_params_.colwidth_sum - 2)) / ((double) MIN_BANK_SIZE));

	const uint32_t origin_col_num = num_columns_;

	uint64_t startIndex = 0;	//for indicating the total order of plans
	std::set<std::vector<uint32_t> > randomPlans;

	/** under each combination, find all massage plans with incremental #rounds **/
	for (i = 0; i < combinations.size(); ++i) {

		uint32_t * targetComb = combinations[i];

		/** Prepare for the RRS plan enumerator **/
		/******************************************************************/
		/**Step 1) get the time for greedy**/
#if 0
		HybridTimer timer_for_greedy;
		timer_for_greedy.Start();

		std::vector<std::vector<uint32_t> > greedyPlans_temp;
		find_greedy_plans(3, greedyPlans_temp, targetComb);
		timer_for_greedy.Stop();

		double cycles_for_greedy = (double)timer_for_greedy.GetNumCycles();
#endif
		double cycles_for_greedy = 1000;

		/**Step 2) Initialize the random enumerator**/
		//uint64_t base_col_card[origin_col_num];

		/**calculate the bit info carried by each bit in each base column**/
		//for (uint32_t colIdx = 0; colIdx < origin_col_num; ++colIdx) {
		//	base_col_card[colIdx] = column_values_[targetComb[colIdx]]->GetCardinality();
		//}

		PlanEnumRRS random_enumerator(num_rows_, compose_params_.colwidth_sum, compose_params_.colwidth_sum,
				std::vector<uint32_t>(targetComb, targetComb+origin_col_num), num_columns_,
				std::vector<uint32_t>(compose_params_.bitwidth, compose_params_.bitwidth+origin_col_num),
				std::vector<uint32_t>(base_col_card, base_col_card+origin_col_num), cycles_for_greedy);

		random_enumerator.Optimize();

		std::set<std::vector<uint32_t> >().swap(randomPlans);
		randomPlans.clear();
		randomPlans = random_enumerator.ObtainPlans();

		assert(!randomPlans.empty());
		/******************************************************************/

		for (round_num = 1; round_num <= max_round_num; ++round_num) {

			std::vector<std::vector<uint32_t> >().swap(split_enum);	//remove the content;
			assert(split_enum.empty());

			/****/
			//std::cout << "Error happens in split_exhaustive_certain_round? round=" << round_num << std::endl;
			split_exhaustive_certain_round(compose_params_.colwidth_sum, 0, round_num, 64, split_enum);
			std::cout << "[INFO ] For column order " << i << ", round_num = " << round_num <<
					", #plans = " << split_enum.size() << std::endl;

			std::set<std::vector<uint32_t> > greedyPlans;
			std::vector<std::vector<uint32_t> > temp_plans;

			/** Greedy plans **/
			if (round_num > 2) {
				find_greedy_plans(round_num, temp_plans, targetComb);
				greedyPlans.insert(temp_plans.begin(), temp_plans.end());

#if 0
				std::cout << "find greedy plan: " << greedyPlans.size() << "==>";
				if (!greedyPlans.empty()) {
					for (uint32_t i = 0; i < greedyPlans[0].size(); ++i) {
						std::cout << greedyPlans[0][i] << "\t";
					}
				}
				std::cout << std::endl;
#endif
			}

			/**Exhaustively estimate all these plans**/
			estimate_plan_collection(targetComb, split_enum, startIndex, greedyPlans, randomPlans);
			startIndex += split_enum.size();
#if 1
			if (startIndex > 1000) {
				return;
			}
#endif
		}
	}

	//destroy combinations
	if (0 == compose_params_.ordered) {
		free(arr);
	}

	for (uint32_t i = 0; i < combinations.size(); ++i) {
		free(combinations[i]);
	}
#endif
}

void StitchComposer::ExhaustiveSearch() {

	HybridTimer timer_exhaustive;
	timer_exhaustive.Start();
#if 1

	/** find all factorial of elements (0, 1, ..., n-1), i.e., n! **/
	uint32_t* arr = (uint32_t*) malloc_aligned(num_columns_ * sizeof(uint32_t));
	uint32_t i;
	for (i = 0; i < num_columns_; ++i) {
		arr[i] = i;
	}

	std::vector<uint32_t *> combinations;
	if (0 == compose_params_.ordered) {	//unordered (GROUP-BY case), need to shuffle the column
		permutations(arr, 0, num_columns_ - 1, combinations);
	} else {
		combinations.push_back(arr);
	}

	assert(combinations.size() > 0);
	std::vector < std::vector<uint32_t> > split_enum;
	uint32_t round_num = 0;

	//calculate method ONE for max_round_num
#if 0
	uint32_t max_round_num = floor(
			((double) (compose_params_.colwidth_sum)) / ((double) MIN_BANK_SIZE + 1))*2+1;
#endif
	//calculate method TWO for max_round_num
#if 0
	uint32_t max_round_num = ceil(
					((double) (compose_params_.colwidth_sum)) / ((double) MIN_BANK_SIZE));
#endif
	//calculate method THREE for max_round_num
	uint32_t max_round_num = ceil(
			((double) (2 * compose_params_.colwidth_sum - 2)) / ((double) MIN_BANK_SIZE));

#if 0
	if (max_round_num > 5) {
		max_round_num = 5;
	}
#endif

	/**warm up the cache by running the baseline**/
	double cost_breakdown_temp[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	std::vector<uint32_t> baseline;
	for (uint32_t idx = 0; idx < num_columns_; ++idx) {
		baseline.push_back(compose_params_.bitwidth[idx]);
	}

	for (int i = 0; i < 10; ++i) {
		run_single_plan(combinations[0], baseline, cost_breakdown_temp);
	}
	/**********************************************/

	uint64_t startIndex = 0;	//for indicating the total order of plans

	/** under each combination, find all massage plans with incremental #rounds **/
	for (i = 0; i < combinations.size(); ++i) {

		uint32_t * targetComb = combinations[i];

		for (round_num = 1; round_num <= max_round_num; ++round_num) {

			std::vector<std::vector<uint32_t> >().swap(split_enum);	//remove the content;
			assert(split_enum.empty());

			split_exhaustive_certain_round(compose_params_.colwidth_sum, 0, round_num, 64, split_enum);
			std::cout << "[INFO ] For column order " << i << ", round_num = " << round_num <<
					", #plans = " << split_enum.size() << std::endl;

			/**Run all these plans to get actual cost**/
			run_plan_collection(targetComb, split_enum, startIndex);
			startIndex += split_enum.size();

#if 0
			//stop if it runs too slow (succeed one hour)
			if (timer_exhaustive.GetAccuSeconds() > 3600.0) {
				return;
			}
#endif

#if 1
			if (startIndex > 1000) {
				return;
			}
#endif
		}
	}
#endif
	timer_exhaustive.Stop();

	//destroy combinations
	if (0 == compose_params_.ordered) {
		free(arr);
	}

	for (uint32_t i = 0; i < combinations.size(); ++i) {
		free(combinations[i]);
	}
}

#if 0
void StitchComposer::optimizer_exhaustive() {

	HybridTimer timer_optimization;
	timer_optimization.Start();
#if 1

	/** find all factorial of elements (0, 1, ..., n-1), i.e., n! **/
	uint32_t* arr = (uint32_t*) malloc_aligned(num_columns_ * sizeof(uint32_t));
	uint32_t i;
	for (i = 0; i < num_columns_; ++i) {
		arr[i] = i;
	}

	std::vector<uint32_t *> combinations;
	if (0 == compose_params_.ordered) {	//unordered (GROUP-BY case), need to shuffle the column
		permutations(arr, 0, num_columns_ - 1, combinations);
	} else {
		combinations.push_back(arr);
	}

	std::vector < std::vector<uint32_t> > split_enum;//enumerate all split methods
	double optimalCost = std::numeric_limits<double>::max();
	double optimalCost_ac = std::numeric_limits<double>::max();
	std::vector < uint32_t > bestSplit;
	std::vector < uint32_t > bestSplit_ac;
	int bestCombIdx = -1;
	int bestCombIdx_ac = -1;
	double baseline_cost = -1.0;
	double baseline_cost_each_enum = -1.0;

	double bestSplit_es2ac = -1.0;
	double bestSplit_es2ac_each_enum = -1.0;

	assert(combinations.size() > 0);
	for (i = 0; i < combinations.size(); ++i) {
		std::vector<std::vector<uint32_t> >().swap(split_enum);	//remove the content;
		assert(split_enum.empty());

		uint32_t * targetComb = combinations[i];
		/** estimate when the average group size hit a threshold **/
		//std::cout << "[INFO ] Estimating the threshold location... " << std::endl;
		//uint32_t threshold_loc = estimateThresholdLoc(targetComb);
		uint32_t threshold_loc = 1;

		assert(
				(threshold_loc >= 0)
						&& (threshold_loc <= compose_params_.colwidth_sum));
		std::cout << "[INFO ] For column order " << i << ", the threshold location is at " << threshold_loc << std::endl;

		//uint32_t maxRounds = 0;

		//disable the prunning rule first
		uint32_t maxRounds = ceil(
				((double) (compose_params_.colwidth_sum)) / ((double) MIN_BANK_SIZE)) + 1;
#if 0
		if (threshold_loc == compose_params_.colwidth_sum) {//no threshold exists
			maxRounds = ceil(
					((double) compose_params_.colwidth_sum)
							/ ((double) MIN_BANK_SIZE));
		} else {
			maxRounds = ceil(
					((double) (threshold_loc + 1)) / ((double) MIN_BANK_SIZE))
					+ ceil(
							((double) (compose_params_.colwidth_sum
									- threshold_loc - 1))
									/ ((double) MIN_BANK_SIZE));
		}
		assert(
				(maxRounds > 0)
						&& maxRounds
								< (ceil(
										((double )compose_params_.colwidth_sum)
												/ ((double )MIN_BANK_SIZE)) + 3));
#endif

		/** find all split methods given an ordering of columns **/
		//split(compose_params_.colwidth_sum, 0, maxRounds, 64, threshold_loc,
		//		compose_params_.colwidth_sum, split_enum);
		split_applyRule1(compose_params_.colwidth_sum, 0, maxRounds, 64, split_enum, false);

		/**sort split_enum according to num of rounds**/
		/** larger rounds first, to warm up the cache **/
		std::vector<plan_rounds_pair_t> plan_pairs;
		plan_rounds_pair_t eachpair;
		for (uint32_t idx = 0; idx < split_enum.size(); ++idx) {
			eachpair.plan = split_enum[idx];
			eachpair.rounds = split_enum[idx].size();
			plan_pairs.push_back(eachpair);
		}
		std::stable_sort(plan_pairs.begin(), plan_pairs.end(), plancompare);
		for (uint32_t idx = 0; idx < split_enum.size(); ++idx) {
			split_enum[idx] = plan_pairs[idx].plan;
		}

		/****************Add the baseline plan to to the END of enumerations*******************/
#if 1
		std::vector<uint32_t> baseline;
		for (uint32_t idx = 0; idx < num_columns_; ++idx) {
			baseline.push_back(compose_params_.bitwidth[targetComb[idx]]);
		}
		split_enum.push_back(baseline);
#endif

		int bestLocalSplitIdx = -1;
		double optimalLocCost = -1.0;
		int bestLocalSplitIdx_ac = -1;
		double optimalLocCost_ac = -1.0;
		/** find the best split strategy for current column ordering **/
		cost_estimation(targetComb, split_enum, bestLocalSplitIdx,
				optimalLocCost, bestLocalSplitIdx_ac, optimalLocCost_ac, baseline_cost_each_enum,
				bestSplit_es2ac_each_enum);

		/**update baseline cost only when it is the first combination**/
		if (0 == i) {
			baseline_cost = baseline_cost_each_enum;
		}

		assert(bestLocalSplitIdx >= 0);
		assert(optimalLocCost > 0);
		assert(bestLocalSplitIdx_ac >= 0);
		assert(optimalLocCost_ac > 0);
		/** update the optimal plan**/
		if (optimalLocCost < optimalCost) {
			bestCombIdx = i;
			optimalCost = optimalLocCost;
			bestSplit = split_enum[bestLocalSplitIdx];

			bestSplit_es2ac = bestSplit_es2ac_each_enum;
		}

		/* Auxiliary*/
		if (optimalLocCost_ac < optimalCost_ac) {
			bestCombIdx_ac = i;
			optimalCost_ac = optimalLocCost_ac;
			bestSplit_ac = split_enum[bestLocalSplitIdx_ac];
		}
	}
	assert(bestCombIdx >= 0);
	assert(optimalCost > 0);
	assert(!bestSplit.empty());
	assert(bestCombIdx_ac >= 0);
	assert(optimalCost_ac > 0);
	assert(!bestSplit_ac.empty());
	assert(baseline_cost > 0);

	/** assign the final decision to the member variables **/
	memcpy(bestOrdering, combinations[bestCombIdx],
			num_columns_ * sizeof(uint32_t));
	bestGlobalSplit = bestSplit;

#if 1	//for statistics
	std::cout << "summary, " << baseline_cost/num_rows_ << ", " << bestSplit_es2ac/num_rows_ << ", "
			<< baseline_cost / bestSplit_es2ac  << ", "
			<< optimalCost_ac/num_rows_ << ", "
			<< baseline_cost / optimalCost_ac << ", ";
	for (i = 0; i < num_columns_; ++i) {	//the baseline
		std::cout << compose_params_.bitwidth[i] << "_";
	}
	std::cout << ", ";
	for (i = 0; i < bestGlobalSplit.size(); ++i) {	//the estimated optimal plan
		std::cout << bestGlobalSplit[i] << "_";
	}
	std::cout << "<";
	for (i = 0; i < num_columns_; ++i) {
		std::cout << bestOrdering[i] << ", ";
	}
	std::cout << ">, ";

	for (i = 0; i < bestSplit_ac.size(); ++i) {	//the actual optimal plan
		std::cout << bestSplit_ac[i] << "_";
	}
	std::cout << "<";
	for (i = 0; i < num_columns_; ++i) {
		std::cout << combinations[bestCombIdx_ac][i] << ", ";
	}
	std::cout << ">, ";

#endif


#if 0
	std::cout << "[INFO ] The best estimated plan is [";
	for (i = 0; i < bestGlobalSplit.size(); ++i) {
		std::cout << bestGlobalSplit[i] << ", ";
	}
	std::cout << "] under the ordering <";
	for (i = 0; i < num_columns_; ++i) {
		std::cout << bestOrdering[i] << ", ";
	}
	std::cout << ">\n";

	std::cout << "[INFO ] The best actual plan is [";
	for (i = 0; i < bestSplit_ac.size(); ++i) {
		std::cout << bestSplit_ac[i] << ", ";
	}
	std::cout << "] under the ordering <";
	for (i = 0; i < num_columns_; ++i) {
		std::cout << combinations[bestCombIdx_ac][i] << ", ";
	}
	std::cout << ">\n";
#endif

	//Destroy
	if (0 == compose_params_.ordered) {
		free(arr);
	}
	for (i = 0; i < combinations.size(); ++i) {
		free(combinations[i]);
	}

#endif

	timer_optimization.Stop();
	setting_.time_optimization += (double) timer_optimization.GetNumCycles();
	std::cout << timer_optimization.GetSeconds() << std::endl;
}
#endif

#if 0
void StitchComposer::optimizer_profiling() {

	HybridTimer timer_optimization;
	timer_optimization.Start();

	/** find all factorial of elements (0, 1, ..., n-1), i.e., n! **/
	uint32_t* arr = (uint32_t*) malloc_aligned(num_columns_ * sizeof(uint32_t));
	uint32_t i;
	for (i = 0; i < num_columns_; ++i) {
		arr[i] = i;
	}

	std::vector<uint32_t *> combinations;
	if (0 == compose_params_.ordered) {	//unordered (GROUP-BY case), need to shuffle the column
		permutations(arr, 0, num_columns_ - 1, combinations);
	} else {
		combinations.push_back(arr);
	}

	std::vector < std::vector<uint32_t> > split_enum;//enumerate all split methods
	double optimalCost = std::numeric_limits<double>::max();
	double optimalCost_ac = std::numeric_limits<double>::max();
	std::vector < uint32_t > bestSplit;
	std::vector < uint32_t > bestSplit_ac;
	int bestCombIdx = -1;
	//int bestCombIdx_ac = -1;
	//double baseline_cost = -1.0;
	double baseline_cost_each_enum = -1.0;

	//double bestSplit_es2ac = -1.0;
	double bestSplit_es2ac_each_enum = -1.0;

	assert(combinations.size() > 0);
	for (i = 0; i < combinations.size(); ++i) {
		std::vector<std::vector<uint32_t> >().swap(split_enum);	//remove the content;
		assert(split_enum.empty());

		uint32_t * targetComb = combinations[i];
		/** estimate when the average group size hit a threshold **/
		//std::cout << "[INFO ] Estimating the threshold location... " << std::endl;
		//uint32_t threshold_loc = estimateThresholdLoc(targetComb);
		uint32_t threshold_loc = 1;

		assert(
				(threshold_loc >= 0)
						&& (threshold_loc <= compose_params_.colwidth_sum));
		//std::cout << "[INFO ] For column order " << i << ", the threshold location is at " << threshold_loc << std::endl;

		//uint32_t maxRounds = 0;

		//disable the prunning rule first
		uint32_t maxRounds = ceil(
				((double) (compose_params_.colwidth_sum)) / ((double) MIN_BANK_SIZE)) + 2;
#if 0
		if (threshold_loc == compose_params_.colwidth_sum) {//no threshold exists
			maxRounds = ceil(
					((double) compose_params_.colwidth_sum)
							/ ((double) MIN_BANK_SIZE));
		} else {
			maxRounds = ceil(
					((double) (threshold_loc + 1)) / ((double) MIN_BANK_SIZE))
					+ ceil(
							((double) (compose_params_.colwidth_sum
									- threshold_loc - 1))
									/ ((double) MIN_BANK_SIZE));
		}
		assert(
				(maxRounds > 0)
						&& maxRounds
								< (ceil(
										((double )compose_params_.colwidth_sum)
												/ ((double )MIN_BANK_SIZE)) + 3));
#endif

		/** find all split methods given an ordering of columns **/
		split_exhaustive_all_round(compose_params_.colwidth_sum, 0, maxRounds, 64, threshold_loc,
				compose_params_.colwidth_sum, split_enum);

#if 0
		/**sort split_enum according to num of rounds**/
		std::vector<plan_rounds_pair_t> plan_pairs;
		plan_rounds_pair_t eachpair;
		for (uint32_t idx = 0; idx < split_enum.size(); ++idx) {
			eachpair.plan = split_enum[idx];
			eachpair.rounds = split_enum[idx].size();
			plan_pairs.push_back(eachpair);
		}
		std::stable_sort(plan_pairs.begin(), plan_pairs.end(), plancompare);
		for (uint32_t idx = 0; idx < split_enum.size(); ++idx) {
			split_enum[idx] = plan_pairs[idx].plan;
		}
#endif

		/****************Add the baseline plan to to the END of enumerations*******************/
#if 1
		std::vector<uint32_t> baseline;
		for (uint32_t idx = 0; idx < num_columns_; ++idx) {
			baseline.push_back(compose_params_.bitwidth[targetComb[idx]]);
		}
		split_enum.push_back(baseline);
#endif

		int bestLocalSplitIdx = -1;
		double optimalLocCost = -1.0;
		int bestLocalSplitIdx_ac = -1;
		double optimalLocCost_ac = -1.0;
		/** find the best split strategy for current column ordering **/
		cost_estimation_profiling(targetComb, split_enum, bestLocalSplitIdx,
				optimalLocCost, bestLocalSplitIdx_ac, optimalLocCost_ac, baseline_cost_each_enum,
				bestSplit_es2ac_each_enum);

		/**update baseline cost only when it is the first combination**/
		//if (0 == i) {
		//	baseline_cost = baseline_cost_each_enum;
		//}

		assert(bestLocalSplitIdx >= 0);
		assert(optimalLocCost > 0);
		assert(bestLocalSplitIdx_ac >= 0);
		assert(optimalLocCost_ac > 0);
		/** update the optimal plan**/
		if (optimalLocCost < optimalCost) {
			bestCombIdx = i;
			optimalCost = optimalLocCost;
			bestSplit = split_enum[bestLocalSplitIdx];

			//bestSplit_es2ac = bestSplit_es2ac_each_enum;
		}

		/* Auxiliary*/
		if (optimalLocCost_ac < optimalCost_ac) {
			//bestCombIdx_ac = i;
			optimalCost_ac = optimalLocCost_ac;
			bestSplit_ac = split_enum[bestLocalSplitIdx_ac];
		}
	}
	assert(bestCombIdx >= 0);
	assert(optimalCost > 0);
	assert(!bestSplit.empty());
	//assert(bestCombIdx_ac >= 0);
	assert(optimalCost_ac > 0);
	assert(!bestSplit_ac.empty());
	//assert(baseline_cost > 0);

	/** assign the final decision to the member variables **/
	memcpy(bestOrdering, combinations[bestCombIdx],
			num_columns_ * sizeof(uint32_t));
	bestGlobalSplit = bestSplit;

#if 0	//for statistics
	std::cout << "summary, " << baseline_cost/num_rows_ << ", " << bestSplit_es2ac / num_rows_ << ", "
			<< baseline_cost / bestSplit_es2ac  << ", "
			<< optimalCost_ac / num_rows_ << ", "
			<< baseline_cost / optimalCost_ac << ", ";
	for (i = 0; i < num_columns_; ++i) {	//the baseline
		std::cout << compose_params_.bitwidth[i] << "_";
	}
	std::cout << ", ";
	for (i = 0; i < bestGlobalSplit.size(); ++i) {	//the estimated optimal plan
		std::cout << bestGlobalSplit[i] << "_";
	}
	std::cout << "<";
	for (i = 0; i < num_columns_; ++i) {
		std::cout << bestOrdering[i] << ", ";
	}
	std::cout << ">, ";

	for (i = 0; i < bestSplit_ac.size(); ++i) {	//the actual optimal plan
		std::cout << bestSplit_ac[i] << "_";
	}
	std::cout << "<";
	for (i = 0; i < num_columns_; ++i) {
		std::cout << combinations[bestCombIdx_ac][i] << ", ";
	}
	std::cout << ">\n";

#endif


#if 0
	std::cout << "[INFO ] The best estimated plan is [";
	for (i = 0; i < bestGlobalSplit.size(); ++i) {
		std::cout << bestGlobalSplit[i] << ", ";
	}
	std::cout << "] under the ordering <";
	for (i = 0; i < num_columns_; ++i) {
		std::cout << bestOrdering[i] << ", ";
	}
	std::cout << ">\n";

	std::cout << "[INFO ] The best actual plan is [";
	for (i = 0; i < bestSplit_ac.size(); ++i) {
		std::cout << bestSplit_ac[i] << ", ";
	}
	std::cout << "] under the ordering <";
	for (i = 0; i < num_columns_; ++i) {
		std::cout << combinations[bestCombIdx_ac][i] << ", ";
	}
	std::cout << ">\n";
#endif

	//Destroy
	if (0 == compose_params_.ordered) {
		free(arr);
	}
	for (i = 0; i < combinations.size(); ++i) {
		free(combinations[i]);
	}

	timer_optimization.Stop();
	setting_.time_optimization += (double) timer_optimization.GetNumCycles();
	std::cout << "summary, " << (double) timer_optimization.GetNumCycles() << std::endl;
}
#endif

}
