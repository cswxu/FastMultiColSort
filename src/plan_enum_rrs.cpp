#include "plan_enum_rrs.h"
#include "common.h"
#include "stitch_composer.h"

namespace multiAttrSort{

double PlanEnumRRS::calculate_single_column_cost(const double groupsz, const uint32_t banksz_type) {

	assert(groupsz > 0);
	assert((banksz_type >= 0) && (banksz_type < 3));
	double group_num = ((double)num_rows_) / groupsz;
	double factor = 1.0;
	double alpha = 1.0;

	if ((groupsz - 1.0) < 0.000001) {
		factor = 0;
	} else {
		factor = 1.0/(1.0+pow(1.0/(groupsz-1.0), alpha));
	}

	switch(banksz_type) {
		case 0:	//banksize 16
			return group_num * factor * StitchComposer::estimate_mergesort_uint16(groupsz);
			break;
		case 1: //banksize 32
			return group_num * factor * StitchComposer::estimate_mergesort_uint32(groupsz);
			break;
		case 2: //banksize 64
			return group_num * factor * StitchComposer::estimate_mergesort_uint64(groupsz);
			break;
		default:
			std::cout << "[ERROR ] No such banksz_type: " << banksz_type << std::endl;
			exit(EXIT_SUCCESS);
	}
}

PlanEnumRRS::PointType PlanEnumRRS::GetRandomSample() {

	//double x = Dice() * 2 * radius_ - radius_;
	//double y = Dice() * 2 * radius_ - radius_;
	//return PointType({x,y});

	/**get a random sample by dicing each round with an integer (<=64)**/
	//col_width_sum, max_round_

	PointType sample;
	uint32_t sample_round_num = (uint32_t)(ceil(Dice() * max_round_));
	assert(sample_round_num > 0);
	assert(sample_round_num <= max_round_);

	uint32_t bit_num_avail = col_width_sum_;
	uint32_t next_width = 0;
	uint32_t cur_round = 0;

	while (true) {
		sample.clear();
		//sample_round_num = min<uint32_t>((uint32_t)(ceil(Dice() * max_round_) + 2), max_round_);
		sample_round_num = min<uint32_t>((uint32_t)(ceil(Dice() * max_round_)), max_round_);

		assert(sample_round_num > 0);
		assert(sample_round_num <= max_round_);

		bit_num_avail = col_width_sum_;
		next_width = 0;
		cur_round = 0;
		bool isValid = false;

		while (cur_round < (sample_round_num - 1)) {
			next_width = 1 + (uint32_t)(floor(Dice() * (double)(bit_num_avail-1)));
			next_width = min<uint32_t>(next_width, 64);
			assert(next_width > 0);
			sample.push_back(next_width);

			bit_num_avail -= next_width;
			cur_round++;

			if (0 == bit_num_avail) {
				break;
			}
		}

		//filling the last round
		isValid = ((bit_num_avail <= 64) && (bit_num_avail > 0));
		if (isValid) {
			sample.push_back(bit_num_avail);
			assert(sample_round_num == sample.size());
			searchedPlans.insert(sample);

			return sample;
		}
	}
}

PlanEnumRRS::CostType PlanEnumRRS::GetCost(const PointType& split) {

	/**
	 * Given a split, return its estimated cost
	 */

	const uint32_t new_column_num = split.size();

	//byte_num_sum = 0;
	uint32_t byte_num_sum_remove_1st = 0;
	uint32_t bytewidths[new_column_num];
	uint32_t idx;
	for (idx = 0; idx < new_column_num; ++idx) {
		assert(split[idx] > 0);
		assert(split[idx] <= 64);
		if (split[idx] <= 16) {

			byte_num_sum_remove_1st += 2;
			bytewidths[idx] = 2;
		} else if (split[idx] <= 32) {

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
	uint64_t cardinality[new_column_num];
	for (idx = 0; idx < new_column_num; ++idx) {
		cardinality[idx] = 1;
	}

	//also need assemble statistics for costing stitching overhead
	uint32_t assembleStatistics[3][3] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

	uint32_t shiftStatistics[3] = {0, 0, 0};

	GetStatistics(split, cardinality, assembleStatistics, shiftStatistics);

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

	double multi_pass_sort_cost = 0.0;
	double avg_group_sz = 0.0;
	//double factor = 1.0;
	//double alpha = 1.0;
	uint32_t colIdx;
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
	}

	/** a) stitching/shifting overhead **/
	double stitch_cost = StitchComposer::estimate_stitch(assembleStatistics, shiftStatistics, num_rows_);

	/** b) tuple reconstruction overhead (i.e., BATproject function)**/
	/** b1. random access. b2. sequential access over new created values **/
	//tuple_reconstruct_cost = ((new_column_num - 1) * num_rows_ +
	//		byte_num_sum_remove_1st * num_rows_ / CACHE_LINE_SIZE) * COST_PER_SEQ_ACCESS;
	double tuple_reconstruct_cost = StitchComposer::estimate_tuple_reconstruct(bytewidths, new_column_num, num_rows_);
	//tuple_reconstruct_cost = 0;

	/** c) group information extraction**/
	double order_group_cost = 0;	//temporary

	/** aggregation for the total estimated cost**/
	double estimatedTotalCost = stitch_cost + tuple_reconstruct_cost + order_group_cost + multi_pass_sort_cost;

	assert(estimatedTotalCost > 0);

	return estimatedTotalCost;
}

bool PlanEnumRRS::StopCritera() {
	//double accuCycles = timer.GetAccuCycles();
	//return accuCycles > time_limit_;
	counter_++;
	return counter_ > 10;
}

PlanEnumRRS::PointType PlanEnumRRS::GetRandomNeighbor(const PointType& point, double ro) {

	std::vector<uint32_t> transformed_split(point.begin(), point.end());

	//uint32_t original_col_num = point.size();

	/**ro is [0, 1]**/
	uint32_t total_steps = ceil(col_width_sum_ * ro);

	uint32_t remain_step = total_steps;

	int column_to_decrease_bits = 0;
	int column_to_increase_bits = 0;
	uint32_t coin = 0;

	while (remain_step > 0) {

		//randomly pick a column to decrease bits
		column_to_decrease_bits = floor(Dice() * transformed_split.size());

		assert(column_to_decrease_bits >= 0);
		assert(column_to_decrease_bits < (int)transformed_split.size());

		//decide the column to increase bit
		coin = ((uint32_t)(100.0 * Dice())) % 2;
		if (coin) {
			column_to_increase_bits = column_to_decrease_bits - 1;
		} else {
			column_to_increase_bits = column_to_decrease_bits + 1;
		}

		//see if we need to add element in the head or tail
		if ((column_to_increase_bits == -1) || (column_to_increase_bits == (int)transformed_split.size())) {
			transformed_split.push_back(0);
			column_to_increase_bits = transformed_split.size() - 1;
		}

		//shift bit only when it is valid
		if ((transformed_split[column_to_decrease_bits] > 0) &&
				(transformed_split[column_to_increase_bits] < 64)) {
			transformed_split[column_to_decrease_bits]--;
			transformed_split[column_to_increase_bits]++;
		}

		remain_step--;
	}

	//remove the element with zero
	std::vector<uint32_t> neighbor_split;
	uint32_t sum = 0;
	for (size_t i = 0; i < transformed_split.size(); ++i) {
		if (0 != transformed_split[i]) {
			neighbor_split.push_back(transformed_split[i]);
			sum += transformed_split[i];
		}
	}
	assert(sum == col_width_sum_);

	searchedPlans.insert(neighbor_split);

	return neighbor_split;
}

void PlanEnumRRS::GetStatistics(const PointType &split,
			uint64_t* cardinality, uint32_t assembleStatistics[][3], uint32_t shiftStatistics[]) {

	/** stitching statistics, for costing stitching/shifting overhead
	 *
	 * for matrix assembleStatistics, x axis is ToColumns, y axis is FromColumns
	 */

	const uint32_t new_col_num = split.size();

	/**
	 * Step a) Find the shifting operations for each new column
	 */
	std::vector < construct_entry_t > constructRecord[new_col_num];
	uint32_t curOldColIdx = 0;
	uint32_t curNewColIdx = 0;
	uint32_t start = 1; //[start, end] is from the MSB to LSB, as the location pointer for current old column
	uint32_t minVal = 0;
	uint32_t curOldColWidth = bw_[columnOrder_[curOldColIdx]];	//important!!!
	uint32_t curNewColWidth = split[curNewColIdx];
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
		entry.finalIdx = columnOrder_[curOldColIdx];	//the index is transferred by bestOrdering
		entry.start = start;
		entry.end = start + minVal - 1;
		constructRecord[curNewColIdx].push_back(entry);
		start += minVal;

		/*********************update stitching/shifting statistics************************/
		if (split[curNewColIdx] <= 16) {
			shiftStatistics[0]++;

			if (bw_[columnOrder_[curOldColIdx]] <= 16) {
				assembleStatistics[0][0]++;
			} else if (bw_[columnOrder_[curOldColIdx]] <= 32) {
				assembleStatistics[1][0]++;
			} else {
				assembleStatistics[2][0]++;
			}

		} else if (split[curNewColIdx] <= 32) {
			shiftStatistics[1]++;

			if (bw_[columnOrder_[curOldColIdx]] <= 16) {
				assembleStatistics[0][1]++;
			} else if (bw_[columnOrder_[curOldColIdx]] <= 32) {
				assembleStatistics[1][1]++;
			} else {
				assembleStatistics[2][1]++;
			}
		} else {	//	<=64
			shiftStatistics[2]++;

			if (bw_[columnOrder_[curOldColIdx]] <= 16) {
				assembleStatistics[0][2]++;
			} else if (bw_[columnOrder_[curOldColIdx]] <= 32) {
				assembleStatistics[1][2]++;
			} else {
				assembleStatistics[2][2]++;
			}
		}
		/*********************************************************************************/

		if (0 == curOldColWidth) {	//an old column is consumed up
			curOldColIdx++;
			if (curOldColIdx < num_columns_) {
				curOldColWidth = bw_[columnOrder_[curOldColIdx]];
			}
			start = 1;
		}

		if (0 == curNewColWidth) {	//a new column is consumed up
			curNewColIdx++;
			if (curNewColIdx < new_col_num) {
				curNewColWidth = split[curNewColIdx];
			}
		}
	}
	assert(curOldColIdx == num_columns_);
	assert(curNewColIdx == new_col_num);

	/**
	 * Step b) Estimate the cardinality info from cardinality of base columns
	 *
	 * use this info: uint32_t *bw = compose_params_.bitwidth;
	 */
	const uint32_t origin_col_num = num_columns_;
	//uint64_t base_col_card[origin_col_num];
	double card_per_bit[origin_col_num];

	/**calculate the bit info carried by each bit in each base column**/
	for (uint32_t colIdx = 0; colIdx < origin_col_num; ++colIdx) {
		//base_col_card[colIdx] = column_values_[columnOrder[colIdx]]->GetCardinality();
		card_per_bit[colIdx] = pow((double)(base_col_card_[colIdx]), 1.0/(double)(bw_[colIdx]));
	}

	/** use this:
	 * std::vector < construct_entry_t > constructRecord[new_col_num];
	 * to *assemble* the cardinality info
	 */
	std::vector < construct_entry_t > records;
	uint32_t finalIdx;
	uint32_t end;
	uint32_t portion_width;
	for (uint32_t new_colIdx = 0; new_colIdx < (new_col_num -1); ++new_colIdx) {
		records = constructRecord[new_colIdx];
		assert(!records.empty());

		for (uint32_t entryIdx = 0; entryIdx < records.size(); ++entryIdx) {
			finalIdx = records[entryIdx].finalIdx;
			start = records[entryIdx].start;
			end = records[entryIdx].end;
			portion_width = end - start + 1;

			cardinality[new_colIdx] *= pow(card_per_bit[finalIdx], (double)portion_width);
		}
	}
}



}	// namespace
