/*******************************************************************************
 * Copyright (c) 2016
 * The Hong Kong Polytechnic University, Database Group
 *
 * Author: Wenjian Xu (cswxu AT comp DOT polyu.edu.hk)
 *
 * See file LICENSE.md for details.
 *******************************************************************************/

#ifndef STITCH_COMPOSER_H
#define STITCH_COMPOSER_H

#include "composer.h"
#include "types.h"
#include "common.h"
#include <cstring>
#include <set>


namespace multiAttrSort{

class StitchComposer: public Composer {
public:

	StitchComposer(setting_t setting, uint32_t ncolumns, uint64_t nrows,
			compose_params_t params, Column **column_values):
		Composer(setting, ncolumns, nrows, params, column_values) {

		//bestOrdering = (uint32_t *) malloc_aligned(num_columns_ * sizeof(uint32_t));

		/**Initialize the cardinality and cardinality_per_bit information for reuse**/
		//const uint32_t origin_col_num = num_columns_;
		base_col_card =(uint64_t *) malloc_aligned(num_columns_ * sizeof(uint64_t));
		card_per_bit = (double *) malloc_aligned(num_columns_ * sizeof(double));

		double all_multiply = 1.0;
		/**calculate the bit info carried by each bit in each base column**/
		for (uint32_t colIdx = 0; colIdx < num_columns_; ++colIdx) {
			base_col_card[colIdx] = column_values_[colIdx]->GetCardinality();
			all_multiply *= base_col_card[colIdx];

			card_per_bit[colIdx] = pow((double)(base_col_card[colIdx]), 1.0/(double)(compose_params_.bitwidth[colIdx]));
#if 0
			std::cout << "base cardi for column " << colIdx << ": " << base_col_card[colIdx] << std::endl;
			std::cout << "card per bit for column " << colIdx << ": " << card_per_bit[colIdx] << std::endl;
#endif
		}

		/**calculate the shrink factor**/
		shrink_factor = pow((double)num_rows_/all_multiply, 1.0/(double)(num_columns_-1));
		//shrink_factor = (double)num_rows_/all_multiply;
	}

	virtual ~StitchComposer(){
		//free(bestOrdering);
		free(base_col_card);
		free(card_per_bit);
	}

	/**
	 * @return the estimated cycles to sort nitems with bank-16 size
	 */
	static double estimate_mergesort_uint16(double nitems);
	/**
	 * @return the estimated cycles to sort nitems with bank-32 size
	 */
	static double estimate_mergesort_uint32(double nitems);
	/**
	 * @return the estimated cycles to sort nitems with bank-64 size
	 */
	static double estimate_mergesort_uint64(double nitems);

	/**
	 * @return the estimated stitching overhead (# of cycles)
	 */
	static double estimate_stitch(uint32_t assembleStatistics[][3], uint32_t shiftStatistics[], uint32_t nitems);

	/**
	 * @return the estimated tuple reconstruction overhead (# of cycles)
	 */
	static double estimate_tuple_reconstruct(const uint32_t *bytewidths, const uint32_t column_num, const uint64_t rows_num);

protected:
	void SortAllColumns_Pack() override;
	void SortAllColumns_Nonpack() override;
	void ExhaustiveSearch() override;
	void TwoRoundsExhaustive() override;
	void RunAnInstance(std::string intance) override;
	void CostModelComparison() override;

	//for hashing
	void RunAnInstance_hash(std::string instance) override;
	void TwoRoundsExhaustive_hash() override;
	void HashAllColumns() override;
private:

	uint64_t * base_col_card;
	double* card_per_bit;
	double shrink_factor;

	//struct plan_rounds_pair_t {
	//	std::vector<uint32_t> plan;
	//	uint32_t rounds;
	//};
    //struct Compare {
    //	bool operator()(plan_rounds_pair_t i, plan_rounds_pair_t j) {return i.rounds > j.rounds;}
    //} plancompare;

	//std::vector<std::vector<uint32_t> > split_enum;	//enumerate all split methods
	uint32_t partitions[64]; 	//auxiliary structure for splitting
	//uint32_t* bestOrdering;		//ordering of columns
	//std::vector<uint32_t> bestGlobalSplit;


	/** select the best split method (including the order of columns in GROUP-BY case) and run it!!**/
	void optimizeAndRun();


	//void optimizer_exhaustive();	/** exhaustive search (finding the actual optimal plan) **/
	//void optimizer_profiling();
	/**
	 * do the stitching based on <bestOrdering> and <bestGlobalSplit>
	 */
	//void stitching();
	void permutations(uint32_t * arr, uint32_t k, uint32_t m, std::vector<uint32_t *> &comb);
	void swap_val(uint32_t * a, uint32_t * b);
	/**
	 *
	 * @params targetComb:a specific order of columns
	 * @return loc: Given a specific order, estimate the location when the average group size reduces to a certain amount (e.g., 256)
	 * (Note: if loc == colwidth_sum, it means that such location does not exist)
	 */
	//uint32_t estimateThresholdLoc(uint32_t *targetComb);

	/**
	 * @return nbits_right_shifted: number of bits removed in LSB
	 */
	//template <class T>
	//uint32_t estimateDetailedLoc(uint64_t cur_accu_ngroups, uint64_t threshold_lb, uint64_t threshold_ub,
	//					T* values, uint32_t width);

	/*********************************************************************************************************************
	 * Functions to enumerate massage plans
	 */

	/**
	 * Given the max round number, Return ALL possible plans.
	 *
	 * @param n 			remained value
	 * @param curRound 		filling the integer for which round
	 * @param maxRound 		the maximum round num
	 * @param maxValue 		the largest value in each round, 64 in our scenario
	 * @param sum 			the sum of all columns
	 * @param split_enum 	plans to be stored
	 */
	void split_exhaustive_all_round(int n, int curRound, const int maxRound,
			const int maxValue, const int sum, std::vector<std::vector<uint32_t> > &split_enum);
	//void split_applyRule1(int n, int curLevel, const int maxLevel,
	//		const int maxValue, std::vector<std::vector<uint32_t> > &split_enum, bool isBroken);

    /** Given a round number, say 3. Return ALL massaging plans with that round number
     *
     * @param n  			remained value
     * @param curRound		filling the integer for which round
     * @param round_num   	how many rounds do we need to fill up
     * @param maxValue		the largest eligible value allowed for one round (i.e., 64 in our scenario)
     * @param split_enum	plans to be stored
     */
	void split_exhaustive_certain_round(int n, int curRound, const int round_num,
			const int maxValue, std::vector<std::vector<uint32_t> > &split_enum);
	/*********************************************************************************************************************/

	/*********************************************************************************************************************
	 * Functions to costing plans
	 */
	//void cost_estimation(uint32_t * columnOrder, std::vector<std::vector<uint32_t> > &split_enum,
	//		int &bestLocalSplitIdx, double &optimalLocCost, int &bestLocalSplitIdx_ac, double &optimalLocCost_ac,
	//		double & baseline_cost, double & bestSplit_es2ac_each_enum);
	//void cost_estimation_profiling(uint32_t * columnOrder, std::vector<std::vector<uint32_t> > &split_enum,
	//		int &bestLocalSplitIdx, double &optimalLocCost, int &bestLocalSplitIdx_ac, double &optimalLocCost_ac,
	//		double & baseline_cost, double & bestSplit_es2ac_each_enum);

	//==>in use of main_tpch_complete_scenario
	//void cost_estimation_plans(uint32_t * columnOrder, std::vector<std::vector<uint32_t> > &split_enum,
	//		int &optimalIdx, double &optimalCost, double& extractCardinalityNumCycles);

	void estimate_plan_collection(uint32_t * columnOrder,
			std::vector<std::vector<uint32_t> > &split_enum, uint64_t startIndex,
			std::set<std::vector<uint32_t> > &greedyPlans,
			std::set<std::vector<uint32_t> > &randomPlans);

	void estimate_plan_collection_greedy(uint32_t * columnOrder,
			std::vector<std::vector<uint32_t> > &greedy_plans,
			int & bestLocalSplitIdx, double & optimalLocCost);

	void find_greedy_plans(uint32_t round_num,
			std::vector<std::vector<uint32_t> > &validPlans, uint32_t * columnOrder);

	//==>in use of cost_estimation_greedy
	bool determineBitRange(const uint32_t banksize, int & bit_range_min,
			int & bit_range_max, const uint32_t bitOffset);

	void determineGroupszRange(const int bit_range_min, const int bit_range_max,
			double & groupsz_min, double & groupsz_max, uint32_t * columnOrder);

	double estimateGroupNum(const uint32_t bit_position, uint32_t * columnOrder);

	void determineBowlGroupsz(double & groupsz_bowl, int & bit_position_bowl,
			const uint32_t banksz_type, const double avg_groupsz, uint32_t * columnOrder);

	int determineOptBitPosition(const int bit_range_min, const int bit_range_max, const int bit_position_bowl,
						const double groupsz_max, const double groupsz_min, const double groupsz_bowl,
						uint32_t banksz_type);

	int determineOptBitPosition_v2(const int bit_range_min, const int bit_range_max,
			uint32_t banksz_type, uint32_t * columnOrder);

	double calculate_single_column_cost(const double groupsz, const uint32_t banksz_type);

	double calculate_bowl_int16(const double avg_groupsz);
	double calculate_bowl_int32(const double avg_groupsz);
	double calculate_bowl_int64(const double avg_groupsz);

	/*********************************************************************************************************************/

	/*********************************************************************************************************************
	 * Functions to running plans
	 */

	void run_plan_collection(uint32_t * columnOrder,
			std::vector<std::vector<uint32_t> > &split_enum, uint64_t startIndex);

	/**
	 * return the actual total cost for the simulation
	 */
	double run_single_plan(uint32_t *columnOrder, std::vector<uint32_t> &curSplit, double *breakdown);

	/**
	 * @return : # of cycles spent on extracting cardinalities
	 */
	void getCardinalityFromStatistics(uint32_t *columnOrder, std::vector<uint32_t> &curSplit,
			double* cardinality, uint32_t assembleStatistics[][3], uint32_t shiftStatistics[]);

	/*********************************************************************************************************************/

	//void getHeuristics(std::vector < std::vector<uint32_t> > & plans, uint32_t colwidth_sum);

    //uint32_t factor(uint32_t n) {
    //	return (1 == n) ? 1 : (n * factor(n-1));
    //}

	/**
	 * Given a round number, say 3. Return ALL possible banksize categories
	 *
	 * for the value in cate_enum, 0 stands for 16-banksize, 1 stands for 32-banksize, 2 stands for 64-banksize
	 */
	void enum_banksize_category(int curRound, const int round_num, std::vector<std::vector<uint32_t> > &cate_enum);
};

}


#endif
