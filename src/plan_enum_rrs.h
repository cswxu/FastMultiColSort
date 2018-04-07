#ifndef PLAN_ENUM_RRS_H
#define PLAN_ENUM_RRS_H
/*
 * implement a plan enumerator with Recursive Random Search
 *
 *  Created on: 30 Oct 2015
 *      Author: xwj
 */

#include	"rrs.h"
#include	<cassert>
#include	<cmath>
#include	<functional>
#include	<vector>
#include	<set>
#include	"hybrid_timer.h"

namespace multiAttrSort{

class PlanEnumRRS : public RRS<std::vector<uint32_t>, double> {

public:
	typedef std::vector<uint32_t> PointType;
	typedef double	CostType;

	PlanEnumRRS(uint64_t num_rows, uint32_t col_width_sum, uint32_t max_round,
			std::vector<uint32_t> column_order, uint32_t num_columns,
			std::vector<uint32_t> bw, std::vector<uint32_t> base_col_card, double time_limit):
		RRS<std::vector<uint32_t>, double>(0.99,0.1,0.99,0.8,0.5,0.001),
			num_rows_(num_rows), col_width_sum_(col_width_sum), max_round_(max_round),
			columnOrder_(column_order), num_columns_(num_columns),
			bw_(bw), base_col_card_(base_col_card), time_limit_(time_limit){

		timer.Start();
	}

	PointType GetRandomSample() override;

	/**
	 * Given a split, return its estimated cost
	 */
	CostType  GetCost(const PointType& point) override;

	bool StopCritera() override;

	PointType GetRandomNeighbor(const PointType& point, double ro) override;

	std::set<std::vector<uint32_t> > ObtainPlans() {return searchedPlans;}

protected:



	size_t counter_ = 0;

	const uint64_t num_rows_;
	const uint32_t col_width_sum_;
	const uint32_t max_round_;

	const std::vector<uint32_t> columnOrder_;
	const uint32_t num_columns_;				//number of base columns
	const std::vector<uint32_t> bw_;				//bit width for base columns
	const std::vector<uint32_t> base_col_card_;	//cardinality for base columns

	const double time_limit_;

	std::set<std::vector<uint32_t> > searchedPlans;	//to store plans
	HybridTimer timer;

private:
	/**functions used in GetCost**/
	/**Copied from getCardinalityFromStatistics()**/
	void GetStatistics(const PointType &split,
			uint64_t* cardinality, uint32_t assembleStatistics[][3], uint32_t shiftStatistics[]);

	/**Copied from calculate_single_column_cost()**/
	double calculate_single_column_cost(const double groupsz, const uint32_t banksz_type);

};




}	// namespace

#endif
