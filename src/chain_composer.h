/*******************************************************************************
 * Copyright (c) 2016
 * The Hong Kong Polytechnic University, Database Group
 *
 * Author: Wenjian Xu (cswxu AT comp DOT polyu.edu.hk)
 *
 * See file LICENSE.md for details.
 *******************************************************************************/

#ifndef CHAIN_COMPOSER_H
#define CHAIN_COMPOSER_H

#include "composer.h"
#include "types.h"
#include "common.h"
#include "mergesort.h"
#include <cstring>



namespace multiAttrSort{

class ChainComposer: public Composer {
public:

	ChainComposer(setting_t setting, uint32_t ncolumns, uint64_t nrows,
			compose_params_t params, Column **column_values):
		Composer(setting, ncolumns, nrows, params, column_values) {
	}

	ChainComposer() {}

	virtual ~ChainComposer(){}

protected:
	void SortAllColumns_Pack() override;
	void SortAllColumns_Nonpack() override;
	void ExhaustiveSearch() override;
	void TwoRoundsExhaustive() override;
	void RunAnInstance(std::string intance) override;
	void CostModelComparison() override;

	//for hashing
	void RunAnInstance_hash(std::string intance) override;
	void TwoRoundsExhaustive_hash() override;
	void HashAllColumns() override;
private:

};

}


#endif
