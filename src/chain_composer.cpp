#include "chain_composer.h"
#include "mergesort.h"
#include <algorithm>


namespace multiAttrSort{

void ChainComposer::SortAllColumns_Pack() {

}

void ChainComposer::SortAllColumns_Nonpack() {
	multiRoundSorting();
}

void ChainComposer::HashAllColumns() {

	if (compose_params_.hash_type == 1) {	//1: use reordering
		multiRoundHashing_reordered();
	} else if (compose_params_.hash_type == 2){	//2L: partitioned reordering
		multiRoundHashing_partitioned_reordered();
	} else {	//non-reordered
		multiRoundHashing_non_reordered();
	}
}

void ChainComposer::ExhaustiveSearch() {

}

void ChainComposer::TwoRoundsExhaustive() {

}
void ChainComposer::RunAnInstance(std::string intance) {

}

void ChainComposer::CostModelComparison() {

}

void ChainComposer::RunAnInstance_hash(std::string intance) {

}

void ChainComposer::TwoRoundsExhaustive_hash() {

}

}
