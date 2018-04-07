#include <iostream>
#include <cstdlib>
#include <climits>
#include <string>
#include <string.h>
#include <vector>
#include <sstream>
#include <cstdio>
#include <algorithm>
#include <fstream>
#include <set>

//using namespace std;

int main(int argc, char* argv[]) {

	std::string scale_factor = std::string(argv[1]);//scale factor

	//baseline time for tpch, tpchskew1, tpchskew2
	std::vector<std::vector<double> > baseline_tpch ( {
		std::vector<double>(18, 0.0),
		std::vector<double>(18, 0.0)//,
		//std::vector<double>(18, 0.0)
	});

	std::vector<std::vector<double> > plan_search_tpch ( {
		std::vector<double>(18, 0.0),
		std::vector<double>(18, 0.0)//,
		//std::vector<double>(18, 0.0)
	});

	std::vector<std::vector<double> > optimized_tpch ( {
		std::vector<double>(18, 0.0),
		std::vector<double>(18, 0.0)//,
		//std::vector<double>(18, 0.0)
	});

	//index for tpch/skew queries
	std::vector<uint32_t> tpch_indices ({
		1, 2, 3, 7, 9, 10, 13, 16, 18
	});

	std::vector<uint32_t> tpch_groupby_indices ({
		1, 3, 7, 9, 10, 16, 18
	});

	std::vector<uint32_t> tpch_orderby_indices ({
		1, 2, 3, 7, 9, 13, 16, 18
	});

	//read all tpch, tpchskew complete results
	std::vector<std::string> tpch_surfix ({
		"",
		"_skew1"//,
		//"_skew2"
	});

	uint32_t dataset_index = 0;
	//iterate through three sets of result
	for (auto dataset: tpch_surfix) {

		//read result about groupby queries
		for (auto groupby_index : tpch_groupby_indices) {

			char result_file[100];
			sprintf(result_file, "../plots/result_Q%d_%s_GROUPBY%s_complete.data",
					groupby_index, scale_factor.c_str(), dataset.c_str());

			std::ifstream infile_groupby(result_file);

			if (!infile_groupby.is_open()) {
				printf("[ERROR ] Cannot open input file:%s\n", result_file);
				exit(EXIT_SUCCESS);
			}

			std::string line;
			//int index;
			char comma;
			double baseline_search;
			double baseline_sort;
			double opt_search;
			double opt_sort;

			//read line one
			std::getline(infile_groupby, line);
			std::stringstream parseline(line);

			parseline >> baseline_search >> comma >> baseline_sort;

			baseline_tpch[dataset_index][groupby_index] += baseline_sort;

			//read line two
			std::getline(infile_groupby, line);
			parseline.str(line);

			parseline >> opt_search >> comma >> opt_sort;

			plan_search_tpch[dataset_index][groupby_index] += opt_search;
			optimized_tpch[dataset_index][groupby_index] += opt_sort;
		}

		//read result about orderby queries
		for (auto orderby_index : tpch_orderby_indices) {

			char result_file[100];
			sprintf(result_file, "../plots/result_Q%d_%s_ORDERBY%s_complete.data",
					orderby_index, scale_factor.c_str(), dataset.c_str());

			std::ifstream infile_orderby(result_file);

			if (!infile_orderby.is_open()) {
				printf("[ERROR ] Cannot open input file:%s\n", result_file);
				exit(EXIT_SUCCESS);
			}

			std::string line;
			//int index;
			char comma;
			double baseline_search;
			double baseline_sort;
			double opt_search;
			double opt_sort;

			//read line one
			std::getline(infile_orderby, line);
			std::stringstream parseline(line);

			parseline >> baseline_search >> comma >> baseline_sort;

			baseline_tpch[dataset_index][orderby_index] += baseline_sort;

			//read line two
			std::getline(infile_orderby, line);
			parseline.str(line);

			parseline >> opt_search >> comma >> opt_sort;

			plan_search_tpch[dataset_index][orderby_index] += opt_search;
			optimized_tpch[dataset_index][orderby_index] += opt_sort;
		}

		dataset_index++;
	}

	//baseline time for tpch-ds
	std::vector<double> baseline_tpcds(4, 0.0);
	std::vector<double> plan_search_tpcds(4, 0.0);
	std::vector<double> optimized_tpcds(4, 0.0);

	//index for tpcds queries
	std::vector<uint32_t> tpcds_indices ({
		36, 67, 70, 86
	});

	std::vector<std::string> tpcds_types ({
		"GROUPBY",
		"PARTITIONBY",
		"ORDERBY"
	});

	uint32_t query_index = 0;

	for (auto tpcds_id : tpcds_indices) {

		for (auto type : tpcds_types) {

			char result_file[100];
			sprintf(result_file, "../plots/result_TPCDS_Q%d_%s_%s_complete.data",
					tpcds_id, scale_factor.c_str(), type.c_str());

			std::ifstream infile_tpcds(result_file);

			if (!infile_tpcds.is_open()) {
				printf("[ERROR ] Cannot open input file:%s\n", result_file);
				exit(EXIT_SUCCESS);
			}

			std::string line;
			//int index;
			char comma;
			double baseline_search;
			double baseline_sort;
			double opt_search;
			double opt_sort;

			//read line one
			std::getline(infile_tpcds, line);
			std::stringstream parseline(line);

			parseline >> baseline_search >> comma >> baseline_sort;

			baseline_tpcds[query_index] += baseline_sort;

			//read line two
			std::getline(infile_tpcds, line);
			parseline.str(line);

			parseline >> opt_search >> comma >> opt_sort;

			plan_search_tpcds[query_index] += opt_search;
			optimized_tpcds[query_index] += opt_sort;

		}
		query_index++;
	}

	//baseline time for real
	std::vector<double> baseline_real(5, 0.0);
	std::vector<double> plan_search_real(5, 0.0);
	std::vector<double> optimized_real(5, 0.0);

	//index for real queries
	std::vector<std::string> real_queries_files ({
		"../plots/result_realdata_q1_ORDERBY_complete.data",
		"../plots/result_realdata_q2_PARTITIONBY_complete.data",
		"../plots/result_realdata_q3_GROUPBY_complete.data",
		"../plots/result_realdata_q4_GROUPBY_complete.data",
		"../plots/result_realdata_q5_PARTITIONBY_complete.data"
	});

	query_index = 0;
	for (auto real_query : real_queries_files) {

		std::ifstream infile_real(real_query.c_str());

		if (!infile_real.is_open()) {
			printf("[ERROR ] Cannot open input file:%s\n", real_query.c_str());
			exit(EXIT_SUCCESS);
		}

		std::string line;
		//int index;
		char comma;
		double baseline_search;
		double baseline_sort;
		double opt_search;
		double opt_sort;

		//read line one
		std::getline(infile_real, line);
		std::stringstream parseline(line);

		parseline >> baseline_search >> comma >> baseline_sort;

		baseline_real[query_index] += baseline_sort;

		//read line two
		std::getline(infile_real, line);
		parseline.str(line);

		parseline >> opt_search >> comma >> opt_sort;

		plan_search_real[query_index] += opt_search;
		optimized_real[query_index] += opt_sort;

		query_index++;
	}

	//output the data file for gnuplot
	std::cout << "B\t";

	//tpch/tpch-skew alternative
	bool is_tpch_skew = true;

	for (auto tpch_query : tpch_indices) {
		is_tpch_skew = !is_tpch_skew;

		if (is_tpch_skew) {
			std::cout << "0.0\t" << baseline_tpch[1][tpch_query] << "\t";
		} else {
			std::cout << "0.0\t" << baseline_tpch[0][tpch_query] << "\t";
		}
	}

	//tpcds
	for (uint32_t tpcds_id = 0; tpcds_id < 4; ++tpcds_id) {
		std::cout << "0.0\t" << baseline_tpcds[tpcds_id] << "\t";
	}

	//real data
	for (uint32_t real_id = 0; real_id < 5; ++real_id) {
		std::cout << "0.0\t" << baseline_real[real_id] << "\t";
	}

	std::cout << "\nOpt\t";

	is_tpch_skew = true;

	for (auto tpch_query : tpch_indices) {
		is_tpch_skew = !is_tpch_skew;

		if (is_tpch_skew) {
			std::cout << plan_search_tpch[1][tpch_query] << "\t"
					<< optimized_tpch[1][tpch_query] << "\t";
		} else {
			std::cout << plan_search_tpch[0][tpch_query] << "\t"
					<< optimized_tpch[0][tpch_query] << "\t";
		}
	}

	//tpcds
	for (uint32_t tpcds_id = 0; tpcds_id < 4; ++tpcds_id) {
		std::cout << plan_search_tpcds[tpcds_id] << "\t"
				<< optimized_tpcds[tpcds_id] << "\t";
	}

	//real
	for (uint32_t real_id = 0; real_id < 5; ++real_id) {
		std::cout << plan_search_real[real_id] << "\t"
				<< optimized_real[real_id] << "\t";
	}

	std::cout << "\n";

	return 0;
}

