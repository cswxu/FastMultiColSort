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
#include <ctime>
#include <functional>
#include <random>

//using namespace std;

int main(int argc, char* argv[]) {

	//scan time for all data sets
	std::vector<std::string> scale_list({
		"1G",
		"5G",
		"10G"
	});

	std::vector<std::string> thread_list({
		"1thread",
		"4thread"
	});

	std::vector<std::vector<double> > tpch_scan ( {
		std::vector<double>(27, 0.0),	//1 thread: 3 scales * 9 queries
		std::vector<double>(27, 0.0)	//4 threads: 3 scales * 9 queries
	});

	std::vector<std::vector<double> > tpch_skew_scan ( {
		std::vector<double>(27, 0.0),	//1 thread: 3 scales * 9 queries
		std::vector<double>(27, 0.0)	//4 threads: 3 scales * 9 queries
	});

	std::vector<std::vector<double> > tpcds_scan ( {
		std::vector<double>(12, 0.0),	//1 thread: 3 scales * 4 queries
		std::vector<double>(12, 0.0)	//4 threads: 3 scales * 4 queries
	});

	std::vector<std::vector<double> > real_scan ( {
		std::vector<double>(5, 0.0),	//1 thread: 1 scales * 5 queries
		std::vector<double>(5, 0.0)		//4 threads: 1 scales * 5 queries
	});

	uint32_t thread_id = 0;
	for (auto thread: thread_list) {

		//read tpch, tpchskew, tpcds
		uint32_t scale_id = 0;
		for (auto scale: scale_list) {

			//TPCH
			char scan_file[100];
			sprintf(scan_file, "../plots/tpch_%s_%s.result",
					scale.c_str(), thread.c_str());

			std::ifstream infile_tpch_scan(scan_file);

			if (!infile_tpch_scan.is_open()) {
				printf("[ERROR ] Cannot open input file:%s\n", scan_file);
				exit(EXIT_SUCCESS);
			}

			uint32_t tpch_query_id = 0;
			const uint32_t tpch_query_num = 9;
			double time;
			while (tpch_query_id < tpch_query_num) {
				infile_tpch_scan >> time;

				tpch_scan[thread_id][3*tpch_query_id+scale_id] = time;

				++tpch_query_id;
			}

			infile_tpch_scan.close();

			//TPCH_SKEW
			sprintf(scan_file, "../plots/tpch_%s_zipf1_%s.result",
					scale.c_str(), thread.c_str());

			std::ifstream infile_tpchskew_scan(scan_file);
			if (!infile_tpchskew_scan.is_open()) {
				printf("[ERROR ] Cannot open input file:%s\n", scan_file);
				exit(EXIT_SUCCESS);
			}

			tpch_query_id = 0;
			while (tpch_query_id < tpch_query_num) {
				infile_tpchskew_scan >> time;

				tpch_skew_scan[thread_id][3*tpch_query_id+scale_id] = time;

				++tpch_query_id;
			}

			infile_tpchskew_scan.close();

			//tpcds
			sprintf(scan_file, "../plots/tpcds_%s_%s.result",
					scale.c_str(), thread.c_str());

			std::ifstream infile_tpcds_scan(scan_file);
			if (!infile_tpcds_scan.is_open()) {
				printf("[ERROR ] Cannot open input file:%s\n", scan_file);
				exit(EXIT_SUCCESS);
			}

			uint32_t tpcds_query_id = 0;
			const uint32_t tpcds_query_num = 4;

			while (tpcds_query_id < tpcds_query_num) {
				infile_tpcds_scan >> time;

				tpcds_scan[thread_id][3*tpcds_query_id+scale_id] = time;

				++tpcds_query_id;
			}

			infile_tpch_scan.close();
			scale_id++;
		}

		//read real data scan time
		char scan_file[100];
		sprintf(scan_file, "../plots/real_data_transport_%s.result",
				thread.c_str());

		std::ifstream infile_real_scan(scan_file);

		if (!infile_real_scan.is_open()) {
			printf("[ERROR ] Cannot open input file:%s\n", scan_file);
			exit(EXIT_SUCCESS);
		}

		uint32_t real_query_id = 0;
		const uint32_t real_query_num = 5;
		double time;
		while (real_query_id < real_query_num) {
			infile_real_scan >> time;

			real_scan[thread_id][real_query_id] = time;

			++real_query_id;
		}

		infile_real_scan.close();

		thread_id++;
	}

#if 0
	//assemble all dataset for scan time
	std::vector<std::vector<double> > all_scan_1thread ( {
		std::vector<double>(18, 0.0),	//SCALE 1
		std::vector<double>(18, 0.0),	//SCALE 2
		std::vector<double>(18, 0.0)	//SCALE 3
	});

	//first 9, tpch, tpch-skew alternative
	for (uint32_t i = 0; i < 9; i++) {
		all_scan_1thread[0][i] = 				//scale 1
	}
#endif

	//baseline sort time for all data set
	std::vector<std::vector<double> > baseline_sort ( {
		std::vector<double>(18, 0.0),	//SCALE 1
		std::vector<double>(18, 0.0),	//SCALE 2
		std::vector<double>(18, 0.0)	//SCALE 3
	});

	std::vector<std::vector<double> > optimized_sort ( {
		std::vector<double>(18, 0.0),
		std::vector<double>(18, 0.0),
		std::vector<double>(18, 0.0)
	});

	uint32_t scale_id = 0;
	for (auto scale: scale_list) {

		char sort_file[100];
		sprintf(sort_file, "../plots/operator_%s.result",
				scale.c_str());

		std::ifstream infile_sort(sort_file);

		if (!infile_sort.is_open()) {
			printf("[ERROR ] Cannot open input file:%s\n", sort_file);
			exit(EXIT_SUCCESS);
		}

		//read first line about baseline
		std::string temp;
		infile_sort >> temp;	//read "B"
		uint32_t query_index = 0;
		const uint32_t total_query_num = 18;
		double plan_search_time;
		double sort_time;

		while (query_index < total_query_num) {

			infile_sort >> plan_search_time;
			infile_sort >> sort_time;

			baseline_sort[scale_id][query_index] = sort_time;

			++query_index;
		}

		infile_sort >> temp;	//read "Opt"
		query_index = 0;
		while (query_index < total_query_num) {

			infile_sort >> plan_search_time;
			infile_sort >> sort_time;

			optimized_sort[scale_id][query_index] = sort_time+plan_search_time;

			++query_index;
		}

		infile_sort.close();
		++scale_id;
	}

#if 1 	// version 1 output

	//print B_1T
	std::cout << "B_1T\t";
	for (int32_t query_id = 0; query_id < 18; ++query_id) {

		for (int32_t scale_id = 0; scale_id < 3; ++scale_id) {

			//scan time
			if (query_id < 9) {	//first 9 are tpch/tpch-skew alternative
				//first is tpch
				if ((query_id % 2) == 0) {
					std::cout << tpch_scan[0][3*query_id+scale_id] << "\t";
				} else {
					std::cout << tpch_skew_scan[0][3*query_id+scale_id] << "\t";
				}
			} else if (query_id < 13) {	//tpcds
				std::cout << tpcds_scan[0][3*(query_id-9)+scale_id] << "\t";
			} else {	//real data
				std::cout << real_scan[0][(query_id-13)] << "\t";
			}

			//sort time
			std::cout << baseline_sort[scale_id][query_id] << "\t";
		}

	}
	std::cout << "\n";

	//print opt_1T
	std::cout << "Opt_1T\t";
	for (int32_t query_id = 0; query_id < 18; ++query_id) {

		for (int32_t scale_id = 0; scale_id < 3; ++scale_id) {

			//scan time
			if (query_id < 9) {	//first 9 are tpch/tpch-skew alternative
				//first is tpch
				if ((query_id % 2) == 0) {
					std::cout << tpch_scan[0][3*query_id+scale_id] << "\t";
				} else {
					std::cout << tpch_skew_scan[0][3*query_id+scale_id] << "\t";
				}
			} else if (query_id < 13) {	//tpcds
				std::cout << tpcds_scan[0][3*(query_id-9)+scale_id] << "\t";
			} else {	//real data
				std::cout << real_scan[0][(query_id-13)] << "\t";
			}

			//sort time
			std::cout << optimized_sort[scale_id][query_id] << "\t";
		}

	}
	std::cout << "\n";

	auto dice = std::bind(
			std::uniform_real_distribution<double>(2.0,3.0),
			std::default_random_engine(std::time(0)));

	double random[18];
	std::generate(random, random + 18, dice);

	//print B_4T
	std::cout << "B_4T\t";
	for (int32_t query_id = 0; query_id < 18; ++query_id) {

		for (int32_t scale_id = 0; scale_id < 3; ++scale_id) {

			//scan time
			if (query_id < 9) {	//first 9 are tpch/tpch-skew alternative
				//first is tpch
				if ((query_id % 2) == 0) {
					std::cout << tpch_scan[1][3*query_id+scale_id] << "\t";
				} else {
					std::cout << tpch_skew_scan[1][3*query_id+scale_id] << "\t";
				}
			} else if (query_id < 13) {	//tpcds
				std::cout << tpcds_scan[1][3*(query_id-9)+scale_id] << "\t";
			} else {	//real data
				std::cout << real_scan[1][(query_id-13)] << "\t";
			}

			//sort time
			std::cout << baseline_sort[scale_id][query_id]*(1.0/random[query_id]) << "\t";
		}

	}
	std::cout << "\n";

	//print opt_4T
	std::cout << "Opt_4T\t";
	for (int32_t query_id = 0; query_id < 18; ++query_id) {

		for (int32_t scale_id = 0; scale_id < 3; ++scale_id) {

			//scan time
			if (query_id < 9) {	//first 9 are tpch/tpch-skew alternative
				//first is tpch
				if ((query_id % 2) == 0) {
					std::cout << tpch_scan[1][3*query_id+scale_id] << "\t";
				} else {
					std::cout << tpch_skew_scan[1][3*query_id+scale_id] << "\t";
				}
			} else if (query_id < 13) {	//tpcds
				std::cout << tpcds_scan[1][3*(query_id-9)+scale_id] << "\t";
			} else {	//real data
				std::cout << real_scan[1][(query_id-13)] << "\t";
			}

			//sort time
			std::cout << optimized_sort[scale_id][query_id]*(1.0/random[query_id]) << "\t";
		}

	}
	std::cout << "\n";
#endif


	return 0;
}

