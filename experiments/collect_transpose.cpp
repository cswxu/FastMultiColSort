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

	//parameter: number of queries involved
	unsigned int query_num = atoi(argv[1]);

	//read original format
	std::string input_filename = std::string(argv[2]);//scale factor

	std::ifstream infile(input_filename.c_str());

	if (!infile.is_open()) {
		printf("[ERROR ] Cannot open input file:%s\n", input_filename.c_str());
		exit(EXIT_SUCCESS);
	}

	const auto num_columns = 6 * query_num;
	double input_matrix[4][num_columns];

	unsigned int method_id = 0;
	std::string temp;
	for (;method_id < 4; method_id++) {
		infile >> temp;	//the begin dummy

		for (unsigned int i = 0; i < num_columns; ++i) {
			infile >> input_matrix[method_id][i];
		}
	}

	//transpose the input
	unsigned int scale_num = 3;
	unsigned int thread_style = 2;

	for (unsigned int style_id = 0; style_id < thread_style; ++style_id) {

		for (unsigned int scale_id = 0; scale_id < scale_num; ++scale_id) {

			std::cout << (style_id * 3 + scale_id + 1) << "\t";

			for (unsigned int query_id = 0; query_id < query_num; ++query_id) {
#if 1
				method_id = 0;
				for (; method_id < 2; ++method_id) {
					std::cout << input_matrix[style_id*2+method_id][6*query_id + 2*scale_id] << "\t"
							<< input_matrix[style_id*2+method_id][6*query_id + 2*scale_id + 1] << "\t";
				}
#endif
			}

			std::cout << std::endl;
		}
	}



	return 0;
}

