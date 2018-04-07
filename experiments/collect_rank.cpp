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

struct pair_t {
	int index;
	double cost;
};

struct Compare {
	bool operator()(pair_t i, pair_t j) {return i.cost < j.cost;}
} mycompare;

int main(int argc, char* argv[]) {

	//parse the input filename
	std::string file_prefix = std::string(argv[1]);

	//read the exhaustive search actual run file name: file_predix.data
	char actual_file[100];
	sprintf(actual_file, "%s.data", file_prefix.c_str());

	std::ifstream infile_actual(actual_file);

	if (!infile_actual.is_open()) {
		printf("[ERROR ] Cannot open input file:%s\n", actual_file);
		exit(EXIT_SUCCESS);
	}

	std::string line;
	//int index;
	char comma;
	//double cost;

	//char* token1;

	std::vector<pair_t> all_actual_costs;
	while(std::getline(infile_actual, line)) {
		std::stringstream parseline(line);
		pair_t each_pair;

		parseline >> each_pair.index;
		parseline >> comma;
		parseline >> each_pair.cost;

		all_actual_costs.push_back(each_pair);
	}

	//read the exhaustive search costmodel file name: file_predix_costmodel.data
	char model_file[100];
	sprintf(model_file, "%s_costmodel.data", file_prefix.c_str());

	std::ifstream infile_model(model_file);

	if (!infile_model.is_open()) {
		printf("[ERROR ] Cannot open input file:%s\n", model_file);
		exit(EXIT_SUCCESS);
	}

	int index = 0;
	double estimated_cost;
	double error;
	double avg_error = 0;


	while(std::getline(infile_model, line)) {
		std::stringstream parseline(line);

		parseline >> index;
		parseline >> comma;
		parseline >> estimated_cost;

		error = abs(estimated_cost - all_actual_costs[index].cost)/all_actual_costs[index].cost;
		avg_error += error;
	}

	std::cout << avg_error/((double)index+1) << std::endl;

	sort(all_actual_costs.begin(), all_actual_costs.end(), mycompare);
#if 1
	//read greedy plans
	char greedy_file[100];
	sprintf(greedy_file, "%s_greedy.data", file_prefix.c_str());

	std::ifstream infile_greedy(greedy_file);

	if (!infile_greedy.is_open()) {
		printf("[ERROR ] Cannot open input file:%s\n", greedy_file);
		exit(EXIT_SUCCESS);
	}

	index = 0;
	std::set<int> all_greedy;
	while(std::getline(infile_greedy, line)) {
		std::stringstream parseline(line);

		parseline >> index;

		all_greedy.insert(index);
	}

	//read random plans
	char random_file[100];
	sprintf(random_file, "%s_random.data", file_prefix.c_str());

	std::ifstream infile_random(random_file);

	if (!infile_random.is_open()) {
		printf("[ERROR ] Cannot open input file:%s\n", random_file);
		exit(EXIT_SUCCESS);
	}

	std::set<int> all_random;
	while(std::getline(infile_random, line)) {
		std::stringstream parseline(line);

		parseline >> index;

		all_random.insert(index);
	}

	//find rank and percentile for greedy and random
	int greedy_rank = 0;
	double greedy_percentile = 0.0;
	for (unsigned int i = 0; i < all_actual_costs.size(); ++i) {
		if (all_greedy.find(all_actual_costs[i].index) != all_greedy.end()) {
			greedy_rank = i + 1;
			greedy_percentile = 1.00 - (double)(greedy_rank-1) / (double)all_actual_costs.size();
			break;
		}
	}

	int random_rank = 0;
	double random_percentile = 0.0;
	for (unsigned int j = 0; j < all_actual_costs.size(); ++j) {
		if (all_random.find(all_actual_costs[j].index) != all_random.end()) {
			random_rank = j + 1;
			random_percentile = 1.00 - (double)(random_rank-1) / (double)all_actual_costs.size();
			break;
		}
	}

	std::cout << file_prefix << "\t" << greedy_rank << "\t" << greedy_percentile
			<< "\t" << random_rank << "\t" << random_percentile << std::endl;
#endif
	return 0;
}

