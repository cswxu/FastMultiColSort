/*******************************************************************************
 * Copyright (c) 2016
 * The Hong Kong Polytechnic University, Database Group
 *
 * Author: Wenjian Xu (cswxu AT comp DOT polyu.edu.hk)
 *
 * See file LICENSE.md for details.
 *******************************************************************************/

#ifndef GENERATOR_H
#define GENERATOR_H

#include "types.h"
#include "zipf.h"
#include "common.h"
#include "column.h"
#include <sstream>
#include <string.h>
#include <string>
#include <fstream>
#include <cassert>
#include <vector>
#include <iostream>
#include <ctime>
#include <algorithm>
#include <cmath>
#include <set>


namespace multiAttrSort {


class DataGenerator{

public:
	void GenAndSaveAsFile(dataset_t *dataset);
	void ReadFromFile(Column **columns, std::string inputfile,
			uint32_t ncolumns, uint64_t nrows, uint32_t *column_widths);
	DataGenerator(){}
private:
	template <class T>
	void ReadValues(Column **columns, uint32_t colIdx, std::ifstream &input, uint64_t nrows);
};	

template <class T>
void DataGenerator::ReadValues(Column **columns, uint32_t colIdx,
			std::ifstream &input, uint64_t nrows) {
	std::set<T> container;
	T v;
	T *values = (T *)malloc_aligned(nrows * sizeof(T));
	for (uint64_t rowIdx = 0; rowIdx < nrows; rowIdx++) {
		input >> v;
		values[rowIdx] = v;

		container.insert(v);
	}
	SpecialiedColumn<T> *col = new SpecialiedColumn<T>();
	col->SetColumn(values);
	col->SetWidth(sizeof(T));
	assert((container.size() > 0) && (container.size() <= nrows));
	col->SetCardinality(container.size());

	columns[colIdx] = col;

	//clear
	std::set<T>().swap(container);
}

void DataGenerator::ReadFromFile(Column **columns, std::string inputfile,
		uint32_t ncolumns, uint64_t nrows, uint32_t *column_widths) {
	
	std::ifstream input(inputfile);
	if (!input.is_open()) {
		printf("[ERROR ] Cannot open input file:%s\n", inputfile.c_str());
		exit(EXIT_SUCCESS);
	}

	uint32_t width;
	for (uint32_t colIdx = 0; colIdx < ncolumns; colIdx++) {
		width = column_widths[colIdx];

		assert(width > 0 && width <= 64);

		//if (width <= 8) {
		//	ReadValues<uint8_t>(columns, colIdx, input, nrows);
		//} else
		if (width <= 16) {
			ReadValues<uint16_t>(columns, colIdx, input, nrows);
		} else if (width > 16 && width <= 32) {
			ReadValues<uint32_t>(columns, colIdx, input, nrows);
		} else {	//width<32 && width<=64
			ReadValues<uint64_t>(columns, colIdx, input, nrows);
		}
	}
	//the memory allocated here are released in main.cpp

	input.close();
}

void DataGenerator::GenAndSaveAsFile(dataset_t *dataset) {
	std::srand(std::time(0));
	uint64_t nrows = dataset->nrows;
	uint32_t ncolumns = dataset->ncolumns;

	double cardinality_perc = ((double) dataset->cardinality) / 100;
	/* convert long to string*/
	std::ostringstream fname_oss;
	std::string bitwidth_str = dataset->column_bitwidth;

	/* the file name to be saved*/
	fname_oss << nrows << "_" << ncolumns << "_" << dataset->column_bitwidth
			<< "_" << dataset->zipf << "_" << dataset->cardinality << "_"
			<< dataset->group_num << ".dat";
	std::string filename = fname_oss.str();
	fname_oss.flush();

	char bw_str_tmp[50];
	strcpy(bw_str_tmp, bitwidth_str.c_str());
	std::ofstream output(filename.c_str());
	char *bw_str = strtok(bw_str_tmp, "B");	//divided by character 'B'
	int bw = 0;
	uint64_t cardinality = 0;
	std::vector<uint64_t> column;

	for (uint32_t colIdx = 0; colIdx < ncolumns; ++colIdx) {
		bw = atoi(bw_str);
		std::cout << "column_width: " << bw << std::endl;
		assert(bw <= 64);
		if (bw < 64) {
			cardinality = ((1ULL << bw) - 1) * cardinality_perc;
		} else { // bw = 64

			cardinality = std::numeric_limits < uint64_t > ::max();
		}

		//random generator
		auto dice = std::bind(
				std::uniform_int_distribution < uint64_t > (0, cardinality),
				std::default_random_engine(std::time(0)));

		uint64_t eachgz = ceil((double) nrows / (double) (dataset->group_num));
		uint64_t val;
		uint32_t ngroup = 0;

		column.assign(nrows, 0);
		for (uint64_t rowIdx = 0; rowIdx < nrows; rowIdx += eachgz) {
			val = dice();
			ngroup++;
			for (uint32_t gidx = 0; (gidx < eachgz) && (rowIdx + gidx < nrows);
					++gidx) {

				column[rowIdx + gidx] = val;
			}
		}

		std::random_shuffle(column.begin(), column.end());

		/* output the column values */
		for (uint64_t rowIdx = 0; rowIdx < nrows; ++rowIdx) {
			output << column[rowIdx] << std::endl;
		}

		std::cout << "#group: " << ngroup << "agv_group_sz: "
				<< (float) nrows / (float) ngroup << std::endl;

		bw_str = strtok(NULL, "B");
	}

	output.close();
}

}

#endif

