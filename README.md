Sorting is a crucial operation that could be used to implement SQL
operators such as GROUP BY, ORDER BY, and SQL:2003 PARTITION BY.
Queries with *multiple attributes* in those clauses are common
in real workloads.
When executing queries of that kind,
state-of-the-art main-memory column-stores require *one round of sorting per input column*.
With the advent of recent fast scans and
denormalization techniques, that kind of 
**multi-column sorting** could become a bottleneck.
In this project, we propose a new technique called *code massaging*, which
manipulates the bits across the columns 
so that the overall sorting time can be reduced by 
eliminating some rounds of sorting and/or 
by improving the degree of SIMD data level parallelism.

# Build from source

## Clone

```bash
git clone --recursive https://github.com/cswxu/FastMultiColSort.git
```

Or this after cloning without `--recursive`:

```bash
git submodule update --init --recursive
```


## Build

You need [CMake](https://cmake.org/) to generate build scripts. Makefile is tested.

To generate debug build:

```bash
mkdir debug
cd debug
cmake -DCMAKE_BUILD_TYPE=debug ..
make -j4
```

To generate release build:

```bash
mkdir release
cd release
cmake -DCMAKE_BUILD_TYPE=release ..
make -j4
```

NOTE: The default build type is `debug`, which may not give optimal
performance.


# Running executables

Executable programs are in 'experiments/' directory.

To see a full list of options:

```bash
experiments/main_run_an_instance --help
```

An example to run microbenchmark:

**Step 1**: generate data for microbenchmark

```bash
./experiments/main_run_an_instance --gen-data --nrows=16777216 --ncolumns=2 --bitwidth=17B33 --ngroups=8192
```
which means generating 16777216 rows of data with two columns, the first column is encoded in 17-bit and the second one is in 33-bit. 
The generated file is `16777216_2_17B33_0_100_8192.dat`

**Step 2**: run two-column-sorting with all cases of bit-borrowing (see the paper) over these two columns:

```bash
./experiments/main_tworound_exhaustive --nrows=16777216 --ncolumns=2 --bitwidth=17B33 --inputfile=16777216_2_17B33_0_100_8192.dat --comptype=st
```



# Running tests

```bash
make check
```

Build tests without running.

```bash
make check-build
```

# File structure

+ `experiments/` - Executables

+ `third-party/` - Third-party libraries

+ `src/` - source files

+ `tests/` - Unit tests written in GoogleTest framework


# Citing this work

Wenjian Xu, Ziqiang Feng and Eric Lo. "**Fast Multi-Column Sorting in Main-Memory Column-Stores.**"
In Proceedings of the 2016 ACM SIGMOD International Conference on
Management of Data, pp. 1263-1278. ACM, 2016.

Download: https://dl.acm.org/citation.cfm?id=2915205

BibTex:
```bash
@inproceedings{Xu:2016:FMS:2882903.2915205,
 author = {Xu, Wenjian and Feng, Ziqiang and Lo, Eric},
 title = {Fast Multi-Column Sorting in Main-Memory Column-Stores},
 booktitle = {Proceedings of the 2016 International Conference on Management of Data},
 series = {SIGMOD '16},
 year = {2016},
 isbn = {978-1-4503-3531-7},
 location = {San Francisco, California, USA},
 pages = {1263--1278},
 numpages = {16},
 url = {http://doi.acm.org/10.1145/2882903.2915205},
 doi = {10.1145/2882903.2915205},
 acmid = {2915205},
 publisher = {ACM},
 address = {New York, NY, USA},
 keywords = {SIMD, column-store, main-memory, multi-column sorting},
} 
```


# Contact

Wenjian Xu ( cswxu at comp dot polyu dot edu dot hk )


# Platform requirements

1. C++ compiler supporting C++11 and AVX2
2. CPU with AVX2 instruction set extension


# Tested platform

This package has been tested with the following configuration:

- Linux 3.13.0-96-generic (64-bit)
- Intel(R) Core(TM) i7-4770 CPU @ 3.40GHz
- g++ 4.9.3





