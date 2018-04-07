#include <iostream>
#include <cstdlib>
#include <climits>
#include <string>
#include <string.h>
#include <vector>
#include <sstream>

using namespace std;

int floyds(int **graph, int **b, int N);

int main() {

	int N, M;

	cin >> N >> M;
	int ** graph = (int **)malloc(N * sizeof(int *));
	int ** copy = (int **)malloc(N * sizeof(int *));
	int i;
	int j;
	for (i = 0; i < N; ++i) {
		graph[i] = (int *)malloc(N * sizeof(int));
		copy[i] = (int *)malloc(N * sizeof(int));
		for (j = 0; j < N; ++j) {
			graph[i][j] = INT_MAX;
			copy[i][j] = INT_MAX;
		}
	}

	int nodeA_idx;
	int nodeB_idx;
	int dist;
	for (i = 0; i < N-1; ++i) {
		cin >> nodeA_idx >> nodeB_idx >> dist;
		graph[nodeA_idx-1][nodeB_idx-1] = dist;
		copy[nodeA_idx-1][nodeB_idx-1] = dist;
	}

	vector<int> result;
	string line;
	string tmp;

	for (i = 0; i < M; ++i) {
		cin >> line;

		if (line[0] == 'Q') {
			result.push_back(floyds(graph, copy, N));
		} else {
			cin >> nodeA_idx >> nodeB_idx >> dist;
			graph[nodeA_idx-1][nodeB_idx-1] = dist;
		}
	}

	for (i = 0; i < result.size(); ++i) {
		cout << result[i] << endl;
	}

	for (i = 0; i < N; ++i) {
		free(graph[i]);
		free(copy[i]);
	}
	free(graph);
	free(copy);

	return 0;
}

int floyds(int **graph, int **b, int N)
{
    int i, j, k;

	for (i = 0; i < N; ++i) {
		for (j = 0; j < N; ++j) {
			if (i == j) {
				b[i][j] = 0;
			} else {
				b[i][j] = graph[i][j];
			}
		}
	}
#if 0
    for (k = 0; k < N; k++)
    {
        for (i = 0; i < N; i++)
        {
            for (j = 0; j < N; j++)
            {
                if ((b[i][k] != INT_MAX) && (b[k][j] != INT_MAX) && (i != j))
                {
                    if ((b[i][k] + b[k][j] < b[i][j]) || (b[i][j] == INT_MAX))
                    {
                        b[i][j] = b[i][k] + b[k][j];
                    }
                }
            }
        }
    }
#endif

#if 1
    for (k = 0; k < N; k++)
    {
        for (i = 0; i < N; i++)
        {
            for (j = 0; j < N; j++)
            {
				if ((b[i][k] + b[k][j] < b[i][j]) && (b[i][k] != INT_MAX) && (b[k][j] != INT_MAX))
				{
					b[i][j] = b[i][k] + b[k][j];
				}
            }
        }
    }
#endif

    int sum = 0;

    for (i = 0; i < N; ++i) {
    	for (j = i+1; j < N; ++j) {
    		sum += b[i][j];
    	}
    }
    return sum;
}
