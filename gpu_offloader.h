#ifndef GPU_OFFLOADER
#define GPU_OFFLOADER

#include <vector>
#include "gpuwork.decl.h"

using namespace std;

void gpu_offload(CProxy_LUBlk block, vector<int> cellIndices, vector<int> iter,
		 vector<vector<float> > rows, int rowlen, int number);

#endif
