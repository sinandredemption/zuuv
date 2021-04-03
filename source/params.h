#ifndef INC_PARAMS_H
#define INC_PARAMS_H

namespace Params {
	const size_t ThresholdUseBasicContProc = 1024;
	const size_t ThresholdCFDivide = 14;
	const size_t DefaultContAllocSize = ThresholdUseBasicContProc >> 6;
	const size_t DefaultConvAllocSize = 512;

	const size_t ThresholdReduc1 = 80;
	const size_t ThresholdReduc2 = 2560;
	const size_t ThresholdReduc3 = 999999;
	const size_t ThresholdReduc4 = 9999999;

	const size_t BenchmarkDefaultSize = 1000000;	// Number of starting digits in CF benchmark
	const size_t BenchmarkIters = 13;	// Number of Iterations to do in benchmark
	const double BenchmarkMultFactor = 2.0;	// By default, double the size every iteration

	const int TruncateConvergents = 10;

	const size_t ContinuantUseParallel = 50000;
	const size_t CFTermsUseParallelMultThreshold = 2500;
}

#endif