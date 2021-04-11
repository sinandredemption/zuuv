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

	const int TruncateConvergents = 10;	// You can probably get away with even 4... just for safety

	const size_t ContinuantUseParallel = 50000;

	const size_t CFTermsUseParallelMultThreshold1 = 2500;
	const size_t CFTermsUseParallelMultThreshold2 = 3000;
	const size_t CFTermsUseDisk = 1000000ULL;
}

#endif