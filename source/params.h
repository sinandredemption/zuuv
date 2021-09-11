#ifndef INC_PARAMS_H
#define INC_PARAMS_H

namespace Params
{
	const size_t ThresholdUseBasicContProc = 1024;
	const size_t DefaultContAllocSize = ThresholdUseBasicContProc >> 6;

	const size_t ThresholdCFUseBasicProcedure = 14;
	const size_t DefaultCFTermsAlloc = 512;
	const size_t ThresholdLowerHalfReductionSize1 = 80;
	const float  LowerHalfReductionFactor1        = 2.0;
	const size_t ThresholdLowerHalfReductionSize2 = 2560;
	const float  LowerHalfReductionFactor2        = 2.8;
	const size_t ThresholdLowerHalfReductionSize3 = 999999;
	const float  LowerHalfReductionFactor3        = 3.0;
	const size_t ThresholdLowerHalfReductionSize4 = 9999999;
	const float  LowerHalfReductionFactor4        = 3.0;

	const int TruncateConvergents = 10;	// You can probably get away with even 4... just for safety

	const size_t ContinuantUseParallel = 50000;

	const size_t ThresholdUseYMPMul = 0;
	const size_t ThresholdAddSubInParallel = 25000;
	const size_t CFTermsUseDisk = 1000000ULL;

	const size_t BaselineRAMUsage = 466;
	const size_t RAMUsagePerMBofFraction = 72;
}

#endif