#ifndef __BC_PLAN__
#define __BC_PLAN__

struct PulseDetection_plan {
	int decimated_timesamples;
	int iteration;
	int nBoxcars;
	int nBlocks;
	int output_shift;
	int shift;
	int startTaps;
	int unprocessed_samples;
	int total_ut;
};

#endif
