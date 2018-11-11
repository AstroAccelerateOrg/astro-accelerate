#ifndef ASTRO_ACCELERATE_DEVICE_BC_PLAN_HPP
#define ASTRO_ACCELERATE_DEVICE_BC_PLAN_HPP

namespace astroaccelerate {

struct PulseDetection_plan {
	int decimated_timesamples;
	int dtm;
	int iteration;
	int nBoxcars;
	int nBlocks;
	int output_shift;
	int shift;
	int startTaps;
	int unprocessed_samples;
	int total_ut;
};

}

#endif
