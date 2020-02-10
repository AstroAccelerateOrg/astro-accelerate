#ifndef ASTRO_ACCELERATE_AA_DEVICE_LOAD_DATA_HPP
#define ASTRO_ACCELERATE_AA_DEVICE_LOAD_DATA_HPP

#include <stdio.h>
#include <math.h>

namespace astroaccelerate {

	/**
	* Loads data from the host memory into the GPU memory.
	* Users should not need to interact with this function directly.
	*/
	void load_data(int i, int *inBin, unsigned short *device_pointer, unsigned short const*const host_pointer, int t_processed, int maxshift, int nchans, float *dmshifts, float *d_dm_shifts, int nbits);
	
	void load_chunk_data(unsigned short *device_pointer, unsigned short const*const host_pointer, int t_processed, int maxshift, int nchans, float *dmshifts, float *d_dm_shifts, int nbits);
	
	void set_dedispersion_constants(int t_processed, int maxshift);

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_DEVICE_LOAD_DATA_HPP
