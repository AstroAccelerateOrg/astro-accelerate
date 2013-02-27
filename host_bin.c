#include <stdio.h>
#include <stdlib.h>

void bin(size_t binsize, float *bin_buffer, size_t inputsize, float *input_buffer, int nchans, int nsamp) {

	int	c;
	int	s;
	int	next_nsamp = nsamp/2;
	
	// NEEDS VECTORISING !!

	#pragma omp parallel for
	for(c = 0; c < nchans; c++) {
		for(s = 0; s < next_nsamp; s++) {
			bin_buffer[((c*next_nsamp) + s)] = (input_buffer[c*nsamp + 2*s] + input_buffer[c*nsamp + 2*s +1])/2;
		}
	}
}

