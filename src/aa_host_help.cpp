#include <stdio.h>
#include <stdlib.h>
#include "aa_host_help.hpp"

namespace astroaccelerate {

void help()
{
	printf("\n **** Astro-Accelerate **** \n");
	printf("\n Astro-accelerate is a program that performs the de-dispersion of time-series filterbank data.");
	printf("\n The code has algorithms for both CPU and NVIDIA GPU acceleration (Fermi and Kepler).");
	printf("\n Currently astro-accelerate works with a single polarisation and 32-bit data (more to come).");
	printf("\n To use GPU acceleration with no limitations on the maximum DM use the flag \"-algorithm GPU-CACHE\" (about 3x slower than smem algorithm).");
	printf("\n");
	printf("\n Please send any bug reports to wes.armour@oerc.ox.ac.uk");
	printf("\n");
	printf("\n Usage: ./astro-accelerate <filename>\t\tAutomatically de-disperse input data file");
	printf("\n Usage: ./astro-accelerate [arguments]\t\tSpecify arguments for dedispersion");
	printf("\n");
	printf("\n Arguments: \n");
	printf("\n");
	printf("\n\t-f <FILE>\t\t\t\t The input file");
	printf("\n");
	printf("\n\t-tdms <integer value>\t\t\t The total number of de-dispersions to perform for each time sample.");
	printf("\n\t\t\t\t\t\t If this isn't specified the code will use a value that is approximately equal to the number of channels.");
	printf("\n");
	printf("\n\t-dm_step <floating point value>\t\t The step size between each trail dispersion curve.");
	printf("\n\t\t\t\t\t\t If this isn't specified the code will compute this automatically.");
	printf("\n");
	printf("\n\t-dm_low <floating point value>\t\t Starting point of the DM search.");
	printf("\n\t\t\t\t\t\t If this isn't specified the code will start the search at 0.0");
	printf("\n");
	printf("\n\t-kernel <type>\t\t\t\t Specifies the type of algorithm to use");
	printf("\n\t\t\t\t\t\t This can be CPU GPU GPU-CACHE.");
	printf("\n\t\t\t\t\t\t CPU will excecute on the CPU.");
	printf("\n\t\t\t\t\t\t GPU is the default and will use the GPU shared memory. This is our fastest algorithm");
	printf("\n\t\t\t\t\t\t GPU-CACHE will excecute on the gpu and use the gpu cache. This is much faster than CPU but about 3x slower than GPU.");
	printf("\n");
	return;
}

} //namespace astroaccelerate
