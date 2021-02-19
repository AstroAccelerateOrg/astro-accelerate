#include <stdio.h>
#include <stdlib.h>
#include "aa_host_help.hpp"

namespace astroaccelerate {

void help()
{
	printf("\n **** AstroAccelerate **** \n");
	printf("\n AstroAccelerate is a GPU accelerated software package for processing time-domain radio astronomy data.");
	printf("\n AstroAccelerate is able to perform de-dispersion, single-pulse search, periodicity search, frequency domain acceleration search.");
	printf("\n The code is optimised for NVIDIA GPUs.");
	printf("\n Currently, AstroAccelerate works with a single polarisation and 4-bit, 8-bit, and 16-bit input data.");
	printf("\n To use AstroAccelerate, please provide a configuration file, an example of which is in the input_files directory, or use provided python interface and python scripts which are located in the python directory.");
	printf("\n Please check our github for more information: https://github.com/AstroAccelerateOrg");
	printf("\n Also, please cite our publications and AstroAccelerate as listed in our https://github.com/AstroAccelerateOrg/astro-accelerate/blob/master/README.md");
	printf("\n");
	printf("\n Please send any bug reports to astroaccelerate@maillist.ox.ac.uk");
	printf("\n");
	printf("\n");
	return;
}

} //namespace astroaccelerate
