#include <stdio.h>
#include <stdlib.h>

#include "aa_params.hpp"

namespace astroaccelerate {

  /** \brief Print run statistics. */
  void statistics(char *string, int i, cudaStream_t stream, clock_t *in_time, clock_t *out_time, int maxshift, int total_ndms, int nchans, int nsamp, float tsamp, float *dm_low, float *dm_high, float *dm_step, int *ndms)
  {

    int num_reg = SNUMREG;
    int divisions_in_t = SDIVINT;
    int divisions_in_dm = SDIVINDM;
    int num_blocks_t = ( nsamp - maxshift ) / ( divisions_in_t * num_reg );
    int num_blocks_dm = total_ndms / divisions_in_dm;

    double time = 0.0;

    if (strcmp(string, "dedisperse in") == 0)
      {
	*in_time = clock();

      }
    else if (strcmp(string, "dedisperse out") == 0)
      {
	*out_time = clock();
	time = double ( *out_time - *in_time ) / CLOCKS_PER_SEC;

	printf("\nPerformed Brute-Force Dedispersion: %f (GPU estimate)", time);

	int telescope_data = ( nsamp - maxshift );
	float telescope_time = (float) telescope_data * tsamp;
	printf("\nTelescope data:\t%d", telescope_data);
	printf("\nTelescope time:\t%f", tsamp);
	printf("\nReal-time speedup factor: %lf", telescope_time / (float) time);

	double time_processed = nsamp - maxshift;
	double dm_t_processed = time_processed * total_ndms;
	double all_processed = dm_t_processed * nchans;

	printf("\nGops based on %d ops per channel per tsamp: %lf", NOPSLOOP, ( ( NOPSLOOP * all_processed ) / ( time ) ) / 1000000000.0);
	printf("\nDevice global memory bandwidth in GB/s: %f", ( ( ( divisions_in_t * divisions_in_dm * num_blocks_t * num_blocks_dm ) * nchans * sizeof(float) ) / ( time ) ) / 1000000000);

	float num_threads = total_ndms * ( nsamp - maxshift ) / ( num_reg );
	float data_size_loaded = ( num_threads * nchans * sizeof(float) ) / 1000000000;
	float bandwidth = data_size_loaded / time;

	printf("\nDevice global memory bandwidth in GB/s: %f", bandwidth);
	printf("\nDevice bandwidth to cache in GB/s: %f", bandwidth * num_reg);

	float size_gb = ( nchans * ( nsamp - maxshift ) ) / 1000000000.0;
	printf("\nTelescope data throughput in Gb/s: %f", size_gb / time);

      }
    else if (strcmp(string, "save in") == 0)
      {

      }
    else if (strcmp(string, "save out") == 0)
      {

      }
    else if (strcmp(string, "analyse in") == 0)
      {

      }
    else if (strcmp(string, "analyse out") == 0)
      {

      }
  }

} //namespace astroaccelerate
