#include <stdio.h>

#include "aa_host_debug.hpp"

namespace astroaccelerate {

  /** \brief Prints the parameters passed to this function. */
  void debug(int test, clock_t start_time, int range, int *outBin, int enable_debug, int analysis, int output_dmt, int multi_file, float sigma_cutoff, float power, int max_ndms, float *user_dm_low, float *user_dm_high, float *user_dm_step, float *dm_low, float *dm_high, float *dm_step, int *ndms, int nchans, int nsamples, int nifs, int nbits, float tsamp, float tstart, float fch1, float foff, int maxshift, float max_dm, int nsamp, size_t gpu_inputsize, size_t gpu_outputsize, size_t inputsize, size_t outputsize) {

    int	i;
    clock_t now;
        
    now = clock();

#ifdef SM_35
    printf("\n Using GPU __ldg() code (version: sm_35)");
#else
    printf("\n Using standard GPU code");
#endif
    if(test == 1) {
      printf("\nrange:\t\t%d", range);
      printf("\ndebug:\t\t%d", enable_debug);
      printf("\nmulti_file:\t%d", multi_file);
      printf("\nanalysis:\t%d", analysis);
      printf("\noutput_dmt:\t%d", output_dmt);
      printf("\nsigma_cutoff:\t%f", sigma_cutoff);
      printf("\npower:\t\t%f", power);
      printf("\nUser requested DM search range:");
      for(i=0; i<range; i++) {
	printf("\n%f\t%f\t%f\t%d", user_dm_low[i], user_dm_high[i], user_dm_step[i], outBin[i]);
      }
      printf("\nRead user input from file, which took:\t\t%.16g(s)\n", (double)(now - start_time) / CLOCKS_PER_SEC);
      fflush(stdout);
    } else if(test == 2) {
      printf("\nInitialised GPU:\t\t%.16g(s)\n", (double)(now - start_time) / CLOCKS_PER_SEC);
      fflush(stdout);
    } else if(test == 3) {
      printf("\ntsamp:\t\t\t%lf", tsamp);
      printf("\ntstart:\t\t\t%lf", tstart);
      printf("\nfch1:\t\t\t%lf", fch1);
      printf("\nfoff:\t\t\t%lf", foff);
      printf("\nnchans:\t\t\t%d", nchans);
      printf("\nnifs:\t\t\t%d", nifs);
      printf("\nnbits:\t\t\t%d", nbits);
      printf("\nnsamples:\t\t%d", nsamples);
      printf("\nnsamp:\t\t\t%d", nsamp);
      printf("\nGot file header info:\t%.16g(s)\n", (double)(now - start_time) / CLOCKS_PER_SEC);
      fflush(stdout);
    } else if(test == 4) {
      printf("\n\nmaximum DM:\t\t%f", max_dm);
      printf("\nmaxshift:\t\t%d", maxshift);
      printf("\nmax_ndms:\t\t%d", max_ndms), fflush(stdout);
      printf("\nActual DM range that will be searched:");
      for(i=0; i<range; i++) {
	printf("\n%f\t%f\t%f\t%d", dm_low[i], dm_high[i], dm_step[i], ndms[i]);
      }
      printf("\nCalculated strategy:\t%.16g(s)\n", (double)(now - start_time) / CLOCKS_PER_SEC);
      fflush(stdout);
    } else if(test == 5) {
      printf("\nMaxshift efficiency:\t\t%.2f%%", 100.0f-((float)maxshift/(float)nsamp)*100.0f); 
      printf("\nHost Input size:\t\t%d MB", (int) (inputsize / 1024 / 1024));
      printf("\nHost Output size:\t\t%d MB", (int) (outputsize / 1024 / 1024));
      printf("\nDevice Input size:\t\t%d MB", (int) (gpu_inputsize / 1024 / 1024));
      printf("\nDevice Output size:\t\t%d MB", (int) (gpu_outputsize / 1024 / 1024));
      printf("\nAllocated memory:\t%.16g(s)\n", (double)(now - start_time) / CLOCKS_PER_SEC);
      fflush(stdout);
    }  else if(test == 6) {
      printf("\nCalculated dm shifts:\t%.16g(s)\n", (double)(now - start_time) / CLOCKS_PER_SEC);
      fflush(stdout);
    }  else if(test == 7) {
      printf("\nGot input filterbank data:\t%.16g(s)\n", (double)(now - start_time) / CLOCKS_PER_SEC);
      fflush(stdout);
    }  else if(test == 8) {
      printf("\nLoaded data onto the GPU:\t%.16g(s)\n", (double)(now - start_time) / CLOCKS_PER_SEC);
      fflush(stdout);
    }
  }

} //namespace astroaccelerate
