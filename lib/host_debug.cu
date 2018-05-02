#include <stdio.h>
#include "headers/headers_mains.h"
//#include <omp.h>

void debug(int test, clock_t start_time, DDTR_InputData *DDTR_data, DDTR_Plan *DDTR_plan, AA_Parameters *AA_params, MSD_Parameters *MSD_params, SPS_Parameters *SPS_params, PRS_Parameters *PRS_params, FDAS_Parameters *FDAS_params){
	int	i;
	clock_t now;
        
	now = clock();

#ifdef SM_35
	printf("\n Using GPU __ldg() code (version: sm_35)");
#else
	printf("\n Using standard GPU code");
#endif
	if(test == 1) {
		printf("\nrange:\t\t%d", DDTR_plan->nRanges);
		printf("\ndebug:\t\t%d", AA_params->enable_debug);
		printf("\nanalysis:\t%d", AA_params->enable_analysis);
		printf("\nacceleration:\t%d", AA_params->enable_acceleration);
		printf("\nperiodicity:\t%d", AA_params->enable_periodicity);
		printf("\nSPS sigma_cutoff:\t%f", SPS_params->sigma_cutoff);
		printf("\nPRS sigma_cutoff:\t%f", PRS_params->sigma_cutoff);
		printf("\nFDAS sigma_cutoff:\t%f", FDAS_params->sigma_cutoff);
		printf("\npower:\t\t%f", DDTR_plan->power);
		printf("\nUser requested DM search range:");
		for(i=0; i<DDTR_plan->nRanges; i++) {
			printf("\n%f\t%f\t%f\t%d", DDTR_plan->user_dm_low[i], DDTR_plan->user_dm_high[i], DDTR_plan->user_dm_step[i], DDTR_plan->inBin[i]);
		}
		printf("\nGot user input:\t\t%.16g(s)\n", (double)(now - start_time) / CLOCKS_PER_SEC);
		fflush(stdout);
	} else if(test == 2) {
		printf("\nInitialised GPU:\t\t%.16g(s)\n", (double)(now - start_time) / CLOCKS_PER_SEC);
		fflush(stdout);
	} else if(test == 3) {
		printf("\ntsamp:\t\t\t%lf", DDTR_data->tsamp);
		printf("\ntstart:\t\t\t%lf", DDTR_data->tstart);
		printf("\nfch1:\t\t\t%lf", DDTR_data->fch1);
		printf("\nfoff:\t\t\t%lf", DDTR_data->foff);
		printf("\nnchans:\t\t\t%d", DDTR_data->nchans);
		printf("\nnifs:\t\t\t%d", DDTR_data->nifs);
		printf("\nnbits:\t\t\t%d", DDTR_data->nbits);
		printf("\nnsamples:\t\t%d", DDTR_data->nsamples);
        printf("\nnsamp:\t\t\t%d", DDTR_data->nsamp);
		printf("\nGot file header info:\t%.16g(s)\n", (double)(now - start_time) / CLOCKS_PER_SEC);
		fflush(stdout);
	} else if(test == 4) {
		printf("\n\nmaximum DM:\t\t%f", DDTR_plan->max_dm);
		printf("\nmaxshift:\t\t%d", DDTR_plan->maxshift);
		printf("\nmax_ndms:\t\t%d", DDTR_plan->max_ndms), fflush(stdout);
		printf("\nActual DM range that will be searched:");
		for(i=0; i<DDTR_plan->nRanges; i++) {
			printf("\n%f\t%f\t%f\t%d", DDTR_plan->dm_low[i], DDTR_plan->dm_high[i], DDTR_plan->dm_step[i], DDTR_plan->ndms[i]);
		}
		printf("\nCalculated strategy:\t%.16g(s)\n", (double)(now - start_time) / CLOCKS_PER_SEC);
		fflush(stdout);
	} else if(test == 5) {
		printf("\nMaxshift efficiency:\t\t%.2f%%", 100.0f-((float)DDTR_plan->maxshift/(float)DDTR_plan->nsamp)*100.0f); 
		printf("\nHost Input size:\t\t%d MB", (int) ((DDTR_plan->nsamp*DDTR_plan->nchans*sizeof(unsigned short)) / 1024 / 1024));
		printf("\nHost Output size:\t\t%d MB", (int) (DDTR_plan->host_outputsize / 1024 / 1024));
		printf("\nDevice Input size:\t\t%d MB", (int) (DDTR_plan->gpu_inputsize / 1024 / 1024));
		printf("\nDevice Output size:\t\t%d MB", (int) (DDTR_plan->gpu_outputsize / 1024 / 1024));
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
	} /*else if(test == 5) {
        	printf("\npos:\t\t\t%ld", pos);
		printf("\nGot input data:\t\t%.16g(s)\n", (double)(now - start_time) / CLOCKS_PER_SECOND);
		fflush(stdout);
	}*/
}
