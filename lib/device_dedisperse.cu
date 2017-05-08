#include <stdio.h>
#include "headers/params.h"
#include "headers/device_bin.h"

//{{{ dedisperse 

void dedisperse(int i, int t_processed, int *inBin, float *dmshifts, unsigned short *d_input, float *d_output, int nchans, int nsamp, int maxshift, float *tsamp, float *dm_low, float *dm_high, float *dm_step, int *ndms, int nbits)
{

	if(nbits == 16) {

		float shift_one = ( SDIVINDM - 1 ) * ( dm_step[i] / ( *tsamp ) );
		int shifta = (int) floorf(shift_one * dmshifts[nchans - 1]) + ( SDIVINT - 1 ) * 2;
		int lineshift = shifta + ( ( SNUMREG - 1 ) * 2 * SDIVINT );
		printf("\n%f", dm_step[i] / ( *tsamp ));
		printf("\n%f", dmshifts[nchans - 1]);
		printf("\n%d", shifta);
		printf("\nlineshift:\t%d", lineshift);

		// Check to see if the threadblock will load a shared memory line that
		// is long enough for the algorithm to run without an out of bounds
		// access...
		if (( ( SDIVINT - 1 ) + ( ( SDIVINDM - 1 ) * SDIVINT ) - 1 ) > lineshift)
		{

			printf("\nUsing fast shared memory kernel");

			//{{{ Dedisperse data on the GPU 
			float startdm = dm_low[i];

			int divisions_in_t = SDIVINT;
			int divisions_in_dm = SDIVINDM;
			int num_blocks_t = t_processed / ( divisions_in_t * 2 * SNUMREG );
			int num_blocks_dm = ndms[i] / divisions_in_dm;


			dim3 threads_per_block(divisions_in_t, divisions_in_dm);
			dim3 num_blocks(num_blocks_t, num_blocks_dm);

			cudaDeviceSetSharedMemConfig (cudaSharedMemBankSizeFourByte);
			cudaFuncSetCacheConfig(shared_dedisperse_kernel_16, cudaFuncCachePreferShared);
			shared_dedisperse_kernel_16<<<num_blocks, threads_per_block>>>(inBin[i], d_input, d_output, (float) ( startdm / ( *tsamp ) ), (float) ( dm_step[i] / ( *tsamp ) ));
		}
		else
		{
			printf("\nERROR: smem line length is too short.\nRun the auto tuner again!\n");
			exit(0);
		}

	}
	else
	{
		// FOR KEPLER SMEM....
		float shift_one = ( SDIVINDM - 1 ) * ( dm_step[i] / ( *tsamp ) );
		int shifta = (int) floorf(shift_one * dmshifts[nchans - 1]) + ( SDIVINT - 1 ) * 2;
		int lineshift = shifta + ( ( SNUMREG - 1 ) * 2 * SDIVINT );
		printf("\n%f", dm_step[i] / ( *tsamp ));
		printf("\n%f", dmshifts[nchans - 1]);
		printf("\n%d", shifta);
		printf("\nlineshift:\t%d", lineshift);

		// Check to see if the threadblock will load a shared memory line that
		// is long enough for the algorithm to run without an out of bounds
		// access...
		if (( ( SDIVINT - 1 ) + ( ( SDIVINDM - 1 ) * SDIVINT ) - 1 ) > lineshift)
		{

			printf("\nUsing fast shared memory kernel");

			//{{{ Dedisperse data on the GPU 
			float startdm = dm_low[i];

			int divisions_in_t = SDIVINT;
			int divisions_in_dm = SDIVINDM;
			int num_blocks_t = t_processed / ( divisions_in_t * 2 * SNUMREG );
			int num_blocks_dm = ndms[i] / divisions_in_dm;


			dim3 threads_per_block(divisions_in_t, divisions_in_dm);
			dim3 num_blocks(num_blocks_t, num_blocks_dm);

			cudaDeviceSetSharedMemConfig (cudaSharedMemBankSizeFourByte);
			cudaFuncSetCacheConfig(shared_dedisperse_kernel, cudaFuncCachePreferShared);
			shared_dedisperse_kernel<<<num_blocks, threads_per_block>>>(inBin[i], d_input, d_output, (float) ( startdm / ( *tsamp ) ), (float) ( dm_step[i] / ( *tsamp ) ));
		}
		else
		{
			printf("\nERROR: smem line length is too short.\nRun the auto tuner again!\n");
			exit(0);
		}

	}
}

//}}}

