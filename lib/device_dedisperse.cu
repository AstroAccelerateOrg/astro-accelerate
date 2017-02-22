#include <stdio.h>
#include "AstroAccelerate/params.h"
#include "AstroAccelerate/device_bin.h"
#include "AstroAccelerate/kernel_params.h"
#include "AstroAccelerate/kernel_functions.h"

#define NUMDMRANGES 8

//{{{ dedisperse 

void dedisperse(int i, int t_processed, int *inBin, float *dmshifts, unsigned short *d_input, float *d_output, int nchans, int nsamp, int maxshift, float *tsamp, float *dm_low, float *dm_high, float *dm_step, int *ndms)
{

	// FOR KEPLER SMEM....
	float shift_one = ( SDIVINDM_ARRAY[i] - 1 ) * ( dm_step[i] / ( *tsamp ) );
	int shifta = (int) floorf(shift_one * dmshifts[nchans - 1]) + ( SDIVINT_ARRAY[i] - 1 ) * 2;
	int lineshift = shifta + ( ( SNUMREG_ARRAY[i] - 1 ) * 2 * SDIVINT_ARRAY[i] );
	printf("\n%f", dm_step[i] / ( *tsamp ));
	printf("\n%f", dmshifts[nchans - 1]);
	printf("\n%d", shifta);
	printf("\nlineshift:\t%d", lineshift);

	// is long enough for the algorithm to run without an out of bounds
	// access...
	if (( ( SDIVINT_ARRAY[i] - 1 ) + ( ( SDIVINDM_ARRAY[i] - 1 ) * SDIVINT_ARRAY[i] ) - 1 ) > lineshift)
	{

		/* FOR FERMI SMEM....
		 float shift_one = (SDIVINDM_ARRAY[i]-1)*(dm_step[i]/(*tsamp));
		 int shift = (int)floorf(shift_one*dmshifts[nchans-1]);
		 int lineshift = shift+((SNUMREG_ARRAY[i]-1)*SDIVINT_ARRAY[i]);
		 //printf("\n%f",dm_step[i]/(*tsamp));
		 //printf("\n%f",dmshifts[nchans-1]);
		 //printf("\n%d",shift);
		 //printf("\nlineshift:\t%d", lineshift);
		 */
		printf("\nUsing fast shared memory kernel");

		//{{{ Dedisperse data on the GPU 
		float startdm = dm_low[i];

		int divisions_in_t = SDIVINT_ARRAY[i];
		int divisions_in_dm = SDIVINDM_ARRAY[i];
		int num_blocks_t = t_processed / ( divisions_in_t * 2 * SNUMREG_ARRAY[i] );
		int num_blocks_dm = ndms[i] / divisions_in_dm;

//		printf("\ntpro:\t%d, numb:\t%d", t_processed, num_blocks_t), fflush(stdout);

//		if(ndms[i]%divisions_in_dm !=0) {
//			printf("\nERROR: dm block size is not a divisor of the dm range!!");
//			exit(0);
//		} 

		dim3 threads_per_block(divisions_in_t, divisions_in_dm);
		dim3 num_blocks(num_blocks_t, num_blocks_dm);

		int array_size = SDIVINT_ARRAY[i]*SDIVINDM_ARRAY[i];
		int shared_mem = UNROLLS_ARRAY[i]*(array_size + 1) * sizeof(ushort2);

		cudaDeviceSetSharedMemConfig (cudaSharedMemBankSizeFourByte);
		cudaFuncSetCacheConfig(shared_dedisperse_kernel_FARRAY[i], cudaFuncCachePreferShared);
		shared_dedisperse_kernel_FARRAY[i]<<<num_blocks, threads_per_block, shared_mem>>>(inBin[i], d_input, d_output, (float) ( startdm / ( *tsamp ) ), (float) ( dm_step[i] / ( *tsamp) ));

	}
	else
	{
		printf("\nERROR: smem line length is too short.\nRun the auto tuner again!\n");
		exit(0);
	}

	//}}}
	//cudaUnbindTexture(inTex);
	//cudaDestroyTextureObject(tex);
}

//}}}
