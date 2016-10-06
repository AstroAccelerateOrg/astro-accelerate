#include <stdio.h>
#include "AstroAccelerate/params.h"
#include "AstroAccelerate/device_bin.h"

<<<<<<< HEAD
void dedisperse(int i, int t_processed, int *inBin, float *dmshifts, unsigned short *d_input, float *d_output, int nchans, int nsamp, int maxshift, float *tsamp, float *dm_low, float *dm_high, float *dm_step, int *ndms) 
{
=======
//{{{ dedisperse 

void dedisperse(int i, int t_processed, int *inBin, float *dmshifts, float *d_input, cudaTextureObject_t tex, float *d_output, int nchans, int nsamp, int maxshift, float *tsamp, float *dm_low, float *dm_high, float *dm_step, int *ndms) {

>>>>>>> 0ec19baf405fa311d6a7ea91dbb146bcccf88229

	// FOR KEPLER SMEM....
	float shift_one = (SDIVINDM-1)*(dm_step[i]/(*tsamp));	
	int shifta = (int)floorf(shift_one*dmshifts[nchans-1]) + (SDIVINT-1)*2;
<<<<<<< HEAD
	int lineshift = shifta+((SNUMREG-1)*2*SDIVINT);
	printf("\n%f",dm_step[i]/(*tsamp));
	printf("\n%f",dmshifts[nchans-1]);
	printf("\n%d",shifta);
	printf("\nlineshift:\t%d", lineshift);
=======
//	int lineshift = shifta+((SNUMREG-1)*2*SDIVINT);
//	printf("\n%f",dm_step[i]/(*tsamp));
//	printf("\n%f",dmshifts[nchans-1]);
//	printf("\n%d",shifta);
//	printf("\nlineshift:\t%d", lineshift);
>>>>>>> 0ec19baf405fa311d6a7ea91dbb146bcccf88229

	// Check to see if the threadblock will load a shared memory line that
	// is long enough for the algorithm to run without an out of bounds
	// access...
<<<<<<< HEAD
	if(((SDIVINT-1)+((SDIVINDM-1)*SDIVINT) - 1) > lineshift)
	{
=======
//	if(((SDIVINT-1)+((SDIVINDM-1)*SDIVINT) - 1) > lineshift) {
>>>>>>> 0ec19baf405fa311d6a7ea91dbb146bcccf88229

	/* FOR FERMI SMEM....
	float shift_one = (SDIVINDM-1)*(dm_step[i]/(*tsamp));	
	int shift = (int)floorf(shift_one*dmshifts[nchans-1]);
	int lineshift = shift+((SNUMREG-1)*SDIVINT);
	//printf("\n%f",dm_step[i]/(*tsamp));
	//printf("\n%f",dmshifts[nchans-1]);
	//printf("\n%d",shift);
	//printf("\nlineshift:\t%d", lineshift);
	*/
		printf("\nUsing fast shared memory kernel");

<<<<<<< HEAD
		// Dedisperse data on the GPU 
=======
		//{{{ Dedisperse data on the GPU 
>>>>>>> 0ec19baf405fa311d6a7ea91dbb146bcccf88229
		float startdm = dm_low[i];

		int divisions_in_t  = SDIVINT;
		int divisions_in_dm = SDIVINDM;
<<<<<<< HEAD
		int num_blocks_t    = t_processed/(divisions_in_t*2*SNUMREG);
=======
		int num_blocks_t    = t_processed/(divisions_in_t*SNUMREG);
>>>>>>> 0ec19baf405fa311d6a7ea91dbb146bcccf88229
		int num_blocks_dm   = ndms[i]/divisions_in_dm;

//		printf("\ntpro:\t%d, numb:\t%d", t_processed, num_blocks_t), fflush(stdout);

//		if(ndms[i]%divisions_in_dm !=0) {
//			printf("\nERROR: dm block size is not a divisor of the dm range!!");
//			exit(0);
//		} 
		
		dim3 threads_per_block(divisions_in_t, divisions_in_dm);
		dim3 num_blocks(num_blocks_t,num_blocks_dm);

<<<<<<< HEAD
		cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeFourByte);
		cudaFuncSetCacheConfig(shared_dedisperse_kernel, cudaFuncCachePreferShared);
		shared_dedisperse_kernel<<< num_blocks, threads_per_block >>>(inBin[i], d_input, d_output, (float)(startdm/(*tsamp)), (float)(dm_step[i]/(*tsamp)));
	} 
	else
	{
		printf("\nERROR: smem line length is too short.\nRun the auto tuner again!\n");
		exit(0);
	}

	//cudaUnbindTexture(inTex);
	//cudaDestroyTextureObject(tex);
}
=======
		cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);
		cudaFuncSetCacheConfig(shared_dedisperse_kernel, cudaFuncCachePreferShared);
		shared_dedisperse_kernel<<< num_blocks, threads_per_block >>>(d_input, d_output, tex, (float)(startdm/(*tsamp)), (float)(dm_step[i]/(*tsamp)));
//	} else {
//		printf("\nERROR: smem line length is too short.\nRun the auto tuner again!\n");
//		exit(0);
//	}

	//}}}
	//cudaUnbindTexture(inTex);
	//cudaDestroyTextureObject(tex);
}

//}}}

>>>>>>> 0ec19baf405fa311d6a7ea91dbb146bcccf88229
