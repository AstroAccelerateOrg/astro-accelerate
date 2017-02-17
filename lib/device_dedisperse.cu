#include <stdio.h>
#include "AstroAccelerate/params.h"
#include "AstroAccelerate/device_bin.h"

#define NUMDMRANGES 8

//{{{ dedisperse 

static const int UNROLLS_DM[] =  {16, 16, 16, 16, 16, 32, 16, 16};
static const int SNUMREG_DM[] =  {14, 12, 12, 14, 12, 6,  10, 8};
static const int SDIVINT_DM[] =  {12, 8,  10, 6,  8,  8,  14, 6};
static const int SDIVINDM_DM[] = {32, 40, 50, 60, 50, 40, 32, 60};

typedef void (*shared_dedisperse_kernel_FTYPE)(int bin, unsigned short *d_input, float *d_output, float mstartdm, float mdmstep);

shared_dedisperse_kernel_FTYPE shared_dedisperse_kernel_FARRAY[NUMDMRANGES] = {
	shared_dedisperse_kernel_range_0,
	shared_dedisperse_kernel_range_1,
	shared_dedisperse_kernel_range_2,
	shared_dedisperse_kernel_range_3,
	shared_dedisperse_kernel_range_4,
	shared_dedisperse_kernel_range_5,
	shared_dedisperse_kernel_range_6,
	shared_dedisperse_kernel_range_7
};

void dedisperse(int i, int t_processed, int *inBin, float *dmshifts, unsigned short *d_input, float *d_output, int nchans, int nsamp, int maxshift, float *tsamp, float *dm_low, float *dm_high, float *dm_step, int *ndms)
{

	// FOR KEPLER SMEM....
	float shift_one = ( SDIVINDM_DM[i] - 1 ) * ( dm_step[i] / ( *tsamp ) );
	int shifta = (int) floorf(shift_one * dmshifts[nchans - 1]) + ( SDIVINT_DM[i] - 1 ) * 2;
	int lineshift = shifta + ( ( SNUMREG_DM[i] - 1 ) * 2 * SDIVINT_DM[i] );
	printf("\n%f", dm_step[i] / ( *tsamp ));
	printf("\n%f", dmshifts[nchans - 1]);
	printf("\n%d", shifta);
	printf("\nlineshift:\t%d", lineshift);

	// is long enough for the algorithm to run without an out of bounds
	// access...
	if (( ( SDIVINT_DM[i] - 1 ) + ( ( SDIVINDM_DM[i] - 1 ) * SDIVINT_DM[i] ) - 1 ) > lineshift)
	{

		/* FOR FERMI SMEM....
		 float shift_one = (SDIVINDM_DM[i]-1)*(dm_step[i]/(*tsamp));
		 int shift = (int)floorf(shift_one*dmshifts[nchans-1]);
		 int lineshift = shift+((SNUMREG_DM[i]-1)*SDIVINT_DM[i]);
		 //printf("\n%f",dm_step[i]/(*tsamp));
		 //printf("\n%f",dmshifts[nchans-1]);
		 //printf("\n%d",shift);
		 //printf("\nlineshift:\t%d", lineshift);
		 */
		printf("\nUsing fast shared memory kernel");

		//{{{ Dedisperse data on the GPU 
		float startdm = dm_low[i];

		int divisions_in_t = SDIVINT_DM[i];
		int divisions_in_dm = SDIVINDM_DM[i];
		int num_blocks_t = t_processed / ( divisions_in_t * 2 * SNUMREG_DM[i] );
		int num_blocks_dm = ndms[i] / divisions_in_dm;

//		printf("\ntpro:\t%d, numb:\t%d", t_processed, num_blocks_t), fflush(stdout);

//		if(ndms[i]%divisions_in_dm !=0) {
//			printf("\nERROR: dm block size is not a divisor of the dm range!!");
//			exit(0);
//		} 

		dim3 threads_per_block(divisions_in_t, divisions_in_dm);
		dim3 num_blocks(num_blocks_t, num_blocks_dm);

		int array_size = SDIVINT_DM[i]*SDIVINDM_DM[i];
		int shared_mem = UNROLLS_DM[i]*(array_size + 1) * sizeof(ushort2);

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

