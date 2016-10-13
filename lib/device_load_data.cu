//#include "AstroAccelerate/device_dedispersion_kernel.h"

//extern "C" void load_data(int i, float *device_pointer, float *host_pointer, size_t size, int nsamp, int maxshift, int nchans, int t_processed_s, int t_processed_c, float *dmshifts);

<<<<<<< HEAD
//{{{ load_data_from_host_to_device

void load_data(int i, int *inBin, unsigned short *device_pointer, unsigned short *host_pointer, int t_processed, int maxshift, int nchans, float *dmshifts) {

	//{{{ Copy data and set up the GPU constants/variables.
	if(i==-1) {
		int length=(t_processed+maxshift);
		size_t size=nchans*length*sizeof(unsigned short);
=======
void load_data(int i, int *inBin, unsigned short *device_pointer, unsigned short *host_pointer, int t_processed, int maxshift, int nchans, float *dmshifts)
{
	// Copy data and set up the GPU constants/variables.
	if(i == -1)
	{
		int length = (t_processed+maxshift);
		size_t size = nchans*length*sizeof(unsigned short);
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
		cudaMemcpyToSymbol(dm_shifts, dmshifts, nchans * sizeof(float));
		cudaMemcpy(device_pointer, host_pointer, size, cudaMemcpyHostToDevice);
		cudaMemcpyToSymbol(i_nchans, &nchans, sizeof(int));
		cudaMemcpyToSymbol(i_nsamp, &length, sizeof(int));
		cudaMemcpyToSymbol(i_t_processed_s, &t_processed, sizeof(int));
<<<<<<< HEAD
	} else if (i > 0) {
		int length=(t_processed+maxshift);
		cudaMemcpyToSymbol(i_nsamp, &length, sizeof(int));
		cudaMemcpyToSymbol(i_t_processed_s, &t_processed, sizeof(int));
	}
	//}}}
=======
	}
	else if (i > 0)
	{
		int length = (t_processed+maxshift);
		cudaMemcpyToSymbol(i_nsamp, &length, sizeof(int));
		cudaMemcpyToSymbol(i_t_processed_s, &t_processed, sizeof(int));
	}
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
	
	float h_sqrt_taps[PD_MAXTAPS+1];
	for(int f=0; f<=PD_MAXTAPS; f++) h_sqrt_taps[f]=(float) sqrt((double) f);
	cudaMemcpyToSymbol(c_sqrt_taps, h_sqrt_taps, (PD_MAXTAPS+1)*sizeof(float));

<<<<<<< HEAD
}

//}}}
=======
}
>>>>>>> fe80b9c735d1c898047cbb64bcf8da05cd6a21da
