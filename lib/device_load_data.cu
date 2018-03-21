//#include "headers/device_dedispersion_kernel.h"

//extern "C" void load_data(int i, float *device_pointer, float *host_pointer, size_t size, int nsamp, int maxshift, int nchans, int t_processed_s, int t_processed_c, float *dmshifts);

//{{{ load_data_from_host_to_device

void load_data(int i, int *inBin, unsigned short *device_pointer, unsigned short *host_pointer, int t_processed, int maxshift, int nchans, float *dmshifts, cudaStream_t stream)
{

	//{{{ Copy data and set up the GPU constants/variables.
	if (i == -1)
	{
		long int length = ( t_processed + maxshift );
		size_t size = nchans * length * sizeof(unsigned short);
		cudaMemcpyToSymbolAsync(dm_shifts, dmshifts, nchans * sizeof(float),0,cudaMemcpyHostToDevice,stream);
//		cudaMemcpy(device_pointer, host_pointer, size, cudaMemcpyHostToDevice);
		cudaMemcpyAsync(device_pointer, host_pointer, size, cudaMemcpyHostToDevice,stream);
		cudaMemcpyToSymbolAsync(i_nchans, &nchans, sizeof(int),0,cudaMemcpyHostToDevice,stream);
		cudaMemcpyToSymbolAsync(i_nsamp, &length, sizeof(int),0,cudaMemcpyHostToDevice,stream);
		cudaMemcpyToSymbolAsync(i_t_processed_s, &t_processed, sizeof(int),0,cudaMemcpyHostToDevice,stream);
	}
	else if (i > 0)
	{
		long int length = ( t_processed + maxshift );
//		cudaHostRegister(i_nsamp,sizeof(int),0);
		cudaMemcpyToSymbolAsync(i_nsamp, &length, sizeof(int),0,cudaMemcpyHostToDevice,stream);
		cudaMemcpyToSymbolAsync(i_t_processed_s, &t_processed, sizeof(int),0,cudaMemcpyHostToDevice,stream);
	}
	//}}}

	float h_sqrt_taps[PD_MAXTAPS + 1];
//	float *h_sqrt_taps = NULL;
//	cudaMallocHost((void **) &h_sqrt_taps,(PD_MAXTAPS+1)*sizeof(float));
	for (int f = 0; f <= PD_MAXTAPS; f++)
		h_sqrt_taps[f] = (float) sqrt((double) f);
	cudaMemcpyToSymbolAsync(c_sqrt_taps, h_sqrt_taps, ( PD_MAXTAPS + 1 ) * sizeof(float),0,cudaMemcpyHostToDevice,stream);
//	cudaFreeHost(h_sqrt_taps);

}

//}}}
