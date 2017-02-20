//
#include <stdio.h>
#include <stdlib.h>
#include <cufft.h>
//#include <omp.h>
//
#include <errno.h>
#include <string.h>
#include <sys/time.h>
//#include <helper_functions.h>
#include <helper_cuda.h>
#include <cuda_profiler_api.h>
//
#include "AstroAccelerate/params.h"
#include "helper_cuda.h"
#include "AstroAccelerate/fdas.h"
#include "AstroAccelerate/fdas_host.h"

void acceleration_fdas(int range,
					   int nsamp,
					   int max_ndms,
					   int processed,
					   int nboots,
					   int num_trial_bins,
					   int navdms,
					   float narrow,
					   float wide,
					   int nsearch,
					   float aggression,
					   float cutoff,
					   float ***output_buffer,
					   int *ndms,
					   int *inBin,
					   float *dm_low,
					   float *dm_high,
					   float *dm_step,
					   float tsamp,
					   int enable_custom_fft,
					   int enable_inbin,
					   int enable_norm)
{

	fdas_params params;
	fdas_new_acc_sig acc_sig;
	cmd_args cmdargs;
	fdas_gpuarrays gpuarrays;
	fdas_cufftplan fftplans;
	float *acc_signal = NULL;
	struct timeval t_start, t_end;
	double t_gpu = 0.0, t_gpu_i = 0.0;

	//set default arguments
	cmdargs.nharms = 1; //
	cmdargs.nsig = 0; //
	cmdargs.duty = 0.10; //
	cmdargs.iter = 1; //
	cmdargs.writef = 1; //
	cmdargs.zval = 4; //
	cmdargs.mul = 1024; //
	cmdargs.wsig = 0; //
	cmdargs.search = 1; //
	cmdargs.thresh = 10.0; //
	cmdargs.freq0 = 100.5; //
	cmdargs.sigamp = 0.1; //
	cmdargs.basic = 0; //
	cmdargs.kfft = 1; //
	if (enable_custom_fft == 1){
		cmdargs.basic = 0; //
		cmdargs.kfft = 1; //
	}
	else{
		cmdargs.basic = 1; //
		cmdargs.kfft  = 0; //
	}
	//
	if (enable_inbin == 1)
		cmdargs.inbin = 1; //
	else
		cmdargs.inbin = 0; //
	//
	if (enable_norm == 1)
		cmdargs.norm = 1; //
	else
		cmdargs.norm = 0; //

	//get signal parameters
	acc_sig.nsamps = cmdargs.mul * 8192;  //
	acc_sig.freq0 = cmdargs.freq0; //
	acc_sig.mul = cmdargs.mul; 	//
	acc_sig.zval = cmdargs.zval; //
	acc_sig.nsig = cmdargs.nsig; //
	acc_sig.nharms = cmdargs.nharms; //
	acc_sig.duty = cmdargs.duty; //
	acc_sig.sigamp = cmdargs.sigamp; //

	int nearest = (int) floorf(log2f((float) processed));
	printf("\nnearest:\t%d", nearest);
	int samps = (int) powf(2.0, nearest);
	processed=samps;
	printf("\nsamps:\t%d", samps);

	params.nsamps = samps;

	/// Print params.h
	fdas_print_params_h();

	// prepare signal

	params.offset = presto_z_resp_halfwidth((double) ZMAX, 0); //array offset when we pick signal points for the overlp-save method
	printf(" Calculated overlap-save offsets: %d\n", params.offset);

	//
	params.sigblock = KERNLEN - 2 * params.offset + 1;
	params.scale = 1.0f / (float) (KERNLEN);
	params.rfftlen = params.nsamps / 2 + 1;
	params.nblocks = params.rfftlen / params.sigblock;
	params.siglen = params.nblocks * params.sigblock;
	params.extlen = params.nblocks * KERNLEN; //signal array extended to array of separate N=KERNLEN segments
	params.ffdotlen = params.siglen * NKERN; // total size of ffdot complex plane in fourier bins
	params.ffdotlen_cpx = params.extlen * NKERN; // total size of ffdot powers plane in fourier bins

	if (cmdargs.inbin)
		params.ffdotlen = params.ffdotlen * 2;

	if (cmdargs.search)
	{
		printf("\nnumber of convolution templates NKERN = %d, template length = %d, acceleration step in bins = %f, zmax = %d, template size power of 2 = %d, scale=%f \n",
				NKERN, KERNLEN, ACCEL_STEP, ZMAX, NEXP, params.scale);

		if (cmdargs.basic)
			printf("\nBasic algorithm:\n---------\n");

		else if (cmdargs.kfft)
			printf("\nCustom algorithm:\n-------\n");

		printf("\nnsamps = %d\ncpx signal length = %d\ntotal length: initial = %d, extended = %d\nconvolution signal segment length = %d\ntemplate length = %d\n# convolution blocks = %d\nffdot length = %u\n",
				params.nsamps, params.rfftlen, params.siglen, params.extlen,
				params.sigblock, KERNLEN, params.nblocks, params.ffdotlen);
		//
		if (cmdargs.basic)
			printf("ffdot length cpx pts = %u\n", params.ffdotlen_cpx);

		//memory required
		size_t mem_ffdot = params.ffdotlen * sizeof(float); // total size of ffdot plane powers in bytes
		size_t mem_ffdot_cpx = params.ffdotlen_cpx * sizeof(float2); // total size of ffdot plane powers in bytes
		size_t mem_kern_array = KERNLEN * NKERN * sizeof(float2);
		size_t mem_signals = (params.nsamps * sizeof(float)) + (params.siglen * sizeof(float2)) + (params.extlen * sizeof(float2));

		//Determining memory usage
		size_t mem_tot_needed = 0;
		double gbyte = 1024.0 * 1024.0 * 1024.0;
		float mbyte = gbyte / 1024.0;
		size_t mfree, mtotal;

		if (cmdargs.basic)
			mem_tot_needed = mem_ffdot + mem_ffdot_cpx + mem_kern_array + mem_signals;
		if (cmdargs.kfft)
			mem_tot_needed = mem_ffdot + mem_kern_array + mem_signals;
		checkCudaErrors(cudaMemGetInfo(&mfree, &mtotal));

		// get available memory info
		printf( "Total memory for this device: %.2f GB\nAvailable memory on this device for data upload: %.2f GB \n", mtotal / gbyte, mfree / gbyte);

		//Allocating gpu arrays
		gpuarrays.mem_insig = params.nsamps * sizeof(float);
		gpuarrays.mem_rfft = params.rfftlen * sizeof(float2);
		gpuarrays.mem_extsig = params.extlen * sizeof(float2);
		gpuarrays.mem_ffdot = mem_ffdot;
		gpuarrays.mem_ffdot_cpx = mem_ffdot_cpx;
		gpuarrays.mem_ipedge = params.nblocks * 2;

		printf("Total memory needed on GPU for arrays to process 1 DM: %.4f GB\nfloat ffdot plane (for power spectrum) = %.4f GB.\nTemplate array %.4f GB\nOne dimensional signals %.4f\n1 GB = %f",
				(double) mem_tot_needed / gbyte, (double) (mem_ffdot) / gbyte, (double) mem_kern_array / gbyte, (double) mem_signals / gbyte, gbyte);
		//
		if (mem_tot_needed >= (mfree - gbyte / 10)) {
			printf("\nNot enough memory available on the device to process this array\nPlease use a shorter signal\nExiting program...\n\n");
			exit(1);
		}
		getLastCudaError("\nCuda Error\n");

		fdas_alloc_gpu_arrays(&gpuarrays, &cmdargs);
		getLastCudaError("\nCuda Error\n");

		// Calculate kernel templates on CPU and upload-fft on GPU
		printf("\nCreating acceleration templates with KERNLEN=%d, NKERN = %d zmax=%d... ",	KERNLEN, NKERN, ZMAX);

		fdas_create_acc_kernels(gpuarrays.d_kernel, &cmdargs);
		printf(" done.\n");
		getLastCudaError("\nCuda Error\n");

		//Create cufft plans
		fdas_cuda_create_fftplans(&fftplans, &params);
		getLastCudaError("\nCuda Error\n");

		// Starting main acceleration search
		//cudaGetLastError(); //reset errors
		printf("\n\nStarting main acceleration search\n\n");

		int iter=cmdargs.iter;
		int titer=1;

		// FFT
		for (int i = 0; i < range; i++)
		{
			for (int dm_count = 0; dm_count < ndms[i] - 1; ++dm_count)
			{

				//first time PCIe transfer and print timing
				gettimeofday(&t_start, NULL); //don't time transfer
				checkCudaErrors( cudaMemcpy(gpuarrays.d_in_signal, output_buffer[i][dm_count], processed*sizeof(float), cudaMemcpyHostToDevice));
				//checkCudaErrors( cudaMemcpy(gpuarrays.d_in_signal, acc_signal, params.nsamps*sizeof(float), cudaMemcpyHostToDevice));

				cudaDeviceSynchronize();
				gettimeofday(&t_end, NULL);
				t_gpu = (double) (t_end.tv_sec + (t_end.tv_usec / 1000000.0)  - t_start.tv_sec - (t_start.tv_usec/ 1000000.0)) * 1000.0;
				t_gpu_i = (t_gpu /(double)titer);
				printf("\n\nAverage vector transfer time of %d float samples (%.2f Mb) from 1000 iterations: %f ms\n\n", params.nsamps, (float)(gpuarrays.mem_insig)/mbyte, t_gpu_i);

				cudaProfilerStart(); //exclude cuda initialization ops
				if(cmdargs.basic)
				{
					gettimeofday(&t_start, NULL); //don't time transfer
				    fdas_cuda_basic(&fftplans, &gpuarrays, &cmdargs, &params );
				    cudaDeviceSynchronize();
				    gettimeofday(&t_end, NULL);
				    t_gpu = (double) (t_end.tv_sec + (t_end.tv_usec / 1000000.0)  - t_start.tv_sec - (t_start.tv_usec/ 1000000.0)) * 1000.0;
				    t_gpu_i = (t_gpu / (double)iter);
				    printf("\n\nConvolution using basic algorithm with cuFFT\nTotal process took: %f ms per iteration \nTotal time %d iterations: %f ms\n", t_gpu_i, iter, t_gpu);
				 }

				 #ifndef NOCUST
				 if (cmdargs.kfft)
				 {
					 printf("\nMain: running FDAS with custom fft\n");
				     gettimeofday(&t_start, NULL); //don't time transfer
				  	fdas_cuda_customfft(&fftplans, &gpuarrays, &cmdargs, &params);

				     cudaDeviceSynchronize();
				     gettimeofday(&t_end, NULL);
				     t_gpu = (double) (t_end.tv_sec + (t_end.tv_usec / 1000000.0)  - t_start.tv_sec - (t_start.tv_usec/ 1000000.0)) * 1000.0;
				     t_gpu_i = (t_gpu / (double)iter);
				     printf("\n\nConvolution using custom FFT:\nTotal process took: %f ms\n per iteration \nTotal time %d iterations: %f ms\n", t_gpu_i, iter, t_gpu);
				 }
				 #endif
				 if(cmdargs.writef)
				 {
					 fdas_write_ffdot(&gpuarrays, &cmdargs, &params, dm_low[i], dm_count, dm_step[i]);
				 }
  				// Call sofias code here pass...
				// output_buffer[i][dm_count],
			}
		}
	}

	if (cmdargs.search)
	{
		cufftDestroy(fftplans.realplan);
	    cufftDestroy(fftplans.forwardplan);
	    // releasing GPU arrays
	    fdas_free_gpu_arrays(&gpuarrays, &cmdargs);
	}

}
