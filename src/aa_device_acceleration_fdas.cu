#include "aa_device_acceleration_fdas.hpp"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <cufft.h>
#include <errno.h>
#include <string.h>
#include <sys/time.h>
#include <cuda_profiler_api.h>

#include "aa_log.hpp"
#include "aa_device_harmonic_summing.hpp"
#include "aa_params.hpp"
#include "aa_fdas_test_parameters.hpp"
#include "aa_fdas.hpp"
#include "aa_fdas_util.hpp"
#include "aa_fdas_host.hpp"
#include "aa_device_MSD.hpp"
#include "aa_device_MSD_plane_profile.hpp"
#include "aa_device_peak_find.hpp"
#include "presto_funcs.hpp"
#include "presto.hpp"

//#define FDAS_DEBUG
//#define FDAS_ACC_SIG_TEST

namespace astroaccelerate {

/**
* \brief Function that performs a fourier domain accelerated search (fdas).
* \brief Users should not interact with this function directly. Instead they should use aa_fdas_plan and aa_fdas_strategy.
*/
void acceleration_fdas(
	int range,
	int nsamp,
	int max_ndms,
	int processed,
	float cutoff,
	float ***output_buffer,
	int const*const ndms,
	int *inBin,
	float *dm_low,
	float *dm_high,
	float *dm_step,
	float tsamp,
	const int fdas_max_nHarmonics,
	const bool enable_custom_fft,
	const bool enable_inbin,
	const bool enable_norm,
	float sigma_constant,
	const bool enable_output_ffdot_plan,
	const bool enable_output_fdas_list,
	const bool enable_harmonic_sum
) {

    astroaccelerate::fdas_params params;
    // fdas_new_acc_sig acc_sig;
    astroaccelerate::cmd_args cmdargs;
    astroaccelerate::fdas_gpuarrays gpuarrays;
    astroaccelerate::fdas_cufftplan fftplans;
    //float *acc_signal = NULL;
    struct timeval t_start, t_end;
    //double t_gpu = 0.0, t_gpu_i = 0.0;

    //set default arguments
    cmdargs.nharms = fdas_max_nHarmonics; //
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
    if (enable_custom_fft){
      cmdargs.basic = 0; //
      cmdargs.kfft = 1; //
    }
    else{
      cmdargs.basic = 1; //
      cmdargs.kfft  = 0; //
    }
    //
    if (enable_inbin)
      cmdargs.inbin = 1; //
    else
      cmdargs.inbin = 0; //
    //
    if (enable_norm)
      cmdargs.norm = 1; //
    else
      cmdargs.norm = 0; //

    //get signal parameters
    /*acc_sig.nsamps = cmdargs.mul * 8192;  //
      acc_sig.freq0 = cmdargs.freq0; //
      acc_sig.mul = cmdargs.mul; 	//
      acc_sig.zval = cmdargs.zval; //
      acc_sig.nsig = cmdargs.nsig; //
      acc_sig.nharms = cmdargs.nharms; //
      acc_sig.duty = cmdargs.duty; //
      acc_sig.sigamp = cmdargs.sigamp; //
    */
    int nearest = (int) floorf(log2f((float) processed));
    printf("\nnearest:\t%d", nearest);
    int samps = (int) powf(2.0, nearest);
    processed=samps;
    printf("\nsamps:\t%d", samps);

    params.nsamps = samps;
    params.tsamp = tsamp;

    /// Print params.h
    astroaccelerate::fdas_print_params_h();

    // prepare signal
    params.offset = presto_z_resp_halfwidth((double) ZMAX, 0); //array offset when we pick signal points for the overlp-save method
    printf(" Calculated overlap-save offsets: %d\n", params.offset);

    //
    params.sigblock = KERNLEN - 2 * params.offset + 1;
    params.scale = sqrt(2) / (float) (KERNLEN);
    params.rfftlen = params.nsamps / 2 + 1;
    params.nblocks = params.rfftlen / params.sigblock;
    params.siglen = params.nblocks * params.sigblock;
    params.extlen = params.nblocks * KERNLEN; //signal array extended to array of separate N=KERNLEN segments
    params.ffdotlen = (unsigned int)params.siglen * (unsigned int)NKERN; // total size of ffdot complex plane in fourier bins
    params.ffdotlen_cpx = params.extlen * NKERN; // total size of ffdot powers plane in fourier bins
    params.max_list_length = params.ffdotlen/4;
	
    if (cmdargs.inbin)
      params.ffdotlen = params.ffdotlen * 2;

    if (cmdargs.search) {
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
		size_t mem_max_list_size = params.max_list_length*4*sizeof(float);
		size_t mem_signals = (params.nsamps * sizeof(float)) + (params.siglen * sizeof(float2)) + (params.extlen * sizeof(float2));

		//Determining memory usage
		size_t mem_tot_needed = 0;
		double gbyte = 1024.0 * 1024.0 * 1024.0;
		float mbyte = gbyte / 1024.0;
		size_t mfree, mtotal;

		if (cmdargs.basic)
		mem_tot_needed = mem_ffdot + mem_ffdot_cpx + mem_kern_array + mem_signals + mem_max_list_size; // KA added + mem_max_list_size
		if (cmdargs.kfft)
		mem_tot_needed = mem_ffdot + mem_kern_array + mem_signals + mem_max_list_size; // KA added + mem_max_list_size
		cudaError_t e = cudaMemGetInfo(&mfree, &mtotal);

		if(e != cudaSuccess) {
		LOG(log_level::error, "Could not cudaMemGetInfo in aa_device_acceleration_fdas.cu (" + std::string(cudaGetErrorString(e)) + ")");
		}

		// get available memory info
		printf( "Total memory for this device: %.2f GB\nAvailable memory on this device for data upload: %.2f GB \n", mtotal / gbyte, mfree / gbyte);

		//Allocating gpu arrays
		gpuarrays.mem_insig = params.nsamps * sizeof(float);
		gpuarrays.mem_rfft = params.rfftlen * sizeof(float2);
		gpuarrays.mem_extsig = params.extlen * sizeof(float2);
		gpuarrays.mem_ffdot = mem_ffdot;
		gpuarrays.mem_ffdot_cpx = mem_ffdot_cpx;
		gpuarrays.mem_ipedge = params.nblocks * 2;
		gpuarrays.mem_max_list_size = mem_max_list_size;

		printf("Total memory needed on GPU for arrays to process 1 DM: %.4f GB\nfloat ffdot plane (for power spectrum) = %.4f GB.\nTemplate array %.4f GB\nOne dimensional signals %.4f\n1 GB = %f",
			(double) mem_tot_needed / gbyte, (double) (mem_ffdot) / gbyte, (double) mem_kern_array / gbyte, (double) mem_signals / gbyte, gbyte);
		//
		if (mem_tot_needed >= (mfree - gbyte / 10)) {
		printf("\nNot enough memory available on the device to process this array\nPlease use a shorter signal\nExiting program...\n\n");
		exit(1);
		}
		//getLastCudaError("\nCuda Error\n");

		fdas_alloc_gpu_arrays(&gpuarrays, &cmdargs);
		//getLastCudaError("\nCuda Error\n");

		// Calculate kernel templates on CPU and upload-fft on GPU
		printf("\nCreating acceleration templates with KERNLEN=%d, NKERN = %d zmax=%d... ",	KERNLEN, NKERN, ZMAX);

		fdas_create_acc_kernels(gpuarrays.d_kernel, &cmdargs);
		printf(" done.\n");
		//getLastCudaError("\nCuda Error\n");

		//Create cufft plans
		fdas_cuda_create_fftplans(&fftplans, &params);
		//getLastCudaError("\nCuda Error\n");

		printf("\n\nStarting main acceleration search\n\n");

		// FFT
		for (int i = 0; i < range; i++) {
			processed=samps/inBin[i];
			printf("Processing range %d containing %d DM steps.\n", range, ndms[i]);
			for (int dm_count = 0; dm_count < ndms[i] - 1; ++dm_count) {

				//first time PCIe transfer and print timing
				gettimeofday(&t_start, NULL); //don't time transfer
						
				//!TEST!: put test signal here
				#ifdef FDAS_CONV_TEST
				printf("\n************** TEST FOR FDAS ***********************\n");
				srand(time(NULL));
				for(int f=0; f<processed; f++) output_buffer[i][dm_count][f]=rand() / (float)RAND_MAX;
						
				if (processed>15000){
					for(int f=15000; f<processed; f++){
						output_buffer[i][dm_count][f] = (f%FDAS_TEST_TOOTH_LENGTH)/500.0;
					}
					
					for(int f=0; f<192; f++){
					output_buffer[i][dm_count][f + 5300] = 10.0;
					}
								
					for(int f=0; f<128; f++){
					output_buffer[i][dm_count][f + 8626] = 10.0;
					}
								
					for(int f=0; f<36; f++){
					output_buffer[i][dm_count][f + 9626] = 10.0;
					}
								
					for(int f=0; f<83; f++){
					output_buffer[i][dm_count][f + 10626] = 10.0;
					}
								
					for(int f=0; f<138; f++){
					output_buffer[i][dm_count][f + 11626] = 10.0;
					}
				}
				
				std::ofstream FILEOUT;
				FILEOUT.open ("acc_conv_test_input_signal.dat", std::ofstream::out);
				for(int f=0; f<processed; f++){
					FILEOUT << output_buffer[i][dm_count][f] << std::endl;
				}
				FILEOUT.close();
				#endif
						
				#ifdef FDAS_ACC_SIG_TEST
				double acc_sig_snr = 1.0;
				fdas_new_acc_sig acc_sig;
						
				acc_sig.freq0 = FDAS_TEST_FREQUENCY;
				acc_sig.nsamps = processed;
				acc_sig.zval = FDAS_TEST_ZVALUE;
				acc_sig.nharms = FDAS_TEST_HAMONICS;
				acc_sig.duty = FDAS_TEST_DUTY_CYCLE/100.0;
				acc_sig.sigamp = FDAS_TEST_SIGNAL_AMPLITUDE;
						
				double t0, tau;
				double omega = 2*M_PI*acc_sig.freq0;
				double accel;
				double tobs;
				double sampling_rate = 0.000064;
				double light_speed = 2.99792458e8;
						
				tobs = (double) (sampling_rate*acc_sig.nsamps);
				accel = ((double)acc_sig.zval * light_speed) / (acc_sig.freq0*tobs*tobs);
				printf("\n\npreparing test signal, observation time = %f s, %d nsamps f0 = %f Hz with %d harmonics\n", tobs, acc_sig.nsamps, acc_sig.freq0, acc_sig.nharms);
				printf("\nz = %d accelereation = %f m/s^2\n", acc_sig.zval, accel);
						
				printf("\nNow creating accelerated signal with fc=%f, accel=%f, harmonics=%d, duty cycle=%.1f, noise=%f signal samples=%d, signal level: %.2f\n", acc_sig.freq0, accel, acc_sig.nharms, acc_sig.duty*100.0, acc_sig_snr, acc_sig.nsamps,acc_sig.sigamp);
						
				for ( int sd=0; sd<acc_sig.nsamps; ++sd){
					t0 = sd*sampling_rate;
					tau = t0 + (t0*(accel*t0) / light_speed /2.0);
					if (acc_sig_snr!=0){
						output_buffer[i][dm_count][sd] = 0;
					}
					for (int j = 1; j <= acc_sig.nharms; ++j){
						output_buffer[i][dm_count][sd] += (2.0/(j*M_PI)*sin(j*M_PI*acc_sig.duty))*acc_sig.sigamp*cos(j*omega*tau); 
					}
				}
						
				#endif
				//!TEST!: put test signal here
						
				e = cudaMemcpy(gpuarrays.d_in_signal, output_buffer[i][dm_count], processed*sizeof(float), cudaMemcpyHostToDevice);

				if(e != cudaSuccess) {
				LOG(log_level::error, "Could not cudaMemcpy in aa_device_acceleration.cu (" + std::string(cudaGetErrorString(e)) + ")");
				}

				cudaDeviceSynchronize();
				gettimeofday(&t_end, NULL);
				//t_gpu = (double) (t_end.tv_sec + (t_end.tv_usec / 1000000.0)  - t_start.tv_sec - (t_start.tv_usec/ 1000000.0)) * 1000.0;
				//t_gpu_i = (t_gpu /(double)titer);
				//printf("\n\nAverage vector transfer time of %d float samples (%.2f Mb) from 1000 iterations: %f ms\n\n", params.nsamps, (float)(gpuarrays.mem_insig)/mbyte, t_gpu_i);

				cudaProfilerStart(); //exclude cuda initialization ops
				if(cmdargs.basic) {
					gettimeofday(&t_start, NULL); //don't time transfer
					fdas_cuda_basic(&fftplans, &gpuarrays, &cmdargs, &params );
					cudaDeviceSynchronize();
					gettimeofday(&t_end, NULL);
					//t_gpu = (double) (t_end.tv_sec + (t_end.tv_usec / 1000000.0)  - t_start.tv_sec - (t_start.tv_usec/ 1000000.0)) * 1000.0;
					//t_gpu_i += (t_gpu / (double)iter);
					//printf("Convolution using basic algorithm with cuFFT: Total process took: %f ms per iteration; Total time %d iterations: %f ms;\n", t_gpu_i, iter, t_gpu);
				}

				if (cmdargs.kfft) {
					gettimeofday(&t_start, NULL); //don't time transfer
					fdas_cuda_customfft(&fftplans, &gpuarrays, &cmdargs, &params);
					cudaDeviceSynchronize();
					gettimeofday(&t_end, NULL);
					//t_gpu = (double) (t_end.tv_sec + (t_end.tv_usec / 1000000.0)  - t_start.tv_sec - (t_start.tv_usec/ 1000000.0)) * 1000.0;
					//t_gpu_i += (t_gpu / (double)iter);
					//printf("Convolution using custom FFT: Total process took: %f ms per iteration; Total time %d iterations: %f ms;\n", t_gpu_i, iter, t_gpu);
				}

				// Calculating base level noise and peak find
				if(cmdargs.basic || cmdargs.kfft){
					int ibin=1;
					if (cmdargs.inbin) ibin=2;

					unsigned int list_size;
					unsigned int *gmem_fdas_peak_pos;
					if ( cudaSuccess != cudaMalloc((void**) &gmem_fdas_peak_pos, 1*sizeof(int))) printf("Allocation error!\n");
					cudaMemset((void*) gmem_fdas_peak_pos, 0, sizeof(int));

					//printf("Dimensions for BLN: ibin:%d; siglen:%d;\n", ibin, params.siglen);

					//checkCudaErrors(cudaGetLastError());

					// Structure of the gpuarrays.d_ffdot_pwr is as follows
					// Max negative acceleration
					// Zero acceleration; start of the zero acceleration is: 
					// Max positive acceleration
					
					// Call the 2D Harmonic Sum kernel
					// TODO: Fix input that is going into harmonic sum and related number of NKERN
					// Output of the harmonic sum is 2d data frequency x acceleration where frequency is fastest changing
					// Size of the output is max_f_idx x max_fdot_idx
					if(enable_harmonic_sum){
						std::vector<int> boxcarwidths{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32};
						float *h_MSD_interpolated = (float*) malloc(boxcarwidths.size()*2*sizeof(float));
						float *d_workarea;
						if (cudaSuccess != cudaMalloc((void**) &d_workarea, sizeof(float)*NKERN*params.siglen*ibin)) LOG(log_level::error, "Allocation error! W " + std::string(cudaGetErrorString(e)));
						float *d_MSD_interpolated;
						if (cudaSuccess != cudaMalloc((void**) &d_MSD_interpolated, boxcarwidths.size()*2*sizeof(float) )  ) printf("Allocation error! N \n");
						
						double total_time = 0;
						double dit_time = 0;
						double MSD_time = 0;
						float *d_MSD_DIT = NULL;
					
						// d_MSD_interpolated memory structure is as follows: 
						// mean[n] and std[n] represent mean and standard deviation for partial sum of n elements
						// Structure: { mean[1], std[1], mean[2], std[2], ... , mean[last], std[last]}
						MSD_plane_profile(d_MSD_interpolated, gpuarrays.d_ffdot_pwr, d_MSD_DIT, d_workarea, false, ibin*params.siglen, NKERN, &boxcarwidths, 0, dm_low[i], dm_high[i], sigma_constant, 0, false, &total_time, &dit_time, &MSD_time);

						e = cudaMemcpy(h_MSD_interpolated, d_MSD_interpolated, boxcarwidths.size()*2*sizeof(float), cudaMemcpyDeviceToHost);
						if(e != cudaSuccess) {
							LOG(log_level::error, "Could not cudaMemcpy d_MSD_interpolated in aa_device_acceleration_fdas.cu (" + std::string(cudaGetErrorString(e)) + ")");
						}
						#ifdef FDAS_DEBUG
						for(int f=0; f<(int)boxcarwidths.size();f++){
							printf("boxcar size: %d; MSD=[%f; %f].\n", f, h_MSD_interpolated[2*f], h_MSD_interpolated[2*f+1]);
						}
						#endif


						size_t zero_position = ((NKERN-1)/2)*ibin*params.siglen;
						size_t nFreqBins = ibin*params.siglen;
						float pos_in_dm = dm_count*dm_step[i] + dm_low[i];
						size_t max_f_idx = ibin*params.siglen/cmdargs.nharms;
						size_t max_fdot_idx = ((NKERN-1)/2) + 1;
						size_t max_half_plane_pos = (((NKERN-1)/2)+1)*max_f_idx;
						
						// Positive half
						periodicity_two_dimensional_greedy_harmonic_summing(
							&gpuarrays.d_ffdot_pwr[zero_position], 
							&gpuarrays.d_ffdot_max[max_half_plane_pos], 
							&gpuarrays.d_ffdot_SNR[max_half_plane_pos], 
							&gpuarrays.d_ffdot_Harmonics[max_half_plane_pos], 
							d_MSD_interpolated, 
							nFreqBins, 
							max_fdot_idx, 
							max_f_idx, 
							max_fdot_idx, 
							cmdargs.nharms
						);
						
						//Negative half
						int nBlocks_x, nBlocks_y, nThreads;
						dim3 gridSize, blockSize;
						nThreads = 256;
						size_t nNegative_acc = ((NKERN-1)/2);
						cudaError_t cuda_error;
						
						cuda_error = cudaMemset(gpuarrays.d_ffdot_max, 0, max_half_plane_pos*sizeof(float));
						cuda_error = cudaMemset(gpuarrays.d_ffdot_SNR, 0, max_half_plane_pos*sizeof(float));
						cuda_error = cudaMemset(gpuarrays.d_ffdot_Harmonics, 0, max_half_plane_pos*sizeof(ushort));
						
						// Negative half plane
						gridSize.x = (nFreqBins + nThreads - 1)/nThreads;
						gridSize.y = (nNegative_acc + 1)/2;
						gridSize.z = 1;
						blockSize.x = nThreads;
						blockSize.y = 1;
						blockSize.z = 1;
						
						call_kernel_flip_negative_ffdot_inplace_float(
							gridSize, 
							blockSize, 
							gpuarrays.d_ffdot_pwr, 
							nNegative_acc + 1, 
							nFreqBins
						);
						
						periodicity_two_dimensional_greedy_harmonic_summing(
							gpuarrays.d_ffdot_pwr, 
							gpuarrays.d_ffdot_max, 
							gpuarrays.d_ffdot_SNR, 
							gpuarrays.d_ffdot_Harmonics, 
							d_MSD_interpolated, 
							nFreqBins, 
							max_fdot_idx, 
							max_f_idx, 
							max_fdot_idx, 
							cmdargs.nharms
						);

						// Correction of the output power, SNR and harm planes
						// So they are compatible with peak-finding.
						gridSize.x = (max_f_idx + nThreads - 1)/nThreads;
						gridSize.y = (max_fdot_idx)/2;
						gridSize.z = 1;
						blockSize.x = nThreads;
						blockSize.y = 1;
						blockSize.z = 1;
						
						call_kernel_flip_negative_ffdot_inplace_float(
							gridSize, 
							blockSize, 
							gpuarrays.d_ffdot_max, 
							max_fdot_idx, 
							max_f_idx
						);
						call_kernel_flip_negative_ffdot_inplace_float(
							gridSize, 
							blockSize, 
							gpuarrays.d_ffdot_SNR, 
							max_fdot_idx, 
							max_f_idx
						);
						call_kernel_flip_negative_ffdot_inplace_ushort(
							gridSize, 
							blockSize, 
							gpuarrays.d_ffdot_Harmonics, 
							max_fdot_idx, 
							max_f_idx
						);
						
						// Peak find for harmonic summed values!
						peak_find_fdas_harm(
							gpuarrays.d_fdas_peak_list,
							gpuarrays.d_ffdot_max, 
							gpuarrays.d_ffdot_SNR, 
							gpuarrays.d_ffdot_Harmonics, 
							max_f_idx, 
							NKERN, 
							((NKERN-1)/2),
							cmdargs.thresh, 
							params.max_list_length,
							gmem_fdas_peak_pos,
							dm_count*dm_step[i] + dm_low[i]
						);
						
						cuda_error = cudaMemcpy(&list_size, gmem_fdas_peak_pos, sizeof(unsigned int), cudaMemcpyDeviceToHost);
						if(cuda_error != cudaSuccess) {
							LOG(log_level::error, "Cannot perform cudaMemcpyDeviceToHost for gmem_fdas_peak_pos (" + std::string(cudaGetErrorString(cuda_error)) + ")");
						}
						
						fdas_write_list_harm(
							&gpuarrays, 
							&cmdargs, 
							&params, 
							h_MSD_interpolated,
							dm_count*dm_step[i] + dm_low[i], 
							list_size
						);
						
						#ifdef FDAS_ACC_SIG_TEST
						{
							char filename[200];
							sprintf(filename, "acc_fdas_2d_harm_whole_plane_%f.dat", dm_count*dm_step[i] + dm_low[i]);
							fdas_write_test_ffdot_harmonic(
								gpuarrays.d_ffdot_max, 
								gpuarrays.d_ffdot_SNR, 
								gpuarrays.d_ffdot_Harmonics, 
								max_f_idx, 
								NKERN, 
								((NKERN-1)/2),
								&params,
								&cmdargs,
								filename
							);
							exit(1);
						}
						#endif
						
						cuda_error = cudaMemset(gpuarrays.d_fdas_peak_list, 0, mem_max_list_size);
						if(cuda_error != cudaSuccess) {
							LOG(log_level::error, "Error setting d_fdas_peak_list to 0 (" + std::string(cudaGetErrorString(cuda_error)) + ")");
						};
						cuda_error = cudaMemset(gmem_fdas_peak_pos, 0, sizeof(unsigned int));
						if(cuda_error != cudaSuccess) {
							LOG(log_level::error, "Error setting gmem_fdas_peak_pos to 0 (" + std::string(cudaGetErrorString(cuda_error)) + ")");
						};

						free(h_MSD_interpolated);
						cudaFree(d_MSD_interpolated);
						cudaFree(d_workarea);
					}
					//============= END OF HARMONIC SUMMING ============
					else {
						//================> No harmonic sum
						float *d_MSD;
						if ( cudaSuccess != cudaMalloc((void**) &d_MSD, sizeof(float)*3)) printf("Allocation error!\n");
						float h_MSD[3];
						
						Find_MSD(d_MSD, gpuarrays.d_ffdot_pwr, params.siglen*ibin, NKERN, 0, sigma_constant, 0);
						
						PEAK_FIND_FOR_FDAS(gpuarrays.d_ffdot_pwr, gpuarrays.d_fdas_peak_list, d_MSD, NKERN, ibin*params.siglen, cmdargs.thresh, params.max_list_length, gmem_fdas_peak_pos, dm_count*dm_step[i] + dm_low[i]);		
						e = cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost);
						if(e != cudaSuccess) {
							LOG(log_level::error, "Could not cudaMemcpy in aa_device_acceleration_fdas.cu (" + std::string(cudaGetErrorString(e)) + ")");
						}
						#ifdef FDAS_DEBUG
							printf("FDAS: Mean and standard deviation h_MSD=[%f; %f; %f]\n", h_MSD[0], h_MSD[1], h_MSD[2]);
						#endif
						
						e = cudaMemcpy(&list_size, gmem_fdas_peak_pos, sizeof(unsigned int), cudaMemcpyDeviceToHost);
						if(e != cudaSuccess) {
							LOG(log_level::error, "Could not cudaMemcpy in aa_device_acceleration_fdas.cu (" + std::string(cudaGetErrorString(e)) + ")");
						}
						#ifdef FDAS_DEBUG
							printf("FDAS: Number of peaks found %u.", list_size);
						#endif
						
						#ifdef FDAS_ACC_SIG_TEST
							fdas_write_list(&gpuarrays, &cmdargs, &params, h_MSD, dm_low[i], dm_count, dm_step[i], list_size);
							fdas_write_ffdot(&gpuarrays, &cmdargs, &params, dm_low[i], dm_count, dm_step[i]);
							exit(1);
						#endif
						
						if (enable_output_fdas_list) {
							if(list_size>0){
								fdas_write_list(&gpuarrays, &cmdargs, &params, h_MSD, dm_low[i], dm_count, dm_step[i], list_size);
							}
						}
						cudaFree(d_MSD);
					}
					
					//!TEST!: do not perform peak find instead export the thing to file.
					#ifdef FDAS_CONV_TEST
					fdas_write_test_ffdot(&gpuarrays, &cmdargs, &params, dm_low[i], dm_count, dm_step[i]);
					exit(1);
					#endif
					//!TEST!: do not perform peak find instead export the thing to file.
					
					
					cudaFree(gmem_fdas_peak_pos);
				}
				
				//Output ffdot planea
				if (enable_output_ffdot_plan) {
					fdas_write_ffdot(&gpuarrays, &cmdargs, &params, dm_low[i], dm_count, dm_step[i]);
				}
			} // for loop -> DMs
		} // for loop -> range
	} // if search

    if (cmdargs.search) {
		cufftDestroy(fftplans.realplan);
		cufftDestroy(fftplans.forwardplan);
		fdas_free_gpu_arrays(&gpuarrays, &cmdargs);
	}

} // acceleration_fdas

} //namespace astroaccelerate
