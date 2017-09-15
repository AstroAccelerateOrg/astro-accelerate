#include <stdio.h>
#include <stdlib.h>
#include <cufft.h>
#include <math.h>
#include "headers/params.h"

#include "headers/device_periodicity_parameters.h"
#include "headers/device_peak_find.h"
#include "headers/device_MSD_BLN_grid.h"
#include "headers/device_MSD_BLN_pw.h"
#include "headers/device_MSD_limited.h"
#include "headers/device_power.h"
#include "headers/device_harmonic_summing.h"

// define to see debug info
#define GPU_PERIODICITY_SEARCH_DEBUG

void Export_data_in_range(float *GPU_data, int nTimesamples, int nDMs, int DM_start, int DM_end, const char *filename, float dm_step, float dm_low) {
	std::ofstream FILEOUT;
	
	float *h_temp;
	h_temp = new float[nTimesamples*nDMs];
	cudaMemcpy(h_temp, GPU_data, nTimesamples*nDMs*sizeof(float), cudaMemcpyDeviceToHost);
	
	FILEOUT.open (filename, std::ofstream::out);
	if(DM_start==DM_end) DM_end++;
	for(int dm=DM_start; dm<DM_end; dm++){
		for(int t=0; t<nTimesamples; t++){
			FILEOUT << t << " " << (dm*dm_step + dm_low) << " " << GPU_data[dm*nTimesamples + t] << std::endl;
		}
		FILEOUT << std::endl;
		FILEOUT << std::endl;
	}
	FILEOUT.close();
	
	delete [] h_temp;
}

void Export_data_in_range(float2 *GPU_data, int nTimesamples, int nDMs, int DM_start, int DM_end, const char *filename, float dm_step, float dm_low) {
	std::ofstream FILEOUT;

	float *h_temp;
	h_temp = new float[nTimesamples*nDMs];
	cudaMemcpy(h_temp, GPU_data, nTimesamples*nDMs*sizeof(float), cudaMemcpyDeviceToHost);
	
	FILEOUT.open (filename, std::ofstream::out);
	if(DM_start==DM_end) DM_end++;
	for(int dm=DM_start; dm<DM_end; dm++){
		for(int t=0; t<nTimesamples; t++){
			FILEOUT << t << " " << (dm*dm_step + dm_low) << " " << GPU_data[dm*nTimesamples + t].x << " " << GPU_data[dm*nTimesamples + t].y << std::endl;
		}
		FILEOUT << std::endl;
		FILEOUT << std::endl;
	}
	FILEOUT.close();
	
	delete [] h_temp;
}

void Export_data_to_file(float *data, int size, const char *filename){
	FILE *fp_out;
	if (( fp_out = fopen(filename, "wb") ) == NULL)	{
		fprintf(stderr, "Error opening output file!\n");
		exit(0);
	}
	fwrite(data, size*sizeof(float), 4, fp_out);
	fclose(fp_out);
}

void GPU_periodicity(int range, int nsamp, int max_ndms, int processed, float sigma_cutoff, float ***output_buffer, int *ndms, int *inBin, float *dm_low, float *dm_high, float *dm_step, float tsamp, int nHarmonics) {
	// processed = maximum number of time-samples through out all ranges
	// nTimesamples = number of time-samples in given range 'i'
	// TODO:
	//     ->Be more clever regarding memory allocations for cuFFT use:
	//			const int NX = 1024;
	//			const int BATCH = 100000;
	//			size_t workSize;
	//			cufftEstimate1d(NX, CUFFT_C2C, BATCH, &workSize);
	// 	     or
	//			cufftHandle plan;
	//			cufftCreate(&plan);
	//			cufftGetSize1d(plan, NX, CUFFT_C2C, BATCH, &workSize);
	//      ->Use callbacks for power calculation
	//      ->Solve peak finding problem from batch to batch (we do not want to find peaks on shared borders)
	//      ->Interbinning which is not performed at the moment
	//      ->max_ndms is possibly the same thing as max_nDMs! Investigate.
	
	//------------------ Temporary for debuging
	int export_data=0;
	
	Periodicity_parameters per_param;
	per_param.assign(sigma_cutoff, nHarmonics);

//--------------------------------------------------------------------
//------> Starting Periodicity from scratch
	printf("\n");
	printf("------------ STARTING PERIODICITY SEARCH ------------\n\n");
	per_param.print_parameters();
	
	//---------> Initial stuff
	int nTimesamples, nDMs, max_nDMs, itemp;
	size_t max_nDMs_in_memory, max_nDMs_in_range;
	GpuTimer timer, periodicity_timer;
	double Total_periodicity_time = 0, Total_calc_time = 0, calc_time_per_range = 0, Total_copy_time = 0, copy_time_per_range = 0;
	
	periodicity_timer.Start();
	
	//---------> Finding nearest lower power of two (because of FFT algorithm)
	int nearest = (int) floorf(log2f((float) processed));
	nTimesamples = (int) powf(2.0, nearest);
	printf("Decreasing number of processed samples to nearest lower power of two, because of FFT algorithm...\n");
	printf("Number of processed timesamples: %d; nearest power of two: %d\n", processed, nTimesamples);
	processed = nTimesamples;
	
	// determining maximum number of DM trials in through all ranges
	if(range>0){
		max_nDMs = ndms[0];
		for(int f=0; f<range; f++){
			if(ndms[f]>max_nDMs) max_nDMs=ndms[f];
		}
	}
	printf("maximum number of DM trials through all ranges is %d\n", max_nDMs);
	
	//--------> Determining maximum number of DM trials we can fit into memory
	size_t free_mem,total_mem;
	cudaMemGetInfo(&free_mem,&total_mem);
	printf("     Memory available:%0.3f MB \n", ((float) free_mem)/(1024.0*1024.0) );
	max_nDMs_in_memory = (free_mem*0.95)/( (processed+2)*(5.5*sizeof(float) + 2*sizeof(ushort))); // 1 for real input real, 2 for complex output, 2 for complex cuFFT, 1 for peaks + 1 ushort
	if( (max_nDMs+32)<max_nDMs_in_memory) { //if we can fit all DM trials from largest range into memory then we need to find nearest higher multiple of PHS_NTHREADS
		itemp = (int) (max_nDMs/PHS_NTHREADS);
		if( (max_nDMs%PHS_NTHREADS)>0 ) itemp++;
		max_nDMs_in_memory = itemp*PHS_NTHREADS;
	}
	itemp = (int) (max_nDMs_in_memory/PHS_NTHREADS); // if we cannot fit all DM trials from largest range into memory we find nearest lower multiple of PHS_NTHREADS
	max_nDMs_in_memory = itemp*PHS_NTHREADS;
	printf("     Maximum number of DM trials which fit into memory is %d; Input plane size: %0.2f MB;\n", max_nDMs_in_memory, (((float) max_nDMs_in_memory*processed*sizeof(float))/(1024.0*1024.0)));
	
	
	//--------> Allocation of GPU memory. We allocate such amount of memory as to accommodate maximum number of DM trials from all ranges.
	unsigned int input_plane_size = (processed+2)*max_nDMs_in_memory;
	float *d_one_A; //for input and interbinned values
	if ( cudaSuccess != cudaMalloc((void **) &d_one_A,  sizeof(float)*input_plane_size )) printf("Periodicity Allocation error! d_one_A\n");
	
	float *d_two_B; //for cuFFT complex output and peaks
	if ( cudaSuccess != cudaMalloc((void **) &d_two_B,  sizeof(float)*2*input_plane_size )) printf("Periodicity Allocation error! d_two_B\n");
	
	float *d_half_C; // for power values
	if ( cudaSuccess != cudaMalloc((void **) &d_half_C,  sizeof(float)*input_plane_size/2 )) printf("Periodicity Allocation error! d_spectra_Real\n");
	
	ushort *d_power_harmonics, *d_interbin_harmonics;
	if ( cudaSuccess != cudaMalloc((void **) &d_power_harmonics, sizeof(ushort)*input_plane_size )) printf("Periodicity Allocation error! d_harmonics\n");
	if ( cudaSuccess != cudaMalloc((void **) &d_interbin_harmonics, sizeof(ushort)*input_plane_size )) printf("Periodicity Allocation error! d_harmonics\n");
	//-----------------------------------------------------------------------------
	cudaMemset((void*) d_power_harmonics, 0, sizeof(ushort)*input_plane_size);
	cudaMemset((void*) d_interbin_harmonics, 0, sizeof(ushort)*input_plane_size);
	//-----------------------------------------------------------------------------
	
	float *d_dedispersed_data, *d_FFT_complex_output, *d_frequency_power, *d_frequency_interbin, *d_frequency_power_CT, *d_frequency_interbin_CT, *d_power_SNR, *d_interbin_SNR, *d_power_list, *d_interbin_list;
	d_dedispersed_data      = d_one_A;
	d_FFT_complex_output    = d_two_B;
	d_frequency_power       = d_half_C;
	d_frequency_interbin    = d_one_A;
	d_frequency_power_CT    = &d_two_B[0];
	d_frequency_interbin_CT = &d_two_B[input_plane_size];
	d_power_SNR             = d_half_C;
	d_interbin_SNR          = d_one_A;
	d_power_list            = &d_two_B[0];
	d_interbin_list         = &d_two_B[input_plane_size];
	
	int *gmem_power_peak_pos;
	if ( cudaSuccess != cudaMalloc((void**) &gmem_power_peak_pos, 1*sizeof(int)) )  printf("Periodicity Allocation error! gmem_power_peak_pos\n");
	int *gmem_interbin_peak_pos;
	if ( cudaSuccess != cudaMalloc((void**) &gmem_interbin_peak_pos, 1*sizeof(int)) )  printf("Periodicity Allocation error! gmem_interbin_peak_pos\n");
	
	float *d_MSD;
	if ( cudaSuccess != cudaMalloc((void**) &d_MSD, sizeof(float)*4)) {printf("Periodicity Allocation error! d_MSD\n");}
	
	checkCudaErrors(cudaGetLastError());
	
	int local_max_list_size = (input_plane_size)/4; //maximum number of peaks per batch
	
	float *h_all_power_peaks, *h_all_interbin_peaks;
	h_all_power_peaks  = (float *)malloc(input_plane_size*sizeof(float));  // this might be too much, but it is very conservative assumption
	h_all_interbin_peaks  = (float *)malloc(input_plane_size*2*sizeof(float));  // this might be too much, but it is very conservative assumption
	size_t max_host_power_peaks = (input_plane_size)/4;
	size_t max_host_interbin_peaks = (input_plane_size*2)/4;
	size_t host_power_peak_pos;
	size_t host_interbin_peak_pos;
	int temp_host_power_peak_pos, temp_host_interbin_peak_pos;

	
	for (int i = 0; i < range; i++) {
		calc_time_per_range = 0; copy_time_per_range = 0;
		host_power_peak_pos = 0; host_interbin_peak_pos = 0;
		
		nTimesamples = processed/inBin[i];
		nDMs = ndms[i];
		max_nDMs_in_range = max_nDMs_in_memory*inBin[i];
		printf("Processing de-dispersion range:%f--%f:%f; inBin:%d; Timesamples:%d; DM trials:%d; max_nDMs:%d;\n", dm_low[i], dm_high[i], dm_step[i], inBin[i], nTimesamples, nDMs, max_nDMs_in_range);
		
		int nRepeats, nRest, DM_shift, DMs_per_cycle;
		std::vector<int> DM_list;
		
		//---------> Setting up batches
		nRepeats = nDMs/max_nDMs_in_range;
		nRest = nDMs - nRepeats*max_nDMs_in_range;
		for(int f=0; f<nRepeats; f++) DM_list.push_back(max_nDMs_in_range);
		if(nRest>0) DM_list.push_back(nRest);
		
		if(nRepeats>0) printf("     Periodicity search will run %d batches each containing %d DM trials. Remainder %d DM trials\n", (int) DM_list.size(), max_nDMs_in_range, nRest);
		else printf("     Periodicity search will run 1 batch containing %d DM trials.\n", nRest);
		
		if(DM_list.size()>0){
			DM_shift = 0;
			for(int dm=0; dm<DM_list.size(); dm++) {
				DMs_per_cycle = DM_list[dm];
				printf("\tBatch %d contains %d DM trials.\n",dm,DMs_per_cycle);
				
				cudaMemset((void*) gmem_power_peak_pos, 0, sizeof(int));
				cudaMemset((void*) gmem_interbin_peak_pos, 0, sizeof(int));
				
				
				//---------> Copy data from the host
				timer.Start();
				for(int ff=0; ff<DMs_per_cycle; ff++){
					checkCudaErrors( cudaMemcpy( &d_one_A[ff*nTimesamples], output_buffer[i][DM_shift + ff], nTimesamples*sizeof(float), cudaMemcpyHostToDevice));
				}
				timer.Stop();
				copy_time_per_range = copy_time_per_range + timer.Elapsed();
				//---------<
				
				
				//---------> cuFFT
				timer.Start();
				cufftHandle plan_input;
				cufftResult cufft_error;
				cufft_error = cufftPlan1d(&plan_input, nTimesamples, CUFFT_R2C, DMs_per_cycle);
				if ( cufft_error != CUFFT_SUCCESS) printf("CUFFT error: %d", cufft_error);
				//cufft_error = cufftExecR2C(plan_input, (cufftReal *)d_one_A, (cufftComplex *)d_two_B);
				cufft_error = cufftExecR2C(plan_input, (cufftReal *)d_dedispersed_data, (cufftComplex *)d_FFT_complex_output);
				if ( cufft_error != CUFFT_SUCCESS) printf("CUFFT error: %d", cufft_error);
				cufftDestroy(plan_input);
				timer.Stop();
				printf("     -> cuFFT took %f ms\n", timer.Elapsed());
				calc_time_per_range = calc_time_per_range + timer.Elapsed();
				//---------<
				
				//-----------------------------------------------------------------------------------
				if(i==0 && dm==0 && export_data) Export_data_in_range(d_one_A, nTimesamples, nDMs, 483, 493, "Input_data.dat", dm_step[i], dm_low[i]);
				//-----------------------------------------------------------------------------------
				
				//-----------------------------------------------------------------------------------
				if(i==0 && dm==0 && export_data) Export_data_in_range(d_two_B, ((nTimesamples>>1)+1), nDMs, 483, 493, "FFT_data.dat", dm_step[i], dm_low[i]);
				//-----------------------------------------------------------------------------------
				
				checkCudaErrors(cudaGetLastError());
				
				//---------> Calculate powers and interbinning
				timer.Start();
				//simple_power_and_interbin( (float2 *) d_two_B, d_half_C, d_one_A, nTimesamples, DMs_per_cycle);
				simple_power_and_interbin( (float2 *) d_FFT_complex_output, d_frequency_power, d_frequency_interbin, nTimesamples, DMs_per_cycle);
				timer.Stop();
				printf("     -> Calculation of powers and interbining took %f ms\n", timer.Elapsed());
				calc_time_per_range = calc_time_per_range + timer.Elapsed();
				//---------<
				
				//-----------------------------------------------------------------------------------
				if(i==0 && dm==0 && export_data) Export_data_in_range(d_half_C, nTimesamples/2, nDMs, 483, 493, "power_data.dat", dm_step[i], dm_low[i]);
				//-----------------------------------------------------------------------------------
				
				//-----------------------------------------------------------------------------------
				if(i==0 && dm==0 && export_data) Export_data_in_range(d_one_A, nTimesamples, nDMs, 483, 493, "Interbin_data.dat", dm_step[i], dm_low[i]);
				//-----------------------------------------------------------------------------------
				
				checkCudaErrors(cudaGetLastError());
				
				//---------> Mean and StDev on powers
				timer.Start();
				//MSD_limited(d_half_C, d_MSD, DMs_per_cycle, (nTimesamples>>1), 0);
				MSD_limited(d_frequency_power, d_MSD, DMs_per_cycle, (nTimesamples>>1), 0);
				//MSD_BLN_pw(d_half_C, d_MSD, DMs_per_cycle, nTimesamples/2, 0, sigma_constant);
				timer.Stop();
				printf("     -> MSD took %f ms\n", timer.Elapsed());
				calc_time_per_range = calc_time_per_range + timer.Elapsed();
				//---------<
				
				
				checkCudaErrors(cudaGetLastError());
				
				//---------> Corner turn
				timer.Start();
				//corner_turn_SM(d_half_C, &d_two_B[0], (nTimesamples>>1), DMs_per_cycle);
				//corner_turn_SM(d_one_A, &d_two_B[input_plane_size], nTimesamples, DMs_per_cycle);
				corner_turn_SM(d_frequency_power, d_frequency_power_CT, (nTimesamples>>1), DMs_per_cycle);
				corner_turn_SM(d_frequency_interbin, d_frequency_interbin_CT, nTimesamples, DMs_per_cycle);
				timer.Stop();
				printf("     -> corner turn took %f ms\n", timer.Elapsed());
				calc_time_per_range = calc_time_per_range + timer.Elapsed();
				//---------<
				
				checkCudaErrors(cudaGetLastError());
				
				//---------> Harmonic summing
				timer.Start();
				//periodicity_simple_harmonic_summing(&d_two_B[0], d_half_C, d_power_harmonics, d_MSD, (nTimesamples>>1), DMs_per_cycle, per_param.nHarmonics);
				//periodicity_simple_harmonic_summing(&d_two_B[input_plane_size], d_one_A, d_interbin_harmonics, d_MSD, nTimesamples, DMs_per_cycle, per_param.nHarmonics);
				periodicity_simple_harmonic_summing(d_frequency_power_CT, d_power_SNR, d_power_harmonics, d_MSD, (nTimesamples>>1), DMs_per_cycle, per_param.nHarmonics);
				periodicity_simple_harmonic_summing(d_frequency_interbin_CT, d_interbin_SNR, d_interbin_harmonics, d_MSD, nTimesamples, DMs_per_cycle, per_param.nHarmonics);
				timer.Stop();
				printf("     -> harmonic summing took %f ms\n", timer.Elapsed());
				calc_time_per_range = calc_time_per_range + timer.Elapsed();
				//---------<
				
				checkCudaErrors(cudaGetLastError());
				
				
				
				//---------> Peak finding
				timer.Start();
				//Peak_find_for_periodicity_search(d_half_C, d_power_harmonics, &d_two_B[0], (nTimesamples>>1), DMs_per_cycle, per_param.sigma_cutoff, local_max_list_size, gmem_power_peak_pos, DM_shift);
				//Peak_find_for_periodicity_search(d_one_A, d_interbin_harmonics, &d_two_B[input_plane_size], nTimesamples, DMs_per_cycle, per_param.sigma_cutoff, local_max_list_size, gmem_interbin_peak_pos, DM_shift);				
				Peak_find_for_periodicity_search(d_power_SNR, d_power_harmonics, d_power_list, (nTimesamples>>1), DMs_per_cycle, per_param.sigma_cutoff, local_max_list_size, gmem_power_peak_pos, DM_shift);
				Peak_find_for_periodicity_search(d_interbin_SNR, d_interbin_harmonics, d_interbin_list, nTimesamples, DMs_per_cycle, per_param.sigma_cutoff, local_max_list_size, gmem_interbin_peak_pos, DM_shift);
				/*
				SNR_limited(d_half_C, &d_two_B[0], d_power_harmonics, d_MSD, 1, DMs_per_cycle, (nTimesamples>>1), 0);
				SNR_limited(d_one_A, &d_two_B[input_plane_size], d_interbin_harmonics, d_MSD, 1, DMs_per_cycle, nTimesamples, 0);
				
				checkCudaErrors(cudaGetLastError());
				
				//-----------------------------------------------------------------------------------
				if(i==0 && dm==0 && export_data) Export_data_in_range(&d_two_B[0], nTimesamples/2, nDMs, 483, 493, "power_SNR_data.dat", dm_step[i], dm_low[i]);
				//-----------------------------------------------------------------------------------
				
				//-----------------------------------------------------------------------------------
				if(i==0 && dm==0 && export_data) Export_data_in_range(&d_two_B[input_plane_size], nTimesamples, nDMs, 483, 493, "Interbin_SNR_data.dat", dm_step[i], dm_low[i]);
				//-----------------------------------------------------------------------------------
				
				Peak_find_for_periodicity_search(&d_two_B[0], d_power_harmonics, d_half_C, DMs_per_cycle, (nTimesamples>>1), per_param.sigma_cutoff, local_max_list_size, gmem_power_peak_pos, DM_shift);
				Peak_find_for_periodicity_search(&d_two_B[input_plane_size], d_interbin_harmonics, d_one_A, DMs_per_cycle, nTimesamples, per_param.sigma_cutoff, local_max_list_size, gmem_interbin_peak_pos, DM_shift);
				*/
				timer.Stop();
				printf("     -> Peak finding took %f ms\n", timer.Elapsed());
				calc_time_per_range = calc_time_per_range + timer.Elapsed();
				//---------<
				
				checkCudaErrors(cudaGetLastError());
				
				//---------> Transferring peaks to the host
				timer.Start();
				
				checkCudaErrors(cudaMemcpy(&temp_host_power_peak_pos, gmem_power_peak_pos, sizeof(int), cudaMemcpyDeviceToHost));
				#ifdef GPU_PERIODICITY_SEARCH_DEBUG
				printf("     -> POWER: Total number of peaks found in this range is %d; maximum number of peaks:%d;\n", temp_host_power_peak_pos, local_max_list_size);
				#endif
				if( (host_power_peak_pos + temp_host_power_peak_pos)<max_host_power_peaks){
					//checkCudaErrors(cudaMemcpy(&h_all_power_peaks[host_power_peak_pos*4], d_half_C, temp_host_power_peak_pos*4*sizeof(float), cudaMemcpyDeviceToHost));
					//checkCudaErrors(cudaMemcpy(&h_all_power_peaks[host_power_peak_pos*4], &d_two_B[0], temp_host_power_peak_pos*4*sizeof(float), cudaMemcpyDeviceToHost));
					checkCudaErrors(cudaMemcpy(&h_all_power_peaks[host_power_peak_pos*4], d_power_list, temp_host_power_peak_pos*4*sizeof(float), cudaMemcpyDeviceToHost));
					host_power_peak_pos = host_power_peak_pos + temp_host_power_peak_pos;
				}
				else printf("     ->      Maximum list size reached! Increase list size or increase sigma cutoff.\n");
				
				checkCudaErrors(cudaMemcpy(&temp_host_interbin_peak_pos, gmem_interbin_peak_pos, sizeof(int), cudaMemcpyDeviceToHost));
				#ifdef GPU_PERIODICITY_SEARCH_DEBUG
				printf("     -> INTERBIN: Total number of peaks found in this range is %d; maximum number of peaks:%d;\n", temp_host_interbin_peak_pos, local_max_list_size);
				#endif
				if( (host_interbin_peak_pos + temp_host_interbin_peak_pos)<max_host_interbin_peaks){
					//checkCudaErrors(cudaMemcpy(&h_all_interbin_peaks[host_interbin_peak_pos*4], d_one_A, temp_host_interbin_peak_pos*4*sizeof(float), cudaMemcpyDeviceToHost));
					//checkCudaErrors(cudaMemcpy(&h_all_interbin_peaks[host_interbin_peak_pos*4], &d_two_B[input_plane_size], temp_host_interbin_peak_pos*4*sizeof(float), cudaMemcpyDeviceToHost));
					checkCudaErrors(cudaMemcpy(&h_all_interbin_peaks[host_interbin_peak_pos*4], d_interbin_list, temp_host_interbin_peak_pos*4*sizeof(float), cudaMemcpyDeviceToHost));
					host_interbin_peak_pos = host_interbin_peak_pos + temp_host_interbin_peak_pos;
				}
				else printf("     ->      Maximum list size reached! Increase list size or increase sigma cutoff.\n");
				
				timer.Stop();
				copy_time_per_range = copy_time_per_range + timer.Elapsed();
				//---------<
				
				DM_shift = DM_shift + DMs_per_cycle;
			} // end of for through batches
			
			//---------> Peak processing on host and export
			char filename[200];
			if(host_power_peak_pos>0){
				#pragma omp parallel for
				for (int count = 0; count < host_power_peak_pos; count++){
					h_all_power_peaks[4*count]     = h_all_power_peaks[4*count]*dm_step[i] + dm_low[i];
					h_all_power_peaks[4*count + 1] = h_all_power_peaks[4*count + 1]*(1.0/(tsamp*nTimesamples));
				}
				
				sprintf(filename, "fourier-dm_%.2f-%.2f.dat", dm_low[i], dm_high[i]);
				Export_data_to_file(h_all_power_peaks, host_power_peak_pos, filename);
			}
			
			if(host_interbin_peak_pos>0){
				#pragma omp parallel for
				for (int count = 0; count < host_interbin_peak_pos; count++){
					h_all_interbin_peaks[4*count]     = h_all_interbin_peaks[4*count]*dm_step[i] + dm_low[i];
					h_all_interbin_peaks[4*count + 1] = h_all_interbin_peaks[4*count + 1]*(1.0/(tsamp*nTimesamples));
				}
				
				sprintf(filename, "fourier_inter-dm_%.2f-%.2f.dat", dm_low[i], dm_high[i]);
				Export_data_to_file(h_all_interbin_peaks, host_interbin_peak_pos, filename);
			}
			//---------<
		}
		
		printf("     -----------------------\n");
		printf("     -> This range calculation time: %f ms\n", calc_time_per_range);
		printf("     -> This range copy time:        %f ms\n", copy_time_per_range);
		printf("\n");
		Total_calc_time = Total_calc_time + calc_time_per_range;
		calc_time_per_range = 0;
		Total_copy_time = Total_copy_time + copy_time_per_range;
		copy_time_per_range = 0;
	}

	periodicity_timer.Stop();
	Total_periodicity_time = periodicity_timer.Elapsed();
	
	printf("\nTimer:\n");
	printf("Total calculation time: %f ms\n", Total_calc_time);
	printf("Total copy time:        %f ms\n", Total_copy_time);
	printf("Total periodicity time: %f ms\n", Total_periodicity_time);

	cudaDeviceSynchronize();
	
	cudaFree(d_MSD);
	cudaFree(d_one_A);
	cudaFree(d_two_B);
	cudaFree(d_half_C);
	cudaFree(d_power_harmonics);
	cudaFree(d_interbin_harmonics);
	cudaFree(gmem_power_peak_pos);
	cudaFree(gmem_interbin_peak_pos);
	
	free(h_all_power_peaks);
	free(h_all_interbin_peaks);
	
	cudaFree(d_dedispersed_data);
	cudaFree(d_FFT_complex_output);
	cudaFree(d_frequency_power);
	cudaFree(d_frequency_interbin);
	cudaFree(d_frequency_power_CT);
	cudaFree(d_frequency_interbin_CT);
	cudaFree(d_power_SNR);
	cudaFree(d_interbin_SNR);
	cudaFree(d_power_list);
	cudaFree(d_interbin_list);
	
	
	
	
	
	
	
}


