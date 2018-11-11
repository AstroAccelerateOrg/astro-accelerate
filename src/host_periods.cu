#include <stdio.h>
#include <stdlib.h>
#include <cufft.h>
#include <math.h>
#include "params.hpp"

namespace astroaccelerate {

// define to see debug info
#define GPU_PERIODICITY_SEARCH_DEBUG

void periodicity(int range, int nsamp, int max_ndms, int processed, int nboots, int num_trial_bins, int navdms, float narrow, float wide, int nsearch, float aggression, float cutoff, float ***output_buffer, int *ndms, int *inBin, float *dm_low, float *dm_high, float *dm_step, float tsamp)
{
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

	/*
//--------------------------------------------------------------------
//------> Starting Periodicity from scratch
	printf("\n");
	printf("------------ STARTING PERIODICITY SEARCH ------------\n\n");
	
	//---------> Initial stuff
	int nTimesamples, nDMs, max_nDMs_in_range, itemp;
	GpuTimer timer;
	float Total_periodicity_time = 0, Time_per_range = 0;
	
	//---------> Finding nearest lower power of two (because of FFT algorithm)
	int nearest = (int) floorf(log2f((float) processed));
	nTimesamples = (int) powf(2.0, nearest);
	printf("Decreasing number of processed samples to nearest lower power of two, because of FFT algorithm...\n");
	printf("Number of processed timesamples: %d; nearest power of two: %d\n", processed, nTimesamples);
	processed = nTimesamples;
	
	//--------> Determining maximum number of DM trials we can fit into memory
	size_t free_mem,total_mem;
	cudaMemGetInfo(&free_mem,&total_mem);
	printf("     Memory available:%0.3f MB \n", ((float) free_mem)/(1024.0*1024.0) );
	size_t max_nDMs_in_memory = (free_mem*0.95)/(processed*(5.5*sizeof(float) + 2*sizeof(ushort))); // 1 for real input real, 2 for complex output, 2 for complex cuFFT, 1 for peaks + 1 ushort
	itemp = (int) (max_nDMs/PER_SEARCH_NTHREADS);
	max_nDMs = itemp*PER_SEARCH_NTHREADS;
	printf("     Maximum number of DM trials which fit into memory is %d; Input plane size: %0.2f MB;\n", max_nDMs_in_memory, (((float) max_nDMs_in_memory*processed*sizeof(float))/(1024.0*1024.0));
	
	
	//--------> Allocation of GPU memory. We allocate such amount of memory as to accommodate maximum number of DM trials from first range.
	unsigned int input_plane_size = (processed+2)*max_nDMs_in_memory; //+2 because of interbinning
	float *d_one_A; //for input and interbinned values
	if ( cudaSuccess != cudaMalloc((void **) &d_one_A,  sizeof(float)*input_plane_size )) printf("Periodicity Allocation error! d_one_A\n");
	
	float *d_two_B; //for cuFFT complex output and peaks
	if ( cudaSuccess != cudaMalloc((void **) &d_two_B,  sizeof(float)*2*input_plane_size )) printf("Periodicity Allocation error! d_two_B\n");
	
	float *d_half_C; // for power values
	if ( cudaSuccess != cudaMalloc((void **) &d_local_peaks,  sizeof(float)*input_plane_size/2 )) printf("Periodicity Allocation error! d_spectra_Real\n");
	
	ushort *d_harmonics;
	if ( cudaSuccess != cudaMalloc((void **) &d_harmonics, sizeof(ushort)*2*input_plane_size )) printf("Periodicity Allocation error! d_harmonics\n");
	
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
		Time_per_range = 0;
		host_power_peak_pos = 0; host_interbin_peak_pos = 0;
		
		nTimesamples = processed/inBin[i];
		nDMs = ndms[i];
		max_nDMs_in_range = max_nDMs_in_memory*inBin[i];
		printf("Processing de-dispersion range:%f--%f:%f; inBin:%d; Timesamples:%d; DM trials:%d; max_nDMs:%d;\n", dm_low[i], dm_high[i], dm_step[i], inBin[i], nTimesamples, nDMs, max_nDMs_in_range);
		
		unsigned int plane_size = nTimesamples*nDMs;
		int nRepeats, nRest, DM_shift, itemp, DMs_per_cycle;
		std::vector<int> DM_list;
		
		//---------> Setting up batches
		nRepeats = nDMs/max_nDMs_in_range;
		nRest = nDMs - nRepeats*max_nDMs_in_range;
		for(int f=0; f<nRepeats; f++) DM_list.push_back(max_nDMs);
		if(nRest>0) DM_list.push_back(nRest);
		
		printf("     Periodicity search will run %d batches each containing %d DM trials. Remainder %d DM trials\n", (int) DM_list.size(), max_nDMs, nRest);
		
		if(DM_list.size()>0){
			DMs_per_cycle = DM_list[0];
			DM_shift = 0;
			for(int f=0; f<DM_list.size(); f++) {
				cudaMemset((void*) gmem_power_peak_pos, 0, sizeof(int));
				cudaMemset((void*) gmem_interbin_peak_pos, 0, sizeof(int));
				
				//---------> Copy data from the host
				checkCudaErrors( cudaMemcpy( d_one_A, output_buffer[i][DM_shift], nTimesamples*DMs_per_cycle*sizeof(float), cudaMemcpyHostToDevice));
				//---------<

				//---------> cuFFT
				timer.Start();
				cufftHandle plan_input;
				cufftResult error;
				if ( cufftPlan1d(&plan_input, processed, CUFFT_R2C, DMs_per_cycle) != CUFFT_SUCCESS) printf("CUFFT error: %d", error);
				cufftExecR2C(plan_input, (cufftReal *)d_one_A, (cufftComplex *)d_two_B, CUFFT_FORWARD);
				cufftDestroy(plan_input);
				timer.Stop();
				printf("     -> cuFFT took %f ms\n", timer.Elapsed());
				Time_per_range = Time_per_range + timer.Elapsed();
				//---------<
				
				//---------> Calculate powers and interbinning
				timer.Start();
				Calculate_power_and_interbin();
				timer.Stop();
				printf("     -> Calculation of powers and interbining took %f ms\n", timer.Elapsed());
				Time_per_range = Time_per_range + timer.Elapsed();
				//---------<
				
				//---------> Mean and StDev on powers
				timer.Start();
				MSD_limited(d_half_C, d_MSD, nDMs, nTimesamples/2, 0);
				//MSD_BLN_pw(d_half_C, d_MSD, nDMs, nTimesamples/2, 0, sigma_constant);
				timer.Stop();
				printf("     -> baselevel MSD took %f ms\n", timer.Elapsed());
				Time_per_range = Time_per_range + timer.Elapsed();
				//---------<
				
				//---------> Harmonic summing
				timer.Start();
				
				timer.Stop();
				printf("     -> harmonic summing took %f ms\n", timer.Elapsed());
				Time_per_range = Time_per_range + timer.Elapsed();
				//---------<
				
				//---------> Peak finding
				timer.Start();
				Peak_find_for_periodicity_search(d_half_C, &d_harmonics[0], &d_two_B[0], nDMs, nTimesamples, cutoff, local_max_list_size, gmem_power_peak_pos, DM_shift);
				Peak_find_for_periodicity_search(d_one_A, &d_harmonics[input_plane_size], &d_two_B[input_plane_size], nDMs, nTimesamples, cutoff, local_max_list_size, gmem_interbin_peak_pos, DM_shift);
				timer.Stop();
				printf("     -> Peak finding took %f ms\n", timer.Elapsed());
				Time_per_range = Time_per_range + timer.Elapsed();
				//---------<
				
				checkCudaErrors(cudaGetLastError());
				
				//---------> Transferring peaks to the host
				checkCudaErrors(cudaMemcpy(&temp_host_power_peak_pos, gmem_power_peak_pos, sizeof(int), cudaMemcpyDeviceToHost));
				#ifdef GPU_PERIODICITY_SEARCH_DEBUG
				printf("POWER: Total number of peaks found in this range is %d; maximum number of peaks:%d;\n", temp_host_power_peak_pos, local_max_list_size);
				#endif
				else if( (host_power_peak_pos + temp_power_peak_pos)<max_power_peak_size){
					checkCudaErrors(cudaMemcpy(&h_all_power_peaks[host_power_peak_pos*4], &d_two_B[0], temp_host_power_peak_pos*4*sizeof(float), cudaMemcpyDeviceToHost));
					host_power_peak_pos = host_power_peak_pos + temp_host_power_peak_pos;
				}
				else printf("     Maximum list size reached! Increase list size or increase sigma cutoff.\n");
				
				checkCudaErrors(cudaMemcpy(&temp_host_interbin_peak_pos, gmem_interbin_peak_pos, sizeof(int), cudaMemcpyDeviceToHost));
				#ifdef GPU_PERIODICITY_SEARCH_DEBUG
				printf("POWER: Total number of peaks found in this range is %d; maximum number of peaks:%d;\n", temp_host_interbin_peak_pos, local_max_list_size);
				#endif
				else if( (host_interbin_peak_pos + temp_interbin_peak_pos)<max_interbin_peak_size){
					checkCudaErrors(cudaMemcpy(&h_all_interbin_peaks[host_interbin_peak_pos*4], &d_two_B[input_plane_size], temp_host_interbin_peak_pos*4*sizeof(float), cudaMemcpyDeviceToHost));
					host_interbin_peak_pos = host_interbin_peak_pos + temp_host_interbin_peak_pos;
				}
				else printf("     Maximum list size reached! Increase list size or increase sigma cutoff.\n");
				//---------<
				
				DM_shift = DM_shift + DM_list[f];
			}
			
			//---------> Peak processing on host and export
			#pragma omp parallel for
			for (int count = 0; count < host_power_peak_pos; count++){
				h_all_power_peaks[4*count]     = h_all_power_peaks[4*count]*dm_step[i] + dm_low[i];
				h_all_power_peaks[4*count + 1] = h_all_power_peaks[4*count + 1]*(1.0/(tsamp*nTimesamples));
			}
			
			#pragma omp parallel for
			for (int count = 0; count < host_interbin_peak_pos; count++){
				h_all_interbin_peaks[4*count]     = h_all_interbin_peaks[4*count]*dm_step[i] + dm_low[i];
				h_all_interbin_peaks[4*count + 1] = h_all_interbin_peaks[4*count + 1]*(1.0/(tsamp*nTimesamples));
			}
			
			FILE *fp_out;
			char filename[200];
			
			if(host_power_peak_pos>0){
				sprintf(filename, "fourier-dm_%.2f-%.2f.dat", dm_low[i], dm_high[i]);
				if (( fp_out = fopen(filename, "wb") ) == NULL)	{
					fprintf(stderr, "Error opening output file!\n");
					exit(0);
				}
				fwrite(h_all_power_peaks, host_power_peak_pos)*sizeof(float), 4, fp_out);
				fclose(fp_out);
			}
			
			if(host_interbin_peak_pos>0){
				sprintf(filename, "fourier_inter-dm_%.2f-%.2f.dat", dm_low[i], dm_high[i]);
				if (( fp_out = fopen(filename, "wb") ) == NULL)	{
					fprintf(stderr, "Error opening output file!\n");
					exit(0);
				}
				fwrite(h_all_interbin_peaks, host_interbin_peak_pos)*sizeof(float), 4, fp_out);
				fclose(fp_out);
			}
			//---------<
		}
				
	}
	
	*/
	
	/*
	// Example FFT....

	printf("\n");

	printf("[1DCUFFT] is starting...\n");

	//FILE	*fp_c, *fp_dm, *fp_harm;
	FILE *fp_c, *fp_dm;
	char filename[200];

	int number_of_candidates = 10;

	float** h_top_list = (float**) malloc(sizeof(float*) * 5);
	for (int a = 0; a < 5; a++)
	{
		h_top_list[a] = (float*) malloc(sizeof(float) * number_of_candidates);
	}
	for (int a = 0; a < 5; a++)
	{
		for (int b = 0; b < number_of_candidates; b++)
		{
			h_top_list[a][b] = 0.0f;
		}
	}

	for (int i = 0; i < range; i++)
	{
		int samps = processed / inBin[i];

		// Allocate memory for signal
		cufftReal* d_signal_in;
		cudaMalloc((void**) &d_signal_in, samps * sizeof(cufftReal));

		cufftComplex* d_signal_out;
		cudaMalloc((void**) &d_signal_out, ( samps / 2 + 1 ) * sizeof(cufftComplex));

		cufftComplex* h_signal = (cufftComplex*) malloc(( samps / 2 + 1 ) * sizeof(cufftComplex));
		float* h_signal_x = (float*) malloc(sizeof(float) * ( samps / 2 + 1 ) * ndms[i]);
		float* h_signal_y = (float*) malloc(sizeof(float) * ( samps / 2 + 1 ) * ndms[i]);
		float* h_signal_p = (float*) malloc(sizeof(float) * ( samps / 2 + 1 ) * ndms[i]);
		float* h_harm = (float*) malloc(sizeof(float) * ( samps / 2 + 1 ) * ndms[i]);
		float* h_signal_inter = (float*) malloc(sizeof(float) * 2 * ( samps / 2 + 1 ) * ndms[i]);

		float** h_candidates = (float**) malloc(sizeof(float*) * ndms[i]);
		for (int a = 0; a < ndms[i]; a++) {
			h_candidates[a] = (float*) malloc(sizeof(float) * ( samps / 2 + 1 ));
		}
		for (int a = 0; a < ndms[i]; a++) {
			for (int b = 0; b < samps / 2 + 1; b++) {
				h_candidates[a][b] = 0.0f;
			}
		}

		sprintf(filename, "fourier-%d.dat", i);
		if (( fp_c = fopen(filename, "w") ) == NULL)
		{
			fprintf(stderr, "Error opening output file!\n");
			exit(0);
		}
		sprintf(filename, "fourier_inter-%d.dat", i);
		if (( fp_dm = fopen(filename, "w") ) == NULL)
		{
			fprintf(stderr, "Error opening output file!\n");
			exit(0);
		}

		
		
		// CUFFT plan
		cufftHandle plan;
		cufftPlan1d(&plan, samps, CUFFT_R2C, 1);

		for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {

			cudaMemcpy(d_signal_in, output_buffer[i][dm_count], samps * sizeof(float), cudaMemcpyHostToDevice);

			// Transform signal 
			//printf("\nTransforming dm: %f using cufftExecR2C\n", dm);
			cufftExecR2C(plan, (cufftReal *) d_signal_in, (cufftComplex *) d_signal_out);

			// Copy device memory to host
			cudaMemcpy(h_signal, d_signal_out, sizeof(cufftComplex) * ( samps / 2 + 1 ), cudaMemcpyDeviceToHost);

			h_signal_p[0 + dm_count * ( samps / 2 )] = 0.0;
			#pragma omp parallel for
			for (int j = 1; j < samps / 2; j++) {
				//	h_signal[j].x = h_signal[j].x-h_signal[0].x;
				//	h_signal[j].y = h_signal[j].y-h_signal[0].y;
				h_signal_x[j + dm_count * ( samps / 2 )] = h_signal[j].x;
				h_signal_y[j + dm_count * ( samps / 2 )] = h_signal[j].y;
				h_signal_p[j + dm_count * ( samps / 2 )] = ( ( h_signal[j].x * h_signal[j].x + h_signal[j].y * h_signal[j].y ) );
				h_signal_inter[2 * j + dm_count * samps] = h_signal_p[j + dm_count * ( samps / 2 )];
				h_signal_inter[2 * j + 1 + dm_count * samps] = 0.616850275 * ( ( h_signal[j].x - h_signal[j + 1].x )*( h_signal[j].x - h_signal[j + 1].x ) + ( h_signal[j].y - h_signal[j + 1].y ) * ( h_signal[j].y - h_signal[j + 1].y ) );
			}
		}

		//Destroy CUFFT context
		cufftDestroy(plan);

		// cleanup memory
		free(h_signal);
		cudaFree(d_signal_in);
		cudaFree(d_signal_out);
		
		
		
		
		//--------------------------------------------------------------------------
		//--------> CPU MSD
		double mean, stddev;
		double total = 0.0;
		
		for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {
			for (int j = 0; j < ( samps / 2 ); j++) {
				total += ( (double) ( h_signal_p[j + dm_count * ( samps / 2 )] ) );
			}
		}
		mean = ( total / (double) ( ( samps / 2 ) * ndms[i] ) ); 
		// Mean for data sample

		// Calculate standard deviation
		total = 0.0;
		for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {
			for (int j = 0; j < ( samps / 2 ); j++) {
				total += (double) ( ( h_signal_p[j + dm_count * ( samps / 2 )] - (float) mean ) * ( h_signal_p[j + dm_count * ( samps / 2 )] - (float) mean ) );
			}
		}
		stddev = sqrt(abs(total) / (double) ( ( samps / 2 ) * ndms[i] )); // Stddev for data sample
		//--------> CPU MSD
		//--------------------------------------------------------------------------
		
		//--------------------------------------------------------------------------
		//--------> Calculation of SNR and export
		for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {
			float dm = dm_low[i] + dm_step[i] * dm_count;
			for (int j = 0; j < ( samps / 2 ); j++) {
				if ((float) ( ( h_signal_p[j + dm_count * ( samps / 2 )] - mean ) / stddev ) > cutoff) {
					fprintf(fp_c, "\n%f\t%f\t%f", dm, j*((1.0f/tsamp)/(samps)), (float) ( ( (double) h_signal_p[j + dm_count * ( samps / 2 )] - mean ) / stddev ));

				}
			}
			fprintf(fp_c, "\n");
			for (int j = 0; j < samps; j++) {
				if ((float) ( ( h_signal_inter[j + dm_count * samps] - mean ) / stddev ) > cutoff) {
					fprintf(fp_dm, "\n%f\t%lf\t%f", dm, j * ( ( 1.0 / tsamp ) / ( 2 * samps ) ), (float) ( ( (double) h_signal_inter[j + dm_count * samps] - mean ) / stddev ));

				}
			}
			fprintf(fp_dm, "\n");
		}
		//--------> Calculation of SNR and export
		//--------------------------------------------------------------------------
		*/
		
		/*
		 int harm_max=32;
		 for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {
		 for(int j=0;j< (samps/2);j++){
		 h_harm[j+dm_count*(samps/2)] = (h_signal_p[j+dm_count*(samps/2)]);
		 }
		 }
		 int harm=1;
		 sprintf(filename, "harmonic-%d-%d.dat", i, harm);
		 if ((fp_harm=fopen(filename, "w")) == NULL) {
		 fprintf(stderr, "Error opening output file!\n");
		 exit(0);
		 }

		 // Calculate the mean
		 total=0.0;
		 for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {
		 for(int j=0;j< (samps/2);j++){
		 total += ((double)(h_harm[j+dm_count*(samps/2)]));
		 }
		 }
		 mean = (total/(double)((samps/2)*ndms[i]));  // Mean for data sample

		 // Calculate standard deviation
		 total = 0.0;
		 for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {
		 for(int j=0;j< (samps/2); j++){
		 total += (double)((h_harm[j+dm_count*(samps/2)]-(float)mean)*(h_harm[j+dm_count*(samps/2)]-(float)mean));
		 }
		 }
		 stddev = sqrt(abs(total) / (double)((samps/2)*ndms[i])); // Stddev for data sample

		 for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {
		 float dm = dm_low[i]+dm_step[i]*dm_count;
		 for(int j=0;j< samps/2; j++){
		 float candidate = (float)(((double)h_harm[j+dm_count*(samps/2)]-mean)/stddev);
		 if(candidate > cutoff) {
		 fprintf(fp_harm, "\n%f\t%f\t%f", dm, j*1*((1.0f/tsamp)/(samps)), candidate);

		 for(int c = 0; c < number_of_candidates; c++) {
		 if(candidate > h_top_list[4][c]) {
		 for(int d = number_of_candidates - 1; d > c; d--) {
		 h_top_list[0][d] = h_top_list[0][d-1];
		 h_top_list[1][d] = h_top_list[1][d-1];
		 h_top_list[2][d] = h_top_list[2][d-1];
		 h_top_list[3][d] = h_top_list[3][d-1];
		 h_top_list[4][d] = h_top_list[4][d-1];
		 }
		 h_top_list[0][c] = dm;
		 h_top_list[1][c] = j*1*((1.0f/tsamp)/(samps));
		 h_top_list[2][c] = harm;
		 h_top_list[3][c] = j;
		 h_top_list[4][c] = candidate;
		 c=number_of_candidates;
		 }
		 }
		 }
		 h_candidates[dm_count][j] = candidate;
		 }
		 fprintf(fp_harm, "\n");
		 }
		 fclose(fp_harm);

		 for(harm = 2; harm <= harm_max; harm=2*harm) {

		 sprintf(filename, "harmonic-%d-%d.dat", i, harm);
		 if ((fp_harm=fopen(filename, "w")) == NULL) {
		 fprintf(stderr, "Error opening output file!\n");
		 exit(0);
		 }

		 for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {
		 for(int j=0;j < (samps/(2*harm))-harm; j++){
		 h_harm[j*harm+dm_count*(samps/2)] += (h_signal_p[j+dm_count*(samps/2)]);
		 for(int lerp = j+1; lerp < j+harm; lerp++) h_harm[lerp+dm_count*(samps/2)] += (h_signal_p[j+dm_count*(samps/2)] +
		 (h_signal_p[j+1+dm_count*(samps/2)]-h_signal_p[j+dm_count*(samps/2)])*((lerp-j)/harm));
		 }
		 }

		 // Calculate the mean
		 total=0.0;
		 for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {
		 for(int j=0;j< (samps/2);j++){
		 total += ((double)(h_harm[j+dm_count*(samps/2)]));
		 }
		 }
		 mean = (total/(double)((samps/2)*ndms[i]));  // Mean for data sample

		 // Calculate standard deviation
		 total = 0.0;
		 for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {
		 for(int j=0;j< (samps/2); j++){
		 total += (double)((h_harm[j+dm_count*(samps/2)]-(float)mean)*(h_harm[j+dm_count*(samps/2)]-(float)mean));
		 }
		 }
		 stddev = sqrt(abs(total) / (double)((samps/2)*ndms[i])); // Stddev for data sample

		 for (int dm_count = 0; dm_count < ndms[i]; dm_count++) {
		 float dm = dm_low[i]+dm_step[i]*dm_count;
		 for(int j=0;j< samps/2; j++){
		 float candidate = (float)(((double)h_harm[j+dm_count*(samps/2)]-mean)/stddev);
		 if(candidate > sqrt(harm)*cutoff) {
		 fprintf(fp_harm, "\n%f\t%f\t%f", dm, j*harm*((1.0f/tsamp)/(samps)), candidate);
		 for(int c = 0; c < number_of_candidates; c++) {
		 if(candidate > h_top_list[4][c]) {
		 for(int d = number_of_candidates - 1; d > c; d--) {
		 h_top_list[0][d] = h_top_list[0][d-1];
		 h_top_list[1][d] = h_top_list[1][d-1];
		 h_top_list[2][d] = h_top_list[2][d-1];
		 h_top_list[3][d] = h_top_list[3][d-1];
		 h_top_list[4][d] = h_top_list[4][d-1];
		 }
		 h_top_list[0][c] = dm;
		 h_top_list[1][c] = j*harm*((1.0f/tsamp)/(samps));
		 h_top_list[2][c] = harm;
		 h_top_list[3][c] = harm;
		 h_top_list[4][c] = candidate;
		 c=number_of_candidates;
		 }
		 }
		 }
		 h_candidates[dm_count][j] = (float)(((double)h_harm[j+dm_count*(samps/2)]-mean)/stddev);
		 }
		 fprintf(fp_harm, "\n");
		 }
		 fclose(fp_harm);
		 }
		 */
}

//	for (int c = 0 ; c < ( number_of_candidates - 1 ); c++) {
//		printf("\nCandidate: %d, DM: %f, PERIOD: %f, HARMONIC: %f, PxH: %f, SNR: %f", c, h_top_list[0][c], 1.0f/h_top_list[1][c], h_top_list[2][c], (1.0f/h_top_list[1][c] * h_top_list[2][c]), h_top_list[4][c]);
//	}


} //namespace astroaccelerate
