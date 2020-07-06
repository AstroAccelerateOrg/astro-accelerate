#include <stdio.h>
#include <stdlib.h>
#include <cufft.h>
#include <cuda_profiler_api.h>
#include <iostream>
#include <vector>

#include "aa_params.hpp"
#include "aa_log.hpp"
#include "aa_gpu_timer.hpp"
#include "aa_timelog.hpp"
#include "presto_funcs.hpp"

#include "aa_jerk_plan.hpp"
#include "aa_jerk_strategy.hpp"


#include "aa_jerk_CandidateList.h"

// MSD
#include "aa_device_MSD_Configuration.hpp"
#include "aa_device_MSD.hpp"
#include "aa_device_MSD_plane_profile.hpp"

// Convolution
#include "aa_device_convolution.hpp"

class GPU_Memory_for_JERK_Search {
	private:
		size_t MSD_interpolated_size_in_bytes;
		size_t MSD_DIT_size_in_bytes;
		
		size_t allocated_DM_trial;
		size_t allocated_DM_trial_ffted;
		size_t allocated_MSD_size;
		size_t allocated_ZW_planes;
		size_t allocated_ZW_candidates;
		size_t allocated_MSD_workarea;
		
	public:
		float *d_DM_trial;
		float2 *d_DM_trial_ffted;
		float *d_interpolated_MSD;
		unsigned int *gmem_peak_pos;
		float *d_ZW_candidates;
		float *d_ZW_planes;
		float *d_MSD_workarea;
};
		

aa_jerk_plan create_plan_from_strategy(aa_jerk_strategy &strategy, size_t nTimesamples, size_t nDMs) {
	bool do_interbinning = (strategy.interbinned_samples()==2?true:false0);
	bool do_high_precision = (strategy.high_precision()==1?true:false);
	aa_jerk_plan jerk_plan(nTimesamples, mDMs, strategy.z_max_search_limit(), strategy.z_search_step(), strategy.w_max_search_limit(), strategy.w_search_step(), do_interbinning, do_high_precision);
	if(strategy.MSD_outlier_rejection()) jerk_plan.enable_MSD_outlier_rejection();
	else jerk_plan.disable_MSD_outlier_rejection();
	jerk_plan.set_outlier_rejection_sigma_cutoff(strategy.OR_sigma_cuttoff());
	jerk_plan.set_candidate_selection_sigma_threshold(strategy.CS_sigma_threshold());
}


int jerk_search_from_ddtr_plan(float ***dedispersed_data, aa_jerk_strategy jerk_strategy, float *dm_low, float *dm_step, int *list_of_ndms, float sampling_time, int *inBin, int nRanges){

	
	//----------> Convolution test
	int JERK_SEARCH_CONVOLUTION_TEST = 0;
	if(JERK_SEARCH_CONVOLUTION_TEST) {
		printf("---------------- JERK Search convolution test --------------------\n");
		std::ofstream FILEOUT;
		char filename[200];
		int convolution_length = 2048;
		
		//--------> Filters
		int nFilters = 56;
		int filter_size_bytes = nFilters*convolution_length*sizeof(float2);
		int filter_halfwidth = nFilters+2;
		
		float2 *h_jerk_filters;	
		float2 *d_jerk_filters;
		
		h_jerk_filters = new float2 [nFilters*convolution_length];
		for (int f = 0; f < nFilters; f++){
			int boxcar_width=(f+1)*2;
			for(int s=0; s<convolution_length; s++){
				h_jerk_filters[f*convolution_length + s].x = 0;
				h_jerk_filters[f*convolution_length + s].y = 0;
	  
				if(s<boxcar_width/2) h_jerk_filters[f*convolution_length + s].x = 1.0;
				if(s>=(convolution_length-boxcar_width/2)) h_jerk_filters[f*convolution_length + s].x = 1.0;
			}
		}
		
		printf("---> Exporting filters:"); fflush(stdout);
		FILEOUT.open("TEST_conv_filters.dat");
		for(int f=0; f<nFilters; f++){
			for(int s=0; s<convolution_length; s++){
				int pos = f*convolution_length + s;
				FILEOUT << sqrt(h_jerk_filters[pos].x*h_jerk_filters[pos].x + h_jerk_filters[pos].y*h_jerk_filters[pos].y) << " " << h_jerk_filters[pos].x << " " << h_jerk_filters[pos].y << std::endl;
			}
			FILEOUT << std::endl;
			FILEOUT << std::endl;
			printf(".");
			fflush(stdout);
		}
		FILEOUT.close();
		printf("\n");
		
		cudaMalloc((void **) &d_jerk_filters,  filter_size_bytes);
		cudaMemcpy(d_jerk_filters, h_jerk_filters, filter_size_bytes, cudaMemcpyHostToDevice);
		forwardCustomFFT(d_jerk_filters, convolution_length, nFilters);
		
		cudaMemcpy(h_jerk_filters, d_jerk_filters, filter_size_bytes, cudaMemcpyDeviceToHost);
		printf("---> Exporting filters:"); fflush(stdout);
		FILEOUT.open("TEST_conv_filters_ffted.dat");
		for(int f=0; f<nFilters; f++){
			for(int s=0; s<convolution_length; s++){
				int pos = f*convolution_length + s;
				FILEOUT << sqrt(h_jerk_filters[pos].x*h_jerk_filters[pos].x + h_jerk_filters[pos].y*h_jerk_filters[pos].y) << " " << h_jerk_filters[pos].x << " " << h_jerk_filters[pos].y << std::endl;
			}
			FILEOUT << std::endl;
			FILEOUT << std::endl;
			printf(".");
			fflush(stdout);
		}
		FILEOUT.close();
		printf("\n");		
		//-----------------------------------<
		
		//--------> Signal
		int signal_length = 524288;
		int signal_size_bytes = signal_length*sizeof(float2);
		float2 *h_signal = NULL;
		float2 *d_signal = NULL;
		
		h_signal = new float2 [signal_length];
		for(int f=15000; f<signal_length; f++){
			h_signal[f].x = (f%4096)/500.0;
			h_signal[f].y = 0;
		}
		
		for(int f=0; f<192; f++){
			h_signal[f + 5300].x = 10.0;
			h_signal[f + 5300].y = 0;
		}
		for(int f=0; f<128; f++){
			h_signal[f + 8626].x = 10.0;
			h_signal[f + 8626].y = 0;
		}
		for(int f=0; f<36; f++){
			h_signal[f + 9626].x = 10.0;
			h_signal[f + 9626].y = 0;
		}
		for(int f=0; f<83; f++){
			h_signal[f + 10626].x = 10.0;
			h_signal[f + 10626].y = 0;
		}
		for(int f=0; f<138; f++){
			h_signal[f + 11626].x = 10.0;
			h_signal[f + 11626].y = 0;
		}
		
		printf("---> Exporting signal.\n");
		FILEOUT.open("TEST_conv_signal.dat");
		for(int s=0; s<signal_length; s++){
			FILEOUT << sqrt(h_signal[s].x*h_signal[s].x + h_signal[s].y*h_signal[s].y) << " " << h_signal[s].x << " " << h_signal[s].y << std::endl;
		}
		FILEOUT.close();
		
		cudaMalloc((void **) &d_signal,  signal_size_bytes);
		cudaMemcpy(d_signal, h_signal, signal_size_bytes, cudaMemcpyHostToDevice);
		//-----------------------------------<
		
		int useful_part_size = convolution_length - 2*filter_halfwidth + 1;
		int nSegments        = (signal_length + useful_part_size - 1)/useful_part_size;
		
		float *d_output = NULL;
		float *h_output = NULL;
		int output_size = nFilters*useful_part_size*nSegments;
		h_output = new float [output_size];
		cudaMalloc((void **) &d_output,  output_size*sizeof(float));
		
		conv_OLS_customFFT(d_signal, d_output, d_jerk_filters, signal_length, convolution_length, useful_part_size, filter_halfwidth, nSegments, nFilters, 1.0/(2048.0f*2048.0f));
		
		cudaMemcpy(h_output, d_output, output_size*sizeof(float), cudaMemcpyDeviceToHost);
		
		printf("---> Exporting output:"); fflush(stdout);
		for(int f=0; f<nFilters; f++){
			sprintf(filename, "TEST_conv_output_filter%d.dat", f);
			FILEOUT.open(filename);
			for(int s=0; s<useful_part_size*nSegments; s++){
				int pos = f*useful_part_size*nSegments + s;
				FILEOUT << h_output[pos] << std::endl;
			}
			FILEOUT.close();
			printf(".");
			fflush(stdout);
		}
		printf("\n");
		
		
		//--------> Clean-up
		delete [] h_signal;
		cudaFree(d_signal);
		
		delete [] h_jerk_filters;
		cudaFree(d_jerk_filters);

		delete [] h_output;
		cudaFree(d_output);
		printf("---> Done!\n");
		printf("---------------- JERK Search convolution test --------------------\n\n");
	}
	//-----------------------------------------------------------<
	
	
	//---------> Time measurements
	GpuTimer timer_total, timer_DM, timer;
	double time_total=0, time_DM=0, time=0;
	timer_total.Start();
	
	jerk_strategy.PrintStrategy();
	
	//---------> Generating filters
	float2 *h_jerk_filters;	
	float2 *d_jerk_filters;
	h_jerk_filters = new float2[jerk_strategy.filter_padded_size()];
	if ( cudaSuccess != cudaMalloc((void **) &d_jerk_filters,  sizeof(float2)*jerk_strategy.filter_padded_size() )) {
		printf("Cannot allocate GPU memory for JERK filters!\n");
		return(1);
	}
	
	jerk_create_acc_filters(h_jerk_filters, &jerk_strategy);
	if ( cudaSuccess != cudaMemcpy(d_jerk_filters, h_jerk_filters, jerk_strategy.filter_padded_size_bytes(), cudaMemcpyHostToDevice) ) {
		printf("Error occured during host -> device transfer!\n");
		return(2);
	}
	
	forwardCustomFFT(d_jerk_filters, jerk_strategy.conv_size(), jerk_strategy.nFilters_total());
	//-----------------------------------------------------------<
	
	
	//---------> Device data
	float *d_DM_trial = NULL;
	if ( cudaSuccess != cudaMalloc((void **) &d_DM_trial,  sizeof(float)*jerk_strategy.nSamples_time_dom() )) {
		printf("Cannot allocate GPU memory for DM trial!\n");
		return(2);
	}
	
	float2 *d_DM_trial_ffted = NULL;
	if ( cudaSuccess != cudaMalloc((void **) &d_DM_trial_ffted,  sizeof(float2)*jerk_strategy.nSamples_freq_dom() )) {
		printf("Cannot allocate GPU memory for FFTed DM trial!\n");
		return(3);
	}
	
	float *d_MSD;
	if ( cudaSuccess != cudaMalloc((void **) &d_MSD,  MSD_PARTIAL_SIZE*jerk_strategy.nHarmonics()*sizeof(float))) {
		printf("Cannot allocate GPU memory for MSD!\n");
		return(4);
	}
	
	unsigned int *gmem_peak_pos;
	if ( cudaSuccess != cudaMalloc((void**) &gmem_peak_pos, sizeof(unsigned int))) {
		printf("Cannot allocate GPU memory for peak position!\n");
	}
	cudaMemset((void*) gmem_peak_pos, 0, sizeof(unsigned int));
	
	float *d_ZW_candidates = NULL;
	float *d_ZW_planes     = NULL;
	float *d_MSD_workarea  = NULL;
	//-----------------------------------------------------------<
	
	size_t default_nTimesamples = jerk_strategy.nTimesamples();
	double MSD_time = 0, Candidate_time = 0, Convolution_time = 0;
	for(int active_range=0; active_range<nRanges; active_range++){
		size_t DM_trial_samples = jerk_strategy.nTimesamples/inBin[active_range];
		size_t nDMs = list_of_ndms[active_range];
		
		//---------> Jerk plan
		JERK_Plan local_plan;
		local_plan = user_plan;
		local_plan.nTimesamples = DM_trial_samples;
		local_plan.nDMs = nDMs;
		
		//---------> Jerk strategy
		cudaMemGetInfo(&free_memory, &total_memory);
		if(free_memory>8589934592) free_memory = 8589934592;
		JERK_Strategy jerk_strategy(local_plan, free_memory);
		jerk_strategy.PrintStrategy();
		
		//---------> Allocation of the output Z-planes and candidates
		size_t ZW_plane_size              = jerk_strategy.output_size_z_plane;
		unsigned int max_nZWCandidates    = (ZW_plane_size/4);
		size_t single_ZW_plane_size_bytes = jerk_strategy.output_size_z_plane*sizeof(float);
		size_t ZW_planes_size_bytes       = jerk_strategy.nZPlanes_per_chunk*single_ZW_plane_size_bytes;
		
		printf("JERK SEARCH -> ZW single plane size: %zu elements = %f MB\n", single_ZW_plane_size_bytes, ((float) single_ZW_plane_size_bytes)/(1024.0*1024.0));
		printf("JERK SEARCH -> ZW plane size: %zu elements = %f MB\n", ZW_planes_size_bytes, ((float) ZW_planes_size_bytes)/(1024.0*1024.0));
		
		if ( cudaSuccess != cudaMalloc((void **) &d_ZW_candidates, single_ZW_plane_size_bytes)) {
			printf("Cannot allocate GPU memory for ZW plane candidates!\n");
		}
		if ( cudaSuccess != cudaMalloc((void **) &d_ZW_planes, ZW_planes_size_bytes)) {
			printf("Cannot allocate GPU memory for ZW planes!\n");
		}
		
		//---------> cuFFT plan
		cufftHandle cuFFT_plan;
		cufftResult cuFFT_error;
		cuFFT_error = cufftPlan1d(&cuFFT_plan, jerk_strategy.nSamples_time_dom, CUFFT_R2C, 1);
		if (cuFFT_error!=CUFFT_SUCCESS) {
			printf("ERROR while constructing cuFFT plan\n");
		}
		
		//---------> MSD configuration
		MSD_Configuration MSD_conf(jerk_strategy.output_size_one_DM, jerk_strategy.nFilters_z, 0, 0);
		MSD_conf.print();
		if ( cudaSuccess != cudaMalloc((void **) &d_MSD_workarea, MSD_conf.nBlocks_total*MSD_PARTIAL_SIZE*sizeof(float))) {
			printf("Cannot allocate GPU memory for MSD workarea!\n");
		}
		

		for(int active_DM = 0; active_DM<nDMs; active_DM++){
			//-------> Copy DM-trial
			if ( cudaSuccess != cudaMemcpy(d_DM_trial, dedispersed_data[active_range][active_DM], jerk_strategy.nSamples_time_dom*sizeof(float), cudaMemcpyHostToDevice) ) {
				printf("ERROR while copying DM-trial to the device!\n");
			}
			
			//-------> FFT of the DM-trial
			cufftExecR2C(cuFFT_plan, d_DM_trial, (cufftComplex *) d_DM_trial_ffted);
			
			checkCudaErrors(cudaGetLastError());
			
			//---------> Candidates
			int ZW_planes_shift = 0; //this is for filters
			//TODO: container for candidates!
			
			JERK_CandidateList allcandidates(jerk_strategy.nFilters_w);
			printf("nSublists: %d\n", allcandidates.getNumberOfSubLists());
			
			for(int f=0; f<(int) jerk_strategy.ZW_chunks.size(); f++){
				int nZW_planes = jerk_strategy.ZW_chunks[f];
				
				// Convolution
				timer.Start();
				printf("JERK DEBUG: Convolution: position of the filters = %d; nTimesamples = %d; conv_size = %d; useful_part_size = %d; offset = %d; nSegments = %d; nFilters = %d;\n", jerk_strategy.nFilters_z*ZW_planes_shift, jerk_strategy.output_size_one_DM, jerk_strategy.conv_size, jerk_strategy.useful_part_size, jerk_strategy.filter_halfwidth, jerk_strategy.nSegments, jerk_strategy.nFilters_z*nZW_planes);
				
				conv_OLS_customFFT(d_DM_trial_ffted, d_ZW_planes, &d_jerk_filters[jerk_strategy.nFilters_z*ZW_planes_shift*jerk_strategy.conv_size], jerk_strategy.nSamples_freq_dom, jerk_strategy.conv_size, jerk_strategy.useful_part_size, jerk_strategy.filter_halfwidth, jerk_strategy.nSegments, jerk_strategy.nFilters_z*nZW_planes, 1.0f);
				timer.Stop(); printf("JERK SEARCH -> Convolution took: %g ms\n", timer.Elapsed());
				Convolution_time += timer.Elapsed();
				
				
				for(int zp=0; zp<nZW_planes; zp++){
					unsigned int nCandidates = 0;
					char filename[200];
					std::ofstream FILEOUT;
					float w = ((float) (ZW_planes_shift + zp))*jerk_strategy.w_search_step - jerk_strategy.w_max_search_limit;
					
					/*
					//-------------------- DIRECT PLANE OUTPUT
					if(active_DM<4){
						float *h_ZW_planes;
						h_ZW_planes = (float *)malloc(ZW_planes_size_bytes);
						checkCudaErrors(cudaMemcpy(h_ZW_planes, d_ZW_planes, ZW_planes_size_bytes, cudaMemcpyDeviceToHost));
						
						printf("---> Exporting plane:"); fflush(stdout);
						sprintf(filename, "TEST_conv_plane_dm%d-w%d.dat", active_DM, (int) ZW_planes_shift + zp);
						FILEOUT.open(filename);
						for(int f=0; f<jerk_strategy.nFilters_z; f++){
							float z_pos = (f-jerk_strategy.nFilters_z_half)*jerk_strategy.z_search_step;
							for(int s=0; s<jerk_strategy.output_size_one_DM; s++){
								int pos = f*jerk_strategy.output_size_one_DM + s;
								FILEOUT << z_pos << " " << (float) (s/(local_plan.nTimesamples*sampling_time*inBin[active_range])) << " " << h_ZW_planes[pos] << std::endl;
							}
							printf(".");
							fflush(stdout);
						}
						FILEOUT.close();
						printf("\n");
						free(h_ZW_planes);
					}
					//-------------------- DIRECT PLANE OUTPUT
					*/
					
					//------------> Calculation of the mean and standard deviation
					timer.Start();
					size_t pos = ((size_t) zp)*((size_t) ZW_plane_size);
					Find_MSD(d_MSD, &d_ZW_planes[pos], d_MSD_workarea, &MSD_conf, user_plan.OR_sigma_cuttoff, user_plan.MSD_outlier_rejection);
					timer.Stop(); printf("JERK SEARCH -> MSD took: %g ms\n", timer.Elapsed());
					MSD_time += timer.Elapsed();
					
					checkCudaErrors(cudaGetLastError());
					
					//------------> Candidate selection
					timer.Start();
					
					cudaMemset((void*) gmem_peak_pos, 0, sizeof(unsigned int));
					PEAK_FIND_FOR_FDAS(&d_ZW_planes[pos], d_ZW_candidates, d_MSD, jerk_strategy.nFilters_z, jerk_strategy.output_size_one_DM, user_plan.CS_sigma_threshold, max_nZWCandidates, gmem_peak_pos, w);
					
					timer.Stop();
					printf("JERK SEARCH -> Candidate selection took: %g ms\n", timer.Elapsed());
					Candidate_time += timer.Elapsed();
					printf("JERK SEARCH -> W Coordinate: %d [%d;%d] \n", w, ZW_planes_shift, zp);
				
					checkCudaErrors(cudaGetLastError());
				
					//------------> Export to host
					checkCudaErrors(cudaMemcpy(&nCandidates, gmem_peak_pos, sizeof(unsigned int), cudaMemcpyDeviceToHost));
					
					
					//-------------------- DIRECT CANDIDATE OUTPUT
					if(nCandidates>0){
						float *h_candidates;
						h_candidates = new float [nCandidates*4];
						cudaMemcpy(h_candidates, d_ZW_candidates, nCandidates*4*sizeof(float), cudaMemcpyDeviceToHost);
						
						for(int f=0; f<nCandidates; f++){
							h_candidates[4*f] = h_candidates[4*f]*jerk_strategy.z_search_step - jerk_strategy.z_max_search_limit;
							h_candidates[4*f + 1] = h_candidates[4*f + 1]/(sampling_time*((float) inBin[active_range])*((double) user_plan.nTimesamples));
						}
						
						FILE *fp_out;
						sprintf(filename, "jerk_search-dm%d-w%d.dat", active_DM, (int) (ZW_planes_shift + zp));
						if (( fp_out = fopen(filename, "wb") ) == NULL)	{
							fprintf(stderr, "Error opening output file!\n");
							exit(0);
						}
						fwrite(h_candidates, nCandidates*sizeof(float), 4, fp_out);
						fclose(fp_out);
						
						delete [] h_candidates;
					}
					//-------------------- DIRECT CANDIDATE OUTPUT
					
					
					float DM = dm_low[active_range] + ((float) active_DM)*dm_step[active_range];
					allcandidates.AddSubListFromGPU(nCandidates, d_ZW_candidates, w, DM, jerk_strategy.nSamples_time_dom, jerk_strategy.nFilters_z, jerk_strategy.z_max_search_limit, jerk_strategy.z_search_step, sampling_time*((float) inBin[active_range]), inBin[active_range]);
					printf("JERK SEARCH -> Number of candidates: %d\n", nCandidates);
					
					checkCudaErrors(cudaGetLastError());
				}
				
				ZW_planes_shift = ZW_planes_shift + nZW_planes;
			}
			//------->
			char str[100];
			sprintf(str, "jerk_results_r%d_dm%d.dat", active_range, active_DM);
			allcandidates.ExportToFile(str);
			
			//Save candidates to disc
		}
		
		//-------> cuFFT
		cufftDestroy(cuFFT_plan);
		
		checkCudaErrors(cudaGetLastError());

		if ( cudaSuccess != cudaFree(d_ZW_candidates)) printf("ERROR while deallocating d_ZW_candidates!\n");
		d_ZW_candidates = NULL;
		if ( cudaSuccess != cudaFree(d_ZW_planes)) printf("ERROR while deallocating d_ZW_planes!\n");
		d_ZW_planes = NULL;
		if ( cudaSuccess != cudaFree(d_MSD_workarea)) printf("ERROR while deallocating d_MSD_workarea!\n");
		d_MSD_workarea = NULL;
		
		checkCudaErrors(cudaGetLastError());
		
		break;
	}
	
	
	
	timer_total.Stop();
	time_total = timer_total.Elapsed()/1000.0;
	printf("Total time for JERK search: %g s \n", time_total);
	printf("Time for convolution: %g s \n", Convolution_time);
	printf("Time for candidate selection: %g s \n", Candidate_time);
	printf("Time for MSD: %g s \n", MSD_time);


	
	delete [] h_jerk_filters;
	cudaFree(d_jerk_filters);
	cudaFree(d_DM_trial);
	cudaFree(d_DM_trial_ffted);
	cudaFree(d_MSD);
	return(0);
}
