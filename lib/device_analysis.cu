//#define GPU_ANALYSIS_DEBUG
//#define MSD_BOXCAR_TEST
//#define GPU_PARTIAL_TIMER
#define GPU_TIMER
#define VERBOSE 1

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include "headers/params.h"

#include "headers/device_MSD_Parameters.h"
#include "headers/device_MSD_plane_profile.h"

#include "headers/device_SPS_DataDescription.h"
#include "headers/device_SPS_Parameters.h"
#include "headers/device_SPS_BC_plan.h"
#include "headers/device_SPS_long.h"
#include "headers/device_peak_find.h"
#include "headers/device_threshold.h"

#include "timer.h"

//TODO:
// Make BC_plan for arbitrary long pulses, by reusing last element in the plane


void Create_list_of_boxcar_widths(std::vector<int> *boxcar_widths, std::vector<int> *BC_widths, int max_boxcar_width){
	int DIT_value, DIT_factor, width;
	DIT_value = 1;
	DIT_factor = 2;
	width = 0;
	for(int f=0; f<(int) BC_widths->size(); f++){
		for(int b=0; b<BC_widths->operator[](f); b++){
			width = width + DIT_value;
			if(width<=max_boxcar_width){
				boxcar_widths->push_back(width);
			}
		}
		DIT_value = DIT_value*DIT_factor;
	}
}


// Extend this to arbitrary size plans
void Create_PD_plan(std::vector<PulseDetection_plan> *PD_plan, std::vector<int> *BC_widths, int nTimesamples){
	int Elements_per_block, itemp, nRest;
	PulseDetection_plan PDmp;
	
	if(BC_widths->size()>0){
		PDmp.shift        = 0;
		PDmp.output_shift = 0;
		PDmp.startTaps    = 0;
		PDmp.iteration    = 0;
		
		PDmp.decimated_timesamples = nTimesamples;
		PDmp.dtm = (nTimesamples>>(PDmp.iteration+1));
		PDmp.dtm = PDmp.dtm - (PDmp.dtm&1);
		
		PDmp.nBoxcars = BC_widths->operator[](0);
		Elements_per_block = PD_NTHREADS*2 - PDmp.nBoxcars;
		itemp = PDmp.decimated_timesamples;
		PDmp.nBlocks = itemp/Elements_per_block;
		nRest = itemp - PDmp.nBlocks*Elements_per_block;
		if(nRest>0) PDmp.nBlocks++;
		PDmp.unprocessed_samples = PDmp.nBoxcars + 6;
		if(PDmp.decimated_timesamples<PDmp.unprocessed_samples) PDmp.nBlocks=0;
		PDmp.total_ut = PDmp.unprocessed_samples;
		
		
		PD_plan->push_back(PDmp);
		
		for(int f=1; f< (int) BC_widths->size(); f++){
			// These are based on previous values of PDmp
			PDmp.shift        = PDmp.nBoxcars/2;
			PDmp.output_shift = PDmp.output_shift + PDmp.decimated_timesamples;
			PDmp.startTaps    = PDmp.startTaps + PDmp.nBoxcars*(1<<PDmp.iteration);
			PDmp.iteration    = PDmp.iteration + 1;
			
			// Definition of new PDmp values
			PDmp.decimated_timesamples = PDmp.dtm;
			PDmp.dtm = (nTimesamples>>(PDmp.iteration+1));
			PDmp.dtm = PDmp.dtm - (PDmp.dtm&1);
			
			PDmp.nBoxcars = BC_widths->operator[](f);
			Elements_per_block=PD_NTHREADS*2 - PDmp.nBoxcars;
			itemp = PDmp.decimated_timesamples;
			PDmp.nBlocks = itemp/Elements_per_block;
			nRest = itemp - PDmp.nBlocks*Elements_per_block;
			if(nRest>0) PDmp.nBlocks++;
			PDmp.unprocessed_samples = PDmp.unprocessed_samples/2 + PDmp.nBoxcars + 6; //
			if(PDmp.decimated_timesamples<PDmp.unprocessed_samples) PDmp.nBlocks=0;
			PDmp.total_ut = PDmp.unprocessed_samples*(1<<PDmp.iteration);
			
			PD_plan->push_back(PDmp);
		}
	}
}


int Get_max_iteration(int max_boxcar_width, std::vector<int> *BC_widths, int *max_width_performed){
	int startTaps, iteration;
	
	startTaps = 0;
	iteration = 0;
	for(int f=0; f<(int) BC_widths->size(); f++){
		startTaps = startTaps + BC_widths->operator[](f)*(1<<f);
		if(startTaps>=max_boxcar_width) {
			iteration = f+1;
			break;
		}
	}
	
	if(max_boxcar_width>startTaps) {
		iteration=(int) BC_widths->size();
	}
	
	*max_width_performed=startTaps;
	return(iteration);
}


void analysis_GPU( float *h_peak_list, size_t *peak_pos, size_t max_peak_size, SPS_DataDescription SPS_data, float *d_SPS_input, SPS_Parameters *SPS_params, MSD_Parameters *MSD_params){
	// Definition of some local variables
	float  local_tsamp  = SPS_data.sampling_time*SPS_data.inBin; // corrected sampling time
	size_t nTimesamples = SPS_data.nTimesamples;
	size_t nDMs         = SPS_data.nDMs;
	if(SPS_params->verbose) {
		printf("\n----------> GPU analysis part\n");
		printf("  Dimensions: nTimesamples:%zu; nDMs:%zu; inBin:%d; sampling time: %f; corrected s. time: %f;\n", nTimesamples, nDMs, SPS_data.inBin, SPS_data.sampling_time, local_tsamp);
	}
	
	//--------> Definition of SPDT boxcar plan
	int max_desired_boxcar_width = (int) (SPS_params->max_boxcar_width_in_sec/local_tsamp);
	int max_width_performed = 0;
	int t_BC_widths[10]={PD_MAXTAPS,16,16,16,8,8,8,8,8,8};
	std::vector<int> BC_widths(t_BC_widths,t_BC_widths+sizeof(t_BC_widths)/sizeof(int));
	std::vector<PulseDetection_plan> PD_plan;
	Create_PD_plan(&PD_plan, &BC_widths, nTimesamples); //PD_plan is independent on maximum boxcar width. which is wrong?
	int max_iteration = Get_max_iteration(max_desired_boxcar_width, &BC_widths, &max_width_performed);
	if(SPS_params->verbose) 
		printf("  Selected iteration:%d; maximum boxcar width requested:%d; maximum boxcar width performed:%d;\n", max_iteration, max_desired_boxcar_width, max_width_performed);
	std::vector<int> h_boxcar_widths;
	Create_list_of_boxcar_widths(&h_boxcar_widths, &BC_widths, max_width_performed);
	
	// It should be like this:
	//   SPS_params should contain BC_widths
	//   SPS_params should also contain function get_maximum_iteration which would give number of iterations required to achieve user defined value in form of max_desired_boxcar_width
	//   Based on maximum_iteration SPS should build PD_plan
	//   Proper error check must be placed so SPS would not die if user chooses wrong maximum search width
	
	//--------> Benchmarking
	double total_time=0, MSD_time=0, SPDT_time=0, PF_time=0;
	

	//---------------------------------------------------------------------------
	//----------> GPU part
	GpuTimer total_timer, timer;
	total_timer.Start();
	
	
	size_t free_mem,total_mem;
	cudaMemGetInfo(&free_mem,&total_mem);
	if(SPS_params->verbose) printf("  Memory required by boxcar filters:%0.3f MB\n",(4.5*nTimesamples*nDMs*sizeof(float) + 2*nTimesamples*nDMs*sizeof(ushort))/(1024.0*1024) );
	if(SPS_params->verbose) printf("  Memory available:%0.3f MB \n", ((float) free_mem)/(1024.0*1024.0) );
	
	
	//-------------------------------------------------------------------------
	//---------> Comparison between interpolated values and computed values
	#ifdef MSD_BOXCAR_TEST
		MSD_plane_profile_boxcars(d_SPS_input, nTimesamples, nDMs, &h_boxcar_widths, MSD_params->OR_sigma_multiplier, SPS_data.d,_low, SPS_data.dm_high, SPS_data.time_start);
	#endif
	//---------> Comparison between interpolated values and computed values
	//-------------------------------------------------------------------------
	
	
	
	//-------------------------------------------------------------------------
	//------------ Using MSD_plane_profile
	size_t MSD_profile_size_in_bytes, MSD_DIT_profile_size_in_bytes, workarea_size_in_bytes;
	cudaMemGetInfo(&free_mem,&total_mem);
	Get_MSD_plane_profile_memory_requirements(&MSD_profile_size_in_bytes, &MSD_DIT_profile_size_in_bytes, &workarea_size_in_bytes, nTimesamples, nDMs, &h_boxcar_widths);
	double dit_time, MSD_only_time;
	float *d_MSD_interpolated;
	float *d_MSD_DIT = NULL;
	float *temporary_workarea;
	cudaMalloc((void **) &d_MSD_interpolated, MSD_profile_size_in_bytes);
	cudaMalloc((void **) &temporary_workarea, workarea_size_in_bytes);
	
	MSD_plane_profile(d_MSD_interpolated, d_SPS_input, d_MSD_DIT, temporary_workarea, false, nTimesamples, nDMs, &h_boxcar_widths, SPS_data.time_start, SPS_data.dm_low, SPS_data.dm_high, MSD_params->OR_sigma_multiplier, MSD_params->enable_outlier_rejection, false, &MSD_time, &dit_time, &MSD_only_time);
	
	#ifdef GPU_PARTIAL_TIMER
		printf("    MSD time: Total: %f ms; DIT: %f ms; MSD: %f ms;\n", MSD_time, dit_time, MSD_only_time);
	#endif
	
	cudaFree(temporary_workarea);
	//------------ Using MSD_plane_profile
	//-------------------------------------------------------------------------	
	
	
	//-------------------------------------------------------------------------
	//------------ Splitting input data into chunks
	std::vector<int> DM_list;
	unsigned long int max_timesamples=(free_mem*0.95)/(5.5*sizeof(float) + 2*sizeof(ushort));
	int DMs_per_cycle = max_timesamples/nTimesamples;
	int nRepeats, nRest, DM_shift, itemp, local_max_list_size;//BC_shift,
	
	itemp = (int) (DMs_per_cycle/THR_WARPS_PER_BLOCK);
	DMs_per_cycle = itemp*THR_WARPS_PER_BLOCK;
	
	nRepeats = nDMs/DMs_per_cycle;
	nRest = nDMs - nRepeats*DMs_per_cycle;
	local_max_list_size = (DMs_per_cycle*nTimesamples)/4;
	
	for(int f=0; f<nRepeats; f++) DM_list.push_back(DMs_per_cycle);
	if(nRest>0) DM_list.push_back(nRest);
	
	if( (int) DM_list.size() > 1 ) 
		printf("  SPS will run %d batches each containing %d DM trials. Remainder %d DM trials\n", (int) DM_list.size(), DMs_per_cycle, nRest);
	else 
		printf("  SPS will run %d batch containing %d DM trials.\n", (int) DM_list.size(), nRest);
	//------------ Splitting input data into chunks
	//-------------------------------------------------------------------------	
	
	
	
	if(DM_list.size()>0){
		DMs_per_cycle = DM_list[0];
		
		float *d_peak_list;
		if ( cudaSuccess != cudaMalloc((void**) &d_peak_list, sizeof(float)*DMs_per_cycle*nTimesamples)) printf("Allocation error! peaks\n");
		
		float *d_decimated;
		if ( cudaSuccess != cudaMalloc((void **) &d_decimated,  sizeof(float)*(((DMs_per_cycle*nTimesamples)/2)+PD_MAXTAPS) )) printf("Allocation error! dedispered\n");
		
		float *d_boxcar_values;
		if ( cudaSuccess != cudaMalloc((void **) &d_boxcar_values,  sizeof(float)*DMs_per_cycle*nTimesamples)) printf("Allocation error! boxcars\n");
		
		float *d_output_SNR;
		if ( cudaSuccess != cudaMalloc((void **) &d_output_SNR, sizeof(float)*2*DMs_per_cycle*nTimesamples)) printf("Allocation error! SNR\n");
		
		ushort *d_output_taps;
		if ( cudaSuccess != cudaMalloc((void **) &d_output_taps, sizeof(ushort)*2*DMs_per_cycle*nTimesamples)) printf("Allocation error! taps\n");
		
		int *gmem_peak_pos;
		cudaMalloc((void**) &gmem_peak_pos, 1*sizeof(int));
		cudaMemset((void*) gmem_peak_pos, 0, sizeof(int));
		
		DM_shift = 0;
		for(int f=0; f<DM_list.size(); f++) {
			//-------------- SPDT
			timer.Start();
			SPDT_search_long_MSD_plane(&d_SPS_input[DM_shift*nTimesamples], d_boxcar_values, d_decimated, d_output_SNR, d_output_taps, d_MSD_interpolated, &PD_plan, max_iteration, nTimesamples, DM_list[f]);
			timer.Stop();
			SPDT_time += timer.Elapsed();
			#ifdef GPU_PARTIAL_TIMER
			printf("    SPDT took:%f ms\n", timer.Elapsed());
			#endif
			//-------------- SPDT
			
			checkCudaErrors(cudaGetLastError());
			
			#ifdef GPU_ANALYSIS_DEBUG
			printf("    BC_shift:%d; DMs_per_cycle:%d; f*DMs_per_cycle:%d; max_iteration:%d;\n", DM_shift*nTimesamples, DM_list[f], DM_shift, max_iteration);
			#endif
			
			if(SPS_params->candidate_algorithm==1){
				//-------------- Thresholding
				timer.Start();
				THRESHOLD(d_output_SNR, d_output_taps, d_peak_list, gmem_peak_pos, SPS_params->sigma_cutoff, DM_list[f], nTimesamples, DM_shift, &PD_plan, max_iteration, local_max_list_size, SPS_data.dm_step, SPS_data.dm_low, local_tsamp, SPS_data.inBin, SPS_data.time_start);
				timer.Stop();
				PF_time += timer.Elapsed();
				#ifdef GPU_PARTIAL_TIMER
				printf("    Thresholding took:%f ms\n", timer.Elapsed());
				#endif
				//-------------- Thresholding
			}
			else {
				//-------------- Peak finding
				timer.Start();
				PEAK_FIND(d_output_SNR, d_output_taps, d_peak_list, DM_list[f], nTimesamples, SPS_params->sigma_cutoff, local_max_list_size, gmem_peak_pos, DM_shift, &PD_plan, max_iteration, SPS_data.dm_step, SPS_data.dm_low, local_tsamp, SPS_data.inBin, SPS_data.time_start);
				timer.Stop();
				PF_time = timer.Elapsed();
				#ifdef GPU_PARTIAL_TIMER
				printf("    Peak finding took:%f ms\n", timer.Elapsed());
				#endif
				//-------------- Peak finding
			}
			
			checkCudaErrors(cudaGetLastError());
			
			int temp_peak_pos = 0;
			checkCudaErrors(cudaMemcpy(&temp_peak_pos, gmem_peak_pos, sizeof(int), cudaMemcpyDeviceToHost));
			#ifdef GPU_ANALYSIS_DEBUG
			printf("    Candidates found:%d; Total #candidates for this chunk:%zu; Maximum #candidates:%zu; Local max. #candidates:%d;\n", temp_peak_pos, (*peak_pos), max_peak_size, local_max_list_size);
			#endif
			if( temp_peak_pos>=local_max_list_size ) {
				printf("    WARNING: Maximum list size reached! Not all candidates will be saved. You can increase sigma cutoff.\n");
				temp_peak_pos = local_max_list_size;
			}
			if( ((*peak_pos) + temp_peak_pos)<max_peak_size){
				checkCudaErrors(cudaMemcpy(&h_peak_list[(*peak_pos)*4], d_peak_list, temp_peak_pos*4*sizeof(float), cudaMemcpyDeviceToHost));
				*peak_pos = (*peak_pos) + temp_peak_pos;
			}
			else printf("    ERROR: Not enough memory to store all candidates on the host!\n");

			DM_shift = DM_shift + DM_list[f];
			cudaMemset((void*) gmem_peak_pos, 0, sizeof(int));
		}
		
		cudaFree(d_peak_list);
		cudaFree(d_boxcar_values);
		cudaFree(d_decimated);
		cudaFree(d_output_SNR);
		cudaFree(d_output_taps);
		cudaFree(gmem_peak_pos);
		cudaFree(d_MSD_DIT);
		cudaFree(d_MSD_interpolated);

	}
	else printf("Error not enough memory to search for pulses\n");

	total_timer.Stop();
	total_time = total_timer.Elapsed();
	#ifdef GPU_TIMER
	printf("\n  TOTAL TIME OF SPS:%f ms\n", total_time);
	printf("  MSD_time: %f ms; SPDT time: %f ms; Candidate selection time: %f ms;\n", MSD_time, SPDT_time, PF_time);
	printf("----------<\n\n");
	#endif
	//----------> GPU part
	//---------------------------------------------------------------------------
	
}
