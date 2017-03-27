//#define OLD_THRESHOLD

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include "headers/params.h"

#include "headers/device_BC_plan.h"
#include "headers/device_peak_find.h"
#include "headers/device_BLN.h"
#include "headers/device_MSD_plane.h"
#include "headers/device_MSD_limited.h"
#include "headers/device_SNR_limited.h"
#include "headers/device_SPS_inplace.h"
#include "headers/device_threshold.h"
#include "headers/device_single_FIR.h"

#include "timer.h"


void Create_PD_plan(std::vector<PulseDetection_plan> *PD_plan, std::vector<int> *BC_widths, int nDMs, int nTimesamples){
	int Elements_per_block, itemp, nRest;
	PulseDetection_plan PDmp;
	
	if(BC_widths->size()>0){
		PDmp.shift        = 0;
		PDmp.output_shift = 0;
		PDmp.startTaps    = 0;
		PDmp.iteration    = 0;
		
		PDmp.decimated_timesamples = nTimesamples;
		PDmp.nBoxcars = BC_widths->operator[](0);
		Elements_per_block = PD_NTHREADS*2 - PDmp.nBoxcars;
		itemp = PDmp.decimated_timesamples;
		PDmp.nBlocks = itemp/Elements_per_block;
		nRest = itemp - PDmp.nBlocks*Elements_per_block;
		if(nRest>0) PDmp.nBlocks++;
		PDmp.unprocessed_samples = PDmp.nBoxcars + 6;
		PDmp.total_ut = PDmp.unprocessed_samples;
		
		PD_plan->push_back(PDmp);
		
		for(int f=1; f< (int) BC_widths->size(); f++){
			PDmp.shift        = PDmp.nBoxcars/2;
			PDmp.output_shift = PDmp.output_shift + nDMs*PDmp.decimated_timesamples;
			PDmp.startTaps    = PDmp.startTaps + PDmp.nBoxcars*(1<<PDmp.iteration);
			PDmp.iteration    = PDmp.iteration + 1;
			
			PDmp.decimated_timesamples = (nTimesamples>>PDmp.iteration);
			PDmp.nBoxcars = BC_widths->operator[](f);
			Elements_per_block=PD_NTHREADS*2 - PDmp.nBoxcars;
			itemp = PDmp.decimated_timesamples;
			PDmp.nBlocks = itemp/Elements_per_block;
			nRest = itemp - PDmp.nBlocks*Elements_per_block;
			if(nRest>0) PDmp.nBlocks++;
			PDmp.unprocessed_samples = PDmp.unprocessed_samples/2 + PDmp.nBoxcars + 6; //
			PDmp.total_ut = PDmp.unprocessed_samples*(1<<PDmp.iteration);
			
			PD_plan->push_back(PDmp);
		}
	}
}


int Get_max_iteration(int max_boxcar_width, std::vector<int> *BC_widths){
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
	
	if(max_boxcar_width>startTaps) iteration=(int) BC_widths->size();
	
	return(iteration);
}

void analysis_GPU(float *h_output_list, size_t *list_pos, size_t max_list_size, float *h_peak_list, size_t *peak_pos, size_t max_peak_size, int i, float tstart, int t_processed, int inBin, int outBin, int *maxshift, int max_ndms, int *ndms, float cutoff, float sigma_constant, float max_boxcar_width_in_sec, float *output_buffer, float *dm_low, float *dm_high, float *dm_step, float tsamp){
	int max_boxcar_width = (int) (max_boxcar_width_in_sec/tsamp);
	//unsigned long int j;
	unsigned long int vals;
	int nTimesamples = t_processed;
	int nDMs = ndms[i];
	int  temp_peak_pos; //temp_list_pos,
	//double total;

	// Calculate the total number of values
	vals = (unsigned long int) ( nDMs*nTimesamples );
	

	double total_time, partial_time;
	float signal_mean_1, signal_sd_1;//, modifier;
	//float max, min, threshold;
	int offset, max_iteration;
	int t_BC_widths[10]={PD_MAXTAPS,16,16,16,8,8,8,8,8,8};
	std::vector<int> BC_widths(t_BC_widths,t_BC_widths+sizeof(t_BC_widths)/sizeof(int));
	std::vector<PulseDetection_plan> PD_plan;

	//---------------------------------------------------------------------------
	//----------> GPU part
	printf("\n GPU analysis part\n\n");
	printf("Dimensions nDMs:%d; nTimesamples:%d; inBin:%d; outBin:%d; maxshift:%d; \n", ndms[i], t_processed, inBin, outBin, *maxshift);
	GpuTimer timer;
	
	float h_MSD[3];
	float *d_MSD;
	if ( cudaSuccess != cudaMalloc((void**) &d_MSD, sizeof(float)*3)) printf("Allocation error!\n");

	total_time = 0;
	
	//-------------- Calculating base level noise using outlier rejection
	timer.Start();
	BLN(output_buffer, d_MSD, 32, 32, nDMs, nTimesamples, 128, sigma_constant); // Those 128 are there because there was a problem with data, I'm not sure if it is still the case.
	timer.Stop();
	partial_time = timer.Elapsed();
	total_time += partial_time;
	printf("MSD limited took:%f ms\n", partial_time);
	
	cudaMemcpy(h_MSD, d_MSD, 3*sizeof(float), cudaMemcpyDeviceToHost);
	signal_mean_1 = h_MSD[0];
	signal_sd_1 = h_MSD[1];
	printf("Bin: %d, Mean: %f, Stddev: %f\n", 1, signal_mean_1, signal_sd_1);
	//-------------- Calculating base level noise using outlier rejection
	
	
	size_t free_mem,total_mem;
	cudaMemGetInfo(&free_mem,&total_mem);
	printf("Memory required by boxcar filters:%0.3f MB\n",(4.5*vals*sizeof(float) + 2*vals*sizeof(ushort))/(1024.0*1024) );
	printf("Memory available:%0.3f MB \n", ((float) free_mem)/(1024.0*1024.0) );
	
	std::vector<int> DM_list;
	unsigned long int max_timesamples=(free_mem*0.95)/(5.5*sizeof(float) + 2*sizeof(ushort));
	#ifdef OLD_THRESHOLD
	max_timesamples=(free_mem*0.95)/(5.5*sizeof(float) + 2*sizeof(ushort));
	#endif
	int DMs_per_cycle = max_timesamples/nTimesamples;
	int nRepeats, nRest, DM_shift, itemp, local_max_list_size;//BC_shift,
	
	itemp = (int) (DMs_per_cycle/THR_WARPS_PER_BLOCK);
	DMs_per_cycle = itemp*THR_WARPS_PER_BLOCK;
	
	nRepeats = nDMs/DMs_per_cycle;
	nRest = nDMs - nRepeats*DMs_per_cycle;
	local_max_list_size = (DMs_per_cycle*nTimesamples)/4;
	
	for(int f=0; f<nRepeats; f++) DM_list.push_back(DMs_per_cycle);
	if(nRest>0) DM_list.push_back(nRest);
	
	
	max_iteration = Get_max_iteration(max_boxcar_width/inBin, &BC_widths);
	printf("Selected iteration:%d; for maximum boxcar width:%d;\n", max_iteration, max_boxcar_width/inBin);
	Create_PD_plan(&PD_plan, &BC_widths, 1, nTimesamples);
	
	if(DM_list.size()>0){
		DMs_per_cycle = DM_list[0];
		
		
		
		//---------> Old thresholding code.
		#ifdef OLD_THRESHOLD
		float *d_list;
		if ( cudaSuccess != cudaMalloc((void**) &d_list, sizeof(float)*DMs_per_cycle*nTimesamples)) printf("Allocation error! list\n");
		#endif
		//---------> Old thresholding code.

		float *d_peak_list;
		if ( cudaSuccess != cudaMalloc((void**) &d_peak_list, sizeof(float)*DMs_per_cycle*nTimesamples)) printf("Allocation error! peaks\n");
		
		float *d_decimated;
		if ( cudaSuccess != cudaMalloc((void **) &d_decimated,  sizeof(float)*((DMs_per_cycle*nTimesamples)/2))) printf("Allocation error! dedispered\n");
		
		float *d_boxcar_values;
		if ( cudaSuccess != cudaMalloc((void **) &d_boxcar_values,  sizeof(float)*DMs_per_cycle*nTimesamples)) printf("Allocation error! boxcars\n");
		
		float *d_output_SNR;
		if ( cudaSuccess != cudaMalloc((void **) &d_output_SNR, sizeof(float)*2*DMs_per_cycle*nTimesamples)) printf("Allocation error! SNR\n");
		
		ushort *d_output_taps;
		if ( cudaSuccess != cudaMalloc((void **) &d_output_taps, sizeof(ushort)*2*DMs_per_cycle*nTimesamples)) printf("Allocation error! taps\n");
		
		//---------> Old thresholding code.
		#ifdef OLD_THRESHOLD
		int *gmem_list_pos;
		cudaMalloc((void**) &gmem_list_pos, 1*sizeof(int));
		cudaMemset((void*) gmem_list_pos, 0, sizeof(int));
		#endif
		//---------> Old thresholding code.
		
		int *gmem_peak_pos;
		cudaMalloc((void**) &gmem_peak_pos, 1*sizeof(int));
		cudaMemset((void*) gmem_peak_pos, 0, sizeof(int));
		
		DM_shift = 0;
		for(int f=0; f<DM_list.size(); f++) {
			//-------------- SPS BLN
			timer.Start();
			//offset=PD_SEARCH_LONG_BLN(&output_buffer[DM_shift*nTimesamples], d_boxcar_values, d_decimated, d_output_SNR, d_output_taps, d_MSD, max_boxcar_width/inBin, DM_list[f], nTimesamples, &max_iteration);
			offset=PD_SEARCH_LONG_BLN_IF(&output_buffer[DM_shift*nTimesamples], d_boxcar_values, d_decimated, d_output_SNR, d_output_taps, d_MSD, &PD_plan, max_iteration, DM_list[f], nTimesamples);
			timer.Stop();
			partial_time = timer.Elapsed();
			total_time += partial_time;
			printf("PD_SEARCH took:%f ms\n", partial_time);
			//-------------- SPS BLN
			
			//-------------- Thresholding
			#ifdef OLD_THRESHOLD
			timer.Start();
			THRESHOLD(d_output_SNR, d_output_taps, d_list, gmem_list_pos, cutoff, DM_list[f], nTimesamples, DM_shift, &PD_plan, max_iteration, local_max_list_size);
			timer.Stop();
			partial_time = timer.Elapsed();
			total_time += partial_time;
			printf("THR_WARP took:%f ms\n", partial_time);
			#endif
			//-------------- Thresholding
			
			printf("BC_shift:%d; DMs_per_cycle:%d; f*DMs_per_cycle:%d; max_iteration:%d; offset:%d;\n", DM_shift*nTimesamples, DM_list[f], DM_shift, max_iteration, offset);
			
			//-------------- Peak finding
			timer.Start();
			PEAK_FIND(d_output_SNR, d_output_taps, d_peak_list, DM_list[f], nTimesamples, cutoff, local_max_list_size, gmem_peak_pos, DM_shift, &PD_plan, max_iteration);
			timer.Stop();
			partial_time = timer.Elapsed();
			total_time += partial_time;
			printf("PEAK_FIND took:%f ms\n", partial_time);
			//-------------- Peak finding
			
			DM_shift = DM_shift + DM_list[f];
			
			//---------> Old thresholding code.
			#ifdef OLD_THRESHOLD
			cudaMemcpy(&temp_list_pos, gmem_list_pos, sizeof(int), cudaMemcpyDeviceToHost);
			printf("temp_peak_pos:%d; host_pos:%d; max:%d;\n", temp_list_pos, (*list_pos), (int) max_list_size);
			if( ((*list_pos) + temp_list_pos)<max_list_size){
				cudaMemcpy(&h_output_list[(*list_pos)*4], d_list, temp_list_pos*4*sizeof(float), cudaMemcpyDeviceToHost);
				*list_pos = (*list_pos) + temp_list_pos;
			}
			else printf("Error list is too small!\n");
			#endif
			//---------> Old thresholding code.
			
			
			cudaMemcpy(&temp_peak_pos, gmem_peak_pos, sizeof(int), cudaMemcpyDeviceToHost);
			printf("temp_peak_pos:%d; host_pos:%d; max:%d;\n", temp_peak_pos, (*peak_pos), (int) max_peak_size);
			if( ((*peak_pos) + temp_peak_pos)<max_peak_size){
				cudaMemcpy(&h_peak_list[(*peak_pos)*4], d_peak_list, temp_peak_pos*4*sizeof(float), cudaMemcpyDeviceToHost);
				*peak_pos = (*peak_pos) + temp_peak_pos;
			}
			else printf("Error peak list is too small!\n");
			
			//---------> Old thresholding code.
			#ifdef OLD_THRESHOLD
			cudaMemset((void*) gmem_list_pos, 0, sizeof(int));
			#endif
			//---------> Old thresholding code.
			cudaMemset((void*) gmem_peak_pos, 0, sizeof(int));
		}
		
		//---------> Old thresholding code.
		#ifdef OLD_THRESHOLD
		cudaFree(d_list);
		#endif
		//---------> Old thresholding code.
		cudaFree(d_peak_list);
		cudaFree(d_boxcar_values);
		cudaFree(d_decimated);
		cudaFree(d_output_SNR);
		cudaFree(d_output_taps);
		//---------> Old thresholding code.
		#ifdef OLD_THRESHOLD
		cudaFree(gmem_list_pos);
		#endif
		//---------> Old thresholding code.
		cudaFree(gmem_peak_pos);

	}
	else printf("Error not enough memory to search for pulses\n");

	
	printf("\n====> TOTAL TIME:%f\n\n", total_time);

	cudaFree(d_MSD);
	//----------> GPU part
	//---------------------------------------------------------------------------
	
}

