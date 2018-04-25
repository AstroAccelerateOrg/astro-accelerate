//#define GPU_ANALYSIS_DEBUG
//#define MSD_BOXCAR_TEST
//#define GPU_PARTIAL_TIMER
#define GPU_TIMER

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include "headers/params.h"

#include "headers/device_BC_plan.h"
#include "headers/device_peak_find.h"
#include "headers/device_MSD_plane_profile.h"
#include "headers/device_SPS_long.h"
#include "headers/device_threshold.h"

#include "timer.h"

//TODO:
// Make BC_plan for arbitrary long pulses, by reusing last element in the plane

void CUDART_CB MyCallback(cudaStream_t stream, cudaError_t status, void *data){
    printf("Inside callback %d\n", (size_t)data);
}


void Create_list_of_boxcar_widths(std::vector<int> *boxcar_widths, std::vector<int> *BC_widths){
	int DIT_value, DIT_factor, width;
	DIT_value = 1;
	DIT_factor = 2;
	width = 0;
	for(int f=0; f<(int) BC_widths->size(); f++){
		for(int b=0; b<BC_widths->operator[](f); b++){
			width = width + DIT_value;
			boxcar_widths->push_back(width);
		}
		DIT_value = DIT_value*DIT_factor;
	}
}

// Extend this to arbitrary size plans
void Create_PD_plan(std::vector<PulseDetection_plan> *PD_plan, std::vector<int> *BC_widths, int nDMs, int nTimesamples){
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


void analysis_GPU(float *h_peak_list, size_t *peak_pos, size_t max_peak_size, int i, float tstart, int t_processed, int inBin, int outBin, int *maxshift, int max_ndms, int *ndms, float cutoff, float OR_sigma_multiplier, float max_boxcar_width_in_sec, float *output_buffer, float *dm_low, float *dm_high, float *dm_step, float tsamp, cudaStream_t streams, int stream_id, int candidate_algorithm, int enable_sps_baselinenoise, float *MSD_workarea, unsigned short *d_output_taps, float *d_MSD_interpolated, float *h_MSD_DIT, float *h_MSD_interpolated, int *gmem_peak_pos, int *temp_peak_pos, unsigned long int maxTimeSamples){
//  void analysis_GPU(float *h_peak_list, size_t *peak_pos, size_t max_peak_size, int i, float tstart, int t_processed, int inBin, int outBin, int *maxshift, int max_ndms, int *ndms, float cutoff, float OR_sigma_multiplier, float max_boxcar_width_in_sec, float *output_buffer, float *dm_low, float *dm_high, float *dm_step, float tsamp, cudaStream_t streams, int candidate_algorithm, int enable_sps_baselinenoise, float *MSD_workarea, unsigned short *d_output_taps, float *d_MSD_interpolated, int *gmem_peak_pos, unsigned long int maxTimeSamples){
	//--------> Task
	int max_boxcar_width = (int) (max_boxcar_width_in_sec/tsamp);
	int max_width_performed=0;
//	size_t vals;
	int nTimesamples = t_processed;
	int nDMs = ndms[i];
//	int *temp_peak_pos;
//	cudaMallocHost((void **) &temp_peak_pos, sizeof(int));
	
	//--------> Benchmarking
	double total_time=0, MSD_time=0, SPDT_time=0, PF_time=0;
	
	//--------> Other
//	size_t free_mem,total_mem;
	char filename[200];

	// Calculate the total number of values
//	vals = ((size_t) nDMs)*((size_t) nTimesamples);


	// Creating events for synchronization on peak list; cudaEventDisableTiming increase 
	// performance and avoid synchronization issues;
	// Default state of Event: occured;
	// cudaEventSynchronize -- blocks host until stream completes all outstanding calls;
	// cudaStreamWaitEvent -- blocks stream until event occurs, blocks launches after this call, does not block host;
	cudaEvent_t waitEvent[NUM_STREAMS];
	for (int s=0; s < NUM_STREAMS; s++)
		checkCudaErrors(cudaEventCreateWithFlags(&waitEvent[s], cudaEventDisableTiming));	
	int second_stream = (stream_id + 1 ) % NUM_STREAMS;
		
	int max_iteration;
	int t_BC_widths[10]={PD_MAXTAPS,16,16,16,8,8,8,8,8,8};
	std::vector<int> BC_widths(t_BC_widths,t_BC_widths+sizeof(t_BC_widths)/sizeof(int));
	std::vector<PulseDetection_plan> PD_plan;

	//---------------------------------------------------------------------------
	//----------> GPU part
	printf("\n----------> GPU analysis part\n");
	printf("  Dimensions nDMs:%d; nTimesamples:%d; inBin:%d; outBin:%d; maxshift:%d; \n", ndms[i], t_processed, inBin, outBin, *maxshift);
//	GpuTimer total_timer, timer;
//	total_timer.Start();
	
	
//	float *d_MSD;
//	checkCudaErrors(cudaGetLastError());
//	if ( cudaSuccess != cudaMalloc((void**) &d_MSD, sizeof(float)*3)) {printf("Allocation error!\n"); exit(201);}
	
	
//	cudaMemGetInfo(&free_mem,&total_mem);
//	printf("  Memory required by boxcar filters:%0.3f MB\n",(4.5*vals*sizeof(float) + 2*vals*sizeof(ushort))/(1024.0*1024) );
//	printf("  Memory available:%0.3f MB \n", ((float) free_mem)/(1024.0*1024.0) );

	
	std::vector<int> DM_list;
//	unsigned long int max_timesamples=(free_mem*0.95)/(5.5*sizeof(float) + 2*sizeof(ushort));
	int DMs_per_cycle = maxTimeSamples/nTimesamples;
//	int DMs_per_cycle = max_timesamples/nTimesamples;
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
	
	
	max_iteration = Get_max_iteration(max_boxcar_width/inBin, &BC_widths, &max_width_performed);
	printf("  Selected iteration:%d; maximum boxcar width requested:%d; maximum boxcar width performed:%d;\n", max_iteration, max_boxcar_width/inBin, max_width_performed);
	Create_PD_plan(&PD_plan, &BC_widths, 1, nTimesamples);
	std::vector<int> h_boxcar_widths;
	Create_list_of_boxcar_widths(&h_boxcar_widths, &BC_widths);
	
	
	//-------------------------------------------------------------------------
	//---------> Comparison between interpolated values and computed values
	#ifdef MSD_BOXCAR_TEST
		MSD_plane_profile_boxcars(output_buffer, nTimesamples, nDMs, &h_boxcar_widths, OR_sigma_multiplier, dm_low[i], dm_high[i], tstart, streams);
	#endif
	//---------> Comparison between interpolated values and computed values
	//-------------------------------------------------------------------------
	
	
	//-------------------------------------------------------------------------
	//------------ Using MSD_plane_profile
	size_t MSD_profile_size_in_bytes, MSD_DIT_profile_size_in_bytes, workarea_size_in_bytes;
	Get_MSD_plane_profile_memory_requirements(&MSD_profile_size_in_bytes, &MSD_DIT_profile_size_in_bytes, &workarea_size_in_bytes, nTimesamples, nDMs, &h_boxcar_widths);
	double dit_time, MSD_only_time;
//	float *d_MSD_interpolated;
	float *d_MSD_DIT = NULL;
//	float *temporary_workarea;
//	cudaMalloc((void **) &d_MSD_interpolated, MSD_profile_size_in_bytes);
//	cudaMalloc((void **) &temporary_workarea, workarea_size_in_bytes);
       	printf("\tSize MSD: %zu \tSize workarea: %zu\n",  MSD_profile_size_in_bytes, workarea_size_in_bytes);

	checkCudaErrors(cudaGetLastError());
     
//	MSD_plane_profile(d_MSD_interpolated, output_buffer, d_MSD_DIT, temporary_workarea, false, nTimesamples, nDMs, &h_boxcar_widths, tstart, dm_low[i], dm_high[i], OR_sigma_multiplier, enable_sps_baselinenoise, false, &MSD_time, &dit_time, &MSD_only_time, streams);
	MSD_plane_profile(d_MSD_interpolated, output_buffer, d_MSD_DIT, h_MSD_DIT, h_MSD_interpolated, MSD_workarea, false, nTimesamples, nDMs, &h_boxcar_widths, tstart, dm_low[i], dm_high[i], OR_sigma_multiplier, enable_sps_baselinenoise, false, &MSD_time, &dit_time, &MSD_only_time, streams);
//	MSD_plane_profile(d_MSD_interpolated, output_buffer, d_MSD_DIT, MSD_workarea, false, nTimesamples, nDMs, &h_boxcar_widths, tstart, dm_low[i], dm_high[i], OR_sigma_multiplier, enable_sps_baselinenoise, false, &MSD_time, &dit_time, &MSD_only_time, streams);


	checkCudaErrors(cudaGetLastError());

//	cudaStreamSynchronize(streams);

      #ifdef GPU_PARTIAL_TIMER
      	printf("    MSD time: Total: %f ms; DIT: %f ms; MSD: %f ms; Size MSD: %zu \n", MSD_time, dit_time, MSD_only_time, MSD_profile_size_in_bytes);
      #endif
      
//	cudaFree(temporary_workarea);
      //------------ Using MSD_plane_profile
      //-------------------------------------------------------------------------	
//     	cudaMemsetAsync((void*) MSD_workarea, 0, sizeof(float)*5.5*maxTimeSamples, streams);

    
	if(DM_list.size()>0){
		DMs_per_cycle = DM_list[0];
//	     	cudaMemsetAsync((void*) MSD_workarea, 0, sizeof(float)*5.5*maxTimeSamples, streams);
	
		int mem_position=0;
		float *d_peak_list;
		d_peak_list = &MSD_workarea[mem_position];
//		d_peak_list = MSD_workarea;
//		if ( cudaSuccess != cudaMalloc((void**) &d_peak_list, sizeof(float)*DMs_per_cycle*nTimesamples)) printf("Allocation error! peaks\n");
		
		float *d_decimated;
		mem_position += nTimesamples*DMs_per_cycle;
		d_decimated = &MSD_workarea[mem_position];
//		if ( cudaSuccess != cudaMalloc((void **) &d_decimated,  sizeof(float)*(((DMs_per_cycle*nTimesamples)/2)+PD_MAXTAPS) )) printf("Allocation error! dedispered\n");
		
		float *d_boxcar_values;
		mem_position += (DMs_per_cycle*nTimesamples)/2+PD_MAXTAPS;
		d_boxcar_values = &MSD_workarea[mem_position];
//		if ( cudaSuccess != cudaMalloc((void **) &d_boxcar_values,  sizeof(float)*DMs_per_cycle*nTimesamples)) printf("Allocation error! boxcars\n");
		
		float *d_output_SNR;
		mem_position += DMs_per_cycle*nTimesamples;
		d_output_SNR = &MSD_workarea[mem_position];
//		if ( cudaSuccess != cudaMalloc((void **) &d_output_SNR, sizeof(float)*2*DMs_per_cycle*nTimesamples)) printf("Allocation error! SNR\n");
		
//		ushort *d_output_taps;
//		if ( cudaSuccess != cudaMalloc((void **) &d_output_taps, sizeof(ushort)*2*DMs_per_cycle*nTimesamples)) printf("Allocation error! taps\n");
		
//		int *gmem_peak_pos;
//		cudaMalloc((void**) &gmem_peak_pos, 1*sizeof(int));
		cudaMemsetAsync((void*) gmem_peak_pos, 0, sizeof(int), streams);
		
		DM_shift = 0;
		for(int f=0; f<DM_list.size(); f++) {
			//-------------- SPDT
//			timer.Start();
			SPDT_search_long_MSD_plane(&output_buffer[DM_shift*nTimesamples], d_boxcar_values, d_decimated, d_output_SNR, d_output_taps, d_MSD_interpolated, &PD_plan, max_iteration, nTimesamples, DM_list[f], streams);
//			timer.Stop();
//			SPDT_time += timer.Elapsed();
			#ifdef GPU_PARTIAL_TIMER
			printf("    SPDT took:%f ms\n", timer.Elapsed());
			#endif
			//-------------- SPDT
//			
			checkCudaErrors(cudaGetLastError());

			
			#ifdef GPU_ANALYSIS_DEBUG
			printf("    BC_shift:%d; DMs_per_cycle:%d; f*DMs_per_cycle:%d; max_iteration:%d;\n", DM_shift*nTimesamples, DM_list[f], DM_shift, max_iteration);
			#endif
	
			if(candidate_algorithm==1){
				//-------------- Thresholding
//				printf("\n\n Threshold ... \n");
//				timer.Start();
				THRESHOLD(d_output_SNR, d_output_taps, d_peak_list, gmem_peak_pos, cutoff, DM_list[f], nTimesamples, DM_shift, &PD_plan, max_iteration, local_max_list_size, streams);
//				cudaDeviceSynchronize();
//				timer.Stop();
//				PF_time += timer.Elapsed();
				#ifdef GPU_PARTIAL_TIMER
				printf("    Thresholding took:%f ms\n", timer.Elapsed());
				#endif
				//-------------- Thresholding
			}
			else {
				//-------------- Peak finding
//				printf("\n\nPeak finding ....\n");
//				timer.StartWithStream(streams);
				PEAK_FIND(d_output_SNR, d_output_taps, d_peak_list, DM_list[f], nTimesamples, cutoff, local_max_list_size, gmem_peak_pos, DM_shift, &PD_plan, max_iteration, streams);
//				timer.StopWithStream(streams);
//				PF_time = timer.ElapsedWithStream(streams);
				#ifdef GPU_PARTIAL_TIMER
				printf("    Peak finding took:%f ms\n", timer.Elapsed());
				#endif
				//-------------- Peak finding
			}
//			cudaEventRecord(waitEvent[stream_id], streams);
//			checkCudaErrors(cudaStreamWaitEvent(streams, waitEvent[second_stream],0));
//			checkCudaErrors(cudaMemcpy(temp_peak_pos, gmem_peak_pos, sizeof(int), cudaMemcpyDeviceToHost));
	//		checkCudaErrors(cudaMemcpyAsync(temp_peak_pos, gmem_peak_pos, sizeof(int), cudaMemcpyDeviceToHost, streams));
//			checkCudaErrors(cudaStreamAddCallback(streams,MyCallback,(void*) stream_id,0));
//			cudaEventRecord(waitEvent[stream_id], streams);
//			checkCudaErrors(cudaGetLastError());
//			bool process_stop = false;
//			while (process_stop == false) {process_stop = cudaEventQuery(waitEvent[stream_id]) == cudaSuccess;}
			checkCudaErrors(cudaStreamSynchronize(streams));
//			checkCudaErrors(cudaMemPrefetchAsync(gmem_peak_pos,2*sizeof(int),cudaCpuDeviceId,streams));
//			checkCudaErrors(cudaEventSynchronize(waitEvent[stream_id]));
//			checkCudaErrors(cudaStreamWaitEvent(streams, waitEvent[stream_id],0));
//			checkCudaErrors(cudaStreamAddCallback(streams,MyCallback,(void*) stream_id,0));
//			checkCudaErrors(cudaMemcpyAsync(temp_peak_pos, gmem_peak_pos, sizeof(int), cudaMemcpyDeviceToHost, streams));
//			checkCudaErrors(cudaMemcpyAsync(&temp_peak_pos, gmem_peak_pos, sizeof(int), cudaMemcpyDeviceToHost, streams));
//			#ifdef GPU_ANALYSIS_DEBUG
//			*gmem_peak_pos = 1;
			printf("    temp_peak_pos:%d; host_pos:%zu; max:%zu; local_max:%d; stream_id:%i;\n", (*gmem_peak_pos), (*peak_pos), max_peak_size, local_max_list_size, stream_id);
//			printf("\nTemp gmem_peak pointer: \t%p\n", gmem_peak_pos);
//			#endif
			if( (*gmem_peak_pos)>=local_max_list_size ) {
				printf("    Maximum list size reached! Increase list size or increase sigma cutoff.\n");
				*temp_peak_pos=local_max_list_size;
			}
			if( ((*peak_pos) + (*gmem_peak_pos))<max_peak_size){
				checkCudaErrors(cudaMemcpyAsync(&h_peak_list[(*peak_pos)*4], d_peak_list, (*gmem_peak_pos)*4*sizeof(float), cudaMemcpyDeviceToHost, streams));
				*peak_pos = (*peak_pos) + (*gmem_peak_pos);
			}
			else printf("Error peak list is too small!\n");

			DM_shift = DM_shift + DM_list[f];
			cudaMemsetAsync((void*) gmem_peak_pos, 0, sizeof(int), streams);
		} //dm_list

		//------------------------> Output
		#pragma omp parallel for
		for (int count = 0; count < (*peak_pos); count++){
			h_peak_list[4*count]     = h_peak_list[4*count]*dm_step[i] + dm_low[i];
			h_peak_list[4*count + 1] = h_peak_list[4*count + 1]*tsamp + tstart;
			h_peak_list[4*count + 2] = h_peak_list[4*count + 2];
			h_peak_list[4*count + 3] = h_peak_list[4*count + 3]*inBin;
		}
        
		FILE *fp_out;
      	
		if(candidate_algorithm==1){
			if((*peak_pos)>0){
				sprintf(filename, "analysed-t_%03.2f-dm_%.2f-%.2f.dat", tstart, dm_low[i], dm_high[i]);
				if (( fp_out = fopen(filename, "wb") ) == NULL)	{
					fprintf(stderr, "Error opening output file!\n");
					exit(0);
				}
				fwrite(h_peak_list, (*peak_pos)*sizeof(float), 4, fp_out);
				fclose(fp_out);
			}
		}
		else {
			if((*peak_pos)>0){
				sprintf(filename, "peak_analysed-t_%3.2f-dm_%.2f-%.2f.dat", tstart, dm_low[i], dm_high[i]);
				if (( fp_out = fopen(filename, "wb") ) == NULL)	{
					fprintf(stderr, "Error opening output file!\n");
					exit(0);
				}
				fwrite(h_peak_list, (*peak_pos)*sizeof(float), 4, fp_out);
				fclose(fp_out);
			}
		}
    	//------------------------> Output
    	
//		cudaFree(d_peak_list);
//		cudaFree(d_boxcar_values);
//		cudaFree(d_decimated);
//		cudaFree(d_output_SNR);
//		cudaFree(d_output_taps);
//		cudaFree(gmem_peak_pos);
//		cudaFree(d_MSD_DIT);
//		cudaFree(d_MSD_interpolated);
      

	}
	else printf("Error not enough memory to search for pulses\n");

//	total_timer.Stop();
//	total_time = total_timer.Elapsed();
	#ifdef GPU_TIMER
	printf("\n  TOTAL TIME OF SPS:%f ms\n", total_time);
	printf("  MSD_time: %f ms; SPDT time: %f ms; Candidate selection time: %f ms;\n", MSD_time, SPDT_time, PF_time);
	printf("----------<\n\n");
	#endif

//	for (int s=0; s < NUM_STREAMS; s++)
//		cudaEventDestroy(waitEvent[s]);
//	cudaFreeHost(&temp_peak_pos);
	//----------> GPU part
	//---------------------------------------------------------------------------
	
}


