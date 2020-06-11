//#define GPU_ANALYSIS_DEBUG
//#define MSD_BOXCAR_TEST
//#define GPU_PARTIAL_TIMER
#define GPU_TIMER

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <utility>
#include "aa_params.hpp"

#include "aa_log.hpp"
#include "aa_timelog.hpp"
#include "aa_device_BC_plan.hpp"
#include "aa_device_peak_find.hpp"
#include "aa_device_MSD_plane_profile.hpp"
#include "aa_device_SPS_long.hpp"
#include "aa_device_threshold.hpp"

#include "aa_gpu_timer.hpp"

#include "aa_device_analysis.hpp"

// \todo Make BC_plan for arbitrary long pulses, by reusing last element in the plane.
// \todo cudaMalloc((void**) &gmem_peak_pos, 1*sizeof(int)); has no corresponding cudaFree.
// \todo cudaMemset((void*) gmem_peak_pos, 0, sizeof(int));  has no corresponding cudaFree.

namespace astroaccelerate {

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


  /**
   * \brief Function that performs analysis on the dedispersed data.
   * \brief Users should not interact with this function. Instead they should use aa_analysis_plan and aa_analysis_strategy.
   * \details Argument int i is the current dm_range.
   */
  bool analysis_GPU(unsigned int *h_peak_list_DM, unsigned int *h_peak_list_TS, float *h_peak_list_SNR, unsigned int *h_peak_list_BW, size_t *peak_pos, size_t max_peak_size, int i, float tstart, int t_processed, int inBin, int *maxshift, int max_ndms, int const*const ndms, float cutoff, float OR_sigma_multiplier, float max_boxcar_width_in_sec, float *output_buffer, float *dm_low, float *dm_high, float *dm_step, float tsamp, int candidate_algorithm, float *d_MSD_workarea, unsigned short *d_output_taps, float *d_MSD_interpolated, unsigned long int maxTimeSamples, int enable_msd_baselinenoise, const bool dump_to_disk, const bool dump_to_user, analysis_output &output){
    //--------> Task
    int max_boxcar_width = (int) (max_boxcar_width_in_sec/tsamp);
    int max_width_performed=0;
    int nTimesamples = t_processed;
    int nDMs = ndms[i];
    int temp_peak_pos;

    TimeLog time_log;
	
    //--------> Benchmarking
    double total_time=0, MSD_time=0, SPDT_time=0, PF_time=0;
	
    //--------> Other
    char filename[200];
	
    //---------------------------------------------------------------------------
    //----------> GPU part
    printf("\n----------> GPU analysis part\n");
    printf("  Dimensions nDMs:%d; nTimesamples:%d; inBin:%d; maxshift:%d; candidate_algorithm:%d \n", ndms[i], t_processed, inBin, *maxshift, candidate_algorithm);
    aa_gpu_timer total_timer, timer;
    total_timer.Start();
	

    //--------> SPDT plan
    int max_iteration;
    int t_BC_widths[10]={PD_MAXTAPS,16,16,16,8,8,8,8,8,8};
    std::vector<int> BC_widths(t_BC_widths,t_BC_widths+sizeof(t_BC_widths)/sizeof(int));
    std::vector<PulseDetection_plan> PD_plan;	
    std::vector<int> DM_list;
    int DMs_per_cycle = maxTimeSamples/nTimesamples;
    int nRepeats, nRest, DM_shift, itemp, local_max_list_size;
	
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
    Create_list_of_boxcar_widths(&h_boxcar_widths, &BC_widths, max_width_performed);
	
	
    //-------------------------------------------------------------------------
    //---------> Comparison between interpolated values and computed values
#ifdef MSD_BOXCAR_TEST
    MSD_plane_profile_boxcars(output_buffer, nTimesamples, nDMs, &h_boxcar_widths, OR_sigma_multiplier, dm_low[i], dm_high[i], tstart);
#endif
    //---------> Comparison between interpolated values and computed values
    //-------------------------------------------------------------------------
	

	
    //-------------------------------------------------------------------------
    //------------ Using MSD_plane_profile
    double dit_time, MSD_only_time;
    float *d_MSD_DIT = NULL; //TODO: make d_MSD_DIT also allocated outside

    MSD_plane_profile(d_MSD_interpolated, output_buffer, d_MSD_DIT, d_MSD_workarea, false, nTimesamples, nDMs, &h_boxcar_widths, tstart, dm_low[i], dm_high[i], OR_sigma_multiplier, enable_msd_baselinenoise, false, &MSD_time, &dit_time, &MSD_only_time);
	
#ifdef GPU_PARTIAL_TIMER
    printf("    MSD time: Total: %f ms; DIT: %f ms; MSD: %f ms;\n", MSD_time, dit_time, MSD_only_time);
#endif
	
    //------------ Using MSD_plane_profile
    //-------------------------------------------------------------------------	
	
	
    if(DM_list.size()>0){
      DMs_per_cycle = DM_list[0];
		
      size_t mem_idx = 0;
		
      unsigned int *d_peak_list_DM;
      d_peak_list_DM  = (unsigned int*) &d_MSD_workarea[mem_idx];
      mem_idx += (nTimesamples*DMs_per_cycle)/4;

      unsigned int *d_peak_list_TS;
      d_peak_list_TS  = (unsigned int*) &d_MSD_workarea[mem_idx];
      mem_idx += (nTimesamples*DMs_per_cycle)/4;
		
      float *d_peak_list_SNR;
      d_peak_list_SNR = &d_MSD_workarea[mem_idx];
      mem_idx += (nTimesamples*DMs_per_cycle)/4;
		
      unsigned int *d_peak_list_BW;
      d_peak_list_BW = (unsigned int*) &d_MSD_workarea[mem_idx];
      mem_idx += (nTimesamples*DMs_per_cycle)/4;
		
      float *d_decimated;
      d_decimated = &d_MSD_workarea[mem_idx];
      mem_idx += (DMs_per_cycle*nTimesamples)/2 + PD_MAXTAPS;
		
      float *d_boxcar_values;
      d_boxcar_values = &d_MSD_workarea[mem_idx];
      mem_idx += DMs_per_cycle*nTimesamples;
		
      float *d_output_SNR;
      d_output_SNR = &d_MSD_workarea[mem_idx];
      mem_idx += DMs_per_cycle*nTimesamples;
		
      // check if we are in limits of allocated memory
		
      int *gmem_peak_pos;
      cudaMalloc((void**) &gmem_peak_pos, 1*sizeof(int));
      cudaMemset((void*) gmem_peak_pos, 0, sizeof(int));

      int *gmem_filteredPeak_pos;
      cudaMalloc((void**) &gmem_filteredPeak_pos, 1*sizeof(int));
      cudaMemset((void*) gmem_filteredPeak_pos, 0, sizeof(int));
		
      DM_shift = 0;
      int DM_list_size = (int)DM_list.size();
      for(int f=0; f<DM_list_size; f++) {
	//-------------- SPDT
	timer.Start();
	SPDT_search_long_MSD_plane(&output_buffer[(size_t)(DM_shift)*(size_t)(nTimesamples)], d_boxcar_values, d_decimated, d_output_SNR, d_output_taps, d_MSD_interpolated, &PD_plan, max_iteration, nTimesamples, DM_list[f]);
	timer.Stop();
	SPDT_time += timer.Elapsed();
	time_log.adding("SPD","SPDT",timer.Elapsed());
#ifdef GPU_PARTIAL_TIMER
	printf("    SPDT took:%f ms\n", timer.Elapsed());
#endif
	//-------------- SPDT
			
	//checkCudaErrors(cudaGetLastError());
			
#ifdef GPU_ANALYSIS_DEBUG
	printf("    BC_shift:%zu; DMs_per_cycle:%d; f*DMs_per_cycle:%d; max_iteration:%d;\n", (size_t)(DM_shift)*(size_t)(nTimesamples), DM_list[f], DM_shift, max_iteration);
#endif
			
	if(candidate_algorithm==1){
	  //-------------- Thresholding
	  timer.Start();
	  SPDT_threshold(d_output_SNR, d_output_taps, d_peak_list_DM, d_peak_list_TS, d_peak_list_SNR, d_peak_list_BW, gmem_peak_pos, cutoff, DM_list[f], nTimesamples, DM_shift, &PD_plan, max_iteration, local_max_list_size);
	  timer.Stop();
	  PF_time += timer.Elapsed();
	  time_log.adding("SPD", "Threshold", timer.Elapsed());
#ifdef GPU_PARTIAL_TIMER
	  printf("    Thresholding took:%f ms\n", timer.Elapsed());
#endif
	  //-------------- Thresholding
	}
	else if(candidate_algorithm==0) {
	  //-------------- Peak finding
	  timer.Start();
	  SPDT_peak_find(d_output_SNR, d_output_taps, d_peak_list_DM, d_peak_list_TS, d_peak_list_SNR, d_peak_list_BW, DM_list[f], nTimesamples, cutoff, local_max_list_size, gmem_peak_pos, DM_shift, &PD_plan, max_iteration);
	  timer.Stop();
	  PF_time = timer.Elapsed();
	  time_log.adding("SPD", "Peak_Find", timer.Elapsed());
#ifdef GPU_PARTIAL_TIMER
	  printf("    Peak finding took:%f ms\n", timer.Elapsed());
#endif
	  //-------------- Peak finding
	}
	else if(candidate_algorithm==2) { //peak filtering
		timer.Start();
		SPDT_peak_find_stencil_7x7(d_output_SNR, d_output_taps, d_peak_list_DM, d_peak_list_TS, d_peak_list_SNR, d_peak_list_BW, DM_list[f], nTimesamples, cutoff, local_max_list_size, gmem_peak_pos, DM_shift, &PD_plan, max_iteration);
		timer.Stop();
		PF_time = timer.Elapsed();
		time_log.adding("SPD", "Stencil_7x7", timer.Elapsed());		
#ifdef GPU_PARTIAL_TIMER
          printf("    Peak finding (stencil 7x7) took: %f ms\n", timer.Elapsed());
#endif
	}
			
	//checkCudaErrors(cudaGetLastError());
			
	cudaError_t e = cudaMemcpy(&temp_peak_pos, gmem_peak_pos, sizeof(int), cudaMemcpyDeviceToHost);

	if(e != cudaSuccess) {
	  LOG(log_level::error, "Could not cudaMemcpy in aa_device_analysis.cu -- temp_peak_pos (" + std::string(cudaGetErrorString(e)) + ")");
		exit(25);
	}
	
#ifdef GPU_ANALYSIS_DEBUG
	printf("    temp_peak_pos:%d; host_pos:%zu; max:%zu; local_max:%d;\n", temp_peak_pos, (*peak_pos), max_peak_size, local_max_list_size);
#endif
	if( temp_peak_pos>=local_max_list_size ) {
	  printf("    Maximum list size reached! Increase list size or increase sigma cutoff.\n");
	  temp_peak_pos=local_max_list_size;
	}
	if( ((*peak_pos) + temp_peak_pos)<max_peak_size){
	  cudaError_t e = cudaMemcpy(&h_peak_list_DM[(*peak_pos)],  d_peak_list_DM,  temp_peak_pos*sizeof(unsigned int), cudaMemcpyDeviceToHost);

	  if(e != cudaSuccess) {
	    LOG(log_level::error, "Could not cudaMemcpy in aa_device_analysis.cu -- peak_list_DM (" + std::string(cudaGetErrorString(e)) + ")");
	  }
	  
	  e = cudaMemcpy(&h_peak_list_TS[(*peak_pos)],  d_peak_list_TS,  temp_peak_pos*sizeof(unsigned int), cudaMemcpyDeviceToHost);

	  if(e != cudaSuccess) {
	    LOG(log_level::error, "Could not cudaMemcpy in aa_device_analysis.cu -- peak_list_TS (" + std::string(cudaGetErrorString(e)) + ")");
	  }
	  
	  e = cudaMemcpy(&h_peak_list_SNR[(*peak_pos)], d_peak_list_SNR, temp_peak_pos*sizeof(float), cudaMemcpyDeviceToHost);

	  if(e != cudaSuccess) {
	    LOG(log_level::error, "Could not cudaMemcpy in aa_device_analysis.cu -- peak_list_SNR (" + std::string(cudaGetErrorString(e)) + ")");
	  }
	  
	  e = cudaMemcpy(&h_peak_list_BW[(*peak_pos)],  d_peak_list_BW,  temp_peak_pos*sizeof(unsigned int), cudaMemcpyDeviceToHost);
	  
	  if(e != cudaSuccess) {
	    LOG(log_level::error, "Could not cudaMemcpy in aa_device_analysis.cu -- peak_list_BW (" + std::string(cudaGetErrorString(e)) + ")");
	  }
	  
	  *peak_pos = (*peak_pos) + temp_peak_pos;
	}
	else printf("Error peak list is too small!\n");

	DM_shift = DM_shift + DM_list[f];
	cudaMemset((void*) gmem_peak_pos, 0, sizeof(int));
      }

	if(candidate_algorithm==2) { //peak filtering
		//------------peak clustering from AA_experimental
		size_t d_output_SNR_size = DMs_per_cycle*nTimesamples;
		size_t d_peak_list_size = DMs_per_cycle*nTimesamples/4;
		unsigned int local_peak_pos = *peak_pos;
		unsigned int *d_peak_list_DM2;
		unsigned int *d_peak_list_TS2;
		unsigned int *d_peak_list_BW2;
		float *d_peak_list_SNR2;
		d_peak_list_DM2  = (unsigned int*) &d_output_SNR[0];
		d_peak_list_TS2  = (unsigned int*) &d_output_SNR[d_peak_list_size];
		d_peak_list_BW2  = (unsigned int*) &d_output_SNR[2*d_peak_list_size];
		d_peak_list_SNR2  = (float*) &d_output_SNR[3*d_peak_list_size];

		timer.Start();
		printf("Number of peaks: %d; which is %f MB; d_output_SNR_size is %f MB;\n", local_peak_pos, (local_peak_pos*4.0*4.0)/(1024.0*1024.0), (d_output_SNR_size*4.0)/(1024.0*1024.0));

		if(d_output_SNR_size > local_peak_pos){
			cudaMemset((void*) d_output_SNR, 0, d_output_SNR_size*sizeof(float));
			cudaMemset((void*) d_peak_list_DM, 0, sizeof(unsigned int)*d_peak_list_size);
			cudaMemset((void*) d_peak_list_TS, 0, sizeof(unsigned int)*d_peak_list_size);
			cudaMemset((void*) d_peak_list_BW, 0, sizeof(unsigned int)*d_peak_list_size);
			cudaMemset((void*) d_peak_list_SNR, 0, sizeof(float)*d_peak_list_size);

			cudaError_t e = cudaMemcpy(d_peak_list_DM2, h_peak_list_DM, sizeof(unsigned int)*local_peak_pos, cudaMemcpyHostToDevice);
			if (e != cudaSuccess){
				LOG(log_level::error, "Could not cudaMemcpy in d_peak_list_DM2 (" + std::string(cudaGetErrorString(e)) + ")");
			}
			e = cudaMemcpy(d_peak_list_TS2, h_peak_list_TS, sizeof(unsigned int)*local_peak_pos, cudaMemcpyHostToDevice);
			if (e != cudaSuccess){
				LOG(log_level::error, "Could not cudaMemcpy in d_peak_list_TS2 (" + std::string(cudaGetErrorString(e)) + ")");
			}		
			e = cudaMemcpy(d_peak_list_BW2, h_peak_list_BW, sizeof(unsigned int)*local_peak_pos, cudaMemcpyHostToDevice);
			if (e != cudaSuccess){
					LOG(log_level::error, "Could not cudaMemcpy in d_peak_list_BW2 (" + std::string(cudaGetErrorString(e)) + ")");
			}
			e = cudaMemcpy(d_peak_list_SNR2, h_peak_list_SNR, sizeof(float)*local_peak_pos, cudaMemcpyHostToDevice);
			if (e != cudaSuccess){
				LOG(log_level::error, "Could not cudaMemcpy in d_peak_list_SNR2 (" + std::string(cudaGetErrorString(e)) + ")");
			}
			
			int filter_size = (int)(PPF_SEARCH_RANGE_IN_MS*0.001/tsamp);
			call_gpu_Filter_peaks(d_peak_list_DM, d_peak_list_TS, d_peak_list_BW, d_peak_list_SNR, d_peak_list_DM2, d_peak_list_TS2, d_peak_list_BW2, d_peak_list_SNR2, local_peak_pos, filter_size, (int)d_peak_list_size, gmem_filteredPeak_pos);

			cudaMemcpy(&temp_peak_pos, gmem_filteredPeak_pos, sizeof(int), cudaMemcpyDeviceToHost);
			local_peak_pos = temp_peak_pos;

			cudaMemcpy(h_peak_list_DM, d_peak_list_DM, local_peak_pos*sizeof(unsigned int), cudaMemcpyDeviceToHost);
			cudaMemcpy(h_peak_list_TS, d_peak_list_TS, local_peak_pos*sizeof(unsigned int), cudaMemcpyDeviceToHost);
			cudaMemcpy(h_peak_list_BW, d_peak_list_BW, local_peak_pos*sizeof(unsigned int), cudaMemcpyDeviceToHost);
			cudaMemcpy(h_peak_list_SNR, d_peak_list_SNR, local_peak_pos*sizeof(float), cudaMemcpyDeviceToHost);
		} else {
			LOG(log_level::error, "Not enough memory for filtering.");
		}

		timer.Stop();
		time_log.adding("SPD", "Clustering", timer.Elapsed());
	#ifdef GPU_PARTIAL_TIMER
		printf("    Clustering took:%f ms\n", timer.Elapsed());
	#endif
		//------------------------------------------------
		*peak_pos = local_peak_pos;
	}
		
      //------------------------> Output
	float *h_peak_list;
	h_peak_list = new float[4*(*peak_pos)];
	int i_peak_pos = (int)(*peak_pos);

	for (int count = 0; count < i_peak_pos; count++){
		h_peak_list[4*count]     = ((double) h_peak_list_DM[count])*dm_step[i] + dm_low[i];
		h_peak_list[4*count + 1] = ((double) h_peak_list_TS[count])*tsamp + tstart;
		h_peak_list[4*count + 2] = ((double) h_peak_list_SNR[count]);
		h_peak_list[4*count + 3] = ((double) h_peak_list_BW[count])*inBin;
	}
        
      FILE *fp_out;
		
      if(candidate_algorithm==1){
	if((*peak_pos)>0){
	  if(dump_to_disk) {
	    sprintf(filename, "analysed-t_%.2f-dm_%.2f-%.2f.dat", tstart, dm_low[i], dm_high[i]);
	    if (( fp_out = fopen(filename, "wb") ) == NULL)	{
	      fprintf(stderr, "Error opening output file!\n");
	      exit(0);
	    }
	    fwrite(h_peak_list, (*peak_pos)*sizeof(float), 4, fp_out);
	    fclose(fp_out);
	  }

	  if(dump_to_user) {
	    output.dm_low  = dm_low [i];
	    output.dm_high = dm_high[i];
	    std::vector<analysis_pulse> pulses;
	    for(auto count = 0; count < i_peak_pos; count++) {
	      analysis_pulse tmp = {h_peak_list[4*count], h_peak_list[4*count + 1], h_peak_list[4*count + 2], h_peak_list[4*count + 3]};
	      pulses.push_back(std::move(tmp));
	    }
	    output.pulses = std::move(pulses);
	  }
	}
      }
      else {
	if((*peak_pos)>0){
	  if(dump_to_disk) {
	    sprintf(filename, "peak_analysed-t_%.2f-dm_%.2f-%.2f.dat", tstart, dm_low[i], dm_high[i]);
	    if (( fp_out = fopen(filename, "wb") ) == NULL)	{
	      fprintf(stderr, "Error opening output file!\n");
	      exit(0);
	    }
	    fwrite(h_peak_list, (*peak_pos)*sizeof(float), 4, fp_out);
	    fclose(fp_out);
	  }

	  if(dump_to_user) {
	    output.dm_low  = dm_low [i];
	    output.dm_high = dm_high[i];
	    std::vector<analysis_pulse> pulses;
            for(auto count = 0; count < i_peak_pos; count++) {
              analysis_pulse tmp = {h_peak_list[4*count], h_peak_list[4*count + 1], h_peak_list[4*count + 2], h_peak_list[4*count + 3]};
              pulses.push_back(std::move(tmp));
            }
            output.pulses = std::move(pulses);
	  }
	}
      }
      delete[] h_peak_list;
      //------------------------> Output
		cudaFree(gmem_peak_pos);
    }
    else printf("Error not enough memory to search for pulses\n");

    total_timer.Stop();
    total_time = total_timer.Elapsed();
    time_log.adding("SPD", "total", total_time);
    time_log.adding("SPD", "MSD", MSD_time);
	time_log.adding("Total", "total", total_time);
#ifdef GPU_TIMER
    printf("\n  TOTAL TIME OF SPS:%f ms\n", total_time);
    printf("  MSD_time: %f ms; SPDT time: %f ms; Candidate selection time: %f ms;\n", MSD_time, SPDT_time, PF_time);
    printf("----------<\n\n");
#endif

    //----------> GPU part
    //---------------------------------------------------------------------------
	
	return true;
  }

} //namespace astroaccelerate

