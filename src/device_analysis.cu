//#define GPU_ANALYSIS_DEBUG
//#define MSD_BOXCAR_TEST
//#define GPU_PARTIAL_TIMER
#define GPU_TIMER

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include "params.hpp"

#include "aa_device_BC_plan.hpp"
#include "device_peak_find.hpp"
#include "device_MSD_plane_profile.hpp"
#include "device_SPS_long.hpp"
#include "device_threshold.hpp"

#include "aa_gpu_timer.hpp"

#include "aa_device_analysis.hpp"

//TODO:
// Make BC_plan for arbitrary long pulses, by reusing last element in the plane

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


  void analysis_GPU(unsigned int *h_peak_list_DM, unsigned int *h_peak_list_TS, float *h_peak_list_SNR, unsigned int *h_peak_list_BW, size_t *peak_pos, size_t max_peak_size, int i, float tstart, int t_processed, int inBin, int *maxshift, int max_ndms, int const*const ndms, float cutoff, float OR_sigma_multiplier, float max_boxcar_width_in_sec, float *output_buffer, float *dm_low, float *dm_high, float *dm_step, float tsamp, int candidate_algorithm, float *d_MSD_workarea, unsigned short *d_output_taps, float *d_MSD_interpolated, unsigned long int maxTimeSamples, int enable_sps_baselinenoise, const bool dump_to_disk, const bool dump_to_user, analysis_output &output){
    //--------> Task
    int max_boxcar_width = (int) (max_boxcar_width_in_sec/tsamp);
    int max_width_performed=0;
    int nTimesamples = t_processed;
    int nDMs = ndms[i];
    int temp_peak_pos;
	
    //--------> Benchmarking
    double total_time=0, MSD_time=0, SPDT_time=0, PF_time=0;
	
    //--------> Other
    char filename[200];
	
    //---------------------------------------------------------------------------
    //----------> GPU part
    printf("\n----------> GPU analysis part\n");
    printf("  Dimensions nDMs:%d; nTimesamples:%d; inBin:%d; maxshift:%d; \n", ndms[i], t_processed, inBin, *maxshift);
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

    MSD_plane_profile(d_MSD_interpolated, output_buffer, d_MSD_DIT, d_MSD_workarea, false, nTimesamples, nDMs, &h_boxcar_widths, tstart, dm_low[i], dm_high[i], OR_sigma_multiplier, enable_sps_baselinenoise, false, &MSD_time, &dit_time, &MSD_only_time);
	
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
		
      DM_shift = 0;
      int DM_list_size = (int)DM_list.size();
      for(int f=0; f<DM_list_size; f++) {
	//-------------- SPDT
	timer.Start();
	SPDT_search_long_MSD_plane(&output_buffer[DM_shift*nTimesamples], d_boxcar_values, d_decimated, d_output_SNR, d_output_taps, d_MSD_interpolated, &PD_plan, max_iteration, nTimesamples, DM_list[f]);
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
			
	if(candidate_algorithm==1){
	  //-------------- Thresholding
	  timer.Start();
	  SPDT_threshold(d_output_SNR, d_output_taps, d_peak_list_DM, d_peak_list_TS, d_peak_list_SNR, d_peak_list_BW, gmem_peak_pos, cutoff, DM_list[f], nTimesamples, DM_shift, &PD_plan, max_iteration, local_max_list_size);
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
	  SPDT_peak_find(d_output_SNR, d_output_taps, d_peak_list_DM, d_peak_list_TS, d_peak_list_SNR, d_peak_list_BW, DM_list[f], nTimesamples, cutoff, local_max_list_size, gmem_peak_pos, DM_shift, &PD_plan, max_iteration);
	  timer.Stop();
	  PF_time = timer.Elapsed();
#ifdef GPU_PARTIAL_TIMER
	  printf("    Peak finding took:%f ms\n", timer.Elapsed());
#endif
	  //-------------- Peak finding
	}
			
	checkCudaErrors(cudaGetLastError());
			
	checkCudaErrors(cudaMemcpy(&temp_peak_pos, gmem_peak_pos, sizeof(int), cudaMemcpyDeviceToHost));
#ifdef GPU_ANALYSIS_DEBUG
	printf("    temp_peak_pos:%d; host_pos:%zu; max:%zu; local_max:%d;\n", temp_peak_pos, (*peak_pos), max_peak_size, local_max_list_size);
#endif
	if( temp_peak_pos>=local_max_list_size ) {
	  printf("    Maximum list size reached! Increase list size or increase sigma cutoff.\n");
	  temp_peak_pos=local_max_list_size;
	}
	if( ((*peak_pos) + temp_peak_pos)<max_peak_size){
	  checkCudaErrors(cudaMemcpy(&h_peak_list_DM[(*peak_pos)],  d_peak_list_DM,  temp_peak_pos*sizeof(unsigned int), cudaMemcpyDeviceToHost));
	  checkCudaErrors(cudaMemcpy(&h_peak_list_TS[(*peak_pos)],  d_peak_list_TS,  temp_peak_pos*sizeof(unsigned int), cudaMemcpyDeviceToHost));
	  checkCudaErrors(cudaMemcpy(&h_peak_list_SNR[(*peak_pos)], d_peak_list_SNR, temp_peak_pos*sizeof(float), cudaMemcpyDeviceToHost));
	  checkCudaErrors(cudaMemcpy(&h_peak_list_BW[(*peak_pos)],  d_peak_list_BW,  temp_peak_pos*sizeof(unsigned int), cudaMemcpyDeviceToHost));
	  *peak_pos = (*peak_pos) + temp_peak_pos;
	}
	else printf("Error peak list is too small!\n");

	DM_shift = DM_shift + DM_list[f];
	cudaMemset((void*) gmem_peak_pos, 0, sizeof(int));
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
	    output.data.insert(output.data.end(), &h_peak_list[0], &h_peak_list[4 * (*peak_pos)]);
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
	    output.data.insert(output.data.end(), &h_peak_list[0], &h_peak_list[4 * (*peak_pos)]);
	  }
	}
      }
      delete[] h_peak_list;
      //------------------------> Output

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

} //namespace astroaccelerate
