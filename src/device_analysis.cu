//#define GPU_ANALYSIS_DEBUG
//#define MSD_BOXCAR_TEST
//#define GPU_PARTIAL_TIMER
#define GPU_TIMER

#include <iostream>
#include <tuple>
#include <vector>

#include <stdio.h>
#include <stdlib.h>

#include "params.hpp"

#include "device_BC_plan.hpp"
#include "device_peak_find.hpp"
#include "device_MSD_plane_profile.hpp"
#include "device_SPS_long.hpp"
#include "device_SPS_plan"
#include "device_threshold.hpp"
#include "device_MSD_parameters.hpp"

#include "gpu_timer.hpp"

void analysis_GPU(bool verbose, float* d_SPS_input, float *h_candidate_list, size_t &number_candidates, SPS_Plan &spsplan){
	// Definition of some local variables
	float local_tsamp  = spsplan.GetCurrentSamplingTime(); // SPS_data.sampling_time*SPS_data.inBin; // corrected sampling time
	size_t nTimesamples = spsplan.GetCurrentTimeSamples();
	size_t nDMs         = spsplan.GetNumberDMs();
	size_t max_candidates = spsplan.GetMaxCandidates();
	if(verbose) {
		std::cout << "----------> Single Pulse GPU Analysis" << std::endl;
		std::cout << "Dimensions: nTimesamples: " << nTimesamples
					<< ", nDMs: " << nDMS
					<< ", original sampling time: " << spsplan.GetOriginalSamplingTime()
					<< ", binning factor: " << spsplan.GetCurrentBinningFactor()
					<< ", binned sampling time: " << local_tsamp << endl;
	}
	
	//--------> Definition of SPDT boxcar plan
	int max_desired_boxcar_width = spsplan.GetCurrentMaxBoxcarWidth();
	int max_width_performed = 0
	int max_iteration = 0;
	
	MSD_Parameters MSD_params = spsplan.GetMSParameters();

	std::tuple<float, float, float> dm_limits = spsplan.GetDMLimits();
	// Old version
	//int t_BC_widths[10]={PD_MAXTAPS,16,16,16,8,8,8,8,8,8};
	//std::vector<int> BC_widths(t_BC_widths,t_BC_widths+sizeof(t_BC_widths)/sizeof(int));
	//std::vector<PulseDetection_plan> PD_plan;
	//Create_PD_plan(&PD_plan, &BC_widths, nTimesamples); //PD_plan is independent on maximum boxcar width. which is wrong?
	//max_iteration = Get_max_iteration(max_desired_boxcar_width, &BC_widths, &max_width_performed);
	//std::vector<int> h_boxcar_widths;
	//Create_list_of_boxcar_widths(&h_boxcar_widths, &BC_widths, max_width_performed);
	
	
	//New version
	// std::vector<PulseDetection_plan> PD_plan;
	// max_iteration = SPS_params->get_max_iteration(&max_width_performed, max_desired_boxcar_width);
	// SPS_params->Create_PD_plan(&PD_plan, &max_width_performed, max_desired_boxcar_width, nTimesamples);
	// if(verbose) 
	// 	printf("  Selected iteration:%d; maximum boxcar width requested:%d; maximum boxcar width performed:%d;\n", max_iteration, max_desired_boxcar_width, max_width_performed);
	// std::vector<int> h_boxcar_widths;
	// SPS_params->Create_list_of_boxcar_widths(&h_boxcar_widths, max_width_performed);
	
	// NOTE: This is the newest version
	//spsplan.CreateSPSPlan();
	max_iteration = spsplan.GetMaxIteration();
	std::vector<int> h_boxcar_widths = spsplan.GetListOfBoxcars();
	//spsplan.UpdateSPSPlan(max_width_performed, max_desired_boxcar_width, nTimesamples)

	if (verbose) {
		spsplan.PrintSPSPlan();
	}
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
	if (verbose) {
		std::cout << "\tMemory required for boxcar filters: " << (4.5 * nTimesamples * nDMs * sizeof(float) + 2 * nTimesamples * nDMs * sizeof(ushort)) / (1024.0*1024)) << "MiB" << std::endl;
		std::cout << "\tMemory available: " << ((float)free_mem)/(1024.0*1024.0)) << "MiB" << std::endl;
	}
	
	//-------------------------------------------------------------------------
	//---------> Comparison between interpolated values and computed values
	#ifdef MSD_BOXCAR_TEST
		MSD_plane_profile_boxcars(d_SPS_input, nTimesamples, nDMs, &h_boxcar_widths, MSD_params.OR_sigma_multiplier, std::get<0>(dm_limits), std::get<1>(dm_limits), spsplan.GetCurrentStartTime());
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
	
	MSD_plane_profile(d_MSD_interpolated, d_SPS_input, d_MSD_DIT, temporary_workarea, false, nTimesamples, nDMs, &h_boxcar_widths, spsplan.GetCurrentStartTime(), std::get<0>(dm_limits), std::get<1>(dm_limits), MSD_params.OR_sigma_multiplier, MSD_params.enable_outlier_rejection, false, &MSD_time, &dit_time, &MSD_only_time);
	
	#ifdef GPU_PARTIAL_TIMER
		std::cout << "\t\tMSD time: Total: " << MSD_time << "ms"
					<< ", DIT: " << dit_time << "ms"
					<< ", MSD: " << MSD_only_time << "ms" << std::endl;
	#endif
	
	cudaFree(temporary_workarea);
	//------------ Using MSD_plane_profile
	//-------------------------------------------------------------------------	
	
	
	//-------------------------------------------------------------------------
	//------------ Splitting input data into chunks
	std::vector<int> DM_list;
	// TODO: What is logic behind this equation
	unsigned long int max_timesamples=(free_mem*0.95)/(5.5*sizeof(float) + 2*sizeof(ushort));
	int DMs_per_cycle = max_timesamples/nTimesamples;
	int nRepeats, nRest, DM_shift, itemp, local_max_list_size;//BC_shift,
	
	itemp = (int) (DMs_per_cycle/THR_WARPS_PER_BLOCK);
	DMs_per_cycle = itemp*THR_WARPS_PER_BLOCK;
	
	nRepeats = nDMs/DMs_per_cycle;
	nRest = nDMs - nRepeats*DMs_per_cycle;
	local_max_list_size = (DMs_per_cycle*nTimesamples)/4;
	
	for (int irepeat = 0; irepeat < nRepeats; ++irepeat) {
		DM_list.push_back(DMs_per_cycle);
	}

	if (nRest > 0) {
		DM_list.push_back(nRest);
	}

	if( (int) DM_list.size() > 1 ) {{
		std::cout << "\tSPS will run " << DM_list.size() << 
					<< " batches each containing " << DMs_per_cycle << " DM trials"
					<< ", remainder  of " << nRest << " DM trials." << std::endl;
	} else { 
		std::cout << "\tSPS will run a single batch containing " << nRest << " DM trials." << std::endl;
	}


	//------------ Splitting input data into chunks
	//-------------------------------------------------------------------------	
	
	if (DM_list.size() > 0) {

		DMs_per_cycle = DM_list[0];
		// TODO: We should really quit on any of these error conditions
		float *d_peak_list;
		if ( cudaSuccess != cudaMalloc((void**) &d_peak_list, sizeof(float) * DMs_per_cycle * nTimesamples)) {
			std::cerr << "ERROR: Device memory allocation for peaks list" << std::endl;
		}
		float *d_decimated;
		if ( cudaSuccess != cudaMalloc((void **) &d_decimated,  sizeof(float)*(((DMs_per_cycle*nTimesamples)/2)+PD_MAXTAPS) )) {
			std::cerr << "ERROR: Device memory allocation for dedispered timeseries" << std::endl;
		}
		float *d_boxcar_values;
		if ( cudaSuccess != cudaMalloc((void **) &d_boxcar_values,  sizeof(float)*DMs_per_cycle*nTimesamples)) {
			std::cerr << "ERROR: Device memory allocation for boxcars list" << std::endl;
		}
		float *d_output_SNR;
		if ( cudaSuccess != cudaMalloc((void **) &d_output_SNR, sizeof(float)*2*DMs_per_cycle*nTimesamples)) {
			std::cerr << "ERROR: Device memory allocation for SNR list" << std::endl;
		}
		ushort *d_output_taps;
		if ( cudaSuccess != cudaMalloc((void **) &d_output_taps, sizeof(ushort)*2*DMs_per_cycle*nTimesamples)) {
			std::cerr << "ERROR: Device memory allocation taps" << std::endl;
		}
		int *gmem_peak_pos;
		cudaMalloc((void**) &gmem_peak_pos, 1*sizeof(int));
		cudaMemset((void*) gmem_peak_pos, 0, sizeof(int));
		
		DM_shift = 0;
		for (int f = 0; f < DM_list.size(); ++f) {
			//-------------- SPDT
			timer.Start();
				SPDT_search_long_MSD_plane(&d_SPS_input[DM_shift*nTimesamples], d_boxcar_values, d_decimated, d_output_SNR, d_output_taps, d_MSD_interpolated, spsplan, max_iteration, nTimesamples, DM_list[f]);
			timer.Stop();
			SPDT_time += timer.Elapsed();
			#ifdef GPU_PARTIAL_TIMER
				std::cout << "\t\tSPDT took: " << timer.Elapsed() << "ms" << std::endl;
			#endif
			//-------------- SPDT
			
			checkCudaErrors(cudaGetLastError());
			
			#ifdef GPU_ANALYSIS_DEBUG
			std::cout << "\t\tBC_shift: " << DM_shift * nTimesamples
						<< ", DMs_per_cycle: " << DM_list[f]
						<< ", f*DMs_per_cycle: " << DM_shift
						<< ", max_iteration: " << max_iteration << std::endl;
			#endif
			
			if (spsplan.GetSPSAlgorithm() == 1) {
				//-------------- Thresholding
				// NOTE: Was THRESHOLD
				timer.Start();
				SPDT_threshold(d_output_SNR, d_output_taps, d_peak_list, gmem_peak_pos, spsplan.GetSigmaCutoff(), DM_list[f], nTimesamples, DM_shift, spsplan, max_iteration, local_max_list_size, dm_limits, local_tsamp, spsplan.GetCurrentBinningFactor(), spsplan.GetCurrentStartTime());
				timer.Stop();
				PF_time += timer.Elapsed();
				#ifdef GPU_PARTIAL_TIMER
					std::cout << "\t\tThresholding took: " << timer.Elapsed() << "ms" << std::endl;
				#endif
				//-------------- Thresholding
			} else if (spsplan.GetSPSAlgorithm() == 0) {
				//-------------- Peak finding
				timer.Start();
				// Note: Was PEAK_FIND
				SPDS_peak_find(d_output_SNR, d_output_taps, d_peak_list, DM_list[f], nTimesamples, spsplan.GetSigmaCutoff(), local_max_list_size, gmem_peak_pos, DM_shift, spsplan, max_iteration, dm_limits, local_tsamp, spsplan.GetCurrentBinningFactor(), spsplan.GetCurrentStartTime());
				timer.Stop();
				PF_time = timer.Elapsed();
				#ifdef GPU_PARTIAL_TIMER
					std::cout << "\t\tPeak finding took: " << timer.Elapsed() << "ms" << std::endl;
				#endif
				//-------------- Peak finding
			}
			
			checkCudaErrors(cudaGetLastError());
			
			int temp_peak_pos = 0;
			checkCudaErrors(cudaMemcpy(&temp_peak_pos, gmem_peak_pos, sizeof(int), cudaMemcpyDeviceToHost));
			#ifdef GPU_ANALYSIS_DEBUG
				std::cout << "\t\tCandidates found:%d; " << temp_peak_pos
							<< ", total candidates for this chunk: " << number_candidates
							<< ", maximum candidates: " << max_candidates
							<< ", local maximum candidates: " << local_max_list_size << std::endl;
			#endif
			if( temp_peak_pos>=local_max_list_size ) {
				std::cout << "\t\tWARNING: Maximum list size reached! Not all candidates will be saved! You should increase detection threshold.\n");
				temp_peak_pos = local_max_list_size;
			}
			if( ((number_candidates) + temp_peak_pos)<max_candidates){
				checkCudaErrors(cudaMemcpy(&h_candidate_list[(number_candidates)*4], d_peak_list, temp_peak_pos*4*sizeof(float), cudaMemcpyDeviceToHost));
				number_candidates = (number_candidates) + temp_peak_pos;
			} else {
				std::cerr << "\t\tERROR: Not enough memory to store all candidates on the host!" << std::endl;
			}

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

	} else {
		// TODO: Is it GPU memory we are talking about
		std::cerr << "ERROR: Single pulse search was not run! Not enough memory to search for pulses" << std::endl;
	}
	total_timer.Stop();
	total_time = total_timer.Elapsed();
	#ifdef GPU_TIMER
	std::cout << std::endl << "TOTAL TIME OF SPS: " << total_time << " ms" << std::cout;
	std::cout << "\tMSD_time: " << MSD_time << "ms"
				<< ", SPDT time: " << SPDT_time << "ms"
				<< ", candidate selection time: " << PF_time << "ms" << std::cout;
	std::cout << "----------<" << std::endl << std::endl;
	#endif
	//----------> GPU part
	//---------------------------------------------------------------------------
	
}