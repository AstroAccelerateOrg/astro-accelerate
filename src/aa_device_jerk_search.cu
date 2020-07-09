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

#include "aa_jerk_CandidateList.hpp"

// MSD
#include "aa_device_MSD_Configuration.hpp"
#include "aa_device_MSD.hpp"
#include "aa_device_MSD_plane_profile.hpp"

// Convolution
#include "aa_device_convolution.hpp"

// Peak find
#include "aa_device_peak_find.hpp"

#define VERBOSE 1

namespace astroaccelerate {

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
		

	void jerk_create_acc_filters(float2 *jerk_filters, aa_jerk_strategy *jerk_strategy){
		int nFilters_z_half  = jerk_strategy->nFilters_z_half();
		int nFilters_z       = jerk_strategy->nFilters_z(); // number of filters must account also for negative accelerations; -1 because z=0;
		int nFilters_w_half  = jerk_strategy->nFilters_w_half();
		int nFilters_w       = jerk_strategy->nFilters_w();
		int convolution_size = jerk_strategy->conv_size();
		int interbinned_samples = jerk_strategy->interbinned_samples();
		float z_step = jerk_strategy->z_search_step();
		float w_step = jerk_strategy->w_search_step();
		//printf("filter generation: nFilters_z_half=%d; nFilters_z=%d; nFilters_w_half=%d; nFilters_w=%d; Total=%d; \n", nFilters_z_half, nFilters_z, nFilters_w_half, nFilters_w, nFilters_z*nFilters_w);
		
		#ifdef EXPORT_FILTERS
		std::ofstream FILEOUT;
		char filename[200];
		#endif
		
		cufftComplex *tempfilter;
		int nFz=0;
		int nFw=0;
		for(int ws=-nFilters_w_half; ws<=nFilters_w_half; ws++){
			for(int zs=-nFilters_z_half; zs<=nFilters_z_half; zs++){
				double z = ((float) zs)*z_step;
				double w = ((float) ws)*w_step;			
				int halfwidth = presto_w_resp_halfwidth(z, w, jerk_strategy->high_precision());
				int filter_size = 2*halfwidth*interbinned_samples;
				int pos = (ws + nFilters_w_half)*nFilters_z*convolution_size + (zs + nFilters_z_half)*convolution_size;
				
				tempfilter = presto_gen_w_response(0.0, interbinned_samples, z, w, filter_size);
				
				#ifdef EXPORT_FILTERS
				sprintf(filename, "filter_z%d_w%d.dat", zs, ws);
				FILEOUT.open(filename);
				for(int c=0; c<filter_size; c++){
					FILEOUT << tempfilter[c].x*tempfilter[c].x + tempfilter[c].y*tempfilter[c].y << " " << tempfilter[c].x << " " << tempfilter[c].y << std::endl;
				}
				FILEOUT.close();
				#endif
				
				presto_place_complex_kernel(tempfilter, filter_size, &jerk_filters[pos], convolution_size);
				free(tempfilter);
				
				if(zs==0) nFw++;
			}
			nFz++;
		}
		
		//printf("nFz=%d; nFw=%d;\n",nFz, nFw);
		
		#ifdef EXPORT_FILTERS
		FILEOUT.open("jerk_filters.dat");
		for(int ws=-nFilters_w_half; ws<nFilters_w_half; ws++){
			for(int zs=-nFilters_z_half; zs<=nFilters_z_half; zs++){
				for(int c=0; c<convolution_size; c++){
					int pos = (zs+nFilters_z_half)*nFilters_w*convolution_size + (ws+nFilters_w_half)*convolution_size + c;
					FILEOUT << jerk_filters[pos].x << " " << jerk_filters[pos].y << std::endl;
				}
				FILEOUT << std::endl;
				FILEOUT << std::endl;
			}
		}
		FILEOUT.close();
		#endif
	}


	int jerk_search_from_ddtr_plan(float ***dedispersed_data, aa_jerk_strategy &jerk_strategy, float *dm_low, float *dm_step, const int *list_of_ndms, float sampling_time, int *inBin, int nRanges){

		
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
		aa_gpu_timer timer_total, timer_DM, timer;
		double time_total=0, time_per_range=0;
		timer_total.Start();
		cudaError_t cudaError;
		if(VERBOSE>0) aa_jerk_strategy::print_info(jerk_strategy);
		
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
			size_t DM_trial_samples = default_nTimesamples/inBin[active_range];
			size_t nDMs = list_of_ndms[active_range];
			
			jerk_strategy.recalculate(DM_trial_samples,nDMs);
			if(VERBOSE>2) aa_jerk_strategy::print_info(jerk_strategy);
			
			//---------> Allocation of the output Z-planes and candidates
			size_t ZW_plane_size              = jerk_strategy.output_size_z_plane();
			unsigned int max_nZWCandidates    = (ZW_plane_size/4);
			size_t single_ZW_plane_size_bytes = jerk_strategy.output_size_z_plane()*sizeof(float);
			std::vector<int> ZW_chunks        = jerk_strategy.ZW_chunks();
			size_t ZW_planes_size_bytes       = ZW_chunks[0]*single_ZW_plane_size_bytes;
			
			
			if(VERBOSE>3) printf("JERK SEARCH -> ZW single plane size: %zu elements = %f MB\n", single_ZW_plane_size_bytes, ((float) single_ZW_plane_size_bytes)/(1024.0*1024.0));
			if(VERBOSE>3) printf("JERK SEARCH -> ZW plane size: %zu elements = %f MB\n", ZW_planes_size_bytes, ((float) ZW_planes_size_bytes)/(1024.0*1024.0));
			
			
			if ( cudaSuccess != cudaMalloc((void **) &d_ZW_candidates, max_nZWCandidates*sizeof(float))) {
				printf("Cannot allocate GPU memory for ZW plane candidates!\n");
			}
			if ( cudaSuccess != cudaMalloc((void **) &d_ZW_planes, ZW_planes_size_bytes)) {
				printf("Cannot allocate GPU memory for ZW planes!\n");
			}
			
			//---------> cuFFT plan
			cufftHandle cuFFT_plan;
			cufftResult cuFFT_error;
			cuFFT_error = cufftPlan1d(&cuFFT_plan, jerk_strategy.nSamples_time_dom(), CUFFT_R2C, 1);
			if (cuFFT_error!=CUFFT_SUCCESS) {
				printf("ERROR while constructing cuFFT plan\n");
			}
			
			//---------> MSD configuration
			MSD_Configuration MSD_conf(jerk_strategy.output_size_one_DM(), jerk_strategy.nFilters_z(), 0, 0);
			if ( cudaSuccess != cudaMalloc((void **) &d_MSD_workarea, MSD_conf.nBlocks_total*MSD_PARTIAL_SIZE*sizeof(float))) {
				printf("Cannot allocate GPU memory for MSD workarea!\n");
			}
			

			for(size_t active_DM = 0; active_DM<nDMs; active_DM++){
				timer_DM.Start();
				
				//-------> Copy DM-trial
				cudaError = cudaMemcpy(d_DM_trial, dedispersed_data[active_range][active_DM], jerk_strategy.nSamples_time_dom()*sizeof(float), cudaMemcpyHostToDevice);
				if ( cudaError != cudaSuccess) {
					printf("ERROR while copying DM-trial to the device!\n");
					printf("Error %s\n", cudaGetErrorString(cudaError));
				}
				
				//-------> FFT of the DM-trial
				cufftExecR2C(cuFFT_plan, d_DM_trial, (cufftComplex *) d_DM_trial_ffted);
				
				//---------> Candidates
				int ZW_planes_shift = 0; //this is for filters
				//TODO: container for candidates!
				
				JERK_CandidateList allcandidates(jerk_strategy.nFilters_w());
				if(VERBOSE>3) printf("nSublists: %zu\n", allcandidates.getNumberOfSubLists());
				
				std::vector<int> ZW_chunks = jerk_strategy.ZW_chunks();
				for(int f=0; f<(int) ZW_chunks.size(); f++){
					int nZW_planes = ZW_chunks[f];
					
					// Convolution
					timer.Start();
					if(VERBOSE>3) printf("JERK DEBUG: Convolution: position of the filters = %d; nTimesamples = %zu; conv_size = %d; useful_part_size = %d; offset = %d; nSegments = %d; nFilters = %d;\n", jerk_strategy.nFilters_z()*ZW_planes_shift, jerk_strategy.output_size_one_DM(), jerk_strategy.conv_size(), jerk_strategy.useful_part_size(), jerk_strategy.filter_halfwidth(), jerk_strategy.nSegments(), jerk_strategy.nFilters_z()*nZW_planes);
					
					conv_OLS_customFFT(d_DM_trial_ffted, d_ZW_planes, &d_jerk_filters[jerk_strategy.nFilters_z()*ZW_planes_shift*jerk_strategy.conv_size()], jerk_strategy.nSamples_freq_dom(), jerk_strategy.conv_size(), jerk_strategy.useful_part_size(), jerk_strategy.filter_halfwidth(), jerk_strategy.nSegments(), jerk_strategy.nFilters_z()*nZW_planes, 1.0f);
					timer.Stop(); 
					if(VERBOSE>2) printf("JERK SEARCH -> Convolution took: %g ms\n", timer.Elapsed());
					Convolution_time += timer.Elapsed();
					
					
					for(int zp=0; zp<nZW_planes; zp++){
						unsigned int nCandidates = 0;
						char filename[200];
						std::ofstream FILEOUT;
						float w = ((float) (ZW_planes_shift + zp))*jerk_strategy.w_search_step() - jerk_strategy.w_max_search_limit();
						
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
						Find_MSD(d_MSD, &d_ZW_planes[pos], d_MSD_workarea, &MSD_conf, jerk_strategy.OR_sigma_cutoff(), jerk_strategy.MSD_outlier_rejection());
						timer.Stop(); 
						if(VERBOSE>2) printf("JERK SEARCH -> MSD took: %g ms\n", timer.Elapsed());
						MSD_time += timer.Elapsed();
						
						//------------> Candidate selection
						timer.Start();
						
						cudaMemset((void*) gmem_peak_pos, 0, sizeof(unsigned int));
						PEAK_FIND_FOR_FDAS(&d_ZW_planes[pos], d_ZW_candidates, d_MSD, jerk_strategy.nFilters_z(), jerk_strategy.output_size_one_DM(), jerk_strategy.CS_sigma_threshold(), max_nZWCandidates, gmem_peak_pos, w);
						
						timer.Stop();
						if(VERBOSE>2) printf("JERK SEARCH -> Candidate selection took: %g ms\n", timer.Elapsed());
						Candidate_time += timer.Elapsed();
						if(VERBOSE>3) printf("JERK SEARCH -> W Coordinate: %e [%d;%d] \n", w, ZW_planes_shift, zp);
					
						//------------> Export to host
						if ( cudaSuccess != cudaMemcpy(&nCandidates, gmem_peak_pos, sizeof(unsigned int), cudaMemcpyDeviceToHost)) {
							printf("ERROR: Cannot copy the number of candidates to the host!\n");
						}
						
						//-------------------- DIRECT CANDIDATE OUTPUT
						if(nCandidates>0){
							float *h_candidates;
							h_candidates = new float [nCandidates*4];
							cudaMemcpy(h_candidates, d_ZW_candidates, nCandidates*4*sizeof(float), cudaMemcpyDeviceToHost);
							
							for(size_t f=0; f<nCandidates; f++){
								h_candidates[4*f] = h_candidates[4*f]*jerk_strategy.z_search_step() - jerk_strategy.z_max_search_limit();
								h_candidates[4*f + 1] = h_candidates[4*f + 1]/(sampling_time*((float) inBin[active_range])*((double) jerk_strategy.nSamples_freq_dom()));
							}
							
							FILE *fp_out;
							sprintf(filename, "jerk_search-dm%d-w%d.dat", (int) active_DM, (int) (ZW_planes_shift + zp));
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
						allcandidates.AddSubListFromGPU(nCandidates, d_ZW_candidates, w, DM, jerk_strategy.nSamples_time_dom(), jerk_strategy.nFilters_z(), jerk_strategy.z_max_search_limit(), jerk_strategy.z_search_step(), sampling_time*((float) inBin[active_range]), inBin[active_range]);
						if(VERBOSE>2) printf("JERK SEARCH -> Number of candidates: %d\n", nCandidates);
					}
					
					ZW_planes_shift = ZW_planes_shift + nZW_planes;
				}
				//------->
				char str[100];
				sprintf(str, "jerk_results_r%d_dm%d.dat", (int) active_range, (int) active_DM);
				size_t nCandidates_for_current_DM = allcandidates.ExportToFile(str);
				
				//Save candidates to disc
				timer_DM.Stop();
				timer_DM.Elapsed();
				time_per_range = time_per_range + timer_DM.Elapsed();
				if(VERBOSE>2) printf("JERK SEARCH -> Time per DM trial %fms\n", timer_DM.Elapsed());
				float DM = dm_low[active_range] + ((float) active_DM)*dm_step[active_range];
				if(VERBOSE>1) printf("JERK search: current DM=%f; Time taken %fms; Candidates found: %zu;\n", DM, timer_DM.Elapsed(), nCandidates_for_current_DM);
			}
			
			//-------> cuFFT
			cufftDestroy(cuFFT_plan);

			if ( cudaSuccess != cudaFree(d_ZW_candidates)) printf("ERROR while deallocating d_ZW_candidates!\n");
			d_ZW_candidates = NULL;
			if ( cudaSuccess != cudaFree(d_ZW_planes)) printf("ERROR while deallocating d_ZW_planes!\n");
			d_ZW_planes = NULL;
			if ( cudaSuccess != cudaFree(d_MSD_workarea)) printf("ERROR while deallocating d_MSD_workarea!\n");
			d_MSD_workarea = NULL;
			if(VERBOSE>0) printf("JERK search: range %d DMs[%f--%f] finished in %fms;\n", active_range, dm_low[active_range], dm_low[active_range] + (nDMs-1)*dm_step[active_range], time_per_range);
			time_per_range = 0;
		}
		
		
		
		timer_total.Stop();
		time_total = timer_total.Elapsed()/1000.0;
		if(VERBOSE>0) printf("Total time for JERK search: %g s \n", time_total);
		if(VERBOSE>0) printf("Time for convolution: %g s \n", Convolution_time);
		if(VERBOSE>0) printf("Time for candidate selection: %g s \n", Candidate_time);
		if(VERBOSE>0) printf("Time for MSD: %g s \n", MSD_time);


		
		delete [] h_jerk_filters;
		if ( cudaSuccess != cudaFree(d_jerk_filters)) printf("ERROR while deallocating d_jerk_filters!\n");
		if ( cudaSuccess != cudaFree(d_DM_trial)) printf("ERROR while deallocating d_DM_trial!\n");
		if ( cudaSuccess != cudaFree(d_DM_trial_ffted)) printf("ERROR while deallocating d_DM_trial_ffted!\n");
		if ( cudaSuccess != cudaFree(d_MSD)) printf("ERROR while deallocating d_MSD!\n");
		return(0);
	}
	
}// namespace
