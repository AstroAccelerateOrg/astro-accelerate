#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <cufft.h>
#include <math.h>
#include <vector>
#include <omp.h>
#include <string.h>

#include "aa_device_periods.hpp"
#include "aa_params.hpp"

#include "aa_periodicity_strategy.hpp"
#include "aa_periodicity_candidates.hpp"
#include "aa_device_MSD_Configuration.hpp"
#include "aa_device_MSD.hpp"
#include "aa_device_MSD_plane_profile.hpp"
#include "aa_device_peak_find.hpp"
#include "aa_device_power.hpp"
#include "aa_device_spectrum_whitening.hpp"
#include "aa_device_harmonic_summing.hpp"
#include "aa_corner_turn.hpp"
#include "aa_device_threshold.hpp"
#include "aa_gpu_timer.hpp"
#include "aa_host_utilities.hpp"
#include "aa_timelog.hpp"
#include "presto_funcs.hpp"
#include "presto.hpp"


namespace astroaccelerate {

  // define to see debug info
  //#define GPU_PERIODICITY_SEARCH_DEBUG
  
  // define to perform CPU spectral whitening debug
  //#define CPU_SPECTRAL_WHITENING_DEBUG
  void checkCudaErrors( cudaError_t CUDA_error){
    if(CUDA_error != cudaSuccess) {
      printf("CUDA error: %d\n", CUDA_error);
    }
  }
  
  void printMem(){
    size_t free = 0;
    size_t total = 0;
    cudaMemGetInfo(&free, &total);
    printf("Free memory: %zu;\n", free);
  }

  /**
   * \class Candidate_List aa_device_periods.cu "src/aa_device_periods.cu"
   * \brief Class for AstroAccelerate to manage the candidate list for periodicity.
   * \brief The user should not use this class for interacting with periodicity.
   * \author -
   * \date -
   */
  class Candidate_List {
  private:


  public:
    static const int el=4; // number of columns in the candidate list
    std::vector<float> list;
    int rangeid;

    /** \brief Constructor for Candidate_List. */
    Candidate_List(int t_rangeid) {
      rangeid = t_rangeid;
      list.clear();
    }

    /** \brief Allocator for list data. */
    void Allocate(int nCandidates) {
      list.resize(nCandidates*el);
    }

    /** \brief Returns size of list data container. */
    size_t size() {
      return((list.size()/el));
    }

    /** \brief Processes the periodicity inBin data group. */
    void Process(float *MSD, aa_periodicity_range *Prange, float mod) {
      float dm_step       = Prange->range.dm_step();
      float dm_low        = Prange->range.dm_low();
      float sampling_time = Prange->range.sampling_time();
      float nTimesamples  = Prange->range.nTimesamples();
      int nPoints_before  = size();
      int harmonics;

      for(size_t c=0; c< size(); c++) {
        harmonics = (int) list[c*el+3];
        list[c*el+0] = list[c*el+0]*dm_step + dm_low;
        list[c*el+1] = list[c*el+1]*(1.0/(sampling_time*nTimesamples*mod));
        list[c*el+2] = (list[c*el+2] - MSD[2*harmonics])/(MSD[2*harmonics+1]);
        list[c*el+3] = list[c*el+3];
      }
    }
    
    /** \brief Rescales and then processes inBin data group. */
    void Rescale_Threshold_and_Process(float *MSD, aa_periodicity_range *Prange, float sigma_cutoff, float mod) {
      float SNR;
      float dm_step       = Prange->range.dm_step();
      float dm_low        = Prange->range.dm_low();
      float sampling_time = Prange->range.sampling_time();
      int nTimesamples    = Prange->range.nTimesamples();
      int nPoints_before = size();
      int nPoints_after;

      std::vector<float> new_list;
      for(int c=0; c<(int)size(); c++) {
        float oldSNR = list[c*el+2];
        int harmonics = (int) list[c*el+3];
        SNR = (list[c*el+2] - MSD[2*harmonics])/(MSD[2*harmonics+1]);

        if(SNR>sigma_cutoff) {
          new_list.push_back(list[c*el+0]*dm_step + dm_low);
          new_list.push_back(list[c*el+1]*(1.0/(sampling_time*nTimesamples*mod)));
          new_list.push_back(SNR);
          new_list.push_back(list[c*el+3]);
        }
      }
      list.clear();
      list = new_list;
      nPoints_after = size();
      printf("   Before: %d; After: %d; sigma_cutoff:%f\n", nPoints_before, nPoints_after, sigma_cutoff);
    }
  };

  /**
   * \class GPU_Memory_for_Periodicity_Search aa_device_periods.cu "src/aa_device_periods.cu"
   * \brief Class for managing GPU memory for periodicity search.
   * \brief It is not clear how much this class interacts with other parts of the codebase to notify of its memory usage.
   * \brief Users should not use this class for interacting with periodicity.
   * \author -
   * \date -
   **/
  class GPU_Memory_for_Periodicity_Search {
  private:
    int MSD_interpolated_size;
    int MSD_DIT_size;

  public:
    float *d_one_A;
    float *d_two_B;
    float *d_half_C;
    
    ushort *d_power_harmonics;
    ushort *d_interbin_harmonics;
    
    // Candidate list
    int *gmem_power_peak_pos;
    int *gmem_interbin_peak_pos;
    
    // MSD
    float *d_MSD;
    float *d_previous_partials;
    float *d_all_blocks;
    
    // cuFFT
    void *cuFFT_workarea;
    
    void Allocate(aa_periodicity_strategy &PSR_strategy){
      int nHarmonics = PSR_strategy.nHarmonics();
      MSD_interpolated_size = nHarmonics;
      MSD_DIT_size = ((int) floorf(log2f((float) nHarmonics))) + 2;
      size_t t_input_plane_size = PSR_strategy.input_plane_size();
        
      if ( cudaSuccess != cudaMalloc((void **) &d_one_A,  sizeof(float)*t_input_plane_size )) printf("Periodicity Allocation error! d_one_A\n");
      if ( cudaSuccess != cudaMalloc((void **) &d_two_B,  sizeof(float)*2*t_input_plane_size )) printf("Periodicity Allocation error! d_two_B\n");
      if ( cudaSuccess != cudaMalloc((void **) &d_half_C,  sizeof(float)*t_input_plane_size/2 )) printf("Periodicity Allocation error! d_spectra_Real\n");
        
      if ( cudaSuccess != cudaMalloc((void **) &d_power_harmonics, sizeof(ushort)*t_input_plane_size )) printf("Periodicity Allocation error! d_harmonics\n");
      if ( cudaSuccess != cudaMalloc((void **) &d_interbin_harmonics, sizeof(ushort)*t_input_plane_size )) printf("Periodicity Allocation error! d_harmonics\n");
        
      if ( cudaSuccess != cudaMalloc((void**) &gmem_power_peak_pos, 1*sizeof(int)) )  printf("Periodicity Allocation error! gmem_power_peak_pos\n");
      if ( cudaSuccess != cudaMalloc((void**) &gmem_interbin_peak_pos, 1*sizeof(int)) )  printf("Periodicity Allocation error! gmem_interbin_peak_pos\n");
        
      if ( cudaSuccess != cudaMalloc((void**) &d_MSD, sizeof(float)*MSD_interpolated_size*2)) {printf("Periodicity Allocation error! d_MSD\n");}
      
      if ( cudaSuccess != cudaMalloc((void**) &d_previous_partials, sizeof(float)*MSD_DIT_size*MSD_PARTIAL_SIZE)) {printf("Periodicity Allocation error! d_previous_partials\n");}
      if ( cudaSuccess != cudaMalloc((void**) &d_all_blocks, sizeof(float)*PSR_strategy.max_total_MSD_blocks()*MSD_PARTIAL_SIZE)) {printf("Periodicity Allocation error! d_MSD\n");}
        
      if ( cudaSuccess != cudaMalloc((void **) &cuFFT_workarea, PSR_strategy.cuFFT_workarea_size()) ) {printf("Periodicity Allocation error! cuFFT_workarea\n");}
    }
    
    void Reset_MSD(){
      cudaMemset(d_MSD, 0, MSD_interpolated_size*2*sizeof(float));
      cudaMemset(d_previous_partials, 0, MSD_DIT_size*MSD_PARTIAL_SIZE*sizeof(float));
    }
    
    void Reset_Candidate_List(){
      cudaMemset(gmem_power_peak_pos, 0, sizeof(int));
      cudaMemset(gmem_interbin_peak_pos, 0, sizeof(int));
    }
    
    int Get_Number_of_Power_Candidates(){
      int temp;
      cudaError_t e = cudaMemcpy(&temp, gmem_power_peak_pos, sizeof(int), cudaMemcpyDeviceToHost);

      if(e != cudaSuccess) {
        LOG(log_level::error, "Could not cudaMemcpy in aa_device_periods.cu (" + std::string(cudaGetErrorString(e)) + ")");
      }
      
      return( (int) temp);
    }
    
    int Get_Number_of_Interbin_Candidates(){
      int temp;
      cudaError_t e = cudaMemcpy(&temp, gmem_interbin_peak_pos, sizeof(int), cudaMemcpyDeviceToHost);
      
      if(e != cudaSuccess) {
        LOG(log_level::error, "Could not cudaMemcpy in aa_device_periods.cu (" + std::string(cudaGetErrorString(e)) + ")");
      }
      
      return( (int) temp);
    }
    
    void Get_MSD(float *h_MSD){
      cudaError_t e = cudaMemcpy(h_MSD, d_MSD, MSD_interpolated_size*2*sizeof(float), cudaMemcpyDeviceToHost);

      if(e != cudaSuccess) {
        LOG(log_level::error, "Could not cudaMemcpy in aa_device_periods.cu (" + std::string(cudaGetErrorString(e)) + ")");
      }
    }
    
    void Get_MSD_partials(float *h_MSD_partials){
      cudaError_t e = cudaMemcpy(h_MSD_partials, d_previous_partials, MSD_DIT_size*MSD_PARTIAL_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
      
      if(e != cudaSuccess) {
        LOG(log_level::error, "Could not cudaMemcpy in aa_device_periods.cu (" + std::string(cudaGetErrorString(e)) + ")");
      }
    }
    
    void Set_MSD_partials(float *h_MSD_partials){
      cudaError_t e = cudaMemcpy(d_previous_partials, h_MSD_partials, MSD_PARTIAL_SIZE*sizeof(float), cudaMemcpyHostToDevice);

      if(e != cudaSuccess) {
        LOG(log_level::error, "Could not cudaMemcpy in aa_device_periods.cu (" + std::string(cudaGetErrorString(e)) + ")");
      }
    }
    
    /** \brief Destructor for GPU_Memory_for_Periodicity_Search. */
    ~GPU_Memory_for_Periodicity_Search(){
      cudaFree(d_one_A);
      cudaFree(d_two_B);
      cudaFree(d_half_C);
      cudaFree(d_power_harmonics);
      cudaFree(d_interbin_harmonics);
      cudaFree(gmem_power_peak_pos);
      cudaFree(gmem_interbin_peak_pos);
      cudaFree(d_MSD);
      cudaFree(d_previous_partials);
      cudaFree(d_all_blocks);
      cudaFree(cuFFT_workarea);
    }
  };


  void Copy_data_for_periodicity_search(float *d_one_A, float **dedispersed_data, aa_periodicity_batch *batch){ //TODO add "cudaStream_t stream1"
    int nStreams = 16;
    cudaStream_t stream_copy[16];
    cudaError_t e;
    float *h_small_dedispersed_data;
    size_t FFT_data_size  = batch->nTimesamples*sizeof(float); // FFT size
    size_t copy_data_size = batch->nTimesamples_to_copy*sizeof(float); // data size to copy
    cudaMallocHost((void **) &h_small_dedispersed_data, nStreams*FFT_data_size);
    memset(h_small_dedispersed_data, 0 , nStreams*FFT_data_size);

    for (int i = 0; i < nStreams; i++){
      e = cudaStreamCreate(&stream_copy[i]);
      if (e != cudaSuccess) {
        LOG(log_level::error, "Could not create streams in periodicity (" + std::string(cudaGetErrorString(e)) + ")");
      }
    }

    size_t stream_offset = batch->nTimesamples;

    #pragma omp parallel for num_threads(nStreams) shared(h_small_dedispersed_data, FFT_data_size, copy_data_size, d_one_A, stream_copy, stream_offset)
    for(int ff=0; ff<(int) batch->nDMs_per_batch; ff++){
      int id_stream = omp_get_thread_num();
      memcpy(h_small_dedispersed_data + id_stream*stream_offset, dedispersed_data[batch->DM_shift + ff], copy_data_size);
      //e = cudaMemcpy( &d_one_A[ff*batch->nTimesamples], dedispersed_data[batch->DM_shift + ff], batch->nTimesamples*sizeof(float), cudaMemcpyHostToDevice);      
      e = cudaMemcpyAsync(&d_one_A[ff*batch->nTimesamples], h_small_dedispersed_data + id_stream*stream_offset, FFT_data_size, cudaMemcpyHostToDevice, stream_copy[id_stream]);      
      cudaStreamSynchronize(stream_copy[id_stream]);
      if(e != cudaSuccess) {
        LOG(log_level::error, "Could not cudaMemcpy in aa_device_periods.cu (" + std::string(cudaGetErrorString(e)) + ")");
      }
    }
    for (int i = 0; i < nStreams; i++){
      e = cudaStreamDestroy(stream_copy[i]);
      if (e != cudaSuccess) {
        LOG(log_level::error, "Could not destroy stream in periodicity (" + std::string(cudaGetErrorString(e)) + ")");
      }
    }

    cudaFreeHost(h_small_dedispersed_data);

  }

  __inline__ float Calculate_frequency(int m, float sampling_time, int nTimesamples){
    return( ((float) m)/(sampling_time*((float) nTimesamples)) );
  }

  void Export_data_in_range(float *GPU_data, int nTimesamples, int nDMs, const char *filename, float dm_step, float dm_low, float sampling_time, int outer_DM_shift, int DMs_per_file=100) {
    char final_filename[100];
    std::ofstream FILEOUT;
    int lower, higher, inner_DM_shift;
    int data_mod = 3;
    if(DMs_per_file<0) DMs_per_file=nDMs;
    
    float *h_data, *h_export;
    size_t data_size = ((size_t) nTimesamples)*((size_t) nDMs);
    size_t export_size = ((size_t) nTimesamples)*((size_t) DMs_per_file)*data_mod;
    h_data = new float[data_size];
    h_export = new float[export_size];
    
    cudaError_t e = cudaMemcpy(h_data, GPU_data, data_size*sizeof(float), cudaMemcpyDeviceToHost);
    
    if(e != cudaSuccess) {
      LOG(log_level::error, "Could not cudaMemcpy in aa_device_periods.cu (" + std::string(cudaGetErrorString(e)) + ")");
    }
    
    int nRepeats = nDMs/DMs_per_file;
    int nRest = nDMs%DMs_per_file;
    std::vector<int> chunk_size;
    for(int f=0; f<nRepeats; f++) chunk_size.push_back(DMs_per_file);
    if(nRest>0) chunk_size.push_back(nRest);
    printf("Data will be exported into %d files\n", (int) chunk_size.size());
    
    inner_DM_shift = 0;
    for(int i=0; i<(int) chunk_size.size(); i++){
      lower = outer_DM_shift + inner_DM_shift;
      higher = outer_DM_shift + inner_DM_shift + chunk_size[i];
      sprintf(final_filename,"%s_%f_%f.dat", filename, lower*dm_step+dm_low, higher*dm_step+dm_low);
      printf("Exporting file %s\n", final_filename);
        
      for(int dm = 0; dm<chunk_size[i]; dm++) {
        for(int t=0; t<nTimesamples; t++){
          int pos = dm*nTimesamples + t;
          h_export[data_mod*pos + 0] = (lower + dm)*dm_step + dm_low;
          h_export[data_mod*pos + 1] = Calculate_frequency(t, sampling_time, nTimesamples);
          h_export[data_mod*pos + 2] = h_data[(inner_DM_shift+dm)*nTimesamples + t];
        }
      }
        
      FILE *fp_out;
      if (( fp_out = fopen(final_filename, "wb") ) == NULL) {
        LOG(log_level::error, "Error opening output file!");
      }
      fwrite(h_export, nTimesamples*chunk_size[i]*sizeof(float), 3, fp_out);
      fclose(fp_out);
        
      inner_DM_shift = inner_DM_shift + chunk_size[i];
    }
    
    delete [] h_data;
    delete [] h_export;
  }

  /**
   * \brief Performs a periodicity search on the GPU.
   * \todo Clarify the difference between Periodicity_search and GPU_periodicity.
   **/
  void Periodicity_search(GPU_Memory_for_Periodicity_Search *gmem, aa_periodicity_strategy PSR_strategy, double *compute_time, size_t input_plane_size, aa_periodicity_range *Prange, aa_periodicity_batch *batch, std::vector<int> *h_boxcar_widths, int harmonic_sum_algorithm, bool enable_scalloping_loss_removal){
    bool transposed_data = true;
    if(harmonic_sum_algorithm != 0){
        transposed_data = false;
    }
    TimeLog time_log;
    
    int local_max_list_size = (input_plane_size)/4;
    if(local_max_list_size > 1073741824) local_max_list_size = 1073741824;
    
    
    float *d_dedispersed_data, *d_FFT_complex_output, *d_frequency_power, *d_frequency_interbin, *d_frequency_power_CT, *d_frequency_interbin_CT, *d_power_SNR, *d_interbin_SNR, *d_power_list, *d_interbin_list, *d_MSD_workarea;
    
    d_dedispersed_data      = gmem->d_one_A;
    d_FFT_complex_output    = gmem->d_two_B;
    d_MSD_workarea          = gmem->d_two_B;
    d_frequency_power       = gmem->d_half_C;
    d_frequency_interbin    = gmem->d_one_A;
    if(transposed_data){
        d_frequency_power_CT    = &gmem->d_two_B[0];
        d_frequency_interbin_CT = &gmem->d_two_B[input_plane_size];
        d_power_SNR             = gmem->d_half_C;
        d_interbin_SNR          = gmem->d_one_A;
        d_power_list            = &gmem->d_two_B[0];
        d_interbin_list         = &gmem->d_two_B[input_plane_size];
        local_max_list_size = local_max_list_size/2;
    }
    else {
        d_frequency_power_CT    = NULL;
        d_frequency_interbin_CT = NULL;
        d_power_SNR             = &gmem->d_two_B[0];
        d_interbin_SNR          = &gmem->d_two_B[input_plane_size];
        d_power_list            = gmem->d_half_C;
        d_interbin_list         = gmem->d_one_A;
    }
    
    int t_nTimesamples      = batch->nTimesamples;
    int t_nTSamplesFFT      = (t_nTimesamples>>1) + 1;
    int t_nDMs_per_batch    = batch->nDMs_per_batch;
    int t_DM_shift          = batch->DM_shift;
    int t_inBin             = 1; // this is because nTimesamples and sampling time is already adjusted for binning
    
    aa_gpu_timer timer;
    
    //---------> cuFFT
    timer.Start();
    cufftHandle plan_input;
    cufftResult cufft_error;
    
    size_t cuFFT_workarea_size;
    
    cufft_error = cufftCreate(&plan_input);
    if (CUFFT_SUCCESS != cufft_error) printf("CUFFT error: %d", cufft_error);
    cufftSetAutoAllocation(plan_input, false);
    
    cufft_error = cufftMakePlan1d(plan_input, t_nTimesamples, CUFFT_R2C, t_nDMs_per_batch, &cuFFT_workarea_size);
    if (CUFFT_SUCCESS != cufft_error) printf("CUFFT error: %d", cufft_error);
    
    cufft_error = cufftSetWorkArea(plan_input, gmem->cuFFT_workarea);
    if (CUFFT_SUCCESS != cufft_error) printf("CUFFT error: %d", cufft_error);

    cufft_error = cufftExecR2C(plan_input, (cufftReal *)d_dedispersed_data, (cufftComplex *)d_FFT_complex_output);
    if ( cufft_error != CUFFT_SUCCESS) printf("CUFFT error: %d\n", cufft_error);
    
    cufft_error = cufftDestroy(plan_input);
    if ( cufft_error != CUFFT_SUCCESS) printf("CUFFT error: %d\n", cufft_error);

    timer.Stop();
    time_log.adding("PSR","cuFFT",timer.Elapsed());
    (*compute_time) = (*compute_time) + timer.Elapsed();
    //---------<
    
    
    #ifdef CPU_SPECTRAL_WHITENING_DEBUG
    CPU_spectral_whitening(d_FFT_complex_output, d_dedispersed_data, Prange->range.dm_step, Prange->range.dm_low, t_nTSamplesFFT, t_nDMs_per_batch, t_DM_shift);
    #endif
    
    //---------> Spectrum whitening
    timer.Start();
    cudaStream_t stream; stream = NULL;
    spectrum_whitening_SGP2((float2 *) d_FFT_complex_output, t_nTSamplesFFT, t_nDMs_per_batch, true, stream);
    timer.Stop();
    time_log.adding("PSR","spectrum whitening",timer.Elapsed());
    (*compute_time) = (*compute_time) + timer.Elapsed();
    //---------<
    
    
    //---------> Calculate powers and interbinning
    timer.Start();
    simple_power_and_interbin( (float2 *) d_FFT_complex_output, d_frequency_power, d_frequency_interbin, t_nTimesamples, t_nDMs_per_batch);
    timer.Stop();
    time_log.adding("PSR","power spectrum",timer.Elapsed());
    (*compute_time) = (*compute_time) + timer.Elapsed();
    //---------<
    
    
    //---------> Mean and StDev on powers
    timer.Start();
    bool perform_continuous = false;

    double total_time, dit_time, MSD_time;
    MSD_plane_profile(gmem->d_MSD, d_frequency_power, gmem->d_previous_partials, d_MSD_workarea, true, (t_nTimesamples>>1), t_nDMs_per_batch, h_boxcar_widths, 0, 0, 0, PSR_strategy.sigma_outlier_rejection_threshold(), PSR_strategy.enable_outlier_rejection(), perform_continuous, &total_time, &dit_time, &MSD_time);
    
    timer.Stop();
    time_log.adding("PSR","MSD",timer.Elapsed());
    (*compute_time) = (*compute_time) + timer.Elapsed();
    //---------<
    
    
    //---------> Harmonic sum
    timer.Start();
    if(harmonic_sum_algorithm == 0){
        //---------> Corner turn
        corner_turn_SM(d_frequency_power, d_frequency_power_CT, (t_nTimesamples>>1), t_nDMs_per_batch);
        corner_turn_SM(d_frequency_interbin, d_frequency_interbin_CT, t_nTimesamples, t_nDMs_per_batch);
        //---------<
        
        //---------> Simple harmonic summing
        periodicity_simple_harmonic_summing(d_frequency_power_CT, d_power_SNR, gmem->d_power_harmonics, gmem->d_MSD, (t_nTimesamples>>1), t_nDMs_per_batch, PSR_strategy.nHarmonics());
        periodicity_simple_harmonic_summing(d_frequency_interbin_CT, d_interbin_SNR, gmem->d_interbin_harmonics, gmem->d_MSD, t_nTimesamples, t_nDMs_per_batch, PSR_strategy.nHarmonics());
        //---------<
    }
    else if(harmonic_sum_algorithm == 1) {
        //---------> Greedy harmonic summing
        periodicity_greedy_harmonic_summing(d_frequency_power, d_power_SNR, gmem->d_power_harmonics, gmem->d_MSD, (t_nTimesamples>>1), t_nDMs_per_batch, PSR_strategy.nHarmonics(), enable_scalloping_loss_removal);
        periodicity_greedy_harmonic_summing(d_frequency_interbin, d_interbin_SNR, gmem->d_interbin_harmonics, gmem->d_MSD, t_nTimesamples, t_nDMs_per_batch, PSR_strategy.nHarmonics(), enable_scalloping_loss_removal);
        //---------<
    }
    else if(harmonic_sum_algorithm == 2) {
        //---------> PRESTO plus harmonic summing
        periodicity_presto_plus_harmonic_summing(d_frequency_power, d_power_SNR, gmem->d_power_harmonics, gmem->d_MSD, (t_nTimesamples>>1), t_nDMs_per_batch, PSR_strategy.nHarmonics(), enable_scalloping_loss_removal);
        periodicity_presto_plus_harmonic_summing(d_frequency_interbin, d_interbin_SNR, gmem->d_interbin_harmonics, gmem->d_MSD, t_nTimesamples, t_nDMs_per_batch, PSR_strategy.nHarmonics(), enable_scalloping_loss_removal);
        //---------<
    }
    else if(harmonic_sum_algorithm == 3) {
        //---------> PRESTO harmonic summing
        periodicity_presto_harmonic_summing(d_frequency_power, d_power_SNR, gmem->d_power_harmonics, gmem->d_MSD, (t_nTimesamples>>1), t_nDMs_per_batch, PSR_strategy.nHarmonics(), enable_scalloping_loss_removal);
        periodicity_presto_harmonic_summing(d_frequency_interbin, d_interbin_SNR, gmem->d_interbin_harmonics, gmem->d_MSD, t_nTimesamples, t_nDMs_per_batch, PSR_strategy.nHarmonics(), enable_scalloping_loss_removal);
        //---------<
    }
    timer.Stop();
    time_log.adding("PSR","harmonic sum",timer.Elapsed());
    (*compute_time) = (*compute_time) + timer.Elapsed();
    //---------<
    
    
    //---------> Peak finding
    int enable_greedy_postprocessing = 0;
    if(harmonic_sum_algorithm == 1) enable_greedy_postprocessing = 1;
    timer.Start();
    if(PSR_strategy.candidate_selection_algorithm()==1){
        //-------------- Thresholding
        if(harmonic_sum_algorithm == 0){
            Threshold_for_periodicity_transposed(d_power_SNR, gmem->d_power_harmonics, d_power_list, gmem->gmem_power_peak_pos, gmem->d_MSD, PSR_strategy.sigma_cutoff(), t_nDMs_per_batch, (t_nTimesamples>>1), t_DM_shift, t_inBin, local_max_list_size);
            Threshold_for_periodicity_transposed(d_interbin_SNR, gmem->d_interbin_harmonics, d_interbin_list, gmem->gmem_interbin_peak_pos, gmem->d_MSD, PSR_strategy.sigma_cutoff(), t_nDMs_per_batch, t_nTimesamples, t_DM_shift, t_inBin, local_max_list_size);
        }
        else {
            Threshold_for_periodicity_normal(d_power_SNR, gmem->d_power_harmonics, d_power_list, gmem->gmem_power_peak_pos, gmem->d_MSD, PSR_strategy.sigma_cutoff(), (t_nTimesamples>>1), t_nDMs_per_batch, t_DM_shift, t_inBin, local_max_list_size, enable_greedy_postprocessing);
            Threshold_for_periodicity_normal(d_interbin_SNR, gmem->d_interbin_harmonics, d_interbin_list, gmem->gmem_interbin_peak_pos, gmem->d_MSD, PSR_strategy.sigma_cutoff(), t_nTimesamples, t_nDMs_per_batch, t_DM_shift, t_inBin, local_max_list_size, enable_greedy_postprocessing);
        }
        //-------------- Thresholding
    }
    else {
        //-------------- Peak finding
        Peak_find_for_periodicity_search(d_power_SNR, gmem->d_power_harmonics, d_power_list, (t_nTimesamples>>1), t_nDMs_per_batch, PSR_strategy.sigma_cutoff(), local_max_list_size, gmem->gmem_power_peak_pos, gmem->d_MSD, t_DM_shift, t_inBin, transposed_data, enable_greedy_postprocessing);
        Peak_find_for_periodicity_search(d_interbin_SNR, gmem->d_interbin_harmonics, d_interbin_list, t_nTimesamples, t_nDMs_per_batch, PSR_strategy.sigma_cutoff(), local_max_list_size, gmem->gmem_interbin_peak_pos, gmem->d_MSD, t_DM_shift, t_inBin, transposed_data, enable_greedy_postprocessing);
        //-------------- Peak finding
    }
    timer.Stop();
    time_log.adding("PSR","candidates",timer.Elapsed());
    (*compute_time) = (*compute_time) + timer.Elapsed();
    //---------<
    
  }

  void Export_Data_To_File(std::vector<Candidate_List> candidates, const char *filename) {
    FILE *fp_out;
    if((fp_out = fopen(filename, "wb")) == NULL) {
      LOG(log_level::error, "Error opening output file!\n");
    }

    for(int f=0; f<(int) candidates.size(); f++) {
      fwrite(&candidates[f].list[0], candidates[f].size()*sizeof(float), Candidate_List::el, fp_out);
    }
    fclose(fp_out);
  }

  int Get_Number_of_Candidates(int *GPU_data){
    int temp;
    cudaError_t e = cudaMemcpy(&temp, GPU_data, sizeof(int), cudaMemcpyDeviceToHost);

    if(e != cudaSuccess) {
      LOG(log_level::error, "Could not cudaMemcpy in aa_device_periods.cu (" + std::string(cudaGetErrorString(e)) + ")");
    }
    
    return(temp);
  }




  /** \brief Function that performs a GPU periodicity search. */
  void GPU_periodicity(aa_periodicity_strategy &PSR_strategy, float ***output_buffer, aa_periodicity_candidates &Power_Candidates, aa_periodicity_candidates &Interbin_Candidates) {
    TimeLog time_log;
    
    LOG(log_level::notice, "------------ STARTING PERIODICITY SEARCH ------------");
    
    
    // Determining available memory (temporary, it should be moved elsewhere) COULD BE DELETED
    //size_t memory_available,total_mem;
    //cudaMemGetInfo(&memory_available,&total_mem);
    
    std::vector<int> h_boxcar_widths; h_boxcar_widths.resize(PSR_strategy.nHarmonics()); 
    for(int f=0; f<PSR_strategy.nHarmonics(); f++) h_boxcar_widths[f]=f+1;
    
    //--------> Allocation of GPU memory
    GPU_Memory_for_Periodicity_Search GPU_memory;
    GPU_memory.Allocate(PSR_strategy);
    
    float h_MSD[PSR_strategy.nHarmonics()*2];
    //----------------------------<
    
    // Timing
    double Total_periodicity_time = 0, Total_calc_time = 0, calc_time_per_range = 0, Total_copy_time = 0, copy_time_per_range = 0;
    aa_gpu_timer timer, periodicity_timer;
    
    periodicity_timer.Start();
    
    LOG(log_level::notice, "------ CONFIGURATION OF PERIODICITY SEARCH DONE -----");
    //---------------------------------------------------------------
    //---------------------------------------------------------------
    
    
    for(int r=0; r<(int) PSR_strategy.nRanges(); r++) {
        aa_periodicity_range current_p_range = PSR_strategy.get_periodicity_range(r);
        //std::vector<Candidate_List> PowerCandidates;
        //std::vector<Candidate_List> InterbinCandidates;

        GPU_memory.Reset_MSD();
        current_p_range.print();
        
        for(int b=0; b<(int) current_p_range.batches.size(); b++) {
            current_p_range.batches[b].print(b);
            
            GPU_memory.Reset_Candidate_List();
            
            //---------> Copy input data to the device
            timer.Start();
            Copy_data_for_periodicity_search(GPU_memory.d_one_A, output_buffer[current_p_range.rangeid], &current_p_range.batches[b]);
            timer.Stop();
            time_log.adding("PSR","Host-To-Device",timer.Elapsed());
            copy_time_per_range = copy_time_per_range + timer.Elapsed();
            //---------<

          // simple harmonic sum 0
          // greedy harmonic sum 1
          // presto plus harmonic sum 2
          // presto harmonic sum 3
          int harmonic_sum_algorithm = PSR_strategy.harmonic_sum_algorithm();
          bool enable_scalloping_loss_removal = PSR_strategy.enable_scalloping_loss_mitigation();
          
          //---------> Periodicity search
          Periodicity_search(&GPU_memory, PSR_strategy, &calc_time_per_range, PSR_strategy.input_plane_size(), &current_p_range, &current_p_range.batches[b], &h_boxcar_widths, harmonic_sum_algorithm, enable_scalloping_loss_removal);
          //---------<
    
          
          GPU_memory.Get_MSD(h_MSD);
          
          //---------> Copy candidates to the host
          timer.Start();
          size_t nPowerCandidates = GPU_memory.Get_Number_of_Power_Candidates();
          size_t nInterbinCandidates = GPU_memory.Get_Number_of_Interbin_Candidates();
          LOG(log_level::debug, " PSR: Total number of candidates found in this range is " + std::to_string(nPowerCandidates) + ";");
          LOG(log_level::debug, " PSR with inter-binning: Total number of candidates found in this range is " + std::to_string(nInterbinCandidates) + ";");
          
          float  range_dm_low  = current_p_range.range.dm_low();
          float  range_dm_step = current_p_range.range.dm_step();
          int    range_id      = current_p_range.rangeid;
          size_t range_nTimesamples = current_p_range.range.nTimesamples(); 
          double range_sampling_time = current_p_range.range.sampling_time();
          float  *pointer_to_candidate_data = NULL;
          int enable_greedy_postprocessing = 0;
          if(harmonic_sum_algorithm == 1) enable_greedy_postprocessing = 1;
          
          // Power candidates
          if(harmonic_sum_algorithm == 0) pointer_to_candidate_data = &GPU_memory.d_two_B[0];
          else pointer_to_candidate_data = GPU_memory.d_half_C;
          Power_Candidates.Add_Candidates(pointer_to_candidate_data, nPowerCandidates, range_id, h_MSD, range_dm_low, range_dm_step, range_sampling_time, range_nTimesamples, 1.0, enable_greedy_postprocessing);
          
          // Interbin candidates
          if(harmonic_sum_algorithm == 0) pointer_to_candidate_data = &GPU_memory.d_two_B[PSR_strategy.input_plane_size()];
          else pointer_to_candidate_data = GPU_memory.d_one_A;
          Interbin_Candidates.Add_Candidates(pointer_to_candidate_data, nInterbinCandidates, range_id, h_MSD, range_dm_low, range_dm_step, range_sampling_time, range_nTimesamples, 2.0, enable_greedy_postprocessing);
          
          timer.Stop();
          time_log.adding("PSR","Device-To-Host",timer.Elapsed());
          copy_time_per_range = copy_time_per_range + timer.Elapsed();
      } //batches
        
        
      Total_calc_time = Total_calc_time + calc_time_per_range;
      calc_time_per_range = 0;
      Total_copy_time = Total_copy_time + copy_time_per_range;
      copy_time_per_range = 0;
    } // inBin ranges
    printf("-----------------------------------------------------------------------------------\n");
    
    periodicity_timer.Stop();
    Total_periodicity_time = periodicity_timer.Elapsed();
    time_log.adding("PSR","total",Total_periodicity_time);

    cudaDeviceSynchronize();
  };

} //namespace astroaccelerate
