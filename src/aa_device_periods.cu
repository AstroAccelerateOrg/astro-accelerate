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

namespace astroaccelerate {

  // define to see debug info
#define GPU_PERIODICITY_SEARCH_DEBUG

  // define to reuse old MSD results to generate a new one (it means new MSD is calculated from more samples) (WORKING WITHIN SAME INBIN VALUE)
#define PS_REUSE_MSD_WITHIN_INBIN

  // experimental and this might not be very useful with real noise BROKEN!
  //#define PS_REUSE_MSD_THROUGH_INBIN

  //#define OLD_PERIODICITY


  // Perhaps it would be best to collapse the structure to either
  // Inherit from dedispersion plan but introduce additional structure which would index dedispersion range, but the values would be different in dedispersion range Problem would be also how to assign ranges into inBin groups and how to link batches to ranges. Batches to ranges could be done vie Id of the range and inBin groups could be linked via list of ranges, but I'm not sure. 


  /**
   * \class Dedispersion_Range aa_device_periods.cu "src/aa_device_periods.cu"
   * \brief Class for AstroAccelerate to manage adding dedispersion measure ranges for processing periodicity.
   * \brief The user should not interact with this class for configuring dedispersion ranges. Instead, the user should use aa_ddtr_plan and aa_ddtr_strategy.
   * \author -
   * \date -
   */

  class Dedispersion_Range {
  public:
    float dm_low;
    float dm_high;
    float dm_step;
    float sampling_time;
    int inBin;
    int nTimesamples;
    int nDMs;
	

    /** \brief Constructor for Dedispersion_Range. */
    Dedispersion_Range() {
      dm_step       = 0;
      dm_low        = 0;
      dm_high       = 0;
      inBin         = 0;
      nTimesamples  = 0;
      nDMs          = 0;
      sampling_time = 0;
    }

    /** \brief Constructor for Dedispersion_Range. */
    Dedispersion_Range(float t_dm_low, float t_dm_high, float t_dm_step, int t_inBin, int t_range_timesamples, int t_ndms, float t_sampling_time) {
      dm_step       = t_dm_step;
      dm_low        = t_dm_low;
      dm_high       = t_dm_high;
      inBin         = t_inBin;
      nTimesamples  = t_range_timesamples;
      nDMs          = t_ndms;
      sampling_time = t_sampling_time;
    }

    /** \brief Method to assign or change values of an instance after construction. */
    void Assign(float t_dm_low, float t_dm_high, float t_dm_step, int t_inBin, int t_range_timesamples, int t_ndms, float t_sampling_time) {
      dm_step       = t_dm_step;
      dm_low        = t_dm_low;
      dm_high       = t_dm_high;
      inBin         = t_inBin;
      nTimesamples  = t_range_timesamples;
      nDMs          = t_ndms;
      sampling_time = t_sampling_time;
    }

    /** \brief Method to print member variables of an instance. */
    void print(void) {
      printf("dedispersion range: %f--%f:%f inBin:%d; nTimesamples:%d; nDMs:%d; sampling time:%f\n", dm_low, dm_high, dm_step, inBin, nTimesamples, nDMs, sampling_time);
    }
  };


  /**
   * \class Dedispersion_Plan aa_device_periods.cu "src/aa_device_periods.cu"
   * \brief Class for AstroAccelerate to manage a dedispersion plan for periodicity.
   * \brief The user should not interact with this class for configuring dedispersion ranges. Instead, the user should use aa_ddtr_plan.
   * \author -
   * \date -
   */
  class Dedispersion_Plan {
  public:
    std::vector<Dedispersion_Range> DD_range;
    float sampling_time;
    int nProcessedTimesamples;
    int max_nDMs;

    /** \brief Method for printing member variable information. */
    void print() {
      printf("De-dispersion ranges:\n");
      for(int f=0; f<(int)DD_range.size(); f++) DD_range[f].print();
    }
  };

  /**
   * \class Periodicity_Batch aa_device_periods.cu "src/aa_device_periods.cu"
   * \brief Class for AstroAccelerate to manage processing periodicity data.
   * \brief The user should not interact with this class for periodicity. Instead, the user should use aa_periodicity_plan and aa_periodicity_strategy.
   * \author -
   * \date -
   */
  class Periodicity_Batch {
  public:
    int nDMs_per_batch;
    int nTimesamples;
    int DM_shift;
    // cuFFT plan? no because cuFFT plan allocates GPU memory
    MSD_Configuration MSD_conf;
    std::vector<float> PowerCandidates;
    std::vector<float> InterbinCandidates;

    /** \brief Constructor for Periodicity_Batch. */
    Periodicity_Batch(int t_nDMs_per_batch, int t_DM_shift, Dedispersion_Range range, int total_blocks) {
      nDMs_per_batch = t_nDMs_per_batch;
      nTimesamples   = range.nTimesamples;
      DM_shift     = t_DM_shift;
		
      MSD_conf = *(new MSD_Configuration((nTimesamples>>1), nDMs_per_batch, 0, total_blocks)); // We divide nTimesamples by 2 because MSD is determined from power spectra (or because of FFT R->C)
    }

    /** \brief Destructor for Periodicity_Batch. */
    ~Periodicity_Batch() {
      //delete MSD_conf;
    }

    /** \brief Method for printing member variable data for an instance. */
    void print() {
      printf("      Batch: nDMs:%d; nTimesamples:%d; DM shift:%d;\n", nDMs_per_batch, nTimesamples, DM_shift);
      MSD_conf.print();
    }
  };


  /**
   * \class Periodicity_Range aa_device_periods.cu "src/aa_device_periods.cu"
   * \brief Class for AstroAccelerate to manage periodicity ranges.
   * \brief The user should not use this class for interacting with periodicity. Instead, the user should use aa_periodicity_plan and aa_periodicity_strategy.
   * \author -
   * \date -
   */
  class Periodicity_Range {
  public:
    Dedispersion_Range range;
    int rangeid;
    std::vector<Periodicity_Batch> batches;
    int total_MSD_blocks;

    void Create_Batches(int max_nDMs_in_range, int address_shift) {
      int nRepeats, nRest;
      total_MSD_blocks = 0;

      nRepeats = range.nDMs/max_nDMs_in_range;
      nRest = range.nDMs - nRepeats*max_nDMs_in_range;
      for(int f=0; f<nRepeats; f++) {
	Periodicity_Batch batch(max_nDMs_in_range, f*max_nDMs_in_range, range, total_MSD_blocks+address_shift);
	batches.push_back(batch);
	total_MSD_blocks = total_MSD_blocks + batch.MSD_conf.nBlocks_total;
	//batches.push_back( *(new Periodicity_Batch(max_nDMs_in_range, f*max_nDMs_in_range, range)) );
      }
      if(nRest>0) {
	Periodicity_Batch batch(nRest, nRepeats*max_nDMs_in_range, range, total_MSD_blocks+address_shift);
	batches.push_back(batch);
	total_MSD_blocks = total_MSD_blocks + batch.MSD_conf.nBlocks_total;
	//batches.push_back( *(new Periodicity_Batch(nRest, nRepeats*max_nDMs_in_range, range)) );
      }
    }

    void Calculate_Range_Properties(Dedispersion_Range t_range) {
      range = t_range;

      // Finding nearest lower power of two (because of FFT algorithm)
      int nearest = (int)floorf(log2f((float)range.nTimesamples));
      range.nTimesamples = (int) powf(2.0, (float) nearest);
      // This should not be necessary
      //range.nTimesamples = range.nTimesamples/range.inBin;
    }

    Periodicity_Range(Dedispersion_Range t_range, int max_nDMs_in_memory, int t_rangeid, int address_shift) {
      int max_nDMs_in_range;
      Calculate_Range_Properties(t_range);
      max_nDMs_in_range = max_nDMs_in_memory*range.inBin;
      rangeid = t_rangeid;
      Create_Batches(max_nDMs_in_range, address_shift);
    }

    /** \brief Destructor for Periodicity_Range. */
    ~Periodicity_Range() {
      //for(int f=0; f<(int)batches.size(); f++) delete batches[f];
      batches.clear();
    }

    /** \brief Method for printing member variable data of an instance. */
    void print(void) {
      int size = (int)batches.size();
      printf("  De-dispersion range: %f--%f:%f inBin:%d; nTimesamples:%d; nDMs:%d;\n", range.dm_low, range.dm_high, range.dm_step, range.inBin, range.nTimesamples, range.nDMs);
      if(size>1) printf("     De-dispersion range will be processed in %d batches each containing %d DM trials with tail of %d DM trials\n", size, batches[0].nDMs_per_batch, batches[size-1].nDMs_per_batch);
      else printf("     Periodicity search will run 1 batch containing %d DM trials.\n", batches[0].nDMs_per_batch);
      printf("     Total number of MSD blocks is %d which is %0.3f MB\n", total_MSD_blocks, (total_MSD_blocks*3.0*sizeof(float))/(1024.0*1024.0));
      printf("--------- Batches ----------\n");
      for(int f=0; f<(int)batches.size(); f++) batches[f].print();
      printf("\n");
    }
  };

  /**
   * \class Periodicity_inBin_Group aa_device_periods.cu "src/aa_device_periods.cu"
   * \brief Class for AstroAccelerate to manage inBin data for periodicity.
   * \brief The user should not use this class for interacting with periodicity. Instead the user should use aa_periodicity_plan and aa_periodicity_strategy.
   * \author -
   * \date -
   */
  class Periodicity_inBin_Group {
  public:
    int total_MSD_blocks;
    std::vector<Periodicity_Range> Prange;
	
    void print(void) {
      for(int f=0; f<(int) Prange.size(); f++){
	printf("  ->Periodicity range: %d\n", f);
	Prange[f].print();
      }		
    }
  };

  /**
   * \class AA_Periodicity_Plan aa_device_periods.cu "src/aa_device_periods.cu"
   * \brief Class for AstroAccelerate to manage periodicity plan.
   * \brief The user should not use this class for interacting with periodicity. Instead the user should use aa_periodicity_plan.
   * \author -
   * \date -
   */
  class AA_Periodicity_Plan {
  public:
    int nHarmonics;
    int max_total_MSD_blocks;
    int max_nTimesamples;
    int max_nDMs;
    int max_nDMs_in_memory;
    size_t input_plane_size;
    size_t cuFFT_workarea_size;
    std::vector<Periodicity_inBin_Group> inBin_group;

    /** \brief Destructor for AA_Periodicity_Plan. */
    ~AA_Periodicity_Plan() {
      inBin_group.clear();
    }
	
    /** \brief Clear inBin data from inBin_group. */
    void clear(){
      inBin_group.clear();
    }
	
    /** \brief Method for printing member variable data of an instance. */
    void print(void) {
      printf("max_total_MSD_blocks: %d MSD blocks;\n", max_total_MSD_blocks);
      printf("max_nTimesamples:     %d time samples;\n", max_nTimesamples);
      printf("max_nDMs:             %d DMs;\n", max_nDMs);
      printf("max_nDMs_in_memory:   %d DMs;\n", max_nDMs_in_memory);
      printf("input_plane_size:     %zu floating point numbers = %0.2f MB;\n", input_plane_size, ((float) input_plane_size*sizeof(float))/(1024.0*1024.0));
      printf("cuFFT_workarea_size:  %zu = %0.2f MB;\n", cuFFT_workarea_size, (float) cuFFT_workarea_size/(1024.0*1024.0));
      for(int f=0; f<(int) inBin_group.size(); f++){
	printf("inBin group: %d\n", f);
	inBin_group[f].print();
      }
    }
  };

  /**
   * \class Candidate_List aa_device_periods.cu "src/aa_device_periods.cu"
   * \brief Class for AstroAccelerate to manage the candidate list for periodicity.
   * \brief The user should not use this class for interacting with periodicity.
   * \author -
   * \date -
   */
  class Candidate_List {
  private:

    /**
     * \brief Calculates the signal-to-noise ratio (SNR) using a linear approximation.
     * \returns A linear approximation of the signal-to-noise ratio (SNR).
     */
    float linear_approximation(float value, int harmonic, float *MSD) {
      float SNR;
      float mean     = MSD[0];
      float sd       = MSD[1];
      float modifier = MSD[2];
      SNR = (value - (harmonic+1)*mean)/(sd + harmonic*modifier);
      return(SNR);
    }

    /**
     * \brief Calculates the signal-to-noise ratio (SNR) using a white-noise approximation.
     * \returns A white noise approximation of the signal-to-noise ratio (SNR).
     */
    float white_noise_approximation(float value, int harmonic, float *MSD) {
      float SNR;
      float mean = MSD[0];
      float sd   = MSD[1];
      SNR = (value - (harmonic+1)*mean)/(sqrt(harmonic+1)*sd);
      return(SNR);
    }

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
    void Process(float *MSD, Periodicity_inBin_Group *inBin_group, float mod) {
      float dm_step       = inBin_group->Prange[rangeid].range.dm_step;
      float dm_low        = inBin_group->Prange[rangeid].range.dm_low;
      float sampling_time = inBin_group->Prange[rangeid].range.sampling_time;
      float nTimesamples    = inBin_group->Prange[rangeid].range.nTimesamples;
      int nPoints_before = size();
      int harmonics;

      for(int c=0; c<(int)size(); c++) {
	harmonics = (int) list[c*el+3];
	list[c*el+0] = list[c*el+0]*dm_step + dm_low;
	list[c*el+1] = list[c*el+1]*(1.0/(sampling_time*nTimesamples*mod));
	//list[c*el+2] = white_noise_approximation(list[c*el+2], list[c*el+3], MSD);
	list[c*el+2] = (list[c*el+2] - MSD[2*harmonics])/(MSD[2*harmonics+1]);
	list[c*el+3] = list[c*el+3];
      }
    }
	
    /** \brief Rescales and then processes inBin data group. */
    void Rescale_Threshold_and_Process(float *MSD, Periodicity_inBin_Group *inBin_group, float sigma_cutoff, float mod) {
      float SNR;
      float dm_step       = inBin_group->Prange[rangeid].range.dm_step;
      float dm_low        = inBin_group->Prange[rangeid].range.dm_low;
      float sampling_time = inBin_group->Prange[rangeid].range.sampling_time;
      int nTimesamples    = inBin_group->Prange[rangeid].range.nTimesamples;
      int nPoints_before = size();
      int nPoints_after;

      std::vector<float> new_list;
      for(int c=0; c<(int)size(); c++) {
	float oldSNR = list[c*el+2];
	int harmonics = (int) list[c*el+3];
	//SNR = white_noise_approximation(list[c*el+2], list[c*el+3], MSD);
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
	
    void Allocate(AA_Periodicity_Plan *P_plan){
      MSD_interpolated_size = P_plan->nHarmonics;
      MSD_DIT_size = ((int) floorf(log2f((float)P_plan->nHarmonics))) + 2;
      size_t t_input_plane_size = P_plan->input_plane_size;
		
      if ( cudaSuccess != cudaMalloc((void **) &d_one_A,  sizeof(float)*t_input_plane_size )) printf("Periodicity Allocation error! d_one_A\n");
      if ( cudaSuccess != cudaMalloc((void **) &d_two_B,  sizeof(float)*2*t_input_plane_size )) printf("Periodicity Allocation error! d_two_B\n");
      if ( cudaSuccess != cudaMalloc((void **) &d_half_C,  sizeof(float)*t_input_plane_size/2 )) printf("Periodicity Allocation error! d_spectra_Real\n");
		
      if ( cudaSuccess != cudaMalloc((void **) &d_power_harmonics, sizeof(ushort)*t_input_plane_size )) printf("Periodicity Allocation error! d_harmonics\n");
      if ( cudaSuccess != cudaMalloc((void **) &d_interbin_harmonics, sizeof(ushort)*t_input_plane_size )) printf("Periodicity Allocation error! d_harmonics\n");
		
      if ( cudaSuccess != cudaMalloc((void**) &gmem_power_peak_pos, 1*sizeof(int)) )  printf("Periodicity Allocation error! gmem_power_peak_pos\n");
      if ( cudaSuccess != cudaMalloc((void**) &gmem_interbin_peak_pos, 1*sizeof(int)) )  printf("Periodicity Allocation error! gmem_interbin_peak_pos\n");
		
      if ( cudaSuccess != cudaMalloc((void**) &d_MSD, sizeof(float)*MSD_interpolated_size*2)) {printf("Periodicity Allocation error! d_MSD\n");}
      printf("DEBUG--------------------------------- MSD_DIT_size %d -------- MSD_PARTIAL_SIZE %d\n", MSD_DIT_size, MSD_PARTIAL_SIZE);
      if ( cudaSuccess != cudaMalloc((void**) &d_previous_partials, sizeof(float)*MSD_DIT_size*MSD_PARTIAL_SIZE)) {printf("Periodicity Allocation error! d_previous_partials\n");}
      if ( cudaSuccess != cudaMalloc((void**) &d_all_blocks, sizeof(float)*P_plan->max_total_MSD_blocks*MSD_PARTIAL_SIZE)) {printf("Periodicity Allocation error! d_MSD\n");}
		
      if ( cudaSuccess != cudaMalloc((void **) &cuFFT_workarea, P_plan->cuFFT_workarea_size) ) {printf("Periodicity Allocation error! cuFFT_workarea\n");}
    }
	
    void Reset_MSD(){
      printf("  Resetting MSD variables\n");
      cudaMemset(d_MSD, 0, MSD_interpolated_size*2*sizeof(float));
      cudaMemset(d_previous_partials, 0, MSD_DIT_size*MSD_PARTIAL_SIZE*sizeof(float));
    }
	
    void Reset_Candidate_List(){
      printf("         Resetting candidate lists\n");
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




  void Create_DD_plan(Dedispersion_Plan *D_plan, int nRanges, float *dm_low, float *dm_high, float *dm_step, int *inBin, int nTimesamples, int const*const ndms, float sampling_time){
    int max_nDMs = 0;
    for(int f=0; f<nRanges; f++){
      Dedispersion_Range trange;
      trange.Assign(dm_low[f], dm_high[f], dm_step[f], inBin[f], nTimesamples/inBin[f], ndms[f], sampling_time);
      D_plan->DD_range.push_back(trange);
      if(ndms[f]>max_nDMs) max_nDMs = ndms[f];
    }
	
    D_plan->nProcessedTimesamples = nTimesamples;
    D_plan->max_nDMs = max_nDMs;
    D_plan->sampling_time = sampling_time;
  }

  int Calculate_max_nDMs_in_memory(size_t max_nTimesamples, size_t max_nDMs, size_t memory_available, float multiple_float, float multiple_ushort) {
    size_t max_nDMs_in_memory, itemp;

    size_t memory_per_DM = ((max_nTimesamples+2)*(multiple_float*sizeof(float) + multiple_ushort*sizeof(ushort)));
    printf("   Memory required for one DM is %0.3f MB available memory is %0.3f MB %zu\n", memory_per_DM/(1024.0*1024.0), (memory_available*0.95)/(1024.0*1024.0), memory_available);
    max_nDMs_in_memory = (memory_available*0.98)/((max_nTimesamples+2)*(multiple_float*sizeof(float) + multiple_ushort*sizeof(ushort))); // 1 for real input real, 2 for complex output, 2 for complex cuFFT, 1 for peaks + 1 ushort
    if((max_nDMs+PHS_NTHREADS)<max_nDMs_in_memory) { //if we can fit all DM trials from largest range into memory then we need to find nearest higher multiple of PHS_NTHREADS
      itemp = (int)(max_nDMs/PHS_NTHREADS);
      if((max_nDMs%PHS_NTHREADS)>0) itemp++;
      max_nDMs_in_memory = itemp*PHS_NTHREADS;
    }
    itemp = (int)(max_nDMs_in_memory/PHS_NTHREADS); // if we cannot fit all DM trials from largest range into memory we find nearest lower multiple of PHS_NTHREADS
    max_nDMs_in_memory = itemp*PHS_NTHREADS;

    return(max_nDMs_in_memory);
  }

  void Create_Periodicity_Plan(AA_Periodicity_Plan *P_plan, Dedispersion_Plan *DD_plan, size_t max_nDMs_in_memory) {
    int oldinBin, last;
    int max_nTimesamples;
    int max_total_MSD_blocks;
	
    int nearest = (int)floorf(log2f((float) DD_plan->nProcessedTimesamples));
    P_plan->max_nTimesamples     = (int) powf(2.0, (float) nearest);
    P_plan->max_nDMs             = DD_plan->max_nDMs;
    P_plan->max_total_MSD_blocks = 0;
	
    oldinBin = 0;
    last = 0;
    max_nTimesamples     = 0;
    for(int f=0; f<(int)DD_plan->DD_range.size(); f++) {
      int localinBin = DD_plan->DD_range[f].inBin;

      if(oldinBin != localinBin) {
	P_plan->inBin_group.push_back( *(new Periodicity_inBin_Group) );
	last = (int) P_plan->inBin_group.size() - 1;
	P_plan->inBin_group[last].total_MSD_blocks = 0;
	oldinBin = localinBin;
      }
		
      Periodicity_Range Prange(DD_plan->DD_range[f], max_nDMs_in_memory, f, P_plan->inBin_group[last].total_MSD_blocks);
      P_plan->inBin_group[last].Prange.push_back(Prange);
      P_plan->inBin_group[last].total_MSD_blocks += Prange.total_MSD_blocks;
      if(max_nTimesamples<Prange.range.nTimesamples) max_nTimesamples = Prange.range.nTimesamples;
    }
	
    max_total_MSD_blocks = 0;
    for(int f=0; f<(int)P_plan->inBin_group.size();f++){
      if(max_total_MSD_blocks<P_plan->inBin_group[f].total_MSD_blocks) max_total_MSD_blocks=P_plan->inBin_group[f].total_MSD_blocks;
    }
	
    P_plan->max_total_MSD_blocks = max_total_MSD_blocks;
    P_plan->max_nTimesamples     = max_nTimesamples;
  }

  void Find_Periodicity_Plan(int *max_nDMs_in_memory, AA_Periodicity_Plan *P_plan, Dedispersion_Plan *DD_plan, size_t memory_available){
    size_t memory_allocated, memory_for_data;
    size_t input_plane_size;
	
    int nearest = (int)floorf(log2f((float) DD_plan->nProcessedTimesamples));
    size_t t_max_nTimesamples = (size_t) powf(2.0, (float) nearest);
    size_t t_max_nDMs = (size_t) DD_plan->max_nDMs;
    size_t t_max_nDMs_in_memory;
	
    float multiple_float  = 5.5;
    float multiple_ushort = 2.0; 
	
    printf("--------------------------------- FINDING PERIODICITY PLAN -------------------------------\n");
    memory_for_data = memory_available;
    do {
      printf("   Memory_for_data: %zu = %0.3f MB\n", memory_for_data, (float) memory_for_data/(1024.0*1024.0));
      t_max_nDMs_in_memory = Calculate_max_nDMs_in_memory(t_max_nTimesamples, t_max_nDMs, memory_for_data, multiple_float, multiple_ushort);
      if(t_max_nDMs_in_memory==0) { printf("Error not enough memory for periodicity search!"); exit(1); }
      input_plane_size = (t_max_nTimesamples+2)*t_max_nDMs_in_memory;
      memory_allocated = input_plane_size*multiple_float*sizeof(float) + multiple_ushort*input_plane_size*sizeof(ushort);
      printf("   Maximum number of DM trials which fit into memory is %zu;\n   Maximum time trials %zu;\n   Input plane size: %0.2f MB;\n", (size_t) t_max_nDMs_in_memory, t_max_nTimesamples, (((float) input_plane_size*sizeof(float))/(1024.0*1024.0)) );
		

      Create_Periodicity_Plan(P_plan, DD_plan, t_max_nDMs_in_memory);
      //max_total_blocks = Find_Max_total_blocks(P_plan);
      //------------------- Additional memory for MSD --------------------------------------
      int nDecimations = ((int) floorf(log2f((float)P_plan->nHarmonics))) + 2;
      size_t additional_data_size = P_plan->max_total_MSD_blocks*MSD_RESULTS_SIZE*sizeof(float) + nDecimations*2*MSD_RESULTS_SIZE*sizeof(float) + P_plan->nHarmonics*2*sizeof(float) + 2*sizeof(int);
      memory_allocated = memory_allocated + additional_data_size;
      printf("   Memory available for the component: %0.3f MB (%zu bytes)\n", (float) memory_available/(1024.0*1024.0), memory_available);
      printf("   Memory allocated by the component: %0.3f MB (%zu bytes)\n", (float) memory_allocated/(1024.0*1024.0), memory_allocated);

      if(t_max_nTimesamples!=(size_t)P_plan->max_nTimesamples) printf("Interesting!\n");
		
      if(memory_allocated>memory_available) {
	printf("Not enough memory for given configuration of periodicity plan. Calculating new plan...\n");
	memory_for_data = memory_for_data - additional_data_size;
	P_plan->clear();
      }
    } while(memory_allocated > memory_available);
    printf("------------------------------------------------------------------------------------------\n");
	
    P_plan->max_nDMs_in_memory = (int) t_max_nDMs_in_memory;
    P_plan->input_plane_size = (P_plan->max_nTimesamples+2)*t_max_nDMs_in_memory;
    //----------------------------------------------
    //----------> Exact calculation of cuFFT workarea size
    cufftHandle plan_input;
    cufftResult cufft_error;
    size_t cuFFT_workarea_size;
    cufftCreate(&plan_input);
    cufftSetAutoAllocation(plan_input, false);
    cufft_error = cufftMakePlan1d(plan_input, P_plan->max_nTimesamples, CUFFT_R2C, t_max_nDMs_in_memory, &cuFFT_workarea_size);
    if (CUFFT_SUCCESS != cufft_error){
      printf("CUFFT error: %d", cufft_error);
    }
    printf("Workarea size:%0.2f MB\n", (float) cuFFT_workarea_size/(1024.0*1024.0));
    cufftDestroy(plan_input);
    P_plan->cuFFT_workarea_size = cuFFT_workarea_size;
    //----------------------------------------------<
    *max_nDMs_in_memory = (int) t_max_nDMs_in_memory;
  }

  void Copy_data_for_periodicity_search(float *d_one_A, float **dedispersed_data, Periodicity_Batch *batch){ //TODO add "cudaStream_t stream1"
	int nStreams = 16;
	cudaStream_t stream_copy[16];
	cudaError_t e;
	float *h_small_dedispersed_data;
	size_t data_size = batch->nTimesamples*sizeof(float);
	cudaMallocHost((void **) &h_small_dedispersed_data, nStreams*data_size);

	for (int i = 0; i < nStreams; i++){
		e = cudaStreamCreate(&stream_copy[i]);
		if (e != cudaSuccess) {
			LOG(log_level::error, "Could not create streams in periodicity (" + std::string(cudaGetErrorString(e)) + ")");
		}
	}

	size_t stream_offset = batch->nTimesamples;

    #pragma omp parallel for num_threads(nStreams) shared(h_small_dedispersed_data, data_size, d_one_A, stream_copy, stream_offset)
    for(int ff=0; ff<batch->nDMs_per_batch; ff++){
      int id_stream = omp_get_thread_num();
      memcpy(h_small_dedispersed_data + id_stream*stream_offset, dedispersed_data[batch->DM_shift + ff], data_size);
//      e = cudaMemcpy( &d_one_A[ff*batch->nTimesamples], dedispersed_data[batch->DM_shift + ff], batch->nTimesamples*sizeof(float), cudaMemcpyHostToDevice);      
      e = cudaMemcpyAsync(&d_one_A[ff*batch->nTimesamples], h_small_dedispersed_data + id_stream*stream_offset, data_size, cudaMemcpyHostToDevice, stream_copy[id_stream]);      
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
	fprintf(stderr, "Error opening output file!\n");
	exit(0);
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
  void Periodicity_search(GPU_Memory_for_Periodicity_Search *gmem, aa_periodicity_strategy per_param, double *compute_time, size_t input_plane_size, int inBin, Periodicity_Batch *batch, std::vector<int> *h_boxcar_widths){ //TODO add "cudaStream_t stream1"
    int local_max_list_size = (input_plane_size)/4;
	
    float *d_dedispersed_data, *d_FFT_complex_output, *d_frequency_power, *d_frequency_interbin, *d_frequency_power_CT, *d_frequency_interbin_CT, *d_power_SNR, *d_interbin_SNR, *d_power_list, *d_interbin_list, *d_MSD_workarea;
    d_dedispersed_data      = gmem->d_one_A;
    d_FFT_complex_output    = gmem->d_two_B;
    d_MSD_workarea          = gmem->d_two_B;
    d_frequency_power       = gmem->d_half_C;
    d_frequency_interbin    = gmem->d_one_A;
    d_frequency_power_CT    = &gmem->d_two_B[0];
    d_frequency_interbin_CT = &gmem->d_two_B[input_plane_size];
    d_power_SNR             = gmem->d_half_C;
    d_interbin_SNR          = gmem->d_one_A;
    d_power_list            = &gmem->d_two_B[0];
    d_interbin_list         = &gmem->d_two_B[input_plane_size];
	
    int t_inBin             = inBin;
    int t_nTimesamples      = batch->nTimesamples;
    int t_nDMs_per_batch    = batch->nDMs_per_batch;
    int t_DM_shift          = batch->DM_shift;
	
    aa_gpu_timer timer;
	
    //checkCudaErrors(cudaGetLastError());
	
    //---------> cuFFT
    timer.Start();
    cufftHandle plan_input;
    cufftResult cufft_error;
	
    size_t cuFFT_workarea_size;
	
    cufftCreate(&plan_input);
    cufftSetAutoAllocation(plan_input, false);
	
    cufft_error = cufftMakePlan1d(plan_input, t_nTimesamples, CUFFT_R2C, t_nDMs_per_batch, &cuFFT_workarea_size);
    if (CUFFT_SUCCESS != cufft_error){
      printf("CUFFT error: %d", cufft_error);
    }
    cufftSetWorkArea(plan_input, gmem->cuFFT_workarea);

    cufft_error = cufftExecR2C(plan_input, (cufftReal *)d_dedispersed_data, (cufftComplex *)d_FFT_complex_output);
    if ( cufft_error != CUFFT_SUCCESS) printf("CUFFT error: %d\n", cufft_error);
	
    cufftDestroy(plan_input);

    timer.Stop();
    printf("         -> cuFFT took %f ms\n", timer.Elapsed());
    (*compute_time) = (*compute_time) + timer.Elapsed();
    //---------<
	
    //checkCudaErrors(cudaGetLastError());
	
    //---------> Calculate powers and interbinning
    timer.Start();
    //simple_power_and_interbin( (float2 *) d_two_B, d_half_C, d_one_A, nTimesamples, t_nDMs_per_batch);
    simple_power_and_interbin( (float2 *) d_FFT_complex_output, d_frequency_power, d_frequency_interbin, t_nTimesamples, t_nDMs_per_batch);
    timer.Stop();
    printf("         -> Calculation of powers and interbining took %f ms\n", timer.Elapsed());
    (*compute_time) = (*compute_time) + timer.Elapsed();
    //---------<
	
	//---------> Spectrum whitening
	timer.Start();
	cudaStream_t stream; stream = NULL;
	spectrum_whitening_SGP1(d_frequency_power, (t_nTimesamples>>1), t_nDMs_per_batch, stream);
	spectrum_whitening_SGP1(d_frequency_interbin, t_nTimesamples, t_nDMs_per_batch, stream);
    timer.Stop();
    printf("         -> Performing spectrum whitening took %f ms\n", timer.Elapsed());
    (*compute_time) = (*compute_time) + timer.Elapsed();	
	//---------<
    //-----------------------------------------------------------------------------------
    //if(i==0 && dm==0 && export_data) Export_data_in_range(d_half_C, nTimesamples/2, t_nDMs_per_batch, "power_data", dm_step[i], dm_low[i], tsamp, DM_shift);
    //-----------------------------------------------------------------------------------
    // Do Exports
    //-----------------------------------------------------------------------------------
    //if(i==0 && dm==0 && export_data) Export_data_in_range(d_one_A, nTimesamples, t_nDMs_per_batch, "Interbin_data", dm_step[i], dm_low[i], tsamp DM_shift);
    //-----------------------------------------------------------------------------------
	
    //checkCudaErrors(cudaGetLastError());
	
    //---------> Mean and StDev on powers
    timer.Start();
    bool perform_continuous = false;
#ifdef PS_REUSE_MSD_WITHIN_INBIN
    perform_continuous = true;
#endif


#ifdef OLD_PERIODICITY
    if(per_param.enable_outlier_rejection()){
#ifdef PS_REUSE_MSD_WITHIN_INBIN
      MSD_outlier_rejection_grid(gmem->d_MSD, d_frequency_power, gmem->d_previous_partials, gmem->d_all_blocks, &batch->MSD_conf, per_param.sigma_constant());
#else
      MSD_outlier_rejection(gmem->d_MSD, d_frequency_power, gmem->d_all_blocks, &batch->MSD_conf, per_param.sigma_constant());
#endif
    }
    else {
#ifdef PS_REUSE_MSD_WITHIN_INBIN
      MSD_normal_continuous(gmem->d_MSD, d_frequency_power, gmem->d_previous_partials, gmem->d_all_blocks, &batch->MSD_conf);
#else
      MSD_normal(gmem->d_MSD, d_frequency_power, gmem->d_all_blocks, &batch->MSD_conf);
#endif
    }
#else
    double total_time, dit_time, MSD_time;
    MSD_plane_profile(gmem->d_MSD, d_frequency_power, gmem->d_previous_partials, d_MSD_workarea, true, (t_nTimesamples>>1), t_nDMs_per_batch, h_boxcar_widths, 0, 0, 0, per_param.sigma_constant(), per_param.enable_msd_baseline_noise(), perform_continuous, &total_time, &dit_time, &MSD_time);
    printf("    MSD time: Total: %f ms; DIT: %f ms; MSD: %f ms;\n", total_time, dit_time, MSD_time);
#endif
	

	
    timer.Stop();
    printf("         -> MSD took %f ms\n", timer.Elapsed());
    (*compute_time) = (*compute_time) + timer.Elapsed();
    //---------<
	
    //checkCudaErrors(cudaGetLastError());
	
    //---------> Corner turn
    timer.Start();
    //corner_turn_SM(d_half_C, &d_two_B[0], (nTimesamples>>1), t_nDMs_per_batch);
    //corner_turn_SM(d_one_A, &d_two_B[input_plane_size], nTimesamples, t_nDMs_per_batch);
    corner_turn_SM(d_frequency_power, d_frequency_power_CT, (t_nTimesamples>>1), t_nDMs_per_batch);
    corner_turn_SM(d_frequency_interbin, d_frequency_interbin_CT, t_nTimesamples, t_nDMs_per_batch);
    timer.Stop();
    printf("         -> corner turn took %f ms\n", timer.Elapsed());
    (*compute_time) = (*compute_time) + timer.Elapsed();
    //---------<
	
    //checkCudaErrors(cudaGetLastError());
	
    //---------> Harmonic summing
    timer.Start();

    /** \todo Re-enable the if statement after debugging. */
    //if(per_param.nHarmonics() > 1) {
      
#ifdef OLD_PERIODICITY
      periodicity_simple_harmonic_summing_old(d_frequency_power_CT, d_power_SNR, gmem->d_power_harmonics, gmem->d_MSD, (t_nTimesamples>>1), t_nDMs_per_batch, per_param.nHarmonics());
      periodicity_simple_harmonic_summing_old(d_frequency_interbin_CT, d_interbin_SNR, gmem->d_interbin_harmonics, gmem->d_MSD, t_nTimesamples, t_nDMs_per_batch, per_param.nHarmonics());
#else
      periodicity_simple_harmonic_summing(d_frequency_power_CT, d_power_SNR, gmem->d_power_harmonics, gmem->d_MSD, (t_nTimesamples>>1), t_nDMs_per_batch, per_param.nHarmonics());
      periodicity_simple_harmonic_summing(d_frequency_interbin_CT, d_interbin_SNR, gmem->d_interbin_harmonics, gmem->d_MSD, t_nTimesamples, t_nDMs_per_batch, per_param.nHarmonics());
#endif
      //}
    timer.Stop();
    printf("         -> harmonic summing took %f ms\n", timer.Elapsed());
    (*compute_time) = (*compute_time) + timer.Elapsed();
    //---------<
	
    //checkCudaErrors(cudaGetLastError());
	
    //---------> Peak finding
    timer.Start();
    if(per_param.candidate_algorithm()){
      //-------------- Thresholding
#ifdef OLD_PERIODICITY
      Threshold_for_periodicity_old(d_power_SNR, gmem->d_power_harmonics, d_power_list, gmem->gmem_power_peak_pos, gmem->d_MSD, per_param.sigma_cutoff(), t_nDMs_per_batch, (t_nTimesamples>>1), t_DM_shift, t_inBin, local_max_list_size);
      Threshold_for_periodicity_old(d_interbin_SNR, gmem->d_interbin_harmonics, d_interbin_list, gmem->gmem_interbin_peak_pos, gmem->d_MSD, per_param.sigma_cutoff(), t_nDMs_per_batch, t_nTimesamples, t_DM_shift, t_inBin, local_max_list_size);
#else
      Threshold_for_periodicity(d_power_SNR, gmem->d_power_harmonics, d_power_list, gmem->gmem_power_peak_pos, gmem->d_MSD, per_param.sigma_cutoff(), t_nDMs_per_batch, (t_nTimesamples>>1), t_DM_shift, t_inBin, local_max_list_size);
      Threshold_for_periodicity(d_interbin_SNR, gmem->d_interbin_harmonics, d_interbin_list, gmem->gmem_interbin_peak_pos, gmem->d_MSD, per_param.sigma_cutoff(), t_nDMs_per_batch, t_nTimesamples, t_DM_shift, t_inBin, local_max_list_size);
#endif
      //-------------- Thresholding
    }
    else {
      //-------------- Peak finding
#ifdef OLD_PERIODICITY
      Peak_find_for_periodicity_search_old(d_power_SNR, gmem->d_power_harmonics, d_power_list, (t_nTimesamples>>1), t_nDMs_per_batch, per_param.sigma_cutoff(), local_max_list_size, gmem->gmem_power_peak_pos, gmem->d_MSD, t_DM_shift, t_inBin);
      Peak_find_for_periodicity_search_old(d_interbin_SNR, gmem->d_interbin_harmonics, d_interbin_list, t_nTimesamples, t_nDMs_per_batch, per_param.sigma_cutoff(), local_max_list_size, gmem->gmem_interbin_peak_pos, gmem->d_MSD, t_DM_shift, t_inBin);
#else
      Peak_find_for_periodicity_search(d_power_SNR, gmem->d_power_harmonics, d_power_list, (t_nTimesamples>>1), t_nDMs_per_batch, per_param.sigma_cutoff(), local_max_list_size, gmem->gmem_power_peak_pos, gmem->d_MSD, t_DM_shift, t_inBin);
      Peak_find_for_periodicity_search(d_interbin_SNR, gmem->d_interbin_harmonics, d_interbin_list, t_nTimesamples, t_nDMs_per_batch, per_param.sigma_cutoff(), local_max_list_size, gmem->gmem_interbin_peak_pos, gmem->d_MSD, t_DM_shift, t_inBin);
#endif
      //-------------- Peak finding
    }
    timer.Stop();
    printf("         -> Peak finding took %f ms\n", timer.Elapsed());
    (*compute_time) = (*compute_time) + timer.Elapsed();
    //---------<
	
    //checkCudaErrors(cudaGetLastError());
  }

  void Export_Data_To_File(std::vector<Candidate_List> candidates, const char *filename) {
    FILE *fp_out;
    if((fp_out = fopen(filename, "wb")) == NULL) {
      fprintf(stderr, "Error opening output file!\n");
      exit(0);
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
  void GPU_periodicity(int range, int nsamp, int max_ndms, int processed, float sigma_cutoff, float ***output_buffer, int const*const ndms, int *inBin, float *dm_low, float *dm_high, float *dm_step, float tsamp, int nHarmonics, bool candidate_algorithm, bool enable_msd_baseline_noise, float OR_sigma_multiplier) {
    // processed = maximum number of time-samples through out all ranges
    // nTimesamples = number of time-samples in given range 'i'
    // TODO:
    //      ->Use callbacks for power calculation
    //      ->Solve peak finding problem from batch to batch (we do not want to find peaks on shared borders)
    //      ->max_ndms is possibly the same thing as max_nDMs! Investigate.
    //      ->check is zero-th element does not mess up statistics
    //      ->prepare data on the host before copying them to device
    //		->There must by one inBin group per stream. We cannot have batches or Pranges per stream because MSD. :(
	
    printf("\n");
    printf("------------ STARTING PERIODICITY SEARCH ------------\n\n");
    // Creating periodicity parameters object (temporary, it should be moved elsewhere)
    aa_periodicity_plan per_param_plan(sigma_cutoff, OR_sigma_multiplier, nHarmonics, 0, candidate_algorithm, enable_msd_baseline_noise); // \warning The periodicity plan uses a hardcoded number (this also used to be the case for the (now deprecated) Periodicity_parameters class.
    aa_periodicity_strategy per_param(per_param_plan);
    per_param.print_parameters();	
	
    std::vector<int> h_boxcar_widths; h_boxcar_widths.resize(nHarmonics); 
    for(int f=0; f<nHarmonics; f++) h_boxcar_widths[f]=f+1;
	
    // Creating DDrange vector (temporary, it should be moved elsewhere)
    Dedispersion_Plan DD_plan;
    Create_DD_plan(&DD_plan, range, dm_low, dm_high, dm_step, inBin, processed, ndms, tsamp);
	
    printf("--------------------------\n");
    DD_plan.print();
    printf("--------------------------\n");
	
    // Determining available memory (temporary, it should be moved elsewhere)
    size_t memory_available,total_mem;
    cudaMemGetInfo(&memory_available,&total_mem);
    printf("   Memory available:%0.3f MB \n", ((float) memory_available)/(1024.0*1024.0) );
	
	
    //----------------------------------------------------------------------------
    //----------- Finding Periodicity plan
    int max_nDMs_in_memory;
    AA_Periodicity_Plan P_plan;
    P_plan.nHarmonics = nHarmonics;
    Find_Periodicity_Plan(&max_nDMs_in_memory, &P_plan, &DD_plan, memory_available);
    size_t input_plane_size = P_plan.input_plane_size;
    P_plan.print();
    printf("\n");
	
    printf("--------------------------------- GPU ALLOCATIONS -------------------------------\n");
    //-----------------------------------------------------------------------------
    //-----------------------------------------------------------------------------
    //--------> Allocation of GPU memory
    GPU_Memory_for_Periodicity_Search GPU_memory;
	
    GPU_memory.Allocate(&P_plan);
	
    float h_MSD[P_plan.nHarmonics*2];
    //----------------------------<
    //-----------------------------------------------------------------------------
    //-----------------------------------------------------------------------------
    printf("---------------------------------------------------------------------------------\n");
    printf("\n");
	
    //checkCudaErrors(cudaGetLastError());
	
    // Timing
    double Total_periodicity_time = 0, Total_calc_time = 0, calc_time_per_range = 0, Total_copy_time = 0, copy_time_per_range = 0;
    aa_gpu_timer timer, periodicity_timer;
	
    periodicity_timer.Start();
	
    printf("--------------------------------- PERIODICITY START -------------------------------\n");
    //---------------------------------------------------------------
    //---------------------------------------------------------------
    for(int p=0; p<(int) P_plan.inBin_group.size(); p++) {
      printf("--------------------------------------------------\n");
      printf("inBin range: %d\n", p);
      std::vector<Candidate_List> PowerCandidates;
      std::vector<Candidate_List> InterbinCandidates;

#ifdef PS_REUSE_MSD_THROUGH_INBIN
      if(p>0){
	float h_MSD_partials[MSD_PARTIAL_SIZE];
	int ratio = P_plan.inBin_group[p].Prange[0].range.inBin/P_plan.inBin_group[p-1].Prange[0].range.inBin;
	ratio = ratio*4;
				
	GPU_memory.Get_MSD_partials(h_MSD_partials);
	printf("This groups inBin: %d; Previous groups inBin: %d; ratio: %d;\n", P_plan.inBin_group[p].Prange[0].range.inBin, P_plan.inBin_group[p-1].Prange[0].range.inBin, P_plan.inBin_group[p].Prange[0].range.inBin/P_plan.inBin_group[p-1].Prange[0].range.inBin);
				
	printf("Old MSD: [%f;%f]\n", h_MSD_partials[0]/h_MSD_partials[2], sqrt(h_MSD_partials[1]/h_MSD_partials[2]));
	h_MSD_partials[0] = h_MSD_partials[0]/ratio;
	h_MSD_partials[1] = h_MSD_partials[1]/(ratio*ratio);
	printf("New MSD: [%f;%f]\n", h_MSD_partials[0]/h_MSD_partials[2], sqrt(h_MSD_partials[1]/h_MSD_partials[2]));
				
	GPU_memory.Set_MSD_partials(h_MSD_partials);
      }
#else
      GPU_memory.Reset_MSD();
#endif
		
      //checkCudaErrors(cudaGetLastError());

      for(int r=0; r<(int) P_plan.inBin_group[p].Prange.size(); r++) {
	printf("  Prange: %d\n", r);
	P_plan.inBin_group[p].Prange[r].print();
			
	for(int b=0; b<(int)P_plan.inBin_group[p].Prange[r].batches.size(); b++) {
	  printf("      Batch: %d\n", b);
	  P_plan.inBin_group[p].Prange[r].batches[b].print();
				
	  GPU_memory.Reset_Candidate_List();
				
	  //checkCudaErrors(cudaGetLastError());

	  //---------> Copy input data to the device
	  timer.Start();
	  Copy_data_for_periodicity_search(GPU_memory.d_one_A, output_buffer[P_plan.inBin_group[p].Prange[r].rangeid], &P_plan.inBin_group[p].Prange[r].batches[b]);
	  timer.Stop();
	  printf("         => Copy of data to device took %f ms\n", timer.Elapsed() );
	  copy_time_per_range = copy_time_per_range + timer.Elapsed();
	  //---------<
				
				
	  //checkCudaErrors(cudaGetLastError());
	  // TODO: cudaStreamSynchronize(stream);
				
				
	  //---------> Periodicity search
	  Periodicity_search(&GPU_memory, per_param, &calc_time_per_range, input_plane_size, P_plan.inBin_group[p].Prange[r].range.inBin, &P_plan.inBin_group[p].Prange[r].batches[b], &h_boxcar_widths);
	  //---------<
				
				
	  //checkCudaErrors(cudaGetLastError());
	  // TODO: cudaStreamSynchronize(stream);
				
				
	  //---------> Copy candidates to the host
	  timer.Start();
				
	  int last_entry;
	  int nPowerCandidates = GPU_memory.Get_Number_of_Power_Candidates();
	  PowerCandidates.push_back(*(new Candidate_List(r)));
	  last_entry = PowerCandidates.size()-1;
	  printf("         -> POWER: Total number of peaks found in this range is %d;\n", nPowerCandidates);
	  PowerCandidates[last_entry].Allocate(nPowerCandidates);
	  cudaError_t e = cudaMemcpy( &PowerCandidates[last_entry].list[0], &GPU_memory.d_two_B[0], nPowerCandidates*Candidate_List::el*sizeof(float), cudaMemcpyDeviceToHost);
	  
	  if(e != cudaSuccess) {
	    LOG(log_level::error, "Could not cudaMemcpy in aa_device_periods.cu (" + std::string(cudaGetErrorString(e)) + ")");
	  }
	  
	  //checkCudaErrors(cudaGetLastError());
				
	  int nInterbinCandidates = GPU_memory.Get_Number_of_Interbin_Candidates();
	  InterbinCandidates.push_back(*(new Candidate_List(r)));
	  last_entry = InterbinCandidates.size()-1;
	  printf("         -> INTERBIN: Total number of peaks found in this range is %d;\n", nInterbinCandidates);
	  InterbinCandidates[last_entry].Allocate(nInterbinCandidates);
	  e = cudaMemcpy( &InterbinCandidates[last_entry].list[0], &GPU_memory.d_two_B[input_plane_size], nInterbinCandidates*Candidate_List::el*sizeof(float), cudaMemcpyDeviceToHost);
	  
	  if(e != cudaSuccess) {
	    LOG(log_level::error, "Could not cudaMemcpy in aa_device_periods.cu (" + std::string(cudaGetErrorString(e)) + ")");
	  }
	  
	  timer.Stop();
	  printf("         => Copy of candidates took %f ms\n", timer.Elapsed() );
	  copy_time_per_range = copy_time_per_range + timer.Elapsed();
	  //---------<
				
#ifndef PS_REUSE_MSD_WITHIN_INBIN
	  GPU_memory.Get_MSD(h_MSD);
					
	  PowerCandidates[last_entry].Process(h_MSD, &P_plan.inBin_group[p],  1.0);
	  InterbinCandidates[last_entry].Process(h_MSD, &P_plan.inBin_group[p], 2.0);
#endif
				
	  printf("\n");
	} //batches

      } // ranges

#ifdef PS_REUSE_MSD_WITHIN_INBIN
      GPU_memory.Get_MSD(h_MSD);

			

      for(int f=0; f<(int)PowerCandidates.size(); f++) {
	//PowerCandidates[f].Rescale_Threshold_and_Process(h_MSD, &P_plan.inBin_group[p], per_param.sigma_cutoff, 1.0);
	PowerCandidates[f].Process(h_MSD, &P_plan.inBin_group[p], 1.0);
      }

      for(int f=0; f<(int)InterbinCandidates.size(); f++) {
	//InterbinCandidates[f].Rescale_Threshold_and_Process(h_MSD, &P_plan.inBin_group[p], per_param.sigma_cutoff, 2.0);
	InterbinCandidates[f].Process(h_MSD, &P_plan.inBin_group[p], 2.0);
      }
#endif

      // Export of the candidate list for inBin range;
      char filename[100];
      float range_dm_low  = P_plan.inBin_group[p].Prange[0].range.dm_low;
      float range_dm_high = P_plan.inBin_group[p].Prange[P_plan.inBin_group[p].Prange.size()-1].range.dm_high; 
      sprintf(filename, "fourier-dm_%.2f-%.2f.dat", range_dm_low, range_dm_high);
      Export_Data_To_File(PowerCandidates, filename);
      sprintf(filename, "fourier_inter-dm_%.2f-%.2f.dat", range_dm_low, range_dm_high);
      Export_Data_To_File(InterbinCandidates, filename);
		
      // Cleanup
      for(int f=0; f<(int) PowerCandidates.size(); f++) PowerCandidates[f].list.clear();
      PowerCandidates.clear();
      for(int f=0; f<(int) InterbinCandidates.size(); f++) InterbinCandidates[f].list.clear();
      InterbinCandidates.clear();
		
		
      printf("     -----------------------\n");
      printf("     -> This range calculation time: %f ms\n", calc_time_per_range);
      printf("     -> This range copy time:        %f ms\n", copy_time_per_range);
      printf("\n");
      Total_calc_time = Total_calc_time + calc_time_per_range;
      calc_time_per_range = 0;
      Total_copy_time = Total_copy_time + copy_time_per_range;
      copy_time_per_range = 0;
    } // inBin ranges
    printf("-----------------------------------------------------------------------------------\n");
	
    periodicity_timer.Stop();
    Total_periodicity_time = periodicity_timer.Elapsed();
    printf("\nTimer:\n");
    printf("Total calculation time: %f ms\n", Total_calc_time);
    printf("Total copy time:        %f ms\n", Total_copy_time);
    printf("Total periodicity time: %f ms\n", Total_periodicity_time);	

    cudaDeviceSynchronize();
  };

} //namespace astroaccelerate
