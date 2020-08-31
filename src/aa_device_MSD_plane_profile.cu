#include <string>
#include <fstream>
#include "aa_device_single_FIR.hpp"
#include "aa_bin_gpu.hpp"
#include "aa_device_MSD.hpp"
#include "aa_gpu_timer.hpp"
#include "aa_log.hpp"


#include "aa_device_MSD_shared_kernel_functions.hpp"
//#define MSD_PLANE_DEBUG
//#define MSD_PLANE_EXPORT

namespace astroaccelerate {

  /**
   * \struct MSD_Data
   * \brief Struct to hold MSD processing parameters.
   */
  struct MSD_Data {
    int width;
    double mean;
    double sd;
  };



  //---------------------------------------------------------------
  //------------- MSD plane profile

  void Do_MSD_normal(float *d_MSD, float *d_input, float *d_MSD_workarea, int nTimesamples, int nDMs, int offset, float OR_sigma_multiplier, int enable_outlier_rejection){
    MSD_Configuration conf(nTimesamples, nDMs, offset, 0);
    if(enable_outlier_rejection){
      MSD_outlier_rejection(d_MSD, d_input, d_MSD_workarea, &conf, OR_sigma_multiplier);
    }
    else {
      MSD_normal(d_MSD, d_input, d_MSD_workarea, &conf);
    }
  }

  void Do_MSD_continuous(float *d_MSD, float *d_input, float *d_previous_partials, float *d_MSD_workarea, int nTimesamples, int nDMs, int offset, float OR_sigma_multiplier, int enable_outlier_rejection){
    MSD_Configuration conf(nTimesamples, nDMs, offset, 0);
    if(enable_outlier_rejection){
      MSD_outlier_rejection_continuous(d_MSD, d_input, d_previous_partials, d_MSD_workarea, &conf, OR_sigma_multiplier);
    }
    else {
      MSD_normal_continuous(d_MSD, d_input, d_previous_partials, d_MSD_workarea, &conf);
    }
  }

  inline void Do_MSD(float *d_MSD, float *d_input, float *d_previous_partials, float *d_MSD_workarea, int nTimesamples, int nDMs, int offset, float OR_sigma_multiplier, int enable_outlier_rejection, bool perform_continuous) {
    if(perform_continuous) Do_MSD_continuous(d_MSD, d_input, d_previous_partials, d_MSD_workarea, nTimesamples, nDMs, offset, OR_sigma_multiplier, enable_outlier_rejection);
    else Do_MSD_normal(d_MSD, d_input, d_MSD_workarea, nTimesamples, nDMs, offset, OR_sigma_multiplier, enable_outlier_rejection);
  }


  void MSD_plane_profile_debug(float *d_MSD, int DIT_value, int nTimesamples, int MSD_pos){
    float h_MSD[MSD_RESULTS_SIZE];
    cudaError_t e = cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
    
    if(e != cudaSuccess) {
      LOG(log_level::error, "Could not cudaMemcpy in aa_device_MSD_plane_profile.cu (" + std::string(cudaGetErrorString(e)) + ")");
	exit(25);
    }
    
    printf("    DiT:%d; nTimesamples:%d; decimated_timesamples:%d; MSD_pos: %d; MSD:[%f; %f; %f]\n", (int) DIT_value, (int) nTimesamples, (int) (nTimesamples>>1), (int) MSD_pos, h_MSD[0], h_MSD[1], h_MSD[2]);
  }


  void MSD_of_input_plane(float *d_MSD_DIT, std::vector<int> *h_MSD_DIT_widths, float *d_input_data, float *d_MSD_DIT_previous, float *d_sudy, float *d_lichy, float *d_MSD_workarea, size_t nTimesamples, size_t nDMs, int nDecimations, int max_width_performed, float OR_sigma_multiplier, int enable_outlier_rejection, bool high_memory, bool perform_continuous, double *total_time, double *dit_time, double *MSD_time){
    aa_gpu_timer timer, total_timer;
    double t_dit_time=0, t_MSD_time=0;
    int nRest;
    size_t decimated_timesamples;
    int DIT_value;

	
    total_timer.Start();
    //----------------------------------------------------------------------------------------
    //-------- DIT = 1
    DIT_value = 1;
	
    timer.Start();
    Do_MSD(d_MSD_DIT, d_input_data, d_MSD_DIT_previous, d_MSD_workarea, nTimesamples, nDMs, 0, OR_sigma_multiplier, enable_outlier_rejection, perform_continuous);
    timer.Stop();	t_MSD_time += timer.Elapsed();
    h_MSD_DIT_widths->push_back(DIT_value);
#ifdef MSD_PLANE_DEBUG
    printf("a------------------ MSD_plane_profile DEBUG - MSD_DIT calculation ------------\n");
    printf("    MSD format: [ mean ; StDev ; nElements ]\n");
    MSD_plane_profile_debug(d_MSD_DIT, DIT_value, nTimesamples, 0);
#endif
    //----------------------------------------------------------------------------------------
	
    //checkCudaErrors(cudaGetLastError());
	
    //----------------------------------------------------------------------------------------
    //-------- DIT = 2
    DIT_value = DIT_value*2;
	
    if(high_memory){
      //printf("High memory: DIT=2 is not split\n");
      timer.Start();
      nRest = GPU_DiT_v2_wrapper(d_input_data, d_lichy, nDMs, nTimesamples);
      decimated_timesamples = (nTimesamples>>1);
      timer.Stop();	t_dit_time += timer.Elapsed();
		
      timer.Start();
      Do_MSD(&d_MSD_DIT[MSD_RESULTS_SIZE], d_lichy, &d_MSD_DIT_previous[MSD_RESULTS_SIZE], d_MSD_workarea, decimated_timesamples, nDMs, nRest, OR_sigma_multiplier, enable_outlier_rejection, perform_continuous);
      timer.Stop();	t_MSD_time += timer.Elapsed();
      h_MSD_DIT_widths->push_back(DIT_value);
		
#ifdef MSD_PLANE_DEBUG
      MSD_plane_profile_debug(&d_MSD_DIT[MSD_RESULTS_SIZE], DIT_value, decimated_timesamples, 1);
#endif
		
      timer.Start();
      nRest = GPU_DiT_v2_wrapper(d_lichy, d_sudy, nDMs, decimated_timesamples);
      timer.Stop();	t_dit_time += timer.Elapsed();
    }
    else {
      //printf("Low memory: DIT=2 is split in two\n");
      // First decimation is split into two parts, that way we can lower the memory requirements for MSD_plane_profile
      // First half of the decimation
      int nDMs_half = (nDMs>>1);
      timer.Start();
      nRest = GPU_DiT_v2_wrapper(d_input_data, d_lichy, nDMs_half, nTimesamples);
      decimated_timesamples = (nTimesamples>>1);
      timer.Stop();	t_dit_time += timer.Elapsed();

      timer.Start();
      Do_MSD_continuous(&d_MSD_DIT[MSD_RESULTS_SIZE], d_lichy, &d_MSD_DIT[2*MSD_RESULTS_SIZE], d_MSD_workarea, decimated_timesamples, nDMs_half, nRest, OR_sigma_multiplier, enable_outlier_rejection);
      timer.Stop();	t_MSD_time += timer.Elapsed();
		
      timer.Start();
      nRest = GPU_DiT_v2_wrapper(d_lichy, d_sudy, nDMs_half, decimated_timesamples);
      timer.Stop();	t_dit_time += timer.Elapsed();
		
      // second half of the decimation
      timer.Start();
      nRest = GPU_DiT_v2_wrapper(&d_input_data[nDMs_half*nTimesamples], d_lichy, nDMs_half, nTimesamples);
      decimated_timesamples = (nTimesamples>>1);
      timer.Stop();	t_dit_time += timer.Elapsed();

      timer.Start();
      Do_MSD_continuous(&d_MSD_DIT[MSD_RESULTS_SIZE], d_lichy, &d_MSD_DIT[2*MSD_RESULTS_SIZE], d_MSD_workarea, decimated_timesamples, nDMs_half, nRest, OR_sigma_multiplier, enable_outlier_rejection);
      timer.Stop();	t_MSD_time += timer.Elapsed();
      h_MSD_DIT_widths->push_back(DIT_value);
		
      timer.Start();
      nRest = GPU_DiT_v2_wrapper(d_lichy, &d_sudy[nDMs_half*(decimated_timesamples>>1)], nDMs_half, decimated_timesamples);
      timer.Stop();	t_dit_time += timer.Elapsed();
		
#ifdef MSD_PLANE_DEBUG
      MSD_plane_profile_debug(&d_MSD_DIT[MSD_RESULTS_SIZE], DIT_value, decimated_timesamples, 1);
#endif
    }
	
    decimated_timesamples = (nTimesamples>>2);
    DIT_value = DIT_value*2;
	
    timer.Start();
    Do_MSD(&d_MSD_DIT[2*MSD_RESULTS_SIZE], d_sudy, &d_MSD_DIT_previous[2*MSD_RESULTS_SIZE], d_MSD_workarea, decimated_timesamples, nDMs, nRest, OR_sigma_multiplier, enable_outlier_rejection, perform_continuous);
    timer.Stop();	t_MSD_time += timer.Elapsed();
    h_MSD_DIT_widths->push_back(DIT_value);	
	
#ifdef MSD_PLANE_DEBUG
    MSD_plane_profile_debug(&d_MSD_DIT[2*MSD_RESULTS_SIZE], DIT_value, decimated_timesamples, 2);
#endif
    //----------------------------------------------------------------------------------------
	
    //checkCudaErrors(cudaGetLastError());
	
    //----------------------------------------------------------------------------------------
    //-------- DIT > 3
    int f = 0;
    size_t last_admited_f = 0;
    size_t last_admited_DIT_value = 0;
    int switch_to_boxcar = 0;
	
    for(f=3; f<=nDecimations; f++){
      timer.Start();
      DIT_value = DIT_value*2;
      if(decimated_timesamples<20) switch_to_boxcar = 1;
      if(DIT_value<=max_width_performed && switch_to_boxcar == 0){
	if(f%2==0){
	  timer.Start();
	  //                            in      out :|
	  nRest = GPU_DiT_v2_wrapper(d_lichy, d_sudy, nDMs, decimated_timesamples);
	  timer.Stop();	t_dit_time += timer.Elapsed();
				
	  if(nRest<0) break;
	  decimated_timesamples = (decimated_timesamples>>1);

	  timer.Start();
	  Do_MSD(&d_MSD_DIT[f*MSD_RESULTS_SIZE], d_sudy, &d_MSD_DIT_previous[f*MSD_RESULTS_SIZE], d_MSD_workarea, decimated_timesamples, nDMs, nRest, OR_sigma_multiplier, enable_outlier_rejection, perform_continuous);
	  timer.Stop();	t_MSD_time += timer.Elapsed();
	}
	else {
	  timer.Start();
	  //                            in      out :|
	  nRest = GPU_DiT_v2_wrapper(d_sudy, d_lichy, nDMs, decimated_timesamples);
	  timer.Stop();	t_dit_time += timer.Elapsed();
				
	  if(nRest<0) break;
	  decimated_timesamples = (decimated_timesamples>>1);

	  timer.Start();
	  Do_MSD(&d_MSD_DIT[f*MSD_RESULTS_SIZE], d_lichy, &d_MSD_DIT_previous[f*MSD_RESULTS_SIZE], d_MSD_workarea, decimated_timesamples, nDMs, nRest, OR_sigma_multiplier, enable_outlier_rejection, perform_continuous);
	  timer.Stop();	t_MSD_time += timer.Elapsed();
	}
	h_MSD_DIT_widths->push_back(DIT_value);
	last_admited_f = f;
	last_admited_DIT_value = DIT_value;
			
#ifdef MSD_PLANE_DEBUG
	MSD_plane_profile_debug(&d_MSD_DIT[f*MSD_RESULTS_SIZE], DIT_value, decimated_timesamples, f);
#endif
      }
      //checkCudaErrors(cudaGetLastError());
    }
    //----------------------------------------------------------------------------------------
	
    //checkCudaErrors(cudaGetLastError());

    //----------------------------------------------------------------------------------------
    //-------- Boxcar for last boxcar width if needed
    if(DIT_value>max_width_performed && switch_to_boxcar == 0){
      f = last_admited_f + 1;
      DIT_value = (DIT_value>>1);
      decimated_timesamples = nTimesamples/DIT_value;
      int nTaps = max_width_performed/DIT_value + 1;
		
      if(f%2==0){
	timer.Start();
	//               in        out   :|
	nRest = PPF_L1(d_lichy, d_sudy, nDMs, decimated_timesamples, nTaps);
	timer.Stop();	t_dit_time += timer.Elapsed();

	//checkCudaErrors(cudaGetLastError());
			
	timer.Start();
			
	Do_MSD(&d_MSD_DIT[f*MSD_RESULTS_SIZE], d_sudy, &d_MSD_DIT_previous[f*MSD_RESULTS_SIZE], d_MSD_workarea, decimated_timesamples, nDMs, nRest, OR_sigma_multiplier, enable_outlier_rejection, perform_continuous);
	timer.Stop();	t_MSD_time += timer.Elapsed();
      }
      else {
	timer.Start();
	//               in       out      :|
	nRest = PPF_L1(d_sudy, d_lichy, nDMs, decimated_timesamples, nTaps);
	timer.Stop();	t_dit_time += timer.Elapsed();
			
	//checkCudaErrors(cudaGetLastError());
			
	timer.Start();
	Do_MSD(&d_MSD_DIT[f*MSD_RESULTS_SIZE], d_lichy, &d_MSD_DIT_previous[f*MSD_RESULTS_SIZE], d_MSD_workarea, decimated_timesamples, nDMs, nRest, OR_sigma_multiplier, enable_outlier_rejection, perform_continuous);
	timer.Stop();	t_MSD_time += timer.Elapsed();
      }
      h_MSD_DIT_widths->push_back(DIT_value*nTaps);

#ifdef MSD_PLANE_DEBUG
      printf("    Performing additional boxcar: nTaps: %d; max_width_performed: %d; DIT_value/2: %d;\n", nTaps, max_width_performed, DIT_value);
      MSD_plane_profile_debug(&d_MSD_DIT[f*MSD_RESULTS_SIZE], DIT_value*nTaps, decimated_timesamples, f);
#endif
    }
    //----------------------------------------------------------------------------------------
	
    //checkCudaErrors(cudaGetLastError());
	
    //----------------------------------------------------------------------------------------
    //-------- Data are decimated too much switching to boxcar filters	
    if(switch_to_boxcar){
#ifdef MSD_PLANE_DEBUG
      printf("    **** Data are decimated too much switching to boxcar filters! ****\n");
#endif
		
      f = last_admited_f + 1; // to determine what array contain input data
      DIT_value = last_admited_DIT_value; // each element of the input array is sum of 'last_admited_DIT_value' elements of initial array
      float *input_pointer, *output_pointer; // setup temporary pointers to input and output based on where is the last output. 
      if(f%2==0){
	input_pointer = d_lichy;
	output_pointer = d_sudy;
      }
      else {
	input_pointer = d_sudy;
	output_pointer = d_lichy;
      }
		
      decimated_timesamples = nTimesamples/DIT_value;
      DIT_value = DIT_value*2;
      while(DIT_value<=max_width_performed){
	int nTaps = DIT_value/last_admited_DIT_value;
	timer.Start();
	//               in       out      :|
	nRest = PPF_L1(input_pointer, output_pointer, nDMs, decimated_timesamples, nTaps);
	timer.Stop();	t_dit_time += timer.Elapsed();
			
	//checkCudaErrors(cudaGetLastError());
			
	timer.Start();
	Do_MSD(&d_MSD_DIT[f*MSD_RESULTS_SIZE], output_pointer, &d_MSD_DIT_previous[f*MSD_RESULTS_SIZE], d_MSD_workarea, decimated_timesamples, nDMs, nRest, OR_sigma_multiplier, enable_outlier_rejection, perform_continuous);
	timer.Stop();	t_MSD_time += timer.Elapsed();
			
#ifdef MSD_PLANE_DEBUG
	printf("    Performing boxcar instead of DIT: nTaps: %d; max_width_performed: %d; DIT_value: %d;\n", nTaps, max_width_performed, DIT_value);
	MSD_plane_profile_debug(&d_MSD_DIT[f*MSD_RESULTS_SIZE], DIT_value, decimated_timesamples, f);
#endif
			
	h_MSD_DIT_widths->push_back(DIT_value);
	f++;
	DIT_value = DIT_value*2;
      }
		
      DIT_value = (DIT_value>>1);
      if(max_width_performed!=DIT_value){
	int nTaps = max_width_performed/last_admited_DIT_value + 1;
	timer.Start();
	//               in       out      :|
	nRest = PPF_L1(input_pointer, output_pointer, nDMs, decimated_timesamples, nTaps);
	timer.Stop();	t_dit_time += timer.Elapsed();
			
	//checkCudaErrors(cudaGetLastError());
			
	timer.Start();
	Do_MSD(&d_MSD_DIT[f*MSD_RESULTS_SIZE], output_pointer, &d_MSD_DIT_previous[f*MSD_RESULTS_SIZE], d_MSD_workarea, decimated_timesamples, nDMs, nRest, OR_sigma_multiplier, enable_outlier_rejection, perform_continuous);
	timer.Stop();	t_MSD_time += timer.Elapsed();
			
	h_MSD_DIT_widths->push_back(last_admited_DIT_value*nTaps);
			
#ifdef MSD_PLANE_DEBUG
	printf("    Performing additional boxcar: nTaps: %d; max_width_performed: %d; DIT_value: %lu;\n", nTaps, max_width_performed, last_admited_DIT_value*nTaps);
	MSD_plane_profile_debug(&d_MSD_DIT[f*MSD_RESULTS_SIZE], last_admited_DIT_value*nTaps, decimated_timesamples, f);
#endif
      }
		
    }
    //----------------------------------------------------------------------------------------
	
#ifdef MSD_PLANE_DEBUG
    printf("------------------------------------------------------<\n");
#endif
	
    //checkCudaErrors(cudaGetLastError());
	
    total_timer.Stop();
    (*total_time) = total_timer.Elapsed();
    (*dit_time) = t_dit_time;
    (*MSD_time) = t_MSD_time;
	
#ifdef GPU_PARTIAL_TIMER
    printf("    MSD of input plane: Total time: %f ms; DiT time: %f ms; MSD time: %f ms;\n", (*total_time), (*dit_time), (*MSD_time));
#endif
  }


  void MSD_Interpolate_linear(float *mean, float *StDev, float desired_width, float *h_MSD_DIT, std::vector<int> *h_MSD_DIT_widths){
    int MSD_DIT_size = h_MSD_DIT_widths->size();
    int position = (int) floorf(log2f((float) desired_width));
	
    float width1 = h_MSD_DIT_widths->operator[](position);
    float mean1 = h_MSD_DIT[(position)*MSD_RESULTS_SIZE];
    float StDev1 = h_MSD_DIT[(position)*MSD_RESULTS_SIZE +1];
	
    if(position == MSD_DIT_size-1 && width1==(int) desired_width) {
      (*mean) = mean1;
      (*StDev) = StDev1;
    }
    else {
      float width2 = h_MSD_DIT_widths->operator[](position+1);
      float distance_in_width = width2 - width1;
		
      float mean2 = h_MSD_DIT[(position+1)*MSD_RESULTS_SIZE];
      float distance_in_mean = mean2 - mean1;
		
      float StDev2 = h_MSD_DIT[(position+1)*MSD_RESULTS_SIZE +1];
      float distance_in_StDev = StDev2 - StDev1;
		
      (*mean) = mean1 + (distance_in_mean/distance_in_width)*((float) desired_width - width1);
      (*StDev) = StDev1 + (distance_in_StDev/distance_in_width)*((float) desired_width - width1);
		
#ifdef MSD_PLANE_DEBUG
      printf("    width:[%f;%f]; mean:[%f;%f]; sd:[%f;%f]\n",width1, width2, mean1, mean2, StDev1, StDev2);
      printf("    distances width %f; mean: %f; StDef: %f\n", distance_in_width, distance_in_mean, distance_in_StDev);
      printf("    desired_width: %f; mean: %f; StDev: %f;\n", desired_width, (*mean), (*StDev));
#endif
    }
  }


  void MSD_Interpolate_square(float *mean, float *StDev, float desired_width, float *h_MSD_DIT, std::vector<int> *h_MSD_DIT_widths){
    int MSD_DIT_size = h_MSD_DIT_widths->size();
    int position = (int) floorf(log2f((float) desired_width));
	
    if(position == MSD_DIT_size-2) position--;
    if(position == MSD_DIT_size-1 && h_MSD_DIT_widths->operator[](position)==(int) desired_width) {
      (*mean)  = h_MSD_DIT[(position)*MSD_RESULTS_SIZE];
      (*StDev) = h_MSD_DIT[(position)*MSD_RESULTS_SIZE +1];
    }
    else {
      float w = desired_width;
		
      float w0 = h_MSD_DIT_widths->operator[](position);
      float mean0  = h_MSD_DIT[(position)*MSD_RESULTS_SIZE];
      float StDev0 = h_MSD_DIT[(position)*MSD_RESULTS_SIZE +1];
		
      float w1 = h_MSD_DIT_widths->operator[](position+1);
      float mean1  = h_MSD_DIT[(position+1)*MSD_RESULTS_SIZE];
      float StDev1 = h_MSD_DIT[(position+1)*MSD_RESULTS_SIZE +1];
		
      float w2 = h_MSD_DIT_widths->operator[](position+2);
      float mean2  = h_MSD_DIT[(position+2)*MSD_RESULTS_SIZE];
      float StDev2 = h_MSD_DIT[(position+2)*MSD_RESULTS_SIZE +1];
		
      float a0 = ((w - w1)*(w - w2))/((w0 - w1)*(w0 - w2));
      float a1 = ((w - w0)*(w - w2))/((w1 - w0)*(w1 - w2));
      float a2 = ((w - w0)*(w - w1))/((w2 - w0)*(w2 - w1));
		
      (*mean)  = a0*mean0 + a1*mean1 + a2*mean2;
      (*StDev) = a0*StDev0 + a1*StDev1 + a2*StDev2;
    }
  }


  void MSD_Export_plane(const char *filename, float *h_MSD_DIT, std::vector<int> *h_MSD_DIT_widths, float *h_MSD_interpolated, std::vector<int> *h_boxcar_widths, int max_width_performed) {
    char str[200];
    std::ofstream FILEOUT;
    int MSD_INTER_SIZE = 2;
	
    sprintf(str,"%s_DIT.dat", filename);
    FILEOUT.open (str, std::ofstream::out);
    size_t h_MSD_DIT_widths_size = h_MSD_DIT_widths->size();
    for(size_t f=0; f<h_MSD_DIT_widths_size; f++){
      FILEOUT << (int) h_MSD_DIT_widths->operator[](f) << " " << h_MSD_DIT[f*MSD_RESULTS_SIZE] << " " << h_MSD_DIT[f*MSD_RESULTS_SIZE + 1] << std::endl;
    }
    FILEOUT.close();
	
    sprintf(str,"%s_Interpolated.dat", filename);
    FILEOUT.open (str, std::ofstream::out);
    size_t h_boxcar_widths_size = h_boxcar_widths->size();
    for(size_t f=0; f< h_boxcar_widths_size; f++){
      if(h_boxcar_widths->operator[](f)<=max_width_performed)
	FILEOUT << (int) h_boxcar_widths->operator[](f) << " " << h_MSD_interpolated[f*MSD_INTER_SIZE] << " " << h_MSD_interpolated[f*MSD_INTER_SIZE + 1] << std::endl;
    }
    FILEOUT.close();
  }


  void MSD_Interpolate_values(float *d_MSD_interpolated, float *d_MSD_DIT, std::vector<int> *h_MSD_DIT_widths, int nMSDs, std::vector<int> *h_boxcar_widths, int max_width_performed, const char *filename){
#ifdef GPU_PARTIAL_TIMER
    aa_gpu_timer timer;
    timer.Start();
#endif
	
    //	float *h_MSD_DIT, *h_MSD_interpolated;
    int nWidths = (int) h_boxcar_widths->size();
    //	h_MSD_DIT = new float[nMSDs*MSD_RESULTS_SIZE];
    //	h_MSD_interpolated = new float[nWidths*MSD_INTER_SIZE];

    // adding memory for the interpolate kernel
    int MSD_DIT_size = h_MSD_DIT_widths->size();
    int *d_MSD_DIT_widths;
    int *d_boxcar;
    cudaError_t e = cudaMalloc((void **) &d_MSD_DIT_widths, sizeof(int)*MSD_DIT_size);

    if(e != cudaSuccess) {
      LOG(log_level::error, "Could not allocate memory in aa_device_MSD_plane_profile.cu (" + std::string(cudaGetErrorString(e)) + ")");
    }
    
    e = cudaMemcpy(d_MSD_DIT_widths, &h_MSD_DIT_widths->operator[](0), sizeof(int)*MSD_DIT_size,cudaMemcpyHostToDevice);

    if(e != cudaSuccess) {
      LOG(log_level::error, "Could not cudaMemcpy (MSD_DIT_width) in aa_device_MSD_plane_profile.cu (" + std::string(cudaGetErrorString(e)) + ")");
    }
    
    e = cudaMalloc((void **) &d_boxcar, sizeof(int)*nWidths);

    if(e != cudaSuccess) {
      LOG(log_level::error, "Could not allocate memory (d_boxcar) in aa_device_MSD_plane_profile.cu (" + std::string(cudaGetErrorString(e)) + ")");
    }
    
    e = cudaMemcpy(d_boxcar, &h_boxcar_widths->operator[](0), sizeof(int)*nWidths,cudaMemcpyHostToDevice);
    
    if(e != cudaSuccess) {
      LOG(log_level::error, "Could not cudaMemcpy (d_boxcar) in aa_device_MSD_plane_profile.cu (" + std::string(cudaGetErrorString(e)) + ")");
    }

	
    //	checkCudaErrors(cudaMemcpy(h_MSD_DIT, d_MSD_DIT, nMSDs*MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost));
	
#ifdef MSD_PLANE_DEBUG
    printf("------------------ MSD_plane_profile DEBUG - Linear interpolation ------------\n");
#endif
    call_kernel_MSD_GPU_Interpolate_linear(1, nWidths,
					   d_MSD_DIT, d_MSD_interpolated, d_MSD_DIT_widths, h_MSD_DIT_widths->size(), d_boxcar, max_width_performed);

    //	for(int f=0; f<nWidths; f++){
    //		if(h_boxcar_widths->operator[](f)<=max_width_performed) {
    //			float mean, StDev;
    //			MSD_Interpolate_linear(&mean, &StDev, (float) h_boxcar_widths->operator[](f), h_MSD_DIT, h_MSD_DIT_widths);
    //			h_MSD_interpolated[f*MSD_INTER_SIZE] = mean;
    //			h_MSD_interpolated[f*MSD_INTER_SIZE+1] = StDev;
    //		}
    //	}
#ifdef MSD_PLANE_DEBUG
    printf("-----------------------------------------------------------------------------<\n");
#endif
	
#ifdef MSD_PLANE_EXPORT
    float *h_MSD_DIT, *h_MSD_interpolated;
    h_MSD_DIT = new float[nMSDs*MSD_RESULTS_SIZE];
    h_MSD_interpolated = new float[nWidths*MSD_RESULTS_SIZE];
    e = cudaMemcpy(h_MSD_DIT, d_MSD_DIT, nMSDs*MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);

    if(e != cudaSuccess) {
      LOG(log_level::error, "Could not cudaMemcpy in aa_device_MSD_plane_profile.cu (" + std::string(cudaGetErrorString(e)) + ")");
    }
    
    e = cudaMemcpy(h_MSD_interpolated, d_MSD_interpolated, nWidths*MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
    
    if(e != cudaSuccess) {
      LOG(log_level::error, "Could not cudaMemcpy in aa_device_MSD_plane_profile.cu (" + std::string(cudaGetErrorString(e)) + ")");
    }
    
    MSD_Export_plane(filename, h_MSD_DIT, h_MSD_DIT_widths, h_MSD_interpolated, h_boxcar_widths, max_width_performed);
    delete[] h_MSD_DIT;
    delete[] h_MSD_interpolated;
#endif
	
    //	checkCudaErrors(cudaMemcpy(d_MSD_interpolated, h_MSD_interpolated, nWidths*MSD_INTER_SIZE*sizeof(float), cudaMemcpyHostToDevice));
	
    //	delete[] h_MSD_DIT;
    //	delete[] h_MSD_interpolated;
	
#ifdef GPU_PARTIAL_TIMER
    timer.Stop();
    printf("    Interpolation step took %f ms;\n", timer.Elapsed());
#endif
  }

  //-------------------------------------------------------------------------<

  void Get_MSD_plane_profile_memory_requirements(size_t *MSD_profile_size_in_bytes, size_t *MSD_DIT_profile_size_in_bytes, size_t *workarea_size_in_bytes, size_t primary_dimension, size_t secondary_dimension, std::vector<int> *boxcar_widths) {
    // temporary work area for decimations. We need 2*1/4 = 1/2.
    size_t t_wsib = (primary_dimension*secondary_dimension*sizeof(float))/2;
	
    // temporary storage for MSD values of decimated input data
    int max_boxcar_width = boxcar_widths->operator[](boxcar_widths->size()-1);
    int nDecimations = ((int) floorf(log2f((float)max_boxcar_width))) + 2;
    t_wsib = t_wsib + nDecimations*MSD_RESULTS_SIZE*sizeof(float);
	
    // temporary storage for calculation of MSD. We have to choose the maximum from all possible variants.
    size_t decimated_pd = primary_dimension;
    int max_nBlocks = 0;
    for(int f=0; f<nDecimations; f++){
      MSD_Configuration conf(decimated_pd, secondary_dimension, 0, 0);
      if(conf.nBlocks_total>max_nBlocks) max_nBlocks = conf.nBlocks_total;
      decimated_pd = (decimated_pd>>1);
    }
    t_wsib = t_wsib + max_nBlocks*MSD_PARTIAL_SIZE*sizeof(float);

#ifdef MSD_PLANE_DEBUG
    printf("Data size primary dim: %zu; secondary dim: %zu;\n", primary_dimension, secondary_dimension);
    printf("Temporary storage for data: %zu bytes = %zu floats;\n", t_wsib, t_wsib/4);
    printf("Size of DIT MSDs: %d elements = %d float = %lu bytes\n", nDecimations, nDecimations*MSD_RESULTS_SIZE, nDecimations*MSD_RESULTS_SIZE*sizeof(float));
    printf("Max MSD blocks: %d blocks = %d float = %lu bytes\n", max_nBlocks, max_nBlocks*MSD_PARTIAL_SIZE, max_nBlocks*MSD_PARTIAL_SIZE*sizeof(float));
#endif
	
    (*workarea_size_in_bytes) = t_wsib;
    (*MSD_profile_size_in_bytes) = boxcar_widths->size()*MSD_RESULTS_SIZE*sizeof(float);
    (*MSD_DIT_profile_size_in_bytes) = nDecimations*MSD_RESULTS_SIZE*sizeof(float);
  }



  void MSD_plane_profile(float *d_MSD_interpolated, float *d_input_data, float *d_MSD_DIT_previous, float *workarea, bool high_memory, size_t primary_dimension, size_t secondary_dimension, std::vector<int> *boxcar_widths, float tstart, float dm_low, float dm_high, float OR_sigma_multiplier, int enable_outlier_rejection, bool perform_continuous, double *total_time, double *dit_time, double *MSD_time){
    int boxcar_widths_size = (int) boxcar_widths->size();
    int max_boxcar_width = boxcar_widths->operator[](boxcar_widths_size-1);
    int nDecimations = ((int) floorf(log2f((float)max_boxcar_width))) + 1;
    int nDIT_widths = nDecimations + 1;
    std::vector<int> h_MSD_DIT_widths;
	
    size_t datasize = primary_dimension*secondary_dimension;
    float *d_sudy, *d_lichy, *d_MSD_DIT, *d_MSD_workarea;
    d_sudy = workarea;
    d_lichy = &workarea[datasize/4];
    if(high_memory) {
      d_MSD_DIT = &workarea[datasize/4 + datasize/2];
      d_MSD_workarea = &workarea[datasize/4 + datasize/2 + (nDecimations+1)*MSD_RESULTS_SIZE];
    }
    else {
      d_MSD_DIT = &workarea[datasize/2];
      d_MSD_workarea = &workarea[datasize/2 + (nDecimations+1)*MSD_RESULTS_SIZE];
    }
	
    cudaMemset((void*) d_MSD_DIT, 0, (nDecimations+1)*MSD_RESULTS_SIZE*sizeof(float));
	
    MSD_of_input_plane(d_MSD_DIT, &h_MSD_DIT_widths, d_input_data, d_MSD_DIT_previous, d_sudy, d_lichy, d_MSD_workarea, primary_dimension, secondary_dimension, nDecimations, max_boxcar_width, OR_sigma_multiplier, enable_outlier_rejection, high_memory, perform_continuous, total_time, dit_time, MSD_time);
	
#ifdef MSD_PLANE_DEBUG
    printf("    Number of calculated MSD values: %d; number of interpolated MSD values: %d;\n",nDIT_widths, boxcar_widths_size);
#endif
	
    char filename[100];
    sprintf(filename,"MSD_plane_profile_i_test-t_%.2f-dm_%.2f-%.2f", tstart, dm_low, dm_high);
    MSD_Interpolate_values(d_MSD_interpolated, d_MSD_DIT, &h_MSD_DIT_widths, nDIT_widths, boxcar_widths, max_boxcar_width, filename);
  }

  //------------- MSD plane profile
  //---------------------------------------------------------------

  //---------------------------------------------------------------
  //------------- MSD plane profile boxcars
/*
  void Create_boxcar_MSD(float *d_data, size_t nTimesamples, size_t nDMs, std::vector<MSD_Data> *boxcar_MSD, std::vector<MSD_Data> *boxcar_MSD_BLN, int max_nTaps, int max_boxcar_width, float OR_sigma_multiplier){
    aa_gpu_timer timer;
    double total_time = 0;
    int nRest;
    MSD_Data mdtemp;
    float *d_boxcar, *d_MSD;
    float h_MSD[MSD_RESULTS_SIZE];
    cudaMalloc((void **) &d_boxcar, nTimesamples*nDMs*sizeof(float));
    cudaMalloc((void **) &d_MSD, MSD_RESULTS_SIZE*sizeof(float));
	
    timer.Start();
	
    MSD_normal(d_MSD, d_data, nTimesamples, nDMs, 0);
    cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
    mdtemp.width = 1; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
    boxcar_MSD->push_back(mdtemp);
	
    MSD_outlier_rejection(d_MSD, d_data, nTimesamples, nDMs, 0, OR_sigma_multiplier);
    cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
    mdtemp.width = 1; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
    boxcar_MSD_BLN->push_back(mdtemp);
	
    timer.Stop();
    total_time = total_time + timer.Elapsed();
    printf("DIT value: %d; took %f ms; Total time %fms\n", 1, timer.Elapsed(), total_time);
  
    for(int f=2; f<=max_nTaps; f++){
      if( ((int)nTimesamples-f+1)>0 ) {
	timer.Start();
			
	nRest = PD_FIR(d_data, d_boxcar, f, nDMs, nTimesamples);
			
	MSD_normal(d_MSD, d_boxcar, nTimesamples, nDMs, nRest);
	cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
	mdtemp.width = f; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
	boxcar_MSD->push_back(mdtemp);
			
	MSD_outlier_rejection(d_MSD, d_boxcar, nTimesamples, nDMs, nRest, OR_sigma_multiplier);
	cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
	mdtemp.width = f; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
	boxcar_MSD_BLN->push_back(mdtemp);
			
	timer.Stop();
	total_time = total_time + timer.Elapsed();
	printf("DIT value: %d; took %f ms; Total time %fms\n", f, timer.Elapsed(), total_time);
      }
    }
	
    //checkCudaErrors(cudaGetLastError());
	
    for(int f=130; f<=256 && f<max_boxcar_width; f+=4){
      printf("nTimesamples: %zu; f: %d; %zu\n", nTimesamples, f, nTimesamples-f+1);
      int itemp = (int) ((int)nTimesamples-f+1);
      if( itemp>0 ) {
	timer.Start();
			
	nRest=PPF_L1(d_data, d_boxcar, nDMs, nTimesamples, f);
			
	if(nRest>0){
	  MSD_normal(d_MSD, d_boxcar, nTimesamples, nDMs, nRest);
	  cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
	  mdtemp.width = f; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
	  boxcar_MSD->push_back(mdtemp);
				
	  MSD_outlier_rejection(d_MSD, d_boxcar, nTimesamples, nDMs, nRest, OR_sigma_multiplier);
	  cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
	  mdtemp.width = f; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
	  boxcar_MSD_BLN->push_back(mdtemp);
	}

	timer.Stop();
	total_time = total_time + timer.Elapsed();
	printf("DIT value: %d; took %f ms; Total time %fms\n", f, timer.Elapsed(), total_time);
      }
      //checkCudaErrors(cudaGetLastError());
    }
	
    //checkCudaErrors(cudaGetLastError());
	
    for(int f=272; f<=512 && f<max_boxcar_width; f+=16){
      printf("nTimesamples: %zu; f: %d; %zu\n", nTimesamples, f, nTimesamples-f+1);
      int itemp = (int) ((int)nTimesamples-f+1);
      if( itemp>0 ) {
	timer.Start();
			
	nRest=PPF_L1(d_data, d_boxcar, nDMs, nTimesamples, f);
			
	if(nRest>0){
	  MSD_normal(d_MSD, d_boxcar, nTimesamples, nDMs, nRest);
	  cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
	  mdtemp.width = f; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
	  boxcar_MSD->push_back(mdtemp);
				
	  MSD_outlier_rejection(d_MSD, d_boxcar, nTimesamples, nDMs, nRest, OR_sigma_multiplier);
	  cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
	  mdtemp.width = f; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
	  boxcar_MSD_BLN->push_back(mdtemp);
	}
			
	timer.Stop();
	total_time = total_time + timer.Elapsed();
	printf("DIT value: %d; took %f ms; Total time %fms\n", f, timer.Elapsed(), total_time);
      }
      //checkCudaErrors(cudaGetLastError());
    }
	
    //checkCudaErrors(cudaGetLastError());
	
    for(int f=544; f<=1024 && f<max_boxcar_width; f+=32){
      printf("nTimesamples: %zu; f: %d; %zu\n", nTimesamples, f, nTimesamples-f+1);
      int itemp = (int) ((int)nTimesamples-f+1);
      if( itemp>0 ) {
	timer.Start();
			
	nRest=PPF_L1(d_data, d_boxcar, nDMs, nTimesamples, f);

	if(nRest>0){
	  MSD_normal(d_MSD, d_boxcar, nTimesamples, nDMs, nRest);
	  cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
	  mdtemp.width = f; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
	  boxcar_MSD->push_back(mdtemp);
				
	  MSD_outlier_rejection(d_MSD, d_boxcar, nTimesamples, nDMs, nRest, OR_sigma_multiplier);
	  cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
	  mdtemp.width = f; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
	  boxcar_MSD_BLN->push_back(mdtemp);
	}
			
	timer.Stop();
	total_time = total_time + timer.Elapsed();
	printf("DIT value: %d; took %f ms; Total time %fms\n", f, timer.Elapsed(), total_time);
      }
      //checkCudaErrors(cudaGetLastError());
    }
	
    //checkCudaErrors(cudaGetLastError());

    for(int f=1088; f<=2048 && f<max_boxcar_width; f+=64){
      printf("nTimesamples: %zu; f: %d; %zu\n", nTimesamples, f, nTimesamples-f+1);
      int itemp = (int) ((int)nTimesamples-f+1);
      if( itemp>0 ) {
	timer.Start();
			
	nRest=PPF_L1(d_data, d_boxcar, nDMs, nTimesamples, f);

	if(nRest>0){
	  MSD_normal(d_MSD, d_boxcar, nTimesamples, nDMs, nRest);
	  cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
	  mdtemp.width = f; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
	  boxcar_MSD->push_back(mdtemp);
				
	  MSD_outlier_rejection(d_MSD, d_boxcar, nTimesamples, nDMs, nRest, OR_sigma_multiplier);
	  cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
	  mdtemp.width = f; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
	  boxcar_MSD_BLN->push_back(mdtemp);
	}
			
	timer.Stop();
	total_time = total_time + timer.Elapsed();
	printf("DIT value: %d; took %f ms; Total time %fms\n", f, timer.Elapsed(), total_time);
      }
      //checkCudaErrors(cudaGetLastError());
    }
	
    //checkCudaErrors(cudaGetLastError());

    for(int f=2176; f<=4096 && f<max_boxcar_width; f+=128){
      printf("nTimesamples: %zu; f: %d; %zu\n", nTimesamples, f, nTimesamples-f+1);
      int itemp = (int) ((int)nTimesamples-f+1);
      if( itemp>0 ) {
	timer.Start();
			
	nRest=PPF_L1(d_data, d_boxcar, nDMs, nTimesamples, f);

	if(nRest>0){		
	  MSD_normal(d_MSD, d_boxcar, nTimesamples, nDMs, nRest);
	  cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
	  mdtemp.width = f; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
	  boxcar_MSD->push_back(mdtemp);
				
	  MSD_outlier_rejection(d_MSD, d_boxcar, nTimesamples, nDMs, nRest, OR_sigma_multiplier);
	  cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
	  mdtemp.width = f; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
	  boxcar_MSD_BLN->push_back(mdtemp);
	}
			
	timer.Stop();
	total_time = total_time + timer.Elapsed();
	printf("DIT value: %d; took %f ms; Total time %fms\n", f, timer.Elapsed(), total_time);
      }
      //checkCudaErrors(cudaGetLastError());
    }
	
    //checkCudaErrors(cudaGetLastError());
	
    for(int f=4352; f<=8192 && f<max_boxcar_width; f+=256){
      printf("nTimesamples: %zu; f: %d; %zu\n", nTimesamples, f, nTimesamples-f+1);
      int itemp = (int) ((int)nTimesamples-f+1);
      if( itemp>0 ) {
	timer.Start();
			
	nRest=PPF_L1(d_data, d_boxcar, nDMs, nTimesamples, f);
			
	if(nRest>0){
	  MSD_normal(d_MSD, d_boxcar, nTimesamples, nDMs, nRest);
	  cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
	  mdtemp.width = f; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
	  boxcar_MSD->push_back(mdtemp);
				
	  MSD_outlier_rejection(d_MSD, d_boxcar, nTimesamples, nDMs, nRest, OR_sigma_multiplier);
	  cudaMemcpy(h_MSD, d_MSD, MSD_RESULTS_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
	  mdtemp.width = f; mdtemp.mean = h_MSD[0]; mdtemp.sd = h_MSD[1];
	  boxcar_MSD_BLN->push_back(mdtemp);
	}
			
	timer.Stop();
	total_time = total_time + timer.Elapsed();
	printf("DIT value: %d; took %f ms; Total time %fms\n", f, timer.Elapsed(), total_time);
      }
      //checkCudaErrors(cudaGetLastError());
    }
	
    //checkCudaErrors(cudaGetLastError());
    cudaError_t e = cudaFree(d_boxcar);

    if(e != cudaSuccess) {
      LOG(log_level::error, "Could not cudaFree in aa_device_MSD_plane_profile.cu (" + std::string(cudaGetErrorString(e)) + ")");
    }
    
    e = cudaFree(d_MSD);

    if(e != cudaSuccess) {
      LOG(log_level::error, "Could not cudaFree in aa_device_MSD_plane_profile.cu (" + std::string(cudaGetErrorString(e)) + ")");
    }
    
  }
*/


/*
  void MSD_plane_profile_boxcars(float *d_input_data, size_t nTimesamples, size_t nDMs, std::vector<int> *boxcar_widths, float OR_sigma_multiplier, float dm_low, float dm_high, float tstart){
    char filename[200];
    std::vector<MSD_Data> h_boxcar_MSD;
    std::vector<MSD_Data> h_boxcar_MSD_OR;
	
    size_t free_mem, total_mem, req_mem;
    cudaMemGetInfo(&free_mem,&total_mem);
    req_mem = nTimesamples*nDMs*sizeof(float);
    printf("Memory available: %f; Memory required: %f;\n", (double) free_mem/(1024.0*1024.0), ((double) req_mem)/(1024.0*1024.0));
    if(free_mem<req_mem) {
      printf("Not enough memory to perform the comparison!\n");
      return;
    }
	
    int boxcar_widths_size = boxcar_widths->size();
    int max_boxcar_width = boxcar_widths->operator[](boxcar_widths_size-1);
	
    Create_boxcar_MSD(d_input_data, nTimesamples, nDMs, &h_boxcar_MSD, &h_boxcar_MSD_OR, 128, max_boxcar_width, OR_sigma_multiplier);
	
    sprintf(filename,"MSD_boxcars_OR%f-t_%.2f-dm_%.2f-%.2f.dat", OR_sigma_multiplier, tstart, dm_low, dm_high);
	
    std::ofstream FILEOUT;
    FILEOUT.open (filename, std::ofstream::out);

    for(size_t f=0; f<h_boxcar_MSD.size(); f++){
      FILEOUT << (int) h_boxcar_MSD[f].width << " " << h_boxcar_MSD[f].mean << " " << h_boxcar_MSD[f].sd << " " << "3" << std::endl;
    }
    FILEOUT << std::endl;
    FILEOUT << std::endl;
    for(size_t f=0; f<h_boxcar_MSD_OR.size(); f++){
      FILEOUT << (int) h_boxcar_MSD_OR[f].width << " " << h_boxcar_MSD_OR[f].mean << " " << h_boxcar_MSD_OR[f].sd << " " << "4" << std::endl;
    }
    FILEOUT << std::endl;
    FILEOUT << std::endl;
	
    FILEOUT.close();
  }
*/

} //namespace astroaccelerate
  
//------------- MSD plane profile boxcars
//---------------------------------------------------------------




