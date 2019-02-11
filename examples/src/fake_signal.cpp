/**
 * Example code for linking against astro-accelerate library.
 *
 * Compile with: g++ -std=c++11 -I/path/to/astro-accelerate/include/ -L/path/to/astro-accelerate/build -Wl,-rpath,/path/to/astro-accelerate/build -lastroaccelerate -I/usr/local/cuda/include/ -L/usr/local/cuda/lib64 -Wl,-rpath,/usr/local/cuda/lib64 -I/usr/local/cuda-8.0/samples/common/inc/ -lcudart dedispersion.cpp -o test
 */

#include <iostream>

#include "aa_ddtr_plan.hpp"
#include "aa_ddtr_strategy.hpp"
#include "aa_filterbank_metadata.hpp"
#include "aa_permitted_pipelines_1.hpp"
#include "aa_device_info.hpp"
#include "aa_permitted_pipelines_2.hpp"
#include "aa_analysis_plan.hpp"
#include "aa_analysis_strategy.hpp"
#include "params.hpp"

#define MAX_VALUE 255

using namespace astroaccelerate;

void fake_generate_data(unsigned short *output, float *scale_factor, int *shifts_index, int width, int signal_start, int signal_max_pos, int nchans){
  for(int i = 0; i < width; i++){
//	std::cout << scale_factor[i] << std::endl;
	int time_pos = (i - signal_max_pos + signal_start)*nchans;
	for(int j = 0; j < nchans; j++){
		output[j + nchans*shifts_index[j] + time_pos] = 255; // (MAX_VALUE*scale_factor[i]);
	}
  }
}

void inverse_gaussian(float* output, int length, float sigma, int *max_pos){
        float wid = 2*sigma*pow(3*log(30),0.5);
        float step = wid*sigma/(length);
        float maximum =0.0f;
        for(int i = 0; i < length; i++){
                float x = (i+1)*step;
                output[i] = pow(1.0*sigma/(2*3.141592*pow(x,3)),0.5)*exp((-sigma*pow(x - 1.0,2))/(2*pow(1.0,2)*x));
                if (maximum <= output[i]) {
                        maximum = output[i];
                        *max_pos = i;
                }
        }
        // normalization part
        for (int i = 0; i < length; i++){
                output[i] = output[i]/maximum;
        }

}

int main() {
  aa_ddtr_plan ddtr_plan;
  ddtr_plan.add_dm(85, 95, 0.1, 1, 1); // Add dm_ranges: dm_low, dm_high, dm_step, inBin, outBin (unused).

  // Filterbank metadata
  // (Data description from "SIGPROC-v3.7 (Pulsar) Signal Processing Programs")
  const double tstart = 50000;
  const double total_bandwidth = 300.0f;
  const double tsamp = 6.4E-5;
  const double nbits = 8;
  unsigned int nsamples = 3.0/tsamp; // 3s of data in current tsamp at minimum; must be more than signal_start + maxshift
  const double fch1 = 1550;
  int nchans = 128;
  const double foff = -total_bandwidth/nchans;
  std::cout << "Nsamples " << nsamples << std::endl;

  // params needed by the fake signal function
  double dm_position = 90.0; // at what dm put the signal
  const int func_width = 1; // width of the signal in terms of # of samples
  int signal_start = 1.5/tsamp; // position of the signal in time
  int signal_max_pos = 0; // position of the peak; usefull for gaussian and inverse gaussian signal
  float* scale_factor;
  scale_factor = (float *)malloc(func_width*sizeof(float));
 
  std::cout << "Signal start at: " << signal_start << std::endl; 
  inverse_gaussian(scale_factor,func_width, 0.5, &signal_max_pos); // generate the scaling factors for the signal with sigma 0.5, also getting position of the max scale

//  for (int n = 0; n < func_width; n++)
//	std::cout << "Vector " << scale_factor[n] << std::endl;
  
  aa_filterbank_metadata metadata(tstart, tsamp, nbits, nsamples, fch1, foff, nchans);

  aa_device_info device_info;
  if(device_info.check_for_devices()) {
    std::cout << "NOTICE: Checked for devices." << std::endl;
  }
  else {
    std::cout << "ERROR: Could not find any devices." << std::endl;
  }

  aa_device_info::CARD_ID selected_card = 1;
  aa_device_info::aa_card_info selected_card_info;
  if(device_info.init_card(selected_card, selected_card_info)) {
    std::cout << "NOTICE: init_card complete." << std::endl;
  }
  else {
    std::cout << "ERROR: init_card incomplete." << std::endl;
  }

  aa_device_info::print_card_info(selected_card_info);
  
  const size_t free_memory = selected_card_info.free_memory; // Free memory on the GPU in bytes
  bool enable_analysis = false; 
  aa_ddtr_strategy strategy(ddtr_plan, metadata, free_memory, enable_analysis);

  if(!(strategy.ready())) {
    std::cout << "There was an error" << std::endl;
    return 0;
  }

  std::vector<float> dm_shifts;
  dm_shifts = strategy.dmshifts();
  int* shifts_index;
  shifts_index = (int*)malloc(sizeof(int)*nchans);
  for (int i = 0; i < nchans; i++){
	shifts_index[i] = floor(dm_shifts[i]*dm_position/tsamp);
  }

  FILE *fp_out3;
  if (( fp_out3 = fopen("shift_data.dat", "wb") ) == NULL) {
    fprintf(stderr, "Error opening output file!\n");
    exit(0);
  }
  for(int i = 0; i < nchans; i++)
	fprintf(fp_out3,"%d\n", shifts_index[i]);
  fclose(fp_out3);


//  std::cout << "dmshift: " << shifts_index[2047] << std::endl;
//  std::cout << "max_shift: " << strategy.maxshift() << std::endl;

  unsigned short* input_data;
  input_data = (unsigned short *)malloc(sizeof(unsigned short)*nsamples*nchans);
  fake_generate_data(input_data, scale_factor, shifts_index, func_width, signal_start, signal_max_pos, nchans);

  FILE *fp_out2;
  if (( fp_out2 = fopen("fake_data.dat", "wb") ) == NULL) {
    fprintf(stderr, "Error opening output file!\n");
    exit(0);
  }

  for (int i = 0; i < nchans; i++)
	for (int j = 0; j < nsamples; j++)
		fprintf(fp_out2,"%d %d %i\n", i, j, input_data[j*nchans + i]);
  fclose(fp_out2);

  FILE *fp_out;
  if (( fp_out = fopen("DD_data.dat", "wb") ) == NULL) {
    fprintf(stderr, "Error opening output file!\n");
    exit(0);
  }
  	
  aa_permitted_pipelines_1<aa_compute::module_option::empty, false> runner(strategy, input_data);
  if(runner.setup()) {
    int chunk_idx = 0;
    std::vector<std::vector<std::vector<float>>> out;
// The user should consume the output vector data
    // upon each iteration of .next(out), since
    // the vector memory is re-allocated for the next chunk.
    while(runner.next(out, chunk_idx)) {
      std::cout << "NOTICE: Pipeline running over next chunk." << std::endl;
    }

    for(size_t i = 0; i < out.at(0).size(); i++)
	for (size_t j = 0; j < out.at(0).at(i).size(); j++)
		fprintf(fp_out, "%zu %zu %.8lf \n", i, j, out.at(0).at(i).at(j));
	
  std::cout << "Value is: " << out.at(0).at(50).at(23437) << std::endl;

  }


  free(scale_factor);
  free(shifts_index);
  free(input_data);
  fclose(fp_out);
  std::cout << "NOTICE: Finished." << std::endl;
  return 0;
}
