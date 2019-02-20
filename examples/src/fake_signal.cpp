#include <iostream>

#include "aa_ddtr_plan.hpp"
#include "aa_ddtr_strategy.hpp"
#include "aa_filterbank_metadata.hpp"
#include "aa_permitted_pipelines_1.hpp"
#include "aa_device_info.hpp"
#include "aa_params.hpp"
#include "aa_log.hpp"
#include "aa_host_fake_signal_generator.hpp"

#define MAX_VALUE 255

using namespace astroaccelerate;

int main() {
  aa_ddtr_plan ddtr_plan;
  ddtr_plan.add_dm(0, 150, 1, 1, 1); // Add dm_ranges: dm_low, dm_high, dm_step, inBin, outBin (unused).


  // Filterbank metadata
  // (Data description from "SIGPROC-v3.7 (Pulsar) Signal Processing Programs")
  const double tstart = 50000;
  const double total_bandwidth = 300.0f;
  const double tsamp = 6.4E-5;
  const double nbits = 8;
  unsigned int nsamples = 2.0/tsamp; // 3s of data in current tsamp at minimum; must be more than signal_start + maxshift
  const double fch1 = 1550;
  int nchans = 128;
  const double foff = -total_bandwidth/nchans;
  std::cout << "Nsamples " << nsamples << std::endl;

  // params needed by the fake signal function
  double dm_position = 90; // at what dm put the signal
  const int func_width = 1/(tsamp*100); // width of the signal in terms of # of samples; now at 1% of samling rate
  int signal_start = 1.0/tsamp; // position of the signal in samples

  // setting the metadata for running fake generator
  aa_fake_signal_metadata f_meta(dm_position, signal_start, func_width);
  // setting the signal metadata
  aa_filterbank_metadata metadata(tstart, tsamp, nbits, nsamples, fch1, foff, nchans);

  // Init the GPU card
  aa_device_info& device_info = aa_device_info::instance();
  if(device_info.check_for_devices()) {
    LOG(log_level::notice, "Checked for devices.");
  }
  else {
    LOG(log_level::error, "Could not find any devices.");
  }

  aa_device_info::CARD_ID selected_card = 0;
  aa_device_info::aa_card_info selected_card_info;
  if(device_info.init_card(selected_card, selected_card_info)) {
    LOG(log_level::notice, "init_card complete. Selected card " + std::to_string(selected_card) + ".");
  }
  else {
    LOG(log_level::error, "init_card incomplete.")
  }

  aa_device_info::print_card_info(selected_card_info);
  
  const size_t free_memory = selected_card_info.free_memory; // Free memory on the GPU in bytes
  bool enable_analysis = false; 

  aa_ddtr_strategy strategy(ddtr_plan, metadata, free_memory, enable_analysis);

  if(!(strategy.ready())) {
    std::cout << "There was an error" << std::endl;
    return 0;
  }

  aa_fake_signal_generator signal;
  signal.produce_mask(0.5, func_width);
  signal.get_fake_signal(strategy, f_meta, metadata);
  std::vector<unsigned short> input_data;
  input_data = signal.signal_data();
//  fake_generate_data(input_data,  scale_factor, shifts_index, func_width, signal_start, signal_max_pos, nchans);

  FILE *fp_out2;
  if (( fp_out2 = fopen("fake_data.dat", "wb") ) == NULL) {
    fprintf(stderr, "Error opening output file!\n");
    exit(0);
  }

  for (int i = 0; i < nchans; i++)
	for (int j = 0; j < nsamples; j++)
		fprintf(fp_out2,"%d %d %i\n", i, j, input_data[j*nchans + i]);
  fclose(fp_out2);

  aa_permitted_pipelines_1<aa_pipeline::component_option::empty, false> runner(strategy, input_data.data());
  
  if(runner.setup()) {
    while(runner.next(true)) {
      LOG(log_level::notice, "Pipeline running over next chunk.");
    }
  }


  float *** ptr = runner.output_buffer();
  const int *ndms = strategy.ndms_data();
  LOG(log_level::dev_debug, "#ndms: " + std::to_string(ndms[0]));
  LOG(log_level::dev_debug, "#t_processed: " + std::to_string( strategy.t_processed()[0][0] ) );

  FILE *fp;
  char filename[200];
  sprintf(filename, "ddtr_data.dat");
    if ((fp=fopen(filename, "wb")) == NULL) {
      fprintf(stderr, "Error opening output file for fake signal!\n");
      exit(0);
    }

  for(int i = 0; i < strategy.range(); i++ ){
	for (int j = 0; j < strategy.ndms(i); j++ ) {
		for (int k = 0; k < strategy.t_processed()[i][0]; k++ ){
			fprintf(fp, "%d %d %lf\n", j, k, ptr[i][j][k]);
		}
	}
  }

  fclose(fp);
  strategy.print_info(strategy);
  
  LOG(log_level::notice, "Finished.");

  return 0;
}
