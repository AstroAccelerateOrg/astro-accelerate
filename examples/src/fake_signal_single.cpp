#include <iostream>

#include "aa_ddtr_plan.hpp"
#include "aa_ddtr_strategy.hpp"
#include "aa_filterbank_metadata.hpp"
#include "aa_permitted_pipelines_1.hpp"
#include "aa_device_info.hpp"
#include "aa_params.hpp"
#include "aa_log.hpp"
#include "aa_host_fake_signal_generator.hpp"

#include "aa_sigproc_input.hpp"


using namespace astroaccelerate;

int main() {
  aa_ddtr_plan ddtr_plan;
  ddtr_plan.add_dm(0, 250, 1, 1, 1); // Add dm_ranges: dm_low, dm_high, dm_step, inBin, outBin (unused).

  // Filterbank metadata
  // (Data description from "SIGPROC-v3.7 (Pulsar) Signal Processing Programs")
  const double tstart = 50000;
  const double total_bandwidth = 300.0f;
  const double tsamp = 6.4E-5;
  const double nbits = 8;
  const unsigned int nsamples = 1.0/tsamp; // 2s of data in current tsamp at minimum; must be more than signal_start + maxshift
  const double fch1 = 1550;
  const int nchans = 128;
  const double foff = -total_bandwidth/nchans;

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

  // params needed by the fake signal function
  double dm_position = 250.0; // at what dm put the signal
  const int func_width = 1/(tsamp*25); // width of the signal in terms of # of samples; now at 1% of sampling rate
  const int signal_start = 0.2/tsamp; // position of the signal in samples
  bool dump_to_disk = true;
  const float sigma = 0.5;

  // setting the metadata for running fake generator
  aa_fake_signal_metadata f_meta(dm_position, signal_start, func_width, sigma);

  // creating the signal 
  aa_fake_signal_generator signal;
  signal.create_signal(strategy, f_meta, metadata, dump_to_disk);

  if(!(signal.ready())) {
	std::cout << "Error in creating fake signal" << std::endl;
	return 0;
  }

  std::vector<unsigned short> input_data;
  input_data = signal.signal_data();

  aa_permitted_pipelines_1<aa_pipeline::component_option::empty, false> runner(strategy, input_data.data());
 
  if(runner.setup()) {
    while(runner.next(true)) {
      LOG(log_level::notice, "Pipeline running over next chunk.");
    }
  }

  float *** ptr = runner.output_buffer();
  FILE *fp;
  char filename[200];
  sprintf(filename, "ddtr_data.dat");
    if ((fp=fopen(filename, "wb")) == NULL) {
      fprintf(stderr, "Error opening output file for fake signal!\n");
      exit(0);
    }

  for(size_t i = 0; i < strategy.range(); i++ ){
	for (int j = 0; j < strategy.ndms(i); j++ ) {
		for (int k = 0; k < strategy.t_processed()[i][0]; k++ ){
			fprintf(fp, "%hu %d %lf\n", j, k, ptr[i][j][k]);
		}
	}
  }

  fclose(fp);

  signal.print_info(f_meta);
  strategy.print_info(strategy);

  LOG(log_level::notice, "Finished.");

  return 0;
}
