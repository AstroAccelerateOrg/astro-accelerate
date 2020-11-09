#include <iostream>

#include "aa_ddtr_plan.hpp"
#include "aa_ddtr_strategy.hpp"
#include "aa_filterbank_metadata.hpp"
#include "aa_device_info.hpp"
#include "aa_params.hpp"
#include "aa_log.hpp"
#include "aa_host_fake_signal_generator.hpp"

#include "aa_permitted_pipelines_generic.hpp"
#include "aa_pipeline_api.hpp"

#include "aa_analysis_plan.hpp"
#include "aa_analysis_strategy.hpp"

#include "aa_periodicity_plan.hpp"
#include "aa_periodicity_strategy.hpp"

#include <time.h>

using namespace astroaccelerate;

int main() {
  aa_ddtr_plan ddtr_plan;
  float dm_search_max = 120.0;
  ddtr_plan.add_dm(0, dm_search_max, 1, 1, 1); // Add dm_ranges: dm_low, dm_high, dm_step, inBin, outBin (unused).

  srand(time(NULL));

  // Filterbank metadata
  // (Data description from "SIGPROC-v3.7 (Pulsar) Signal Processing Programs")
  const double tstart = 50000;
  const double total_bandwidth = 300.0f;
  const double tsamp = 6.4E-5;
  const double nbits = 8;
  const unsigned int nsamples = 1.5/tsamp; // 2s of data in current tsamp at minimum; must be more than signal_start + maxshift
  const double fch1 = 1550;
  const int nchans = 128;
  const double foff = -total_bandwidth/nchans;

  // setting the signal metadata
  aa_filterbank_metadata metadata(tstart, tsamp, nbits, nsamples, fch1, foff, nchans);

  // Init the GPU card
  int device = 0;
  aa_device_info selected_device(device);

  selected_device.print_card_info();
  
  const size_t free_memory = selected_device.free_memory(); // Free memory on the GPU in bytes
  bool enable_analysis = false; 

  aa_ddtr_strategy strategy(ddtr_plan, metadata, free_memory, enable_analysis, &selected_device);

  if(!(strategy.ready())) {
    std::cout << "There was an error" << std::endl;
    return 0;
  }

  // params needed by the fake signal function
  double dm_position = rand()%100; // at what dm put the signal
  const int func_width = 1/(tsamp*10); 
  const int period = rand()%100 + 100; // pulsar period with ms
  bool dump_to_disk = true;
  const float sigma = 0.5;

  // setting the metadata for running fake generator
  aa_fake_signal_metadata f_meta(dm_position, func_width, sigma, period);

  aa_fake_signal_generator signal;
  signal.create_signal(strategy, f_meta, metadata, dump_to_disk);

  if(!(signal.ready())) {
	std::cout << "Error in creating fake signal" << std::endl;
	return 0;
  }

  std::vector<unsigned short> input_data;
  input_data = signal.signal_data();

        aa_pipeline::pipeline pipeline_components;
        pipeline_components.insert(aa_pipeline::component::dedispersion);

  aa_pipeline::pipeline_option pipeline_options;
  pipeline_options.insert(aa_pipeline::component_option::copy_ddtr_data_to_host);

	aa_pipeline_api<unsigned short> runner(pipeline_components, pipeline_options, metadata, input_data.data(), selected_device);
        runner.bind(ddtr_plan);

	//test if the runner is ready.
        if (runner.ready()) {
                LOG(log_level::notice, "Pipeline is ready.");
        }
        else {
                LOG(log_level::notice, "Pipeline is not ready.");
        }

        //------------- Run the pipeline
                aa_pipeline_runner::status status_code;
                while(runner.run(status_code)){
                }
        //-------------<

                //------------  Get data 
                // Please note that the ddtr output is not scaled to true units
                float ***ptr = runner.output_buffer();
                int dm_index_position = (int)dm_position/(int)strategy.dm(0).step;
		int signal_start = -signal.max_pos() - 2*strategy.maxshift(); //f_meta.get_signal_start();
		int repeats = (nsamples + 2*strategy.maxshift())/(period) + 1;
		int sig_add = period/1000.0/tsamp;
		int r_pos = 0;
		while (signal_start < 0) {
			signal_start += sig_add;
			r_pos++;
		}

		int count = 0;
		int count_success = 0;
		for(int r = r_pos; r < repeats; r++){
			if ( (unsigned int)signal_start + (unsigned int)strategy.maxshift() < nsamples){
				count++;
		                if ( ptr[0][dm_index_position][signal_start] == 255 ){
		                        LOG(log_level::notice, "Peak #" + std::to_string(r-r_pos) + " found at position. [" + std::to_string(dm_index_position) + ", " + std::to_string(signal_start) + ", " + std::to_string(ptr[0][dm_index_position][signal_start]) + "], period: " + std::to_string(period) + " ms.");
					count_success++;
		                }
		                else {
		                        LOG(log_level::notice, "Peak #" + std::to_string(r_pos) + " not found. [" + std::to_string(dm_index_position) + ", " + std::to_string(signal_start) + ", " + std::to_string(ptr[0][dm_index_position][signal_start]) + "], period:" + std::to_string(period) + " ms.");
		                }
				signal_start+=sig_add;
			} 
			else {
				break;
			}
		}

		if ( (count == count_success) && (count_success > 0)) {
			LOG(log_level::notice, "Test passed. " + std::to_string(count_success) + "/" + std::to_string(count) + ".");
		}
		else {
			LOG(log_level::notice, "Test failed." + std::to_string(count_success) + "/" + std::to_string(count) + ".");
		}

//  signal.print_info(f_meta);
//  strategy.print_info(strategy);

  LOG(log_level::notice, "Finished.");

  return 0;
}
