#include <iostream>
#include <time.h>

#include "aa_ddtr_plan.hpp"
#include "aa_ddtr_strategy.hpp"
#include "aa_filterbank_metadata.hpp"
#include "aa_permitted_pipelines_generic.hpp"
#include "aa_pipeline_api.hpp"
#include "aa_device_info.hpp"
#include "aa_params.hpp"
#include "aa_log.hpp"
#include "aa_host_fake_signal_generator.hpp"

#include "aa_sigproc_input.hpp"

using namespace astroaccelerate;

int main() {
	//-----------------------  Init the GPU card
	int device = 0;
	aa_device_info selected_device(device);
	
	selected_device.print_card_info();
	//-------------------------------------------

	//-------- Define user DM plan 
	aa_ddtr_plan ddtr_plan;
	float dm_search_max = 250.0;
	int number_of_tests = 10;
	int count_success_test = 0;
	ddtr_plan.add_dm(0, dm_search_max, 1, 1, 1); // Add dm_ranges: dm_low, dm_high, dm_step, inBin, outBin (unused).
	//----------------------------
	
	//initialize the random seed
	srand(time(NULL));

	for(int t=0; t < number_of_tests; t++){	
		// Filterbank metadata
		// (Data description from "SIGPROC-v3.7 (Pulsar) Signal Processing Programs")
		const double tstart = 50000;
		const double total_bandwidth = 300.0f;
		const double tsamp = 6.4E-5; // ----> 15625 samples per second
		const double nbits = 8;
		const unsigned int nsamples = 2/tsamp;  // set to 2s of observation 
		const double fch1 = 1550; // central frequency --> 1400 MHz
		const int nchans = 128;
		const double foff = -total_bandwidth/nchans;
		// params needed by the fake signal function
		double dm_position = rand()%(int)(0.8*dm_search_max); // at what dm put the signal
		const int func_width = 1/(tsamp*25); // width of the signal in terms of # of samples;
		const int signal_start = (rand()%10)*1000 + (rand()%10)*100 + (rand()%10)*10;// 0.2/tsamp; // position of the signal in samples; mean the position of the peak;
		bool dump_signal_to_disk = false; // this option writes the generated signal to a file 'fake_signal.dat'
		const float sigma = 0.5;
		//---------------------------------------------------------------------------
		
		// setting the signal metadata
		aa_filterbank_metadata metadata(tstart, tsamp, nbits, nsamples, fch1, foff, nchans);
		// setting the metadata for running fake generator
		aa_fake_signal_metadata f_meta(dm_position, signal_start, func_width, sigma);	
	
		const size_t free_memory = selected_device.free_memory(); // Free memory on the GPU in bytes
		bool enable_analysis = false; 
		
		aa_ddtr_strategy strategy(ddtr_plan, metadata, free_memory, enable_analysis, &selected_device);
		if(!(strategy.ready())) {
			std::cout << "There was an error" << std::endl;
			return 0;
		}
	
		// creating the signal -------------------------
		aa_fake_signal_generator signal;
		signal.create_signal(strategy, f_meta, metadata, dump_signal_to_disk);
		if(!(signal.ready())) {
		      std::cout << "Error in creating fake signal" << std::endl;
		      return 0;
		}
	
		std::vector<unsigned short> input_data;
		input_data = signal.signal_data();
		//-----------------------------------------------
	
		aa_pipeline::pipeline pipeline_components;
	        pipeline_components.insert(aa_pipeline::component::dedispersion);
	
		aa_pipeline::pipeline_option pipeline_options;
	        //insert option to copy the DDTR output data from GPU memory to the host memory
	        //do not insert this option if the output is not needed
	        pipeline_options.insert(aa_pipeline::component_option::copy_ddtr_data_to_host);
	
		aa_pipeline_api<unsigned short> runner(pipeline_components, pipeline_options, metadata, input_data.data(), selected_device);
		runner.bind(ddtr_plan);
	 
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
	
		if ( ptr[0][dm_index_position][signal_start] == 255 ){
			LOG(log_level::notice, "Peak found at position. [" + std::to_string(dm_index_position) + ", " + std::to_string(signal_start) + ", " + std::to_string(ptr[0][dm_index_position][signal_start]) + "].");
			count_success_test++;
		}
		else {
			LOG(log_level::notice, "Peak not found. [" + std::to_string(dm_index_position) + ", " + std::to_string(signal_start) + ", " + std::to_string(ptr[0][dm_index_position][signal_start]) + "].");
		}
	} // for loop for certain number of tests.

	if (count_success_test == number_of_tests){
		LOG(log_level::notice, "Test OK. Passed: " + std::to_string(count_success_test) + "/" + std::to_string(number_of_tests) + ".");
	}
	else {
		LOG(log_level::notice, "Test not passed: " + std::to_string(count_success_test) + "/" + std::to_string(number_of_tests) + ".");
	}

//	signal.print_info(f_meta);	// print info about signal
//	strategy.print_info(strategy);  // print info about strategy
	return 0;
}
