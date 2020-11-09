#include <iostream>

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
	//-------------------------------------------

	//-------- Define user DM plan 
	aa_ddtr_plan ddtr_plan;
	ddtr_plan.add_dm(0, 250, 1, 1, 1); // Add dm_ranges: dm_low, dm_high, dm_step, inBin, outBin (unused).
	//----------------------------
	
	// Filterbank metadata
	// (Data description from "SIGPROC-v3.7 (Pulsar) Signal Processing Programs")
	const double tstart = 50000;
	const double total_bandwidth = 300.0f;
	const double tsamp = 6.4E-5;
	const double nbits = 8;
	const unsigned int nsamples = 1.0/tsamp;
	const double fch1 = 1550;
	const int nchans = 128;
	const double foff = -total_bandwidth/nchans;
	// params needed by the fake signal function
	double dm_position = 250.0; // at what dm put the signal
	const int func_width = 1/(tsamp*25); // width of the signal in terms of # of samples;
	const int signal_start = 0.2/tsamp; // position of the signal in samples; mean the position of the peak;
	bool dump_signal_to_disk = true; // this option writes the generated signal to a file 'fake_signal.dat'
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
	

	//------------  Get data and write to a file 'ddtr_data.dat' 
	// Please note that the ddtr output is not scaled to true units
	float ***ptr = runner.output_buffer();
	
	FILE *fp;
	char filename[200];
	sprintf(filename, "ddtr_data.dat");
	if ((fp=fopen(filename, "wb")) == NULL) {
		fprintf(stderr, "Error opening output file for fake signal!\n");
		exit(0);
	}
	
	for(size_t i = 0; i < strategy.get_nRanges(); i++ ){
	      for (int j = 0; j < strategy.ndms(i); j++ ) {
	      	for (int k = 0; k < strategy.t_processed()[i][0]; k++ ){
	      		fprintf(fp, "%hu %d %lf\n", j, k, ptr[i][j][k]);
	      	}
	      }
	}
	
	fclose(fp);
	//----------------------------------------------------------
	
	signal.print_info(f_meta);
	strategy.print_info(strategy);
	
	LOG(log_level::notice, "Finished.");
	
	return 0;
}
