/**
 * Example code for linking against astro-accelerate library.
 *
 * Compile with: g++ -std=c++11 -I/path/to/astro-accelerate/include/ -L/path/to/astro-accelerate/build -Wl,-rpath,/path/to/astro-accelerate/build -lastroaccelerate -I/usr/local/cuda/include/ -L/usr/local/cuda/lib64 -Wl,-rpath,/usr/local/cuda/lib64 -I/usr/local/cuda-8.0/samples/common/inc/ -lcudart dedispersion.cpp -o test
 */

#include <iostream>

#include "aa_ddtr_plan.hpp"
#include "aa_ddtr_strategy.hpp"
#include "aa_filterbank_metadata.hpp"
#include "aa_sigproc_input.hpp"
#include "aa_permitted_pipelines_generic.hpp"
#include "aa_pipeline_api.hpp"
#include "aa_device_info.hpp"

using namespace astroaccelerate;

int main(int argc, char *argv[]) {
	if (argc != 2) {
		LOG(log_level::notice, "Not enough arguments. To run type: ./examples_dedispersion <path to the fil file>");
		return 0;
	}
	//-------------- Select de-dispersion plan
	aa_ddtr_plan ddtr_plan;
	ddtr_plan.add_dm(0, 370, 0.307, 1, 1); // Add dm_ranges: dm_low, dm_high, dm_step, inBin, outBin (unused).
	ddtr_plan.add_dm(370, 740, 0.652, 1, 1);
	ddtr_plan.add_dm(740, 1480, 1.266, 4, 4);
	ddtr_plan.add_dm(1480, 2950, 25.12, 8, 8);
	ddtr_plan.add_dm(2950, 5000, 4.000, 16, 16);
	//--------------<
	
	// Filterbank metadata
	aa_sigproc_input filterbank_datafile(argv[1]);
	aa_filterbank_metadata metadata = filterbank_datafile.read_metadata();
	filterbank_datafile.read_signal();
	
	//Select desired device and initialize it by creating aa_device_info
	int device = 0;
	aa_device_info selected_device(device);

	//-------------- Configure pipeline. Select components and their options
	aa_pipeline::pipeline pipeline_components;
	pipeline_components.insert(aa_pipeline::component::dedispersion); // pipeline must always contain dedispersion step
	//pipeline_components.insert(aa_pipeline::component::analysis); //optional
	//pipeline_components.insert(aa_pipeline::component::periodicity); // optional
	//pipeline_components.insert(aa_pipeline::component::fdas); // optional
	
	aa_pipeline::pipeline_option pipeline_options;
	pipeline_options.insert(aa_pipeline::component_option::zero_dm);
	//insert option to copy the DDTR output data from GPU memory to the host memory
	//do not insert this option if the output is not needed
	pipeline_options.insert(aa_pipeline::component_option::copy_ddtr_data_to_host);
	//--------------<
	
	aa_pipeline_api<unsigned short> pipeline_runner(pipeline_components, pipeline_options, metadata, filterbank_datafile.input_buffer().data(), selected_device);
	pipeline_runner.bind(ddtr_plan);

        if (pipeline_runner.ready()) {
                LOG(log_level::notice, "Pipeline is ready.");
        }
        else {
                LOG(log_level::notice, "Pipeline is not ready.");
        }

	//------------- Run the pipeline
	aa_pipeline_runner::status status_code;
	while(pipeline_runner.run(status_code)){
	}
	//-------------<
	
	float ***ddtr_output;
	ddtr_output = pipeline_runner.output_buffer();
	aa_ddtr_strategy plan = pipeline_runner.ddtr_strategy();
	
	//ddtr_output[pos_range][n_dms][samples_pos]*plan.dm(pos_range).inBin
	LOG(log_level::notice, "DDDTR output[0][0][0]: " + std::to_string(ddtr_output[0][0][0]*plan.dm(0).inBin));

	std::cout << "NOTICE: Finished." << std::endl;

	return 0;
}
