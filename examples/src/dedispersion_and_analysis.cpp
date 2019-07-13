/**
 * Example code for linking against astro-accelerate library.
 *
 * Compile with: g++ -std=c++11 -I/path/to/astro-accelerate/include/ -L/path/to/astro-accelerate/build -Wl,-rpath,/path/to/astro-accelerate/build -lastroaccelerate -I/usr/local/cuda/include/ -L/usr/local/cuda/lib64 -Wl,-rpath,/usr/local/cuda/lib64 -I/usr/local/cuda-8.0/samples/common/inc/ -lcudart dedispersion.cpp -o test
 */

#include <iostream>

#include "aa_ddtr_plan.hpp"
#include "aa_ddtr_strategy.hpp"
#include "aa_filterbank_metadata.hpp"
#include "aa_permitted_pipelines_generic.hpp"

#include "aa_analysis_plan.hpp"
#include "aa_analysis_strategy.hpp"

#include "aa_periodicity_plan.hpp"
#include "aa_periodicity_strategy.hpp"

#include "aa_fdas_plan.hpp"
#include "aa_fdas_strategy.hpp"

#include "aa_log.hpp"
#include "aa_sigproc_input.hpp"

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>

using namespace astroaccelerate;

int initialise_device(int device_id) {
	int deviceCount;
	cudaError_t error_id;
	error_id = cudaGetDeviceCount(&deviceCount);
	if (error_id != cudaSuccess) {
		printf("CUDA ERROR: %s\n", cudaGetErrorString(error_id));
		return(1);
	}
	if (device_id>=deviceCount) {
		printf("Selected device is not available! Device id is %d;\n", device_id);
		return(1);
	}
	if (cudaSetDevice(device_id) != cudaSuccess) {
		printf("ERROR! unable to set the device with id %d.\n", device_id);
		return(1);
	}
	return(0);
}

int main() {
	//------------- Initialise device and get available memory
	if(initialise_device(0)!=0) {
		return(1);
	}
	size_t free_memory,total_memory;
	cudaMemGetInfo(&free_memory,&total_memory);
	//--------------<
	
	//-------------- Select de-dispersion plan
	aa_ddtr_plan ddtr_plan;
	ddtr_plan.add_dm(0, 370, 0.307, 1, 1); // Add dm_ranges: dm_low, dm_high, dm_step, inBin, outBin (unused).
	ddtr_plan.add_dm(370, 740, 0.652, 2, 2);
	ddtr_plan.add_dm(740, 1480, 1.266, 4, 4);
	//--------------<

	//-------------- Read filterbank metadata and data
	aa_sigproc_input filterbank_datafile("/mnt/data/AstroAccelerate/filterbank/BenMeerKAT.fil");
	aa_filterbank_metadata filterbank_metadata = filterbank_datafile.read_metadata();

	if (!filterbank_datafile.read_signal()) {
		std::cout << "ERROR: Could not read telescope data." << std::endl;
		return 0;
	}
	//--------------<
	
	//-------------- Configure pipeline. Select components and their options.
	aa_pipeline::pipeline pipeline_components;
	pipeline_components.insert(aa_pipeline::component::dedispersion); // pipeline must always contain dedispersion step
	pipeline_components.insert(aa_pipeline::component::analysis); //optional
	//pipeline.insert(aa_pipeline::component::periodicity); // optional
	//pipeline.insert(aa_pipeline::component::fdas); // optional
	
	aa_pipeline::pipeline_option pipeline_options;
	pipeline_options.insert(aa_pipeline::component_option::zero_dm);
	//--------------<


	//-------------- Calculate dedispersion strategy
	bool enable_analysis = true;       // The strategy will be optimised to run just dedispersion
	aa_ddtr_strategy ddtr_strategy(ddtr_plan, filterbank_metadata, free_memory, enable_analysis);	
	if (!(ddtr_strategy.ready())) {
		LOG(log_level::error, "ddtr_strategy not ready.");
		return 0;
	}
	//--------------<
	
	//-------------- Configure single pulse detection plan and calculate strategy
	const float sigma_cutoff = 6.0;
	const float sigma_constant = 4.0;
	const float max_boxcar_width_in_sec = 0.5;
	const bool  enable_MSD_outlier_rejection = true;
	const aa_analysis_plan::selectable_candidate_algorithm algo = aa_analysis_plan::selectable_candidate_algorithm::off;

	aa_analysis_plan analysis_plan(ddtr_strategy, sigma_cutoff, sigma_constant, max_boxcar_width_in_sec, algo, enable_MSD_outlier_rejection);
	aa_analysis_strategy analysis_strategy(analysis_plan);

	if (!(analysis_strategy.ready())) {
		LOG(log_level::error, "analysis_strategy not ready.");
		return 0;
	}
	//--------------<
	
	//-------------- Create empty strategy object for unused components
	aa_fdas_strategy empty_fdas_strategy;
	aa_periodicity_strategy empty_periodicity_strategy;
	//--------------<
	
	
	aa_permitted_pipelines_generic pipeline_runner(pipeline_components, pipeline_options, ddtr_strategy, analysis_strategy, empty_periodicity_strategy, empty_fdas_strategy, false, false, false, false, false, filterbank_datafile.input_buffer().data());
	
	if (pipeline_runner.setup()) {
		while (pipeline_runner.next()) {
			LOG(log_level::notice, "Pipeline running over next chunk.");
		}
	}

	
	//aa_pipeline_runner::status &status_code
	
	
	
	
	
	LOG(log_level::notice, "Finished.");
	return 0;
}
