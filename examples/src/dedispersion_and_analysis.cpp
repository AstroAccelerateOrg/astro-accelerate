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

int main() {
	//-------------- Select de-dispersion plan
	aa_ddtr_plan ddtr_plan;
	ddtr_plan.add_dm(0, 370, 0.307, 1, 1); // Add dm_ranges: dm_low, dm_high, dm_step, inBin, outBin (unused).
	ddtr_plan.add_dm(370, 740, 0.652, 1, 1);
	ddtr_plan.add_dm(740, 1480, 1.266, 4, 4);
	ddtr_plan.add_dm(1480, 2950, 25.12, 8, 8);
	ddtr_plan.add_dm(2950, 5000, 4.000, 16, 16);
	//--------------<
	
	// Filterbank metadata
	aa_sigproc_input filterbank_datafile("/home/novotny/filterbank/aa_test_file_dm90_snr10_w064_tobs30.fil");
	aa_filterbank_metadata metadata = filterbank_datafile.read_metadata();
	filterbank_datafile.read_signal();

	aa_device_info& device_info = aa_device_info::instance();
	aa_device_info::CARD_ID selected_card_number = 0;
	aa_device_info::aa_card_info selected_card_info; 
        device_info.init_card(selected_card_number, selected_card_info);

	//-------------- Configure pipeline. Select components and their options
	aa_pipeline::pipeline pipeline_components;
	pipeline_components.insert(aa_pipeline::component::dedispersion); // pipeline must always contain dedispersion step
        pipeline_components.insert(aa_pipeline::component::analysis); //optional
        //pipeline_components.insert(aa_pipeline::component::periodicity); // optional
        //pipeline_components.insert(aa_pipeline::component::fdas); // optional
	
        aa_pipeline::pipeline_option pipeline_options;
	pipeline_options.insert(aa_pipeline::component_option::msd_baseline_noise);
	//--------------<
	
        //-------------- Configure single pulse detection plan and calculate strategy
        const float sigma_cutoff = 6.0;
        const float sigma_constant = 4.0;
        const float max_boxcar_width_in_sec = 0.5;
        const bool  enable_MSD_outlier_rejection = true;
	aa_analysis_plan::selectable_candidate_algorithm candidate_algorithm = aa_analysis_plan::selectable_candidate_algorithm::off;

	aa_pipeline_api<unsigned short> pipeline_runner(pipeline_components, pipeline_options, metadata, filterbank_datafile.input_buffer().data(), selected_card_info);

	pipeline_runner.bind(ddtr_plan);
       
	aa_analysis_plan analysis_plan(pipeline_runner.ddtr_strategy(), sigma_cutoff, sigma_constant, max_boxcar_width_in_sec, candidate_algorithm, enable_MSD_outlier_rejection);
	pipeline_runner.bind(analysis_plan);

        if (pipeline_runner.ready()) {
                LOG(log_level::notice, "Pipeline is ready.");
        }
        else {
                LOG(log_level::notice, "Pipeline is not ready.");
        }

	//------------- Run the pipeline
	size_t nCandidates;
	unsigned int* dm;
	unsigned int* ts;
	unsigned int* width;
	float* snr;
	int c_range, c_tchunk;
	long ts_inc;
	aa_pipeline_runner::status status_code;
	while(pipeline_runner.run(status_code)){
		if ((int)status_code == 1){
			nCandidates = pipeline_runner.SPD_nCandidates();
			dm = pipeline_runner.h_SPD_dm();
			ts = pipeline_runner.h_SPD_width();
			snr = pipeline_runner.h_SPD_snr();
			width = pipeline_runner. h_SPD_width();
			c_range = pipeline_runner.get_current_range();
			c_tchunk = pipeline_runner.get_current_tchunk();
			ts_inc = pipeline_runner.get_current_inc();
		}
	}
	//-------------<
	
	std::cout << "NOTICE: Finished." << std::endl;

	return 0;
}
