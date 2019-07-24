/**
 * Example code for linking against astro-accelerate library.
 *
 * Compile with: g++ -std=c++11 -I/path/to/astro-accelerate/include/ -L/path/to/astro-accelerate/build -Wl,-rpath,/path/to/astro-accelerate/build -lastroaccelerate -I/usr/local/cuda/include/ -L/usr/local/cuda/lib64 -Wl,-rpath,/usr/local/cuda/lib64 -I/usr/local/cuda-8.0/samples/common/inc/ -lcudart dedispersion.cpp -o test
 */

#include <iostream>
#include <fstream>

#include "aa_ddtr_plan.hpp"
#include "aa_ddtr_strategy.hpp"
#include "aa_filterbank_metadata.hpp"
#include "aa_sigproc_input.hpp"
#include "aa_permitted_pipelines_generic.hpp"
#include "aa_pipeline_api.hpp"
#include "aa_device_info.hpp"

using namespace astroaccelerate;

void write_scale_candidates(aa_filterbank_metadata metadata, aa_pipeline_api<unsigned short> &pipeline, long int tprocessed, size_t nCandidates, unsigned int* dm, unsigned int* time_samples, float* snr, unsigned int* width, int current_range, int current_tchunk){
	float scaled_dm, scaled_time_sample, scaled_time;
	int scaled_width;
	char filename[100];
	std::ofstream output_file;

	const int *list_ndms = pipeline.get_ndms_array();
	aa_ddtr_strategy plan = pipeline.ddtr_strategy();
	float dm_low =  plan.dm(current_range).low;
	float dm_high = pipeline.dm_low(current_range) + list_ndms[current_range]*plan.dm(current_range).step;
	sprintf(filename, "results_t-%d_dm-%.3f-%.3f.txt", current_tchunk, dm_low, dm_high);
	output_file.open(filename);
//	printf("DM: %lf %lf %lf %d", dm_low, dm_high, plan.dm(current_range).step,list_ndms[current_range]);
	for (int i = 0; i < (int)nCandidates; i++){
		scaled_dm = dm[i]*plan.dm(current_range).step + dm_low;
                scaled_time_sample = time_samples[i]*plan.dm(current_range).inBin + tprocessed;
                scaled_time = time_samples[i]*metadata.tsamp()*plan.dm(current_range).inBin + tprocessed*metadata.tsamp();
                scaled_width = width[i];

		output_file << scaled_dm << "\t" << snr[i] << "\t" << scaled_time_sample << "\t" << scaled_time << "\t" << scaled_width << "\n";
	}

	output_file.close();
}

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
	aa_sigproc_input filterbank_datafile("/home/novotny/filterbank/ska-mid-b2-small.fil");
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
	size_t SPD_nCandidates; // number of candidates found in single pulse detection
	unsigned int* SPD_candidates_dm;
	unsigned int* SPD_candidates_timesample;
	unsigned int* SPD_candidates_width;
	float* SPD_candidates_snr;
	int c_range, c_tchunk;
	long int timesamples_processed_sofar;
	aa_pipeline_runner::status status_code;
	while(pipeline_runner.run(status_code)){
		if ((int)status_code == 1){
			SPD_nCandidates = pipeline_runner.SPD_nCandidates();
			SPD_candidates_dm = pipeline_runner.h_SPD_dm();
			SPD_candidates_timesample = pipeline_runner.h_SPD_ts();
			SPD_candidates_snr = pipeline_runner.h_SPD_snr();
			SPD_candidates_width = pipeline_runner. h_SPD_width();
			c_range = pipeline_runner.get_current_range();
			c_tchunk = pipeline_runner.get_current_tchunk();
			timesamples_processed_sofar = pipeline_runner.get_current_inc();
			write_scale_candidates(metadata, pipeline_runner, timesamples_processed_sofar, SPD_nCandidates, SPD_candidates_dm, SPD_candidates_timesample, SPD_candidates_snr, SPD_candidates_width, c_range, c_tchunk);		
			printf("Current range:%d; Current time chunk:%d; Time samples proceesed by pipeline so far:%zu;\n", c_range, c_tchunk, timesamples_processed_sofar);
		}
	}
	//-------------<
	
	std::cout << "NOTICE: Finished." << std::endl;

	return 0;
}
