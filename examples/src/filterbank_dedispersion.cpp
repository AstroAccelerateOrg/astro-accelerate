/**
 * Example code for linking against astro-accelerate library.
 * 
 * Compile with: g++ -std=c++11 -I/path/to/astro-accelerate/include/ -L/path/to/astro-accelerate/build -Wl,-rpath,/path/to/astro-accelerate/build -lastroaccelerate -I/usr/local/cuda/include/ -L/usr/local/cuda/lib64 -Wl,-rpath,/usr/local/cuda/lib64 -I/usr/local/cuda/samples/common/inc/ -lcudart dedispersion.cpp -o test
*/

#include "aa_sigproc_input.hpp"
#include "aa_pipeline_wrapper_functions.hpp"

using namespace astroaccelerate;

int main() {
	//------------ Read the filterbank data (please enter the location of existing filterbank)
	aa_sigproc_input       filterbank_datafile("...<enter_desired_filterbank_file>...");
	aa_filterbank_metadata filterbank_metadata = filterbank_datafile.read_metadata();

	
	if(!filterbank_datafile.read_signal()) {
		std::cout << "ERROR: Could not read telescope data." << std::endl;
		return 0;
	}
	//------------<
	
	//------------ Create desired ddtr plan
	std::vector<aa_ddtr_plan::dm> dm_ranges;
	aa_ddtr_plan::dm range1 = {0, 150, 0.1, 1, 1};
	aa_ddtr_plan::dm range2 = {150, 300, 0.2, 1, 1};
	aa_ddtr_plan::dm range3 = {300, 500, 0.25, 1, 1};
	aa_ddtr_plan::dm range4 = {500, 900, 0.4, 2, 2};
	aa_ddtr_plan::dm range5 = {900, 1200, 0.6, 4, 4};
	aa_ddtr_plan::dm range6 = {1200, 1500, 0.8, 4, 4};
	aa_ddtr_plan::dm range7 = {1500, 2000, 1.0, 4, 4};
	aa_ddtr_plan::dm range8 = {2000, 3000, 2.0, 8, 8};
	
	dm_ranges.push_back(range1);
	dm_ranges.push_back(range2);
	dm_ranges.push_back(range3);
	dm_ranges.push_back(range4);
	dm_ranges.push_back(range5);
	dm_ranges.push_back(range6);
	dm_ranges.push_back(range7);
	dm_ranges.push_back(range8);
	//------------<
	
	const aa_pipeline::pipeline_option pipeline_options = {aa_pipeline::component_option::zero_dm};

	//------------ Run the pipeline with computing only ddtr on the provided filterbank file and plan
	dedisperse_telescope_data(filterbank_metadata, pipeline_options, dm_ranges, filterbank_datafile.input_buffer());
	//------------<

	std::cout << "NOTICE: Finished." << std::endl;
	return 0;
}
