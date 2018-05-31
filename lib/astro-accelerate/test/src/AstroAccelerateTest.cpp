#include "AstroAccelerateTest.h"
#include "../../AstroAccelerateParameters.h"
#include "../../DedispersionStrategy.h"
#include "../../DmTime.h"
#include "../../AstroAccelerate.h"

#include <vector>

namespace astroaccelerate {
namespace test {

AstroAccelerateTest::AstroAccelerateTest()
    : ::testing::Test()
{
}

AstroAccelerateTest::~AstroAccelerateTest()
{
}

void AstroAccelerateTest::SetUp()
{
}

void AstroAccelerateTest::TearDown()
{
}

class TestParams : public AstroAccelerateParameters<TestParams> {};

// test AstroAccelerate class
TEST_F(AstroAccelerateTest, AstroAccelerate_call)
{
/*
	// Internal code variables
	// File pointers
	FILE *fp = NULL;
	// Counters and flags
	int i, t, dm_range;
	int range = 0;
	int nb_selected_dm = 0;
	int enable_debug = 0;
	int enable_analysis = 0;
	int enable_acceleration = 0;
	int enable_periodicity = 0;
	int output_dmt = 0;
	int enable_zero_dm = 0;
	int enable_zero_dm_with_outliers = 0;
	int enable_rfi = 0;
	int enable_sps_baselinenoise=0;
	int enable_fdas_custom_fft = 0;
	int enable_fdas_inbin = 0;
	int enable_fdas_norm = 0;
	int enable_output_ffdot_plan = 0;
	int enable_output_fdas_list = 0;
	int analysis_debug = 0;
    int *inBin = NULL;
	int *outBin = NULL;
	int *ndms = NULL;
	int maxshift = 0;
	int max_ndms = 0;
	int max_samps = 0;
	int num_tchunks = 0;
	int total_ndms = 0;
	int multi_file = 1;
	float max_dm = 0.0f;
	int candidate_algorithm=0;
	int failsafe = 0;
	// Memory sizes and pointers
	size_t outputsize = 0;
	size_t gpu_inputsize = 0;
	size_t gpu_outputsize = 0;
	unsigned short *d_input = NULL;
	float *d_output = NULL;
	float *dmshifts = NULL;
	float *user_dm_low = NULL;
	float *user_dm_high = NULL;
	float *user_dm_step = NULL;
	float *dm_low = NULL;
	float *dm_high = NULL;
	float *dm_step = NULL;
	float *selected_dm_low = NULL;
	float *selected_dm_high = NULL;
	// Telescope parameters
	int nchans = 0;
	int nsamp = 0;
	int nbits = 0;
	int nsamples = 0;
	int nifs = 0;
	int nboots = -1;
	int ntrial_bins;
	int navdms = 1;
	int nsearch = 3;
	float aggression = 2.5;
	float narrow = 0.001f;
	float wide = 0.1f;
	float tstart = 0.0f;
	float tsamp = 0.0f;
	float fch1 = 0.0f;
	float foff = 0.0f;
	// Analysis variables
	float power = 2.0f;
	float sigma_cutoff = 6.0f;
	float sigma_constant = 4.0f;
	float max_boxcar_width_in_sec = 0.5f;

		// Users desired de-dispersion strategy. Pick up user defined values from the CLI.
		get_user_input(&fp, my_argc, my_argv, &multi_file, &enable_debug, &enable_analysis,
		    &enable_periodicity, &enable_acceleration, &enable_output_ffdot_plan,
		    &enable_output_fdas_list, &output_dmt, &enable_zero_dm,
		    &enable_zero_dm_with_outliers, &enable_rfi, &enable_fdas_custom_fft,
		    &enable_fdas_inbin, &enable_fdas_norm, &nboots, &ntrial_bins, &navdms,
		    &narrow, &wide, &aggression, &nsearch, &inBin, &outBin, &power, &sigma_cutoff,
		    &sigma_constant, &max_boxcar_width_in_sec, &range, &user_dm_low, &user_dm_high,
		    &user_dm_step, &candidate_algorithm, &enable_sps_baselinenoise, &selected_dm_low, &selected_dm_high, &nb_selected_dm, &analysis_debug, &failsafe);
		// Reads telescope parameters from the header of the input file and then counts the number of samples in the input data file.
		get_file_data(&fp, &nchans, &nsamples, &nsamp, &nifs, &nbits, &tsamp, &tstart,
		&fch1, &foff);

		std::vector<float> vec_user_dm_low(range);
		std::vector<float> vec_user_dm_high(range);
		std::vector<float> vec_user_dm_step(range);
		std::vector<int> vec_inBin(range);
		std::vector<int> vec_outBin(range);

		for(int i = 0; i < range; ++i)
		{
			vec_user_dm_low[i] = user_dm_low[i];
			vec_user_dm_high[i] = user_dm_high[i];
			vec_user_dm_step[i] = user_dm_step[i];
			vec_inBin[i] = inBin[i];
			vec_outBin[i] = outBin[i];
		}

		// Initialise the GPU.
		int device_id = 0;
		size_t gpu_memory = 0;
		cudaSetDevice(device_id);
		size_t mem_free, total;
		cudaMemGetInfo(&mem_free, &total);
		gpu_memory = ( mem_free );

		DedispersionStrategy dedispersion_strategy
											 (vec_user_dm_low
											 ,vec_user_dm_high
											 ,vec_user_dm_step
											 ,vec_inBin
											 ,vec_outBin
											 ,gpu_memory
											 ,power
											 ,range
											 ,nchans
											 ,nsamp
											 ,nifs
											 ,nbits
											 ,tsamp
											 ,fch1
											 ,foff
											 ,sigma_cutoff
											 ,sigma_constant
											 ,max_boxcar_width_in_sec
											 ,narrow
											 ,wide
											 ,nboots
											 ,navdms
											 ,ntrial_bins
											 ,nsearch
											 ,aggression);

		// input buffer
		unsigned short *input_buffer = nullptr;
		size_t inputsize = 0;
		allocate_memory_cpu_input(dedispersion_strategy.get_nsamp(), dedispersion_strategy.get_nchans(), &input_buffer,&inputsize);
		get_recorded_data(&fp, dedispersion_strategy.get_nsamp(), dedispersion_strategy.get_nchans(), dedispersion_strategy.get_nbits(),
						  &input_buffer, &inputsize);

		printf("\nDD + SPS starting\n");
		// dedispersed data
		DmTime<float> output_buffer(dedispersion_strategy);
		// output of sps - assume it's a quarter of the output size
		std::vector<float> output_sps;

		astroaccelerate::AstroAccelerate<TestParams> astroaccelerate(dedispersion_strategy);
		astroaccelerate.run_dedispersion_sps(device_id
						,input_buffer
						,output_buffer
						,output_sps
						);


		/*
		 * Note: this function was written with a buggy fdas in the master, so it's not running.
		 * A fix has been done on the fdas code in the master, merging it to the interface should allow this function to run.
		 * It's basically filling a vector of candidates in the same way as in run_dedispersion_sps
		 *
		 * Todo: find a way to unit test it. Currently just write to file the desired dm and check if the results are
		 * the same as in the master branch
		 *
		 * 1 candidate is 4 elements of the vector: acceleration, frequency derivative, power, DM
		 */
	/*
		printf("\nFDAS is starting\n");
		std::vector<float> output_fdas;
		astroaccelerate.run_fdas(device_id, output_buffer, output_fdas);
  

		fclose(fp);
		free(inBin);
		free(outBin);
		free(user_dm_low);
		free(user_dm_high);
		free(user_dm_step);
		free(input_buffer);
*/
}
//}

} // namespace test
} // namespace astroaccelerate

