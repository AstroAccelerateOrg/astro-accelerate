#include "AstroAccelerateTest.h"
#include "../../AstroAccelerateParameters.h"
#include "../../DedispersionStrategy.h"
#include "../../DedispersionStrategyFile.h"
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
/*
// test dedispersion strategy class
TEST_F(AstroAccelerateTest, test_dedispersion_strategy)
{

	// Following is ok for ska_karel.txt
	char* filename = my_argv[1] + strlen(my_argv[1]) - 13;
	if(strcmp(filename, "ska_karel.txt") == 0)
	{
		// Internal code variables
		// File pointers
		FILE *fp = nullptr;
		// Counters and flags
		int range = 0;
		int enable_debug = 0;
		int enable_analysis = 0;
		int enable_acceleration = 0;
		int enable_periodicity = 0;
		int output_dmt = 0;
		int enable_zero_dm = 0;
		int enable_zero_dm_with_outliers = 0;
		int enable_rfi = 0;
		int enable_fdas_custom_fft = 0;
		int enable_fdas_inbin = 0;
		int enable_fdas_norm = 0;
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
		// Memory sizes and pointers
		float *user_dm_low = nullptr;
		float *user_dm_high = nullptr;
		float *user_dm_step = nullptr;
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
		int maxshift_original;
		double tsamp_original;
		long int inc = 0;
		float tstart = 0.0f;
		float tstart_local = 0.0f;
		float tsamp = 0.0f;
		float fch1 = 0.0f;
		float foff = 0.0f;
		// Analysis variables
		float power = 2.0f;
		float sigma_cutoff = 6.0f;
		float sigma_constant = 4.0f;
		float max_boxcar_width_in_sec = 0.5f;
		// Initialise the GPU.
		int device_id = 0; // hard-coded, would be a parameter
		size_t gpu_memory = 0;
		cudaSetDevice(device_id);
		size_t mem_free, total;
		cudaMemGetInfo(&mem_free, &total);
		gpu_memory = ( mem_free/4 );

		// Users desired de-dispersion strategy. Pick up user defined values from the CLI.
		get_user_input(&fp, my_argc, my_argv, &multi_file, &enable_debug, &enable_analysis,
		&enable_periodicity, &enable_acceleration, &output_dmt, &enable_zero_dm,
		&enable_zero_dm_with_outliers, &enable_rfi, &enable_fdas_custom_fft,
		&enable_fdas_inbin, &enable_fdas_norm, &nboots, &ntrial_bins, &navdms,
		&narrow, &wide, &aggression, &nsearch, &inBin, &outBin, &power, &sigma_cutoff, &sigma_constant, &max_boxcar_width_in_sec,
		&range, &user_dm_low, &user_dm_high, &user_dm_step, &candidate_algorithm);
		// Reads telescope parameters from the header of the input file and then counts the number of samples in the input data file.
		get_file_data(&fp, &nchans, &nsamples, &nsamp, &nifs, &nbits, &tsamp, &tstart,
		&fch1, &foff);

		// dedispersion
		DedispersionStrategy dedispersion_strategy;
		DedispersionStrategyFile(&fp, my_argc, my_argv, dedispersion_strategy, gpu_memory);

		//
		EXPECT_EQ(-1, dedispersion_strategy.get_nboots());
		EXPECT_EQ(0, dedispersion_strategy.get_ntrial_bins());
		EXPECT_EQ(1, dedispersion_strategy.get_navdms());
		EXPECT_FLOAT_EQ(0.001f, dedispersion_strategy.get_narrow());
		EXPECT_FLOAT_EQ(2.5, dedispersion_strategy.get_aggression());
		EXPECT_EQ(3, dedispersion_strategy.get_nsearch());
		EXPECT_FLOAT_EQ(2.0f, dedispersion_strategy.get_power());
		EXPECT_FLOAT_EQ(7.0f, dedispersion_strategy.get_sigma_cutoff());
		EXPECT_FLOAT_EQ(0.1f, dedispersion_strategy.get_wide());
		EXPECT_EQ(6, dedispersion_strategy.get_range());
		// dm:
		EXPECT_FLOAT_EQ(0.0, dedispersion_strategy.get_user_dm_low()[0]);
		EXPECT_FLOAT_EQ(150.0, dedispersion_strategy.get_user_dm_low()[1]);
		EXPECT_FLOAT_EQ(300.0, dedispersion_strategy.get_user_dm_low()[2]);
		EXPECT_FLOAT_EQ(500.0, dedispersion_strategy.get_user_dm_low()[3]);
		EXPECT_FLOAT_EQ(900.0, dedispersion_strategy.get_user_dm_low()[4]);
		EXPECT_FLOAT_EQ(1200.0, dedispersion_strategy.get_user_dm_low()[5]);
		EXPECT_FLOAT_EQ(150.0, dedispersion_strategy.get_user_dm_high()[0]);
		EXPECT_FLOAT_EQ(300.0, dedispersion_strategy.get_user_dm_high()[1]);
		EXPECT_FLOAT_EQ(500.0, dedispersion_strategy.get_user_dm_high()[2]);
		EXPECT_FLOAT_EQ(900.0, dedispersion_strategy.get_user_dm_high()[3]);
		EXPECT_FLOAT_EQ(1200.0, dedispersion_strategy.get_user_dm_high()[4]);
		EXPECT_FLOAT_EQ(1500.0, dedispersion_strategy.get_user_dm_high()[5]);
		EXPECT_FLOAT_EQ(0.1, dedispersion_strategy.get_user_dm_step()[0]);
		EXPECT_FLOAT_EQ(0.2, dedispersion_strategy.get_user_dm_step()[1]);
		EXPECT_FLOAT_EQ(0.25, dedispersion_strategy.get_user_dm_step()[2]);
		EXPECT_FLOAT_EQ(0.4, dedispersion_strategy.get_user_dm_step()[3]);
		EXPECT_FLOAT_EQ(0.6, dedispersion_strategy.get_user_dm_step()[4]);
		EXPECT_FLOAT_EQ(0.8, dedispersion_strategy.get_user_dm_step()[5]);
		//
		EXPECT_EQ(1, dedispersion_strategy.get_in_bin()[0]);
		EXPECT_EQ(1, dedispersion_strategy.get_in_bin()[1]);
		EXPECT_EQ(2, dedispersion_strategy.get_in_bin()[2]);
		EXPECT_EQ(2, dedispersion_strategy.get_in_bin()[3]);
		EXPECT_EQ(4, dedispersion_strategy.get_in_bin()[4]);
		EXPECT_EQ(4, dedispersion_strategy.get_in_bin()[5]);
		EXPECT_EQ(1, dedispersion_strategy.get_out_bin()[0]);
		EXPECT_EQ(1, dedispersion_strategy.get_out_bin()[1]);
		EXPECT_EQ(2, dedispersion_strategy.get_out_bin()[2]);
		EXPECT_EQ(2, dedispersion_strategy.get_out_bin()[3]);
		EXPECT_EQ(4, dedispersion_strategy.get_out_bin()[4]);
		EXPECT_EQ(4, dedispersion_strategy.get_out_bin()[5]);
		//
		EXPECT_EQ(4096, dedispersion_strategy.get_nchans());
		EXPECT_EQ(0, dedispersion_strategy.get_nsamples());
		EXPECT_EQ(234496, dedispersion_strategy.get_nsamp());
		EXPECT_EQ(1, dedispersion_strategy.get_nifs());
		EXPECT_EQ(8, dedispersion_strategy.get_nbits());
		EXPECT_FLOAT_EQ(0.000064, dedispersion_strategy.get_tsamp());
		EXPECT_FLOAT_EQ(50000, dedispersion_strategy.get_tstart());
		EXPECT_FLOAT_EQ(1550, dedispersion_strategy.get_fch1());
		EXPECT_NEAR(-0.073242, dedispersion_strategy.get_foff(), 0.000001);
		//
		EXPECT_NEAR(0.000000, dedispersion_strategy.get_dmshifts()[0], 0.000001);
		EXPECT_NEAR(0.000002, dedispersion_strategy.get_dmshifts()[10], 0.000001);
		EXPECT_NEAR(0.000003, dedispersion_strategy.get_dmshifts()[20], 0.000001);
		EXPECT_NEAR(0.000005, dedispersion_strategy.get_dmshifts()[30], 0.000001);
		EXPECT_NEAR(0.000007, dedispersion_strategy.get_dmshifts()[40], 0.000001);
		EXPECT_NEAR(0.000008, dedispersion_strategy.get_dmshifts()[50], 0.000001);
		//
		EXPECT_EQ(1520, dedispersion_strategy.get_max_ndms());
		EXPECT_FLOAT_EQ(1562, dedispersion_strategy.get_max_dm());
		EXPECT_EQ(5080, dedispersion_strategy.get_total_ndms());
		EXPECT_EQ(22752, dedispersion_strategy.get_maxshift());
		EXPECT_EQ(6, dedispersion_strategy.get_num_tchunks());
		//
		EXPECT_FLOAT_EQ(0.0, dedispersion_strategy.get_dm_low()[0]);
		EXPECT_FLOAT_EQ(152.0, dedispersion_strategy.get_dm_high()[0]);
		EXPECT_FLOAT_EQ(0.1, dedispersion_strategy.get_dm_step()[0]);
		EXPECT_FLOAT_EQ(152.0, dedispersion_strategy.get_dm_low()[1]);
		EXPECT_FLOAT_EQ(304.0, dedispersion_strategy.get_dm_high()[1]);
		EXPECT_FLOAT_EQ(0.2, dedispersion_strategy.get_dm_step()[1]);
		EXPECT_FLOAT_EQ(304.0, dedispersion_strategy.get_dm_low()[2]);
		EXPECT_FLOAT_EQ(514.0, dedispersion_strategy.get_dm_high()[2]);
		EXPECT_FLOAT_EQ(0.25, dedispersion_strategy.get_dm_step()[2]);
		EXPECT_FLOAT_EQ(514.0, dedispersion_strategy.get_dm_low()[3]);
		EXPECT_FLOAT_EQ(930.0, dedispersion_strategy.get_dm_high()[3]);
		EXPECT_FLOAT_EQ(0.4, dedispersion_strategy.get_dm_step()[3]);
		EXPECT_FLOAT_EQ(930.0, dedispersion_strategy.get_dm_low()[4]);
		EXPECT_FLOAT_EQ(1242.0, dedispersion_strategy.get_dm_high()[4]);
		EXPECT_FLOAT_EQ(0.6, dedispersion_strategy.get_dm_step()[4]);
		EXPECT_FLOAT_EQ(1242.0, dedispersion_strategy.get_dm_low()[5]);
		EXPECT_FLOAT_EQ(1562.0, dedispersion_strategy.get_dm_high()[5]);
		EXPECT_FLOAT_EQ(0.8, dedispersion_strategy.get_dm_step()[5]);
		//
		EXPECT_EQ(40896, dedispersion_strategy.get_t_processed()[0][0]);
		EXPECT_EQ(40896, dedispersion_strategy.get_t_processed()[1][0]);
		EXPECT_EQ(20448, dedispersion_strategy.get_t_processed()[2][0]);
		EXPECT_EQ(20448, dedispersion_strategy.get_t_processed()[3][0]);
		EXPECT_EQ(10224, dedispersion_strategy.get_t_processed()[4][0]);
		EXPECT_EQ(10224, dedispersion_strategy.get_t_processed()[5][0]);
		//
		fclose(fp);
		free(inBin);
		free(outBin);
		free(user_dm_low);
		free(user_dm_high);
		free(user_dm_step);
	}
}
*/
// test AstroAccelerate class
TEST_F(AstroAccelerateTest, AstroAccelerate_call)
{
	// Following is ok for ska_karel.txt
	char* filename = my_argv[1] + strlen(my_argv[1]) - 13;
	if(strcmp(filename, "ska_karel.txt") == 0)
	{
		// Internal code variables
		// File pointers
		FILE *fp = nullptr;
		// Counters and flags
		int range = 0;
		int enable_debug = 0;
		int enable_analysis = 0;
		int enable_acceleration = 0;
		int enable_periodicity = 0;
		int output_dmt = 0;
		int enable_zero_dm = 0;
		int enable_zero_dm_with_outliers = 0;
		int enable_rfi = 0;
		int enable_fdas_custom_fft = 0;
		int enable_fdas_inbin = 0;
		int enable_fdas_norm = 0;
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
		// Memory sizes and pointers
		float *user_dm_low = nullptr;
		float *user_dm_high = nullptr;
		float *user_dm_step = nullptr;
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
		int maxshift_original;
		double tsamp_original;
		long int inc = 0;
		float tstart = 0.0f;
		float tstart_local = 0.0f;
		float tsamp = 0.0f;
		float fch1 = 0.0f;
		float foff = 0.0f;
		// Analysis variables
		float power = 2.0f;
		float sigma_cutoff = 6.0f;
		float sigma_constant = 4.0f;
		float max_boxcar_width_in_sec = 0.5f;
		// Initialise the GPU.
		int device_id = 0; // hard-coded, would be a parameter
		size_t gpu_memory = 0;
		cudaSetDevice(device_id);
		size_t mem_free, total;
		cudaMemGetInfo(&mem_free, &total);
		gpu_memory = ( mem_free/2 );

		// Users desired de-dispersion strategy. Pick up user defined values from the CLI.
		get_user_input(&fp, my_argc, my_argv, &multi_file, &enable_debug, &enable_analysis,
		&enable_periodicity, &enable_acceleration, &output_dmt, &enable_zero_dm,
		&enable_zero_dm_with_outliers, &enable_rfi, &enable_fdas_custom_fft,
		&enable_fdas_inbin, &enable_fdas_norm, &nboots, &ntrial_bins, &navdms,
		&narrow, &wide, &aggression, &nsearch, &inBin, &outBin, &power, &sigma_cutoff, &sigma_constant, &max_boxcar_width_in_sec,
		&range, &user_dm_low, &user_dm_high, &user_dm_step, &candidate_algorithm);
		// Reads telescope parameters from the header of the input file and then counts the number of samples in the input data file.
		get_file_data(&fp, &nchans, &nsamples, &nsamp, &nifs, &nbits, &tsamp, &tstart,
		&fch1, &foff);

		// dedispersion
		//DedispersionStrategy dedispersion_strategy;
		//DedispersionStrategyFile(&fp, my_argc, my_argv, dedispersion_strategy, gpu_memory);
		DedispersionStrategy dedispersion_strategy
									 (user_dm_low
									 ,user_dm_high
									 ,user_dm_step
									 ,inBin
									 ,outBin
									 ,gpu_memory
									 ,power
									 ,range
									 ,nchans
									 ,nsamples
									 ,nsamp
									 ,nifs
									 ,nbits
									 ,tsamp
									 ,tstart
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

		printf("\nAA is starting\n");
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
		//*/
		fclose(fp);
		free(inBin);
		free(outBin);
		free(user_dm_low);
		free(user_dm_high);
		free(user_dm_step);
		free(input_buffer);
	}
}

} // namespace test
} // namespace astroaccelerate

