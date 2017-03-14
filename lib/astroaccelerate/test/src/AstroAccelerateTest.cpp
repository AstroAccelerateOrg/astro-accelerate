#include "AstroAccelerateTest.h"
#include "../../AstroAccelerateParameters.h"
#include "../../DedispersionStrategy.h"
#include "../../DmTime.h"
#include "../../AstroAccelerate.h"

#include <vector>

namespace ska {
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

// test dedispersion strategy class
TEST_F(AstroAccelerateTest, test_dedispersion_strategy)
{
	// Following is ok for ska_karel.txt
	char* filename = my_argv[1] + strlen(my_argv[1]) - 13;
	if(strcmp(filename, "ska_karel.txt") == 0)
	{
		// declare objects
		DedispersionStrategy dedispersion_strategy;
		// read user input
		FILE *fp = nullptr;
		dedispersion_strategy.get_user_input(&fp, my_argc, my_argv);
		// check class members values after run
		EXPECT_EQ(1, dedispersion_strategy.get_multi_file());
		EXPECT_EQ(1, dedispersion_strategy.get_enable_debug());
		EXPECT_EQ(1, dedispersion_strategy.get_enable_analysis());
		EXPECT_EQ(0, dedispersion_strategy.get_enable_periodicity());
		//EXPECT_EQ(1, dedispersion_strategy.get_enable_acceleration());
		EXPECT_EQ(0, dedispersion_strategy.get_output_dmt());
		EXPECT_EQ(0, dedispersion_strategy.get_enable_zero_dm());
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
		// get file data
		dedispersion_strategy.get_file_data(&fp);
		// check if it updates correctly
		EXPECT_EQ(4096, dedispersion_strategy.get_nchans());
		EXPECT_EQ(0, dedispersion_strategy.get_nsamples());
		EXPECT_EQ(234496, dedispersion_strategy.get_nsamp());
		EXPECT_EQ(1, dedispersion_strategy.get_nifs());
		EXPECT_EQ(8, dedispersion_strategy.get_nbits());
		EXPECT_FLOAT_EQ(0.000064, dedispersion_strategy.get_tsamp());
		EXPECT_FLOAT_EQ(50000, dedispersion_strategy.get_tstart());
		EXPECT_FLOAT_EQ(1550, dedispersion_strategy.get_fch1());
		EXPECT_NEAR(-0.073242, dedispersion_strategy.get_foff(), 0.000001);
		//EXPECT_FLOAT_EQ(-0.073242, dedispersion_strategy.get_foff());
		// it fails the test:
		// Value of: dedispersion_strategy.get_foff()
		// Actual: -0.073242188
		// Expected: -0.073242
		// Which is: -0.073242001
		//
		// Initialise the GPU.
		int device_id = 0; // hard-coded, would be a parameter
		size_t gpu_memory = 0;
		cudaSetDevice(device_id);
		size_t mem_free, total;
		cudaMemGetInfo(&mem_free, &total);
		gpu_memory = ( mem_free/4 );
		// Call the strategy method
		dedispersion_strategy.make_strategy(gpu_memory);
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
		//*/
		fclose(fp);
	}
}

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
		//
		astroaccelerate::DedispersionStrategy dedispersion_strategy;
		// Input buffer
		size_t inputsize = 0;
		unsigned short *input_buffer = NULL;
		// read user input
		dedispersion_strategy.get_user_input(&fp, my_argc, my_argv);
		// get file data
		dedispersion_strategy.get_file_data(&fp);
		//
		// Allocate memory on host.
		allocate_memory_cpu_input(dedispersion_strategy.get_nsamp()
								  ,dedispersion_strategy.get_nchans()
								  ,&input_buffer
								  ,&inputsize);
		// Store the recorded telescope data contained in the input filterbank file
		// in the allocated memory.
		get_recorded_data(&fp
						  ,dedispersion_strategy.get_nsamp()
						  ,dedispersion_strategy.get_nchans()
						  ,dedispersion_strategy.get_nbits()
						  ,&input_buffer
						  ,&inputsize);
		//
		//
		int device_id = 0; // hard-coded, would be a parameter
		size_t gpu_memory = 0;
		cudaSetDevice(device_id);
		size_t mem_free, total;
		cudaMemGetInfo(&mem_free, &total);
		gpu_memory = ( mem_free/2 );
		// Call the strategy method
		dedispersion_strategy.make_strategy(gpu_memory);
		//
		std::vector<float> output_sps;

		// dedispersed data
		DmTime<float> output_buffer(dedispersion_strategy);

		// call AstroAccelerate main method here
		astroaccelerate::AstroAccelerate<TestParams> astroaccelerate(dedispersion_strategy);
		astroaccelerate.run_dedispersion_sps_fdas(device_id
												  ,dedispersion_strategy
												  ,input_buffer
												  ,output_buffer
												  ,output_sps
												  );

		// close file
		fclose(fp);

		// write output here, not in the library
	}
}

} // namespace test
} // namespace astroaccelerate
}
