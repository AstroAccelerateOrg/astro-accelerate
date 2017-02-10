#include "SpsTest.h"
#include "../../SpsParameters.h"
#include "../../InputData.h"
#include "../../UserInput.h"
#include "../../DedispersionPlan.h"
#include "../../Sps.h"

#include <vector>


namespace ska {
namespace astroaccelerate {
namespace sps {
namespace test {

SpsTest::SpsTest()
    : ::testing::Test()
{
}

SpsTest::~SpsTest()
{
}

void SpsTest::SetUp()
{
}

void SpsTest::TearDown()
{
}

class TestParams : public SpsParameters<TestParams> {};

// test class user input class
TEST_F(SpsTest, test_user_input)
{
	// Following is ok for ska_karel.txt
	char* filename = my_argv[1] + strlen(my_argv[1]) - 13;
	if(strcmp(filename, "ska_karel.txt") == 0)
	{
		// create object
		sps::UserInput user_input;

		// first, check constructor:
		EXPECT_EQ(1, user_input.get_multi_file());
		EXPECT_EQ(0, user_input.get_enable_debug());
		EXPECT_EQ(0, user_input.get_enable_analysis());
		EXPECT_EQ(0, user_input.get_enable_periodicity());
		EXPECT_EQ(0, user_input.get_enable_acceleration());
		EXPECT_EQ(0, user_input.get_output_dmt());
		EXPECT_EQ(0, user_input.get_enable_zero_dm());
		EXPECT_EQ(-1, user_input.get_nboots());
		EXPECT_EQ(0, user_input.get_ntrial_bins());
		EXPECT_EQ(1, user_input.get_navdms());
		EXPECT_EQ(0.001f, user_input.get_narrow());
		EXPECT_EQ(2.5, user_input.get_aggression());
		EXPECT_EQ(3, user_input.get_nsearch());
		EXPECT_EQ(2.0f, user_input.get_power());
		EXPECT_EQ(6.0f, user_input.get_sigma_cutoff());
		EXPECT_EQ(0.1f, user_input.get_wide());
		EXPECT_EQ(0, user_input.get_range());
		EXPECT_EQ(NULL, user_input.get_user_dm_low());
		EXPECT_EQ(NULL, user_input.get_user_dm_high());
		EXPECT_EQ(NULL, user_input.get_user_dm_step());

		// read user input
		FILE *fp = NULL;
		user_input.get_user_input(&fp, my_argc, my_argv);

		// check class members values after run
		EXPECT_EQ(1, user_input.get_multi_file());
		EXPECT_EQ(1, user_input.get_enable_debug());
		EXPECT_EQ(1, user_input.get_enable_analysis());
		EXPECT_EQ(0, user_input.get_enable_periodicity());
		EXPECT_EQ(0, user_input.get_enable_acceleration());
		EXPECT_EQ(0, user_input.get_output_dmt());
		EXPECT_EQ(0, user_input.get_enable_zero_dm());
		EXPECT_EQ(-1, user_input.get_nboots());
		EXPECT_EQ(0, user_input.get_ntrial_bins());
		EXPECT_EQ(1, user_input.get_navdms());
		EXPECT_FLOAT_EQ(0.001f, user_input.get_narrow());
		EXPECT_FLOAT_EQ(2.5, user_input.get_aggression());
		EXPECT_EQ(3, user_input.get_nsearch());
		EXPECT_FLOAT_EQ(2.0f, user_input.get_power());
		EXPECT_FLOAT_EQ(10.0f, user_input.get_sigma_cutoff());
		EXPECT_FLOAT_EQ(0.1f, user_input.get_wide());
		EXPECT_EQ(6, user_input.get_range());
		// dm:
		EXPECT_FLOAT_EQ(0.0, user_input.get_user_dm_low()[0]);
		EXPECT_FLOAT_EQ(150.0, user_input.get_user_dm_low()[1]);
		EXPECT_FLOAT_EQ(300.0, user_input.get_user_dm_low()[2]);
		EXPECT_FLOAT_EQ(500.0, user_input.get_user_dm_low()[3]);
		EXPECT_FLOAT_EQ(900.0, user_input.get_user_dm_low()[4]);
		EXPECT_FLOAT_EQ(1200.0, user_input.get_user_dm_low()[5]);
		EXPECT_FLOAT_EQ(150.0, user_input.get_user_dm_high()[0]);
		EXPECT_FLOAT_EQ(300.0, user_input.get_user_dm_high()[1]);
		EXPECT_FLOAT_EQ(500.0, user_input.get_user_dm_high()[2]);
		EXPECT_FLOAT_EQ(900.0, user_input.get_user_dm_high()[3]);
		EXPECT_FLOAT_EQ(1200.0, user_input.get_user_dm_high()[4]);
		EXPECT_FLOAT_EQ(1500.0, user_input.get_user_dm_high()[5]);
		EXPECT_FLOAT_EQ(0.1, user_input.get_user_dm_step()[0]);
		EXPECT_FLOAT_EQ(0.2, user_input.get_user_dm_step()[1]);
		EXPECT_FLOAT_EQ(0.25, user_input.get_user_dm_step()[2]);
		EXPECT_FLOAT_EQ(0.4, user_input.get_user_dm_step()[3]);
		EXPECT_FLOAT_EQ(0.6, user_input.get_user_dm_step()[4]);
		EXPECT_FLOAT_EQ(0.8, user_input.get_user_dm_step()[5]);
		//
		EXPECT_EQ(1, user_input.get_in_bin()[0]);
		EXPECT_EQ(1, user_input.get_in_bin()[1]);
		EXPECT_EQ(2, user_input.get_in_bin()[2]);
		EXPECT_EQ(2, user_input.get_in_bin()[3]);
		EXPECT_EQ(4, user_input.get_in_bin()[4]);
		EXPECT_EQ(4, user_input.get_in_bin()[5]);
		EXPECT_EQ(1, user_input.get_out_bin()[0]);
		EXPECT_EQ(1, user_input.get_out_bin()[1]);
		EXPECT_EQ(2, user_input.get_out_bin()[2]);
		EXPECT_EQ(2, user_input.get_out_bin()[3]);
		EXPECT_EQ(4, user_input.get_out_bin()[4]);
		EXPECT_EQ(4, user_input.get_out_bin()[5]);


		fclose(fp);
	}
}

// test dedispersion plan class
TEST_F(SpsTest, test_dedispersion_plan)
{
	// Following is ok for ska_karel.txt
	char* filename = my_argv[1] + strlen(my_argv[1]) - 13;
	if(strcmp(filename, "ska_karel.txt") == 0)
	{
		// declare objects
		sps::UserInput user_input;
		sps::DedispersionPlan dedispersion_plan;
		// first, check constructor
		EXPECT_EQ(0   , dedispersion_plan.get_maxshift());
		EXPECT_EQ(NULL, dedispersion_plan.get_dm_low());
		EXPECT_EQ(NULL, dedispersion_plan.get_dm_high());
		EXPECT_EQ(NULL, dedispersion_plan.get_dm_step());
		EXPECT_EQ(NULL, dedispersion_plan.get_dmshifts());
		EXPECT_EQ(NULL, dedispersion_plan.get_ndms());
		EXPECT_EQ(0,    dedispersion_plan.get_max_ndms());
		EXPECT_EQ(0,    dedispersion_plan.get_total_ndms());
		EXPECT_FLOAT_EQ(0.0f, dedispersion_plan.get_max_dm());
		EXPECT_EQ(0, dedispersion_plan.get_range());
		EXPECT_EQ(NULL, dedispersion_plan.get_t_processed());
		EXPECT_EQ(0, dedispersion_plan.get_nbits());
		EXPECT_EQ(0, dedispersion_plan.get_nifs());
		EXPECT_FLOAT_EQ(0.0f, dedispersion_plan.get_tstart());
		EXPECT_FLOAT_EQ(0.0f, dedispersion_plan.get_tsamp());
		EXPECT_EQ(0, dedispersion_plan.get_nsamp());
		EXPECT_EQ(0, dedispersion_plan.get_nsamples());
		EXPECT_EQ(0, dedispersion_plan.get_max_samps());
		EXPECT_EQ(0, dedispersion_plan.get_nchans());
		EXPECT_FLOAT_EQ(0.0f, dedispersion_plan.get_fch1());
		EXPECT_FLOAT_EQ(0.0f, dedispersion_plan.get_foff());
		EXPECT_EQ(0, dedispersion_plan.get_num_tchunks());
		// read user input
		FILE *fp = NULL;
		user_input.get_user_input(&fp, my_argc, my_argv);
		// set dedispersion plan
		dedispersion_plan.set_power(user_input.get_power());
		dedispersion_plan.set_range(user_input.get_range());
		// chech if it updates correctly
		EXPECT_FLOAT_EQ(2.0f, dedispersion_plan.get_power());
		EXPECT_EQ(6, dedispersion_plan.get_range());
		// get file data
		dedispersion_plan.get_file_data(&fp);
		// check if it updates correctly
		EXPECT_EQ(4096, dedispersion_plan.get_nchans());
		EXPECT_EQ(0, dedispersion_plan.get_nsamples());
		EXPECT_EQ(234496, dedispersion_plan.get_nsamp());
		EXPECT_EQ(1, dedispersion_plan.get_nifs());
		EXPECT_EQ(8, dedispersion_plan.get_nbits());
		EXPECT_FLOAT_EQ(0.000064, dedispersion_plan.get_tsamp());
		EXPECT_FLOAT_EQ(50000, dedispersion_plan.get_tstart());
		EXPECT_FLOAT_EQ(1550, dedispersion_plan.get_fch1());
		EXPECT_NEAR(-0.073242, dedispersion_plan.get_foff(), 0.000001);
		//EXPECT_FLOAT_EQ(-0.073242, dedispersion_plan.get_foff());
		// it fails the test:
		// Value of: dedispersion_plan.get_foff()
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
		dedispersion_plan.make_strategy(user_input.get_user_dm_low(),
										user_input.get_user_dm_high(),
										user_input.get_user_dm_step(),
										user_input.get_in_bin(),
										gpu_memory
										);
		//
		EXPECT_NEAR(0.000000, dedispersion_plan.get_dmshifts()[0], 0.000001);
		EXPECT_NEAR(0.000002, dedispersion_plan.get_dmshifts()[10], 0.000001);
		EXPECT_NEAR(0.000003, dedispersion_plan.get_dmshifts()[20], 0.000001);
		EXPECT_NEAR(0.000005, dedispersion_plan.get_dmshifts()[30], 0.000001);
		EXPECT_NEAR(0.000007, dedispersion_plan.get_dmshifts()[40], 0.000001);
		EXPECT_NEAR(0.000008, dedispersion_plan.get_dmshifts()[50], 0.000001);
		//
		EXPECT_EQ(1520, dedispersion_plan.get_max_ndms());
		EXPECT_FLOAT_EQ(1562, dedispersion_plan.get_max_dm());
		EXPECT_EQ(5080, dedispersion_plan.get_total_ndms());
		EXPECT_EQ(22880, dedispersion_plan.get_maxshift());
		EXPECT_EQ(6, dedispersion_plan.get_num_tchunks());
		//
		EXPECT_FLOAT_EQ(0.0, dedispersion_plan.get_dm_low()[0]);
		EXPECT_FLOAT_EQ(152.0, dedispersion_plan.get_dm_high()[0]);
		EXPECT_FLOAT_EQ(0.1, dedispersion_plan.get_dm_step()[0]);
		EXPECT_FLOAT_EQ(152.0, dedispersion_plan.get_dm_low()[1]);
		EXPECT_FLOAT_EQ(304.0, dedispersion_plan.get_dm_high()[1]);
		EXPECT_FLOAT_EQ(0.2, dedispersion_plan.get_dm_step()[1]);
		EXPECT_FLOAT_EQ(304.0, dedispersion_plan.get_dm_low()[2]);
		EXPECT_FLOAT_EQ(514.0, dedispersion_plan.get_dm_high()[2]);
		EXPECT_FLOAT_EQ(0.25, dedispersion_plan.get_dm_step()[2]);
		EXPECT_FLOAT_EQ(514.0, dedispersion_plan.get_dm_low()[3]);
		EXPECT_FLOAT_EQ(930.0, dedispersion_plan.get_dm_high()[3]);
		EXPECT_FLOAT_EQ(0.4, dedispersion_plan.get_dm_step()[3]);
		EXPECT_FLOAT_EQ(930.0, dedispersion_plan.get_dm_low()[4]);
		EXPECT_FLOAT_EQ(1242.0, dedispersion_plan.get_dm_high()[4]);
		EXPECT_FLOAT_EQ(0.6, dedispersion_plan.get_dm_step()[4]);
		EXPECT_FLOAT_EQ(1242.0, dedispersion_plan.get_dm_low()[5]);
		EXPECT_FLOAT_EQ(1562.0, dedispersion_plan.get_dm_high()[5]);
		EXPECT_FLOAT_EQ(0.8, dedispersion_plan.get_dm_step()[5]);
		//
		EXPECT_EQ(40960, dedispersion_plan.get_t_processed()[0][0]);
		EXPECT_EQ(40960, dedispersion_plan.get_t_processed()[1][0]);
		EXPECT_EQ(20480, dedispersion_plan.get_t_processed()[2][0]);
		EXPECT_EQ(20480, dedispersion_plan.get_t_processed()[3][0]);
		EXPECT_EQ(10240, dedispersion_plan.get_t_processed()[4][0]);
		EXPECT_EQ(10240, dedispersion_plan.get_t_processed()[5][0]);
		//
		fclose(fp);
	}
}

// test input/output data class
TEST_F(SpsTest, test_InputData)
{
	// Following is ok for ska_karel.txt
	char* filename = my_argv[1] + strlen(my_argv[1]) - 13;
	if(strcmp(filename, "ska_karel.txt") == 0)
	{
		// declare objects
		sps::UserInput user_input;
		sps::DedispersionPlan dedispersion_plan;
		sps::InputData input_data;
		// read user input
		FILE *fp = NULL;
		user_input.get_user_input(&fp, my_argc, my_argv);
		// set dedispersion plan
		dedispersion_plan.set_power(user_input.get_power());
		dedispersion_plan.set_range(user_input.get_range());
		// get file data
		dedispersion_plan.get_file_data(&fp);
		// Initialise the GPU.
		int device_id = 0; // hard-coded, would be a parameter
		size_t gpu_memory = 0;
		cudaSetDevice(device_id);
		size_t mem_free, total;
		cudaMemGetInfo(&mem_free, &total);
		gpu_memory = ( mem_free/4 );
		// Call the strategy method
		dedispersion_plan.make_strategy(user_input.get_user_dm_low(),
										user_input.get_user_dm_high(),
										user_input.get_user_dm_step(),
										user_input.get_in_bin(),
										gpu_memory
										);
		// allocate memory cpu input
		input_data.allocate_memory_cpu_input(dedispersion_plan);
		EXPECT_EQ(1832, (int)(input_data.get_input_size() / 1024 / 1024));
		// get recorded data
		input_data.get_recorded_data(&fp, dedispersion_plan);
		//
		EXPECT_EQ(113, input_data.get_input_buffer()[0]);
		EXPECT_EQ(129, input_data.get_input_buffer()[500]);
		EXPECT_EQ(130, input_data.get_input_buffer()[1000]);
		EXPECT_EQ(150, input_data.get_input_buffer()[1500]);
		EXPECT_EQ(145, input_data.get_input_buffer()[2000]);
		EXPECT_EQ(93,  input_data.get_input_buffer()[2500]);
		EXPECT_EQ(150, input_data.get_input_buffer()[3000]);
		EXPECT_EQ(186, input_data.get_input_buffer()[3500]);
		EXPECT_EQ(181, input_data.get_input_buffer()[4000]);
		//
		fclose(fp);
	}
}

// test sps class
TEST_F(SpsTest, sps_call)
{
	// Following is ok for ska_karel.txt
	char* filename = my_argv[1] + strlen(my_argv[1]) - 13;
	if(strcmp(filename, "ska_karel.txt") == 0)
	{
 
		// Internal code variables
		// File pointers
		FILE *fp = NULL;
		// declare objects
		ska::astroaccelerate::sps::UserInput user_input;
		// read user input
		user_input.get_user_input(&fp, my_argc, my_argv);

		ska::astroaccelerate::sps::DedispersionPlan dedispersion_plan;
		// set dedispersion plan
		dedispersion_plan.set_power(user_input.get_power());
		dedispersion_plan.set_range(user_input.get_range());
		// get file data
		dedispersion_plan.get_file_data(&fp);

		ska::astroaccelerate::sps::InputData input_data;
		// allocate memory cpu input
		input_data.allocate_memory_cpu_input(dedispersion_plan);
		// get recorded data
		input_data.get_recorded_data(&fp, dedispersion_plan);

		//
		int device_id = 0; // hard-coded, would be a parameter
		size_t gpu_memory = 0;
		cudaSetDevice(device_id);
		size_t mem_free, total;
		cudaMemGetInfo(&mem_free, &total);
		gpu_memory = ( mem_free/2 );

		// Call the strategy method
		dedispersion_plan.make_strategy(user_input.get_user_dm_low(),
										user_input.get_user_dm_high(),
										user_input.get_user_dm_step(),
										user_input.get_in_bin(),
										gpu_memory
										);

		std::vector<float> output_sps;

		// call sps main method here
		ska::astroaccelerate::sps::Sps<TestParams> sps(input_data, dedispersion_plan, user_input);
		sps(device_id, input_data, dedispersion_plan, user_input, gpu_memory, output_sps);


		printf("\nHost DD  output size:\t\t%d MB", (int) (sps.get_output_size() / 1024 / 1024 / sizeof(float)));
		printf("\nHost SPS output size:\t\t%d MB", (int) (output_sps.size() / 1024 / 1024));

		for (int i = 0; i < 100; i+=4)
		{
			printf("\n%f %f %f %f", output_sps[i], output_sps[i+1], output_sps[i+2], output_sps[i+3]);
		}

		// close file
		fclose(fp);

		// write output here, not in the library


	}
}


} // namespace test
} // namespace sps
} // namespace astroaccelerate
} // namespace ska
