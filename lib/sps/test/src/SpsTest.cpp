#include "SpsTest.h"
#include "../../SpsParameters.h"
#include "../../IOData.h"
#include "../../UserInput.h"
#include "../../DedispersionPlan.h"
#include "../../Sps.h"

#include "../../../AstroAccelerate/host_get_file_data.h"

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

// example on how gtest works
/*TEST_F(SpsTest, test_handlers)
{
    unsigned device_id = 0;
    sps::DedispersionPlan dedispersion_plan;
    sps::IOData io_data;
    bool sps_handler_called = false;
    bool dm_handler_called = false;
    {
        sps::Sps<TestParams> sps;
        //sps(device_id, io_data, dedispersion_plan, [&]() { sps_handler_called = true; }, [&]() { dm_handler_called =true;} );
        //sps(device_id, io_data, plan);
    }
    ASSERT_FALSE(dm_handler_called);
    ASSERT_FALSE(sps_handler_called);
}
*/

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
		EXPECT_FLOAT_EQ(7.0f, user_input.get_sigma_cutoff());
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
		int nchans = 0;
		int nsamp = 0;
		int nbits = 0;
		int nsamples = 0;
		int nifs = 0;
		float tsamp = 0.0f;
		float tstart = 0.0f;
		float fch1 = 0.0f;
		float foff = 0.0f;
		printf("\n\n ============ output of get_file_data =========\n\n");
		get_file_data(&fp, &nchans, &nsamples, &nsamp, &nifs, &nbits, &tsamp,
					  &tstart, &fch1, &foff);
		printf("\n\n ==============================================\n\n");
		// set dedispersion plan
		dedispersion_plan.set_nchans(nchans);
		dedispersion_plan.set_nsamples(nsamples);
		dedispersion_plan.set_nsamp(nsamp);
		dedispersion_plan.set_nifs(nifs);
		dedispersion_plan.set_nbits(nbits);
		dedispersion_plan.set_tsamp(tsamp);
		dedispersion_plan.set_tstart(tstart);
		dedispersion_plan.set_fch1(fch1);
		dedispersion_plan.set_foff(foff);
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
		fclose(fp);
	}
}

// test input/output data class
TEST_F(SpsTest, test_iodata)
{
	// Following is ok for ska_karel.txt
	char* filename = my_argv[1] + strlen(my_argv[1]) - 13;
	if(strcmp(filename, "ska_karel.txt") == 0)
	{
		// declare objects
		sps::UserInput user_input;
		sps::DedispersionPlan dedispersion_plan;
		sps::IOData io_data;
		// read user input
		FILE *fp = NULL;
		user_input.get_user_input(&fp, my_argc, my_argv);
		// set dedispersion plan
		dedispersion_plan.set_power(user_input.get_power());
		dedispersion_plan.set_range(user_input.get_range());
		// get file data
		int nchans = 0;
		int nsamp = 0;
		int nbits = 0;
		int nsamples = 0;
		int nifs = 0;
		float tsamp = 0.0f;
		float tstart = 0.0f;
		float fch1 = 0.0f;
		float foff = 0.0f;
		printf("\n\n ============ output of get_file_data =========\n\n");
		get_file_data(&fp, &nchans, &nsamples, &nsamp, &nifs, &nbits, &tsamp,
					  &tstart, &fch1, &foff);
		printf("\n\n ==============================================\n\n");
		// set dedispersion plan
		dedispersion_plan.set_nchans(nchans);
		dedispersion_plan.set_nsamples(nsamples);
		dedispersion_plan.set_nsamp(nsamp);
		dedispersion_plan.set_nifs(nifs);
		dedispersion_plan.set_nbits(nbits);
		dedispersion_plan.set_tsamp(tsamp);
		dedispersion_plan.set_tstart(tstart);
		dedispersion_plan.set_fch1(fch1);
		dedispersion_plan.set_foff(foff);
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
		io_data.allocate_memory_cpu_input(dedispersion_plan);
		EXPECT_EQ(1832, (int)(io_data.get_input_size() / 1024 / 1024));
		// get recorded data
		// check this one
		io_data.get_recorded_data(&fp, dedispersion_plan.get_nchans(),dedispersion_plan.get_nbits());

		// allocate memory cpu output
		io_data.allocate_memory_cpu_output(dedispersion_plan);
		EXPECT_EQ(2784, (int)(io_data.get_output_size() / 1024 / 1024));
		// allocate memory gpu
		io_data.allocate_memory_gpu(dedispersion_plan);
		EXPECT_EQ(498, (int)(io_data.get_gpu_input_size() / 1024 / 1024));
		EXPECT_EQ(997, (int)(io_data.get_gpu_output_size() / 1024 / 1024));
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
		/*
		// declare objects
		sps::UserInput user_input;
		sps::DedispersionPlan dedispersion_plan;
		sps::IOData io_data;
		//sps::Sps<TestParams> sps_object;
		sps::Sps<TestParams> sps_object;
		// read user input
		FILE *fp = NULL;
		user_input.get_user_input(&fp, my_argc, my_argv);
		// set dedispersion plan
		dedispersion_plan.set_power(user_input.get_power());
		dedispersion_plan.set_range(user_input.get_range());
		// get file data
		int nchans = 0;
		int nsamp = 0;
		int nbits = 0;
		int nsamples = 0;
		int nifs = 0;
		float tsamp = 0.0f;
		float tstart = 0.0f;
		float fch1 = 0.0f;
		float foff = 0.0f;
		printf("\n\n ============ output of get_file_data =========\n\n");
		get_file_data(&fp, &nchans, &nsamples, &nsamp, &nifs, &nbits, &tsamp,
					  &tstart, &fch1, &foff);
		printf("\n\n ==============================================\n\n");
		// set dedispersion plan
		dedispersion_plan.set_nchans(nchans);
		dedispersion_plan.set_nsamples(nsamples);
		dedispersion_plan.set_nsamp(nsamp);
		dedispersion_plan.set_nifs(nifs);
		dedispersion_plan.set_nbits(nbits);
		dedispersion_plan.set_tsamp(tsamp);
		dedispersion_plan.set_tstart(tstart);
		dedispersion_plan.set_fch1(fch1);
		dedispersion_plan.set_foff(foff);
		// Initialise the GPU.
		int device_id = 0; // hard-coded, would be a parameter
		size_t gpu_memory = 0;
		cudaSetDevice(device_id);
		size_t mem_free, total;
		cudaMemGetInfo(&mem_free, &total);
		gpu_memory = ( mem_free/4 );
		// allocate memory cpu input
		io_data.allocate_memory_cpu_input(dedispersion_plan);
		// get recorded data
		io_data.get_recorded_data(&fp, dedispersion_plan.get_nchans(),dedispersion_plan.get_nbits());
		//
		// call sps_object main method

		//
		fclose(fp);
		*/
	}
}


} // namespace test
} // namespace sps
} // namespace astroaccelerate
} // namespace ska
