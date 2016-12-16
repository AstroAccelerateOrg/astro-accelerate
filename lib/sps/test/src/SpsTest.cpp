#include "SpsTest.h"
#include "sps/Sps.h"
#include "../../SpsParameters.h"
#include "../../UserInput.h"

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

TEST_F(SpsTest, test_handlers)
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


TEST_F(SpsTest, test_user_input)
{
	// Following is ok for ska_karel.txt
	char* filename = my_argv[1] + strlen(my_argv[1]) - 13;
	if(strcmp(filename, "ska_karel.txt") == 0)
	{
		// create objects
		sps::UserInput user_input;
		sps::DedispersionPlan dedispersion_plan;
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
		user_input.get_user_input(&fp, my_argc, my_argv, dedispersion_plan);

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
		fclose(fp);
	}


}

TEST_F(SpsTest, test_dedispersion_plan)
{
	// declare objects
	sps::UserInput user_input;
	sps::DedispersionPlan dedispersion_plan;
	// first, check constructor
	EXPECT_EQ(NULL, dedispersion_plan.get_in_bin());
	EXPECT_EQ(NULL, dedispersion_plan.get_out_bin());
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
	EXPECT_FLOAT_EQ(2.0f, dedispersion_plan.get_power());
	// get user input
	FILE *fp = NULL;
	user_input.get_user_input(&fp, my_argc, my_argv, dedispersion_plan);
	//

}

} // namespace test
} // namespace sps
} // namespace astroaccelerate
} // namespace ska
