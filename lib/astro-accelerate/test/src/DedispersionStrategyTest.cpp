#include "DedispersionStrategyTest.h"
#include "../../DedispersionStrategy.h"
#include <memory>
#include <exception>



namespace astroaccelerate {
namespace test {


DedispersionStrategyTest::DedispersionStrategyTest()
    : ::testing::Test()
{
}

DedispersionStrategyTest::~DedispersionStrategyTest()
{
}

void DedispersionStrategyTest::SetUp()
{
}

void DedispersionStrategyTest::TearDown()
{
}

TEST_F(DedispersionStrategyTest, test_channel_freq_order_consistency)
{
    // Use Case: frequencies can be provided in either high to low or low to high order
    // Expected Should be consistent values
    std::vector<float> user_dm_low({0});
    std::vector<float> user_dm_high({60});
    std::vector<float> user_dm_step({10});
    std::vector<int> in_bin({1});
    std::vector<int> out_bin({1});

    unsigned nsamples(1<<18);

    // add some frequencies
    std::vector<float> chan_frequencies_desc;
    std::vector<float> chan_frequencies_asc;
    float freq = 1000;
    float low_freq = freq-10*30;;
    for(unsigned i=0; i < 10; ++i )
    {
        chan_frequencies_desc.push_back(freq);
        chan_frequencies_asc.push_back(low_freq);
        freq -= 30;
        low_freq += 30;
    }

    astroaccelerate::DedispersionStrategy ds_1(
        user_dm_low.data()
        , user_dm_high.data()
        , user_dm_step.data()
        , in_bin.data()
        , out_bin.data()
        , 0.5e9 //gpu_memory
        , 8 //power
        , (int) user_dm_low.size()
        , nsamples 
        , 1 //nifs -- not used by DedispersionStrategy
        , 8 //nbits -- not used by DedispersionStrategy
        , 150e-3 //tsamp (assumed to be in seconds)
        , 55000.0000//tf_data.start_time().value() // is this used?
        , 6 //sigma_cutoff
        , 6 //sigma_constant
        , 0.5 //max_boxcar_width_in_sec
        , 0 //narrow -- not used
        , 0 //wide -- not used
        , 0 //nboots -- not used
        , 0 //navdms -- not used
        , 0 //ntrial_bins -- not used
        , 0 //nsearch -- not used
        , 0 //aggression -- not used
        , chan_frequencies_desc //chan frequency array
    );

    astroaccelerate::DedispersionStrategy ds_2(
        user_dm_low.data()
        , user_dm_high.data()
        , user_dm_step.data()
        , in_bin.data()
        , out_bin.data()
        , 0.5e9 //gpu_memory
        , 8 //power
        , (int) user_dm_low.size()
        , nsamples 
        , 1 //nifs -- not used by DedispersionStrategy
        , 8 //nbits -- not used by DedispersionStrategy
        , 150e-3 //tsamp (assumed to be in seconds)
        , 55000.0000//tf_data.start_time().value() // is this used?
        , 6 //sigma_cutoff
        , 6 //sigma_constant
        , 0.5 //max_boxcar_width_in_sec
        , 0 //narrow -- not used
        , 0 //wide -- not used
        , 0 //nboots -- not used
        , 0 //navdms -- not used
        , 0 //ntrial_bins -- not used
        , 0 //nsearch -- not used
        , 0 //aggression -- not used
        , chan_frequencies_asc //chan frequency array
    );

    ASSERT_EQ(ds_1.get_maxshift(), ds_2.get_maxshift());
    ASSERT_EQ(ds_1.get_dedispersed_time_samples(), ds_2.get_dedispersed_time_samples());
    // TODO check other params are thesame where appropriate
}

TEST_F(DedispersionStrategyTest, test_medium_number_of_samples_single_dm_range)
{
    std::vector<float> user_dm_low({0});
    std::vector<float> user_dm_high({60});
    std::vector<float> user_dm_step({10});
    std::vector<int> in_bin({1});
    std::vector<int> out_bin({1});

    // add some frequencies
    std::vector<float> chan_frequencies;
    float freq = 1000;
    for(unsigned i=0; i < 10; ++i )
    {
        chan_frequencies.push_back(freq);
        freq -= 30;
    }

    unsigned nsamples(1<<16);
    astroaccelerate::DedispersionStrategy ds(
        user_dm_low.data()
        , user_dm_high.data()
        , user_dm_step.data()
        , in_bin.data()
        , out_bin.data()
        , 0.5e9 //gpu_memory
        , 8 //power
        , (int) user_dm_low.size()
        , nsamples 
        , 1 //nifs -- not used by DedispersionStrategy
        , 8 //nbits -- not used by DedispersionStrategy
        , 150e-3 //tsamp (assumed to be in seconds)
        , 55000.0000//tf_data.start_time().value() // is this used?
        , 6 //sigma_cutoff
        , 6 //sigma_constant
        , 0.5 //max_boxcar_width_in_sec
        , 0 //narrow -- not used
        , 0 //wide -- not used
        , 0 //nboots -- not used
        , 0 //navdms -- not used
        , 0 //ntrial_bins -- not used
        , 0 //nsearch -- not used
        , 0 //aggression -- not used
        , chan_frequencies //chan frequency array
    );

    ASSERT_EQ(nsamples, ds.get_nsamp());
    
    // TODO verify the calculation of the number of variables we expect is sane
    ASSERT_NE(0, ds.get_dedispersed_time_samples());
    ASSERT_NE(0, ds.get_maxshift());
}

TEST_F(DedispersionStrategyTest, test_overlapping_dm_ranges)
{
    // Use case: multiple DM ranges are provided that overlap
    // Expected reaction: throw an error?
    std::vector<float> user_dm_low({0, 300, 600});
    std::vector<float> user_dm_high({340, 600, 1000});
    std::vector<float> user_dm_step({20, 20, 50});

    std::vector<int> in_bin({1,1,1});
    std::vector<int> out_bin({1,1,1});

    // add some frequencies
    std::vector<float> chan_frequencies;
    float freq = 1000;
    for(unsigned i=0; i < 10; ++i )
    {
        chan_frequencies.push_back(freq);
        freq -= 30;
    }

    unsigned nsamples(1<<16);

    std::unique_ptr<astroaccelerate::DedispersionStrategy> ds_ptr;
    ASSERT_THROW(ds_ptr.reset(new astroaccelerate::DedispersionStrategy(
        user_dm_low.data()
        , user_dm_high.data()
        , user_dm_step.data()
        , in_bin.data()
        , out_bin.data()
        , 0.5e9 //gpu_memory
        , 8 //power
        , (int) user_dm_low.size()
        , nsamples 
        , 1 //nifs -- not used by DedispersionStrategy
        , 8 //nbits -- not used by DedispersionStrategy
        , 150e-3 //tsamp (assumed to be in seconds)
        , 55000.0000//tf_data.start_time().value() // is this used?
        , 6 //sigma_cutoff
        , 6 //sigma_constant
        , 0.5 //max_boxcar_width_in_sec
        , 0 //narrow -- not used
        , 0 //wide -- not used
        , 0 //nboots -- not used
        , 0 //navdms -- not used
        , 0 //ntrial_bins -- not used
        , 0 //nsearch -- not used
        , 0 //aggression -- not used
        , chan_frequencies //chan frequency array
    )), std::runtime_error);

}

} // namespace test
} // namespace astroaccelerate
