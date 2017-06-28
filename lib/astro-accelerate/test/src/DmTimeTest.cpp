#include "DmTimeTest.h"
#include "../../DmTime.h"


namespace astroaccelerate {
namespace test {


DmTimeTest::DmTimeTest()
    : ::testing::Test()
{
}

DmTimeTest::~DmTimeTest()
{
}

void DmTimeTest::SetUp()
{
}

void DmTimeTest::TearDown()
{
}

TEST_F(DmTimeTest, test_dmtime)
{
    std::vector<float> user_dm_low({10, 300, 600});
    std::vector<float> user_dm_high({100, 580, 1000});
    std::vector<float> user_dm_step({20, 40, 50});
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
    DmTime<float> dm_time(ds);
    ASSERT_EQ(ds.get_maxshift(), nsamples - dm_time.max_nsamples());
}

} // namespace test
} // namespace astroaccelerate
