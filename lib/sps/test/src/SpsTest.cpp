#include "SpsTest.h"
#include "sps/Sps.h"
#include "../../SpsParameters.h"

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
    ASSERT_TRUE(dm_handler_called);
    ASSERT_TRUE(sps_handler_called);
}

} // namespace test
} // namespace sps
} // namespace astroaccelerate
} // namespace ska
