#ifndef SKA_ASTROACCELERATE_SPS_TEST_SPSTEST_H
#define SKA_ASTROACCELERATE_SPS_TEST_SPSTEST_H

#include <gtest/gtest.h>

namespace ska {
namespace astroaccelerate {
namespace sps {
namespace test {

/**
 * @brief
 * 
 * @details
 * 
 */

class SpsTest : public ::testing::Test
{
    protected:
        void SetUp() override;
        void TearDown() override;

    public:
        SpsTest();

        ~SpsTest();

    private:
};


} // namespace test
} // namespace sps
} // namespace astroaccelerate
} // namespace ska

#endif // SKA_ASTROACCELERATE_SPS_TEST_SPSTEST_H
