#ifndef ASTROACCELERATE_ASTROACCELERATE_TEST_DMTIMETEST_H
#define ASTROACCELERATE_ASTROACCELERATE_TEST_DMTIMETEST_H

#include <gtest/gtest.h>

namespace astroaccelerate {
namespace test {

/**
 * @brief
 * 
 * @details
 * 
 */

class DmTimeTest : public ::testing::Test
{
    protected:
        void SetUp() override;
        void TearDown() override;

    public:
        DmTimeTest();

        ~DmTimeTest();

    private:
};


} // namespace test
} // namespace astroaccelerate

#endif // ASTROACCELERATE_ASTROACCELERATE_TEST_DMTIMETEST_H
