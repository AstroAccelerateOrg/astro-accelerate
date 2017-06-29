#ifndef ASTROACCELERATE_ASTROACCELERATE_TEST_DEDISPERSIONSTRATEGYTEST_H
#define ASTROACCELERATE_ASTROACCELERATE_TEST_DEDISPERSIONSTRATEGYTEST_H

#include <gtest/gtest.h>

namespace astroaccelerate {
namespace test {

/**
 * @brief
 *    Unit tests for the desiseprsin strategy
 * @details
 * 
 */

class DedispersionStrategyTest : public ::testing::Test
{
    protected:
        void SetUp() override;
        void TearDown() override;

    public:
        DedispersionStrategyTest();
        ~DedispersionStrategyTest();

    private:
};


} // namespace test
} // namespace astroaccelerate

#endif // ASTROACCELERATE_ASTROACCELERATE_TEST_DEDISPERSIONSTRATEGYTEST_H
