#ifndef SKA_ASTROACCELERATE_SPS_TEST_SPSTEST_H
#define SKA_ASTROACCELERATE_SPS_TEST_SPSTEST_H

#include <gtest/gtest.h>

// global variables are bad. avoid them when possible
// used here because: http://stackoverflow.com/questions/4818785/how-to-pass-parameters-to-the-gtest
// see accepted answer
extern int my_argc;
extern char** my_argv;

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
