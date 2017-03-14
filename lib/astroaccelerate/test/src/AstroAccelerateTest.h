#ifndef ASTROACCELERATE__TEST_ASTROACCELERATETEST_H
#define ASTROACCELERATE__TEST_ASTROACCELERATETEST_H

#include <gtest/gtest.h>

// global variables are bad. avoid them when possible
// used here because: http://stackoverflow.com/questions/4818785/how-to-pass-parameters-to-the-gtest
// see accepted answer
extern int my_argc;
extern char** my_argv;

namespace ska {
namespace astroaccelerate {
namespace test {

/**
 * @brief
 * 
 * @details
 * 
 */

class AstroAccelerateTest : public ::testing::Test
{
    protected:
        void SetUp() override;
        void TearDown() override;

    public:
        AstroAccelerateTest();

        ~AstroAccelerateTest();

    private:
};


} // namespace test
} // namespace astroaccelerate
}
#endif // ASTROACCELERATE__TEST_ASTROACCELERATETEST_H
