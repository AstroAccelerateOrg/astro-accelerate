#include <gtest/gtest.h>

int main(int argc, char* argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
    //return ska::panda::test::panda_gtest_main(argc, argv);
}
