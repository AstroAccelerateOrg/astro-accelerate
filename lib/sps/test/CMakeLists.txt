include_directories(${GTEST_INCLUDE_DIR})

link_directories(${GTEST_LIBRARY_DIR})

set(
    gtest_Sps_src
    src/gtest_Sps.cpp
    src/SpsTest.cpp
)
add_executable(gtest_Sps ${gtest_Sps_src} )
target_link_libraries(gtest_Sps ${ASTRO_ACCELERATE_TEST_UTILS} ${ASTRO_ACCELERATE_LIBRARIES} ${GTEST_LIBRARIES})
add_test(gtest_Sps gtest_Sps --test_data "${CMAKE_CURRENT_LIST_DIR}/data")
