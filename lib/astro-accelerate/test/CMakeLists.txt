include_directories(${GTEST_INCLUDE_DIR})

link_directories(${GTEST_LIBRARY_DIR})

set(
    gtest_AstroAccelerate_src
    src/gtest_AstroAccelerate.cpp
    src/AstroAccelerateTest.cpp
)
add_executable(gtest_AstroAccelerate ${gtest_AstroAccelerate_src} )
target_link_libraries(gtest_AstroAccelerate ${ASTRO_ACCELERATE_TEST_UTILS} ${ASTRO_ACCELERATE_LIBRARIES} ${GTEST_LIBRARIES})
add_test(gtest_AstroAccelerate gtest_AstroAccelerate --test_data "${CMAKE_CURRENT_LIST_DIR}/data")
