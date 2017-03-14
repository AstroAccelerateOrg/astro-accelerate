include_directories(${GTEST_INCLUDE_DIR})

link_directories(${GTEST_LIBRARY_DIR})

set(
   gtest_utils_src
    src/gtest_utils.cpp
)
add_executable(gtest_utils ${gtest_utils_src} )
target_link_libraries(gtest_utils ${ASTRO_ACCELERATE_TEST_UTILS} ${ASTRO_ACCELERATE_LIBRARIES} ${GTEST_LIBRARIES})
add_test(gtest_utils gtest_utils --test_data "${CMAKE_CURRENT_LIST_DIR}/data")
