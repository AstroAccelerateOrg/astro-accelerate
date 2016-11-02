# == googletest build - must be build with the same compiler flags
if(APPLE)
add_definitions(-DGTEST_USE_OWN_TR1_TUPLE=1)
else(APPLE)
add_definitions(-DGTEST_USE_OWN_TR1_TUPLE=0)
endif(APPLE)

add_subdirectory("thirdparty/googletest")
set(GTEST_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/thirdparty/googletest/include)
set(GTEST_LIBRARY_DIR ${CMAKE_BINARY_DIR}/thirdparty/googletest)
set(GTEST_LIBRARIES gtest)
