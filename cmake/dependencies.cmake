# thirdparty dependencies
include(cmake/googletest.cmake)

# other dependencies
find_package(Boost COMPONENTS filesystem system program_options REQUIRED)

include_directories(SYSTEM ${BOOST_INCLUDE_DIR})
set(DEPENDENCY_LIBRARIES
    ${Boost_LIBRARIES}
    )

include(compiler_settings)
include(cmake/cuda.cmake)
include(cmake/doxygen.cmake)
