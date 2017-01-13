# thirdparty dependencies
include(cmake/googletest.cmake)

# other dependencies
set(Boost_ADDITIONAL_VERSIONS "1.50 1.50.0")
find_package(Boost COMPONENTS filesystem system program_options REQUIRED)

include_directories(SYSTEM ${BOOST_INCLUDE_DIR})
set(DEPENDENCY_LIBRARIES
#    ${Boost_LIBRARIES}
    )

include(compiler_settings)
include(cmake/cuda.cmake)
include(cmake/doxygen.cmake)
