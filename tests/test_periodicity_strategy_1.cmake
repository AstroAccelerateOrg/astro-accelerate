set(TEST_NAME test_periodicity_strategy_1)
message(STATUS "INFO: Add ${TEST_NAME}")

add_executable(${TEST_NAME} "${CMAKE_CURRENT_LIST_DIR}/${TEST_NAME}")
include_directories(${PROJECT_BASE_DIR})
include_directories($ENV{CUDA_INSTALL_PATH}/include/)
include_directories($ENV{CUDA_INSTALL_PATH}/samples/common/inc/)

set_target_properties(${TEST_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests")
target_sources(${TEST_NAME}
               PUBLIC
               "include/aa_periodicity_strategy.hpp"
)
target_link_libraries(${TEST_NAME} PRIVATE ${CUDA_LIBRARIES} ${CUDA_CUFFT_LIBRARIES} ${CUDA_curand_LIBRARY} astroaccelerate)
add_test(NAME ${TEST_NAME} COMMAND tests/${TEST_NAME})
set_tests_properties(${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "Runs")