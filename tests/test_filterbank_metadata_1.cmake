set(TEST_NAME test_filterbank_metadata_1)
message(STATUS "INFO: Add ${TEST_NAME}")
add_executable(${TEST_NAME} "${CMAKE_CURRENT_LIST_DIR}/${TEST_NAME}")
set_target_properties(${TEST_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests")
target_sources(${TEST_NAME}
               PRIVATE
               PUBLIC
               "include/aa_filterbank_metadata.hpp"
)

add_test(NAME ${TEST_NAME} COMMAND tests/${TEST_NAME})
set_tests_properties(${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "test_filterbank_metadata_1 pass")
