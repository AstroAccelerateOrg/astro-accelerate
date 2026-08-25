set(TEST_NAME test_dedispersion_index)
message(STATUS "INFO: Add ${TEST_NAME}")

add_executable(${TEST_NAME} "${CMAKE_CURRENT_LIST_DIR}/${TEST_NAME}.cpp")
target_include_directories(${TEST_NAME} PRIVATE "${PROJECT_BASE_DIR}/include")
set_target_properties(${TEST_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests")

add_test(NAME ${TEST_NAME} COMMAND tests/${TEST_NAME})
set_tests_properties(${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "Runs")
