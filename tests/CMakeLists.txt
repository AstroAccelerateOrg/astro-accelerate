# CMakeLists.txt for tests
include(CTest)

add_executable(test_ddtr_strategy_1)
set_target_properties(test_ddtr_strategy_1 PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests")
target_sources(test_ddtr_strategy_1
	       PRIVATE
	       "src/aa_ddtr_strategy.cpp"
	       "${CMAKE_CURRENT_LIST_DIR}/test_ddtr_strategy_1.cpp"
	       PUBLIC
	       "include/aa_ddtr_strategy.hpp"
)

add_test(test_ddtr_strategy_1 tests/test_ddtr_strategy_1 COMMAND 1)
set_tests_properties(test_ddtr_strategy_1 PROPERTIES PASS_REGULAR_EXPRESSION "Runs")
