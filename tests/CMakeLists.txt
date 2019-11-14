# CMakeLists.txt for tests
include(CTest)

include(tests/test_filterbank_metadata_1.cmake)
include(tests/test_ddtr_strategy_1.cmake)
include(tests/test_ddtr_strategy_2.cmake)
include(tests/test_periodicity_strategy_1.cmake)
include(tests/test_device_memory_request.cmake)
include(tests/test_single_pulse_search_scan_1.cmake)
include(tests/test_ddtr_fakesignal_single.cmake)
include(tests/test_ddtr_fakesignal_period.cmake)

