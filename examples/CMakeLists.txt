# Examples

# Include directory
include_directories(../include/)

# Example dedispersion
add_executable(examples_dedispersion src/dedispersion.cpp)
target_link_libraries(examples_dedispersion astroaccelerate)

# Example dedispersion and analysis
add_executable(examples_dedispersion_and_analysis src/dedispersion_and_analysis.cpp)
target_link_libraries(examples_dedispersion_and_analysis astroaccelerate)

# Example periodicity
add_executable(examples_periodicity src/periodicity.cpp)
target_link_libraries(examples_periodicity astroaccelerate)

# Examples filterbank processing
add_executable(examples_filterbank_dedispersion src/filterbank_dedispersion.cpp)
target_link_libraries(examples_filterbank_dedispersion astroaccelerate)

# Examples fake signal single
add_executable(examples_fake_single_signal src/fake_signal_single.cpp)
target_link_libraries(examples_fake_single_signal astroaccelerate)

# Examples fake signal periodic
add_executable(examples_fake_periodic_signal src/fake_signal_periodic.cpp)
target_link_libraries(examples_fake_periodic_signal astroaccelerate)
