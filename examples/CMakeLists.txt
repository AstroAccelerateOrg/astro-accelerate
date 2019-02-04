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
