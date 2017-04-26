set(module_astro-accelerate_lib_src_cpu
    src/DedispersionStrategy.cpp
    src/DedispersionStrategyFile.cpp
    PARENT_SCOPE
)

#TEST_UTILS()

add_subdirectory(test)
