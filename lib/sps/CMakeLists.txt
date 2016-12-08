set(module_sps_lib_src_cpu
    src/DedispersionPlan.cpp
    src/IOData.cpp
    src/UserInput.cpp
    PARENT_SCOPE
)

#TEST_UTILS()

add_subdirectory(test)
