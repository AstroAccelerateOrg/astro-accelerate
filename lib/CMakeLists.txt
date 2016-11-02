include(subpackage)

# define the astro_accelerate libraries
set(ASTRO_ACCELERATE_LIBRARIES ${CMAKE_PROJECT_NAME} ${DEPENDENCY_LIBRARIES})

# define the astro_accelerate test_utils libraries
set(ASTRO_ACCELERATE_TEST_UTILS ${CMAKE_PROJECT_NAME}_test_utils)
#list(APPEND ASTRO_ACCELERATE_TEST_UTILS ${PANDA_TEST_LIBRARIES})

# ensure we can access our headers with #include "AstroAccelerate/...."
include_directories(.)
include_directories(${PROJECT_BINARY_DIR})

# subpackages to define ${lib_src} & ${test_utils_src} lists
# in dependency order
include(Backward.cmake)
#SUBPACKAGE(utils) # not dependent on any other subpackage

# should come after all SUBPACKGE directives
include_subpackage_files()

# -- the main cuda library
if(ENABLE_CUDA)
  cuda_compile(cuda_objects ${lib_src_cuda}) 
  list(APPEND cuda_objects ${lib_obj_cuda})
  cuda_generate_link_file(cuda_link_file ${CMAKE_PROJECT_NAME} ${cuda_objects})
  list(APPEND cuda_objects ${cuda_link_file})
endif(ENABLE_CUDA)

# -- the main library target
add_library(${CMAKE_PROJECT_NAME} ${lib_src_cpu} ${cuda_objects})
set_target_properties(${CMAKE_PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)

if(ENABLE_CUDA)
  cuda_add_cufft_to_target(${CMAKE_PROJECT_NAME})
  target_link_libraries(${CMAKE_PROJECT_NAME}
        ${CUDA_LIBRARIES}
  )
  # We need to set the linker language based on what the expected generated file
  # would be. CUDA_C_OR_CXX is computed based on CUDA_HOST_COMPILATION_CPP.
  set_target_properties(${cuda_target}
          PROPERTIES
          LINKER_LANGUAGE ${CUDA_C_OR_CXX}
          )
endif(ENABLE_CUDA)

# test_utils library target
#add_library(${CMAKE_PROJECT_NAME}_test_utils ${test_utils_src_cpu})
#target_include_directories(${CMAKE_PROJECT_NAME}_test_utils PUBLIC ${GTEST_INCLUDE_DIR})

