# Check if cmake has the required version
CMAKE_MINIMUM_REQUIRED(VERSION 2.6.0 FATAL_ERROR)
set(PROJECT_NAME_STR astro-accelerate)
set(PROJECT_LIB_NAME ASTRO_ACCELERATE)
list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_SOURCE_DIR}/cmake)
find_package(CUDA REQUIRED)

# Default installation locations (note these will be prefixed by ${CMAKE_INSTALL_PREFIX})
set(include_dir_name "astroaccelerate")
if(NOT LIBRARY_INSTALL_DIR)
    set(LIBRARY_INSTALL_DIR "lib")
endif(NOT LIBRARY_INSTALL_DIR)

if(NOT INCLUDE_INSTALL_DIR)
    set(INCLUDE_INSTALL_DIR "include/${include_dir_name}")
endif(NOT INCLUDE_INSTALL_DIR)

if(NOT MODULES_INSTALL_DIR)
    set (MODULES_INSTALL_DIR "share/${include_dir_name}")
endif(NOT MODULES_INSTALL_DIR)

if(NOT DOC_INSTALL_DIR)
    set (DOC_INSTALL_DIR ${MODULES_INSTALL_DIR}/doc)
endif(NOT DOC_INSTALL_DIR)

if(NOT BINARY_INSTALL_DIR)
    set(BINARY_INSTALL_DIR "bin/")
endif(NOT BINARY_INSTALL_DIR)

# Project Version
include(git_version)
configure_file(${PROJECT_SOURCE_DIR}/cmake/version.h.in
    ${PROJECT_BINARY_DIR}/${include_dir_name}/version.h)
install(FILES ${PROJECT_BINARY_DIR}/${include_dir_name}/version.h
    DESTINATION ${INCLUDE_INSTALL_DIR})


# Set NVCC flags
set(PROJECT_BASE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(CUDA_LINK_LIBRARIES_KEYWORD PUBLIC)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS})
set(CUDA_PROPAGATE_HOST_FLAGS OFF)

list(APPEND CUDA_NVCC_FLAGS --ptxas-options= -Xptxas -dlcm=cg;)
list(APPEND CUDA_NVCC_FLAGS "-gencode arch=compute_50,code=sm_50")
list(APPEND CUDA_NVCC_FLAGS -g;)
list(APPEND CUDA_NVCC_FLAGS -lineinfo;)
list(APPEND CUDA_NVCC_FLAGS -Xcompiler;-O3;)
list(APPEND CUDA_NVCC_FLAGS -Xcompiler;"-fopenmp";)
list(APPEND CMAKE_CXX_FLAGS "-std=c++11 -O3 -fopenmp -lm -Wall -ffast-math")

# Status information
message(STATUS "Using: ${CMAKE_CXX_COMPILER} with compiler ID ${CMAKE_CXX_COMPILER_ID} and compiler version ${CXX_VERSION}")
message(STATUS "Using CXX compilation flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "Using CUDA NVCC flags ${CUDA_NVCC_FLAGS}")

# Include and linker directories
link_directories(${CUDA_LIBRARY_DIRS})
include_directories(include/)
include_directories(${PROJECT_BASE_DIR})
include_directories($ENV{CUDA_INSTALL_PATH}/include/)
include_directories($ENV{CUDA_INSTALL_PATH}/samples/common/inc/)
include_directories(${CUDA_LIBRARY_DIRS})

# CUDA library object
file(GLOB_RECURSE GPU_SOURCE_FILES "src/*.cu")
cuda_add_library(astroaccelerate SHARED ${GPU_SOURCE_FILES})
target_link_libraries(astroaccelerate PRIVATE ${CUDA_LIBRARY_DIRS} ${CUDA_CUFFT_LIBRARIES} ${CUDA_curand_LIBRARY})
set_target_properties(astroaccelerate PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties(astroaccelerate PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(astroaccelerate PROPERTIES BUILD_SHARED_LIBS ON)

# Standalone executable to link against CUDA library
file(GLOB_RECURSE SOURCE_FILES "src/*.cpp")
add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} astroaccelerate)
