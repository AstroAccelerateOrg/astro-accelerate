# Check if cmake has the required version
CMAKE_MINIMUM_REQUIRED(VERSION 3.2.0 FATAL_ERROR)
set(PROJECT_NAME_STR ASTRO_ACCELERATE)
set(PROJECT_LIB_NAME ASTRO_ACCELERATE)
PROJECT(${PROJECT_NAME_STR} LANGUAGES CXX CUDA)
find_package(CUDA REQUIRED)

# Set NVCC flags
set(PROJECT_BASE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(CUDA_LINK_LIBRARIES_KEYWORD PUBLIC)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS})

list(APPEND CUDA_NVCC_FLAGS --ptxas-options= -Xptxas -dlcm=cg;)
list(APPEND CUDA_NVCC_FLAGS -arch=sm_50;)
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
target_link_libraries(astroaccelerate ${CUDA_LIBRARY_DIRS})
set_target_properties(astroaccelerate PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties(astroaccelerate PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(astroaccelerate PROPERTIES BUILD_SHARED_LIBS ON)

# Standalone executable to link against CUDA library
file(GLOB_RECURSE SOURCE_FILES "src/*.cpp")
add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} astroaccelerate)
