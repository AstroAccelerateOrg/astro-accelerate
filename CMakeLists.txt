# Check if cmake has the required version
CMAKE_MINIMUM_REQUIRED(VERSION 3.9.0 FATAL_ERROR)

project(astro-accelerate)
set(PROJECT_NAME astro-accelerate)
set(PROJECT_VERSION 1.3.1)	# Tagged release version
set(PROJECT_LIB_NAME ASTRO_ACCELERATE)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(JoinPaths)

if(POLICY CMP0076)
    cmake_policy(SET CMP0076 NEW)
endif()

if(POLICY CMP0074)
    cmake_policy(SET CMP0074 NEW)
endif()

if(POLICY CMP0126)
    cmake_policy(SET CMP0126 NEW)
endif()


# Set date and time of build
string(TIMESTAMP ASTRO_ACCELERATE_BUILD_DATE_UTC "%d/%m/%Y" UTC)
string(TIMESTAMP ASTRO_ACCELERATE_BUILD_TIME_UTC "%H:%M:%S" UTC)
if(NOT ASTRO_ACCELERATE_BUILD_DATE_UTC)
set(ASTRO_ACCELERATE_BUILD_DATE_UTC "EMPTY")
endif()

# Setup Doxygen
cmake_policy(SET CMP0012 NEW)
option(ENABLE_DOCS "Enable docs" OFF)
if(${ENABLE_DOCS})
	message("-- INFO: DOCS enabled.")
	set(doc_all_target "ALL")
	if(NOT DOC_INSTALL_DIR)
	       set (DOC_INSTALL_DIR ${MODULES_INSTALL_DIR}/doc)
	       message("-- INFO: DOC_INSTALL_DIR ${DOC_INSTALL_DIR}")
	endif(NOT DOC_INSTALL_DIR)
	include(cmake/doxygen.cmake)
endif()

# Find Git and Git info
find_package(Git)
if(GIT_FOUND)
  message("-- INFO: Git found: ${GIT_EXECUTABLE}")
  execute_process(COMMAND git log -1 --format=%H
  	          WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
		  OUTPUT_VARIABLE GIT_COMMIT_HASH
		  OUTPUT_STRIP_TRAILING_WHITESPACE)
		  
endif()
add_definitions("-DGIT_COMMIT_HASH=${GIT_COMMIT_HASH}")

if(NOT ASTRO_ACCELERATE_BUILD_TIME_UTC)
set(ASTRO_ACCELERATE_BUILD_TIME_UTC "EMPTY")
endif()

find_package(CUDA REQUIRED)

# Set NVCC flags
set(PROJECT_BASE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(CUDA_LINK_LIBRARIES_KEYWORD PUBLIC)
set(CUDA_PROPAGATE_HOST_FLAGS OFF)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS})

list(APPEND CUDA_NVCC_FLAGS --use_fast_math)
list(APPEND CUDA_NVCC_FLAGS -g;)
list(APPEND CUDA_NVCC_FLAGS -Xptxas -O3 -std=c++11;)
list(APPEND CUDA_NVCC_FLAGS -lineinfo;)
list(APPEND CUDA_NVCC_FLAGS -Xcompiler;-O3;)
list(APPEND CUDA_NVCC_FLAGS -Xcompiler;-fopenmp;)
list(APPEND CUDA_NVCC_FLAGS -Xcompiler;-lm;)
list(APPEND CUDA_NVCC_FLAGS -Xcompiler;-Wall;)
list(APPEND CMAKE_CXX_FLAGS "-std=c++11 -O3 -lm -Wall -Wpedantic -Wextra -fopenmp")

if(NOT DEFINED CUDA_ARCH)  
	set(CUDA_ARCH "ALL")
	message("-- INFO: Setting CUDA_ARCH to ALL.")
	message("-- INFO: The target CUDA architecture can be specified using:")
	message("-- INFO:   -DCUDA_ARCH=\"<arch>\"")
	message("-- INFO: where <arch> is one or more of:")
	message("-- INFO:   3.5, 3.7, 5.0, 5.2, 6.0, 6.1, 6.2, 7.0, 7.5, 8.6, 9.0 or ALL.")
	message("-- INFO: Separate multiple architectures with semicolons.")
endif()

foreach(ARCH ${CUDA_ARCH})
	if(ARCH MATCHES ALL|[Aa]ll)
		message("-- INFO: Building CUDA device code for Kepler,")
		message("-- INFO: Maxwell and Pascal architectures")
	   	list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_35,code=sm_35)
		list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_37,code=sm_37)
	   	list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_50,code=sm_50)
	   	list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_52,code=sm_52)
	   	list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_60,code=sm_60)
	   	list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_61,code=sm_61)
	   	list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_62,code=sm_62)
		list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_70,code=sm_70)
		list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_75,code=sm_75)
		list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_80,code=sm_80)
		list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_86,code=sm_86)
		list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_90,code=sm_90)
		set(ASTRO_ACCELERATE_CUDA_ARCH_VERSION "35,37,50,52,60,61,62,70,75,80,86,90")
		set(ASTRO_ACCELERATE_CUDA_SM_VERSION "35,37,50,52,60,61,62,70,75,80,86,90")
	elseif(ARCH MATCHES 3.5)
		message("-- INFO: Building CUDA device code for architecture 3.5")
		list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_35,code=sm_35)
		set(ASTRO_ACCELERATE_CUDA_ARCH_VERSION "35")
		set(ASTRO_ACCELERATE_CUDA_SM_VERSION "35")
	elseif(ARCH MATCHES 3.7)
		message("-- INFO: Building CUDA device code for architecture 3.7")
		list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_37,code=sm_37)
		set(ASTRO_ACCELERATE_CUDA_ARCH_VERSION "37")
		set(ASTRO_ACCELERATE_CUDA_SM_VERSION "37")
	elseif(ARCH MATCHES 5.0)
		message("-- INFO: Building CUDA device code for architecture 5.0")
		list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_50,code=sm_50)
		set(ASTRO_ACCELERATE_CUDA_ARCH_VERSION "50")
		set(ASTRO_ACCELERATE_CUDA_SM_VERSION "50")
	elseif(ARCH MATCHES 5.2)
		message("-- INFO: Building CUDA device code for architecture 5.2")
		list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_52,code=sm_52)
		set(ASTRO_ACCELERATE_CUDA_ARCH_VERSION "52")
		set(ASTRO_ACCELERATE_CUDA_SM_VERSION "52")
	elseif(ARCH MATCHES 6.0)
		message("-- INFO: Building CUDA device code for architecture 6.0")
		list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_60,code=sm_60)
		set(ASTRO_ACCELERATE_CUDA_ARCH_VERSION "60")
		set(ASTRO_ACCELERATE_CUDA_SM_VERSION "60")
	elseif(ARCH MATCHES 6.1)
		message("-- INFO: Building CUDA device code for architecture 6.1")
		list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_61,code=sm_61)
		set(ASTRO_ACCELERATE_CUDA_ARCH_VERSION "61")
		set(ASTRO_ACCELERATE_CUDA_SM_VERSION "61")
	elseif(ARCH MATCHES 6.2)
		message("-- INFO: Building CUDA device code for architecture 6.2")
		list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_62,code=sm_62)
		set(ASTRO_ACCELERATE_CUDA_ARCH_VERSION "62")
		set(ASTRO_ACCELERATE_CUDA_SM_VERSION "62")
	elseif(ARCH MATCHES 7.0)
		message("-- INFO: Building CUDA device code for architecture 7.0")
		list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_70,code=sm_70)
		set(ASTRO_ACCELERATE_CUDA_ARCH_VERSION "70")
		set(ASTRO_ACCELERATE_CUDA_SM_VERSION "70")
	elseif(ARCH MATCHES 7.5)
		message("-- INFO: Building CUDA device code for architecture 7.5")
		list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_75,code=sm_75)
		set(ASTRO_ACCELERATE_CUDA_ARCH_VERSION "75")
		set(ASTRO_ACCELERATE_CUDA_SM_VERSION "75")
	elseif(ARCH MATCHES 8.0)
		message("-- INFO: Building CUDA device code for architecture 8.0")
		list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_80,code=sm_80)
		set(ASTRO_ACCELERATE_CUDA_ARCH_VERSION "80")
		set(ASTRO_ACCELERATE_CUDA_SM_VERSION "80")
	elseif(ARCH MATCHES 8.6)
		message("-- INFO: Building CUDA device code for architecture 8.6")
		list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_86,code=sm_86)
		set(ASTRO_ACCELERATE_CUDA_ARCH_VERSION "86")
		set(ASTRO_ACCELERATE_CUDA_SM_VERSION "86")
	elseif(ARCH MATCHES 9.0)
		message("-- INFO: Building CUDA device code for architecture 9.0")
		list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_90,code=sm_90)
		set(ASTRO_ACCELERATE_CUDA_ARCH_VERSION "90")
		set(ASTRO_ACCELERATE_CUDA_SM_VERSION "90")
	else()
	message(FATAL_ERROR "-- CUDA_ARCH ${ARCH} not recognised or not defined")
	endif()
endforeach()

option(LOGGING_FACILITY "Enable logging facility" ON)
if(${LOGGING_FACILITY})
        message("-- INFO: LOGGING_FACILITY ON")
	add_definitions(-DASTRO_ACCELERATE_LOGGING_FACILITY_ENABLE=1)
else()
	message("-- INFO: LOGGING_FACILITY OFF")
	add_definitions(-DASTRO_ACCELERATE_LOGGING_FACILITY_ENABLE=0)
endif()

# Configure file
list(APPEND CMAKE_PREFIX_PATH "${PROJECT_SOURCE_DIR}/cmake")
configure_file("${PROJECT_SOURCE_DIR}/cmake/version.h.in" "${PROJECT_SOURCE_DIR}/include/aa_version.hpp")
add_definitions("-DASTRO_ACCELERATE_VERSION_H_DEFINED=1")

# Python configure file
configure_file("${PROJECT_SOURCE_DIR}/cmake/py_astro_accelerate.py.in" "${PROJECT_SOURCE_DIR}/python/py_astro_accelerate.py")


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
file(GLOB_RECURSE GPU_SOURCE_FILES "src/*.cu" "src/*.cpp")
cuda_add_library(astroaccelerate SHARED ${GPU_SOURCE_FILES})
target_link_libraries(astroaccelerate PRIVATE ${CUDA_LIBRARIES} ${CUDA_CUFFT_LIBRARIES} ${CUDA_curand_LIBRARY})
set_target_properties(astroaccelerate PROPERTIES CUDA_SEPARABLE_COMPILATION OFF)
set_target_properties(astroaccelerate PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(astroaccelerate PROPERTIES BUILD_SHARED_LIBS ON)


# Standalone executable to link against CUDA library
file(GLOB_RECURSE SOURCE_FILES "src/aa_main.cpp")
add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} PRIVATE astroaccelerate)

install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME}Config
	ARCHIVE DESTINATION lib
	LIBRARY DESTINATION lib
	RUNTIME DESTINATION bin
)

# copy header files
install(DIRECTORY include/ DESTINATION include/aa
FILES_MATCHING PATTERN "*.hpp")

# Create/export CMAKE config file for dependency
install(EXPORT ${PROJECT_NAME}Config DESTINATION cmake)
export(TARGETS ${PROJECT_NAME} FILE ${PROJECT_NAME}Config.cmake)

# Create pkg_config
join_paths(CMAKE_INSTALL_LIBDIR "lib/pkgconfig" "${prefix}")

set(DEST_DIR "${CMAKE_INSTALL_LIBDIR}")
configure_file(
  "${PROJECT_SOURCE_DIR}/astroaccelerate.pc.in"
  "${PROJECT_BINARY_DIR}/astroaccelerate.pc"
  @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/astroaccelerate.pc" DESTINATION "${CMAKE_INSTALL_LIBDIR}")

# Add tests/ directory
option(ENABLE_TESTS "Enable tests" OFF)
if(${ENABLE_TESTS})
	message(STATUS "INFO: ENABLE_TESTS ON")
	enable_testing()
	include(tests/CMakeLists.txt)
endif()

add_subdirectory(examples)
