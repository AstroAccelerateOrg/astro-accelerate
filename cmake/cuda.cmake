option(ENABLE_CUDA "Enable CUDA algorithms. Requires CUDA toolkit to be installed" ON)
if(ENABLE_CUDA)

  cmake_minimum_required(VERSION 2.8)
  find_package(CUDA REQUIRED)
  include_directories(${CUDA_TOOLKIT_INCLUDE})

  # add sdk samples useful headerfiles like cuda_helpers.h
  if(CUDA_SMP_INC)
	include_directories(${CUDA_SMP_INC})
  endif(CUDA_SMP_INC)

  set(CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
  set(CUDA_PROPAGATE_HOST_FLAGS OFF)

  # Pass options to NVCC ( -ccbin /path  --compiler-options -lfftw3f --compiler-options -lm --verbose)
  list(APPEND CUDA_NVCC_FLAGS -DENABLE_CUDA -g -std=c++11 -lineinfo -Xcompiler -fopenmp -O3 -lm -arch=sm_61 --use_fast_math --ptxas-options=-v -Xptxas -dlcm=cg)
  list(APPEND CUDA_NVCC_FLAGS_DEBUG --debug; --device-debug; --generate-line-info -Xcompiler)
  #list(APPEND CUDA_NVCC_FLAGS_DEBUG --debug; --device-debug; --generate-line-info -Xcompiler "-Werror")

  #list(APPEND CUDA_NVCC_FLAGS -arch compute_35) # minumum compute level (Sps restriction)
  #list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_52,code=sm_52) # TitanX
  #list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_50,code=sm_50) # Maxwell
  #list(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_37,code=sm_37) # K80

  set(CMAKE_CXX_FLAGS "-DENABLE_CUDA ${CMAKE_CXX_FLAGS}")

  get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
  foreach(dir ${dirs})
    list(APPEND CUDA_NVCC_INCLUDE_ARGS_USER ${CMAKE_INCLUDE_SYSTEM_FLAG_CXX} ${dir})
  endforeach()

  macro(CUDA_GENERATE_LINK_FILE cuda_link_file cuda_target)

      # Compute the file name of the intermedate link file used for separable
      # compilation.
      CUDA_COMPUTE_SEPARABLE_COMPILATION_OBJECT_FILE_NAME(cuda_link_file ${cuda_target} "${${cuda_target}_SEPARABLE_COMPILATION_OBJECTS}")

      # Add a link phase for the separable compilation if it has been enabled.  If
      # it has been enabled then the ${cuda_target}_SEPARABLE_COMPILATION_OBJECTS
      # variable will have been defined.
      CUDA_LINK_SEPARABLE_COMPILATION_OBJECTS("${cuda_link_file}" ${cuda_target} "${_options}" "${${cuda_target}_SEPARABLE_COMPILATION_OBJECTS}")

  endmacro()

  #
  # @macro CUDA_SUBPACKAGE_COMPILE
  # Use to specify that certain objects should be compiled with 
  # alternative flags (e.g. a specific architecture)
  # @example
  # cuda_subpackage_compile(${my_cuda_files} OPTIONS "-arch compute_35")
  #
  macro(CUDA_SUBPACKAGE_COMPILE)
    FILE(APPEND ${SUBPACKAGE_FILENAME}
      "CUDA_ADD_CUDA_INCLUDE_ONCE()\n"
      "cuda_compile(_cuda_objects "
    )
    foreach(arg ${ARGN})
        if(EXISTS "${arg}")
            set(_file "${arg}")
        else()
            set(_file "${CMAKE_CURRENT_SOURCE_DIR}/${arg}")
            if(NOT EXISTS "${_file}")
                set(_file ${arg})
            endif()
        endif()
        FILE(APPEND ${SUBPACKAGE_FILENAME}
            "${_file}\n"
        )
    endforeach(arg ${ARGV})
    FILE(APPEND ${SUBPACKAGE_FILENAME}
      ")\n"
      "list(APPEND lib_obj_cuda \${_cuda_objects})\n"
    )
  endmacro(CUDA_SUBPACKAGE_COMPILE)

else(ENABLE_CUDA)
  # dummies
  macro(CUDA_SUBPACKAGE_COMPILE)
  endmacro(CUDA_SUBPACKAGE_COMPILE)
endif(ENABLE_CUDA)
