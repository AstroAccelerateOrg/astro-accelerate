#
# Wrapper to FindPackageHandleStandardArgs required for message printing
# in find_package macros to fix compatibility with cmake before 2.5
#
if ("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER 2.5)
    include(FindPackageHandleStandardArgs)
else ("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER 2.5)
    MACRO(FIND_PACKAGE_HANDLE_STANDARD_ARGS package)
        if (${package}_LIBRARY AND ${package}_INCLUDE_DIR)
            set(${package}_FOUND TRUE)
        endif (${package}_LIBRARY AND ${package}_INCLUDE_DIR)
    ENDMACRO(FIND_PACKAGE_HANDLE_STANDARD_ARGS package)
endif ("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER 2.5)
