#    The MIT License (MIT)
#    Copyright (c) 2016 The SKA organisation
#    Permission is hereby granted, free of charge, to any person obtaining a copy
#    of this software and associated documentation files (the "Software"), to deal
#    in the Software without restriction, including without limitation the rights
#    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#    copies of the Software, and to permit persons to whom the Software is
#    furnished to do so, subject to the following conditions:
#    The above copyright notice and this permission notice shall be included in all
#    copies or substantial portions of the Software.
#    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#    SOFTWARE.

# adapted from https://github.com/bilke/cmake-modules/blob/master/CodeCoverage.cmake

if(CMAKE_BUILD_TYPE MATCHES profile)
    #  macro to generate a coverage report of the specified target
    #  the target will be executed (e.g. 'make targetname' and coverage report generated for that launch.
    macro(COVERAGE_TARGET _targetname )
        set(_outputname coverage_${_targetname})
        set(_outputdir coverage)
        set(coverage_target coverage_report_${_targetname})
        ADD_CUSTOM_TARGET(${coverage_target}

            # Cleanup lcov
            COMMAND ${LCOV_PATH} --directory ${_outputdir} --zerocounters

            # run the target
            COMMAND ${CMAKE_MAKE_PROGRAM} ${_targetname} ${ARGN}

            # Capturing lcov counters and generating report
            COMMAND ${LCOV_PATH} --directory . --capture --output-file ${_outputdir}/${_outputname}.info
            COMMAND ${LCOV_PATH} --remove ${_outputdir}/${_outputname}.info 'panda/*' 'test/*' '/usr/*' --output-file ${_outputdir}/${_outputname}.info.cleaned
            COMMAND ${GENHTML_PATH} -o ${_outputdir}/${_outputname} ${_outputdir}/${_outputname}.info.cleaned
            COMMAND ${CMAKE_COMMAND} -E remove ${_outputdir}/${_outputname}.info ${_outputdir}/${_outputname}.info.cleaned
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        )

        # Show info where to find the report
        ADD_CUSTOM_COMMAND(TARGET coverage_report_${_targetname} POST_BUILD
            COMMAND ;
            COMMENT "see ${_outputdir}/${_outputname}/index.html for the coverage report."
        )

    endmacro(COVERAGE_TARGET)

    # -- find the coverage generators
    find_program( GCOV_PATH gcov )
    find_program( LCOV_PATH lcov )
    find_program( GENHTML_PATH genhtml )

    # -- sanity checking
    if(NOT GCOV_PATH) 
        MESSAGE(FATAL_ERROR "gcov not found. specify with GCOV_PATH")
    endif(NOT GCOV_PATH) 
    
    if(NOT LCOV_PATH) 
        MESSAGE(FATAL_ERROR "lcov not found. specify with LCOV_PATH")
    endif(NOT LCOV_PATH) 

    if(NOT GENHTML_PATH)
        MESSAGE(FATAL_ERROR "genhtml not found. specify with GENHTML_PATH")
    endif(NOT GENHTML_PATH)

    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/coverage)

else(CMAKE_BUILD_TYPE MATCHES profile)

    # -- if coverage is not available, output a suitable message if someone tries to build the target 
    macro(COVERAGE_TARGET _targetname )
        set(coverage_target coverage_report_${_targetname})
        ADD_CUSTOM_TARGET(${coverage_target}
            COMMENT "target ${coverage_target} requires a profile build. Did you forget to set CMAKE_BUILD_TYPE=profile?"
        )
    endmacro(COVERAGE_TARGET)

endif()

# -- targets
coverage_target(test)

