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

