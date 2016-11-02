# -- find valgrind
find_program( VALGRIND_PATH valgrind )

if(NOT VALGRIND_PATH) 
    # add a valgrind not installed message for the target
    macro(VALGRIND_TARGET _targetname )
        set(valgrind_target valgrind_report_${_targetname})
        ADD_CUSTOM_TARGET(${valgrind_target}
            COMMENT "target ${valgrind_target} requires valgrind to be available"
        )
    endmacro(VALGRIND_TARGET)
else(NOT VALGRIND_PATH) 
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/valgrind)
    macro(VALGRIND_TARGET _targetname )
        set(_outputdir valgrind)
        set(_outputname valgrind_${_targetname})
        set(valgrind_target valgrind_report_${_targetname})
        ADD_CUSTOM_TARGET(${valgrind_target}
            COMMAND ${VALGRIND_PATH} --xml=yes --xml-file=${_outputdir}/${_outputname} --trace-children=yes ${CMAKE_MAKE_PROGRAM} ${_targetname}
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        )
        ADD_CUSTOM_COMMAND(TARGET valgrind_report_${_targetname} POST_BUILD
            COMMAND ;
            COMMENT "see ${_outputdir}/${_outputname} for the valgrind report."
        )
    endmacro(VALGRIND_TARGET)
endif(NOT VALGRIND_PATH) 


# -- targets
valgrind_target(test)
