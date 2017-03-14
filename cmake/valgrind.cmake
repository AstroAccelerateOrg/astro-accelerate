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
