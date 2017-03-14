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

# thirdparty dependencies
include(cmake/googletest.cmake)

# other dependencies
set(Boost_ADDITIONAL_VERSIONS "1.50 1.50.0")
find_package(Boost COMPONENTS filesystem system program_options REQUIRED)

include_directories(SYSTEM ${BOOST_INCLUDE_DIR})
set(DEPENDENCY_LIBRARIES
#    ${Boost_LIBRARIES}
    )

include(compiler_settings)
include(cmake/cuda.cmake)
include(cmake/doxygen.cmake)
