set(ENABLE_DOC "true")
set(doc_all_target "ALL")

if(ENABLE_DOC)
find_package(Doxygen REQUIRED)

if(DOXYGEN_FOUND)
set(DOXYGEN_BUILD_DIR ${CMAKE_BINARY_DIR}/doc)
set(CMAKE_DOC_OUTPUT_DIR ${CMAKE_BINARY_DIR}/doc)
file(MAKE_DIRECTORY ${DOXYGEN_BUILD_DIR})
set(DOXYFILE_IN ${PROJECT_SOURCE_DIR}/cmake/DoxyfileAPI.in)
configure_file(${DOXYFILE_IN} ${CMAKE_BINARY_DIR}/DoxyfileAPI @ONLY)
add_custom_target(doc ${doc_all_target}
      ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/DoxyfileAPI
      WORKING_DIRECTORY ${DOXYGEN_BUILD_DIR}
      COMMENT "Generating documentation with Doxygen" VERBATIM)
install(DIRECTORY ${DOXYGEN_BUILD_DIR}/ DESTINATION ${DOC_INSTALL_DIR} PATTERN "${DOXYGEN_BUILD_DIR}/*")
else(DOXYGEN_FOUND)
    add_custom_target(doc ${doc_all_target} COMMAND ${CMAKE_COMMAND} -E echo COMMENT "No doc target configured. Doxygen not found" VERBATIM)
   endif(DOXYGEN_FOUND)
else(ENABLE_DOC)
    add_custom_target(doc ${doc_all_target} COMMAND ${CMAKE_COMMAND} -E echo COMMENT "No doc target configured. Rebuild with -DENABLE_DOC=true" VERBATIM)
endif(ENABLE_DOC)