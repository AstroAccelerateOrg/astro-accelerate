if(ENABLE_DOC)
  # add a target to generate API documentation with Doxygen
  find_package(Doxygen)
  if(DOXYGEN_FOUND)
    set(DOXYGEN_BUILD_DIR ${CMAKE_BINARY_DIR}/doc)
    file(MAKE_DIRECTORY ${DOXYGEN_BUILD_DIR})
    configure_file(${CMAKE_SOURCE_DIR}/cmake/Doxyfile.in ${CMAKE_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doc
      ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/Doxyfile
      WORKING_DIRECTORY ${DOXYGEN_BUILD_DIR}
      COMMENT "Generating API documentation with Doxygen" VERBATIM)

    install(DIRECTORY ${DOXYGEN_BUILD_DIR}/
            DESTINATION ${DOC_INSTALL_DIR}
            PATTERN "${DOXYGEN_BUILD_DIR}/*"
           )
   endif(DOXYGEN_FOUND)
else(ENABLE_DOC)
    add_custom_target(doc
      COMMENT "No doc traget configured. Rebuild with -DENABLE_DOC=true" VERBATIM
    )
endif(ENABLE_DOC)
