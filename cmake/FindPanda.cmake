# Find the native PANDA includes and library
#
#  PANDA_INSTALL_DIR - top where panda framwork has been installed (lib and include dirs included)
#  PANDA_LIBRARY_DIR - explicitly define directory where to find panda libraries
#  PANDA_INCLUDE_DIR - where to find panda includes
#  PANDA_LIBRARIES   - List of libraries when using panda.
#  PANDA_TEST_LIBRARIES - List of panda test support libraries
#  PANDA_FOUND       - True if panda found.

# Already in cache, be silent
message("searching ${PANDA_INSTALL_DIR}")

IF (PANDA_INCLUDE_DIR)
    SET(PANDA_INC_DIR ${PANDA_INCLUDE_DIR})
    UNSET(PANDA_INCLUDE_DIR)
ENDIF (PANDA_INCLUDE_DIR)

FIND_PATH(PANDA_INCLUDE_DIR panda/Version.h 
    PATHS ${PANDA_INC_DIR}
	  ${PANDA_INSTALL_DIR}/include
	  /usr/local/include 
	  /usr/include )
message("Found ${PANDA_INCLUDE_DIR} : ${PANDA_INSTALL_DIR}")

SET(PANDA_NAMES panda)
FOREACH( lib ${PANDA_NAMES} )
    FIND_LIBRARY(PANDA_LIBRARY_${lib} 
	NAMES ${lib}
	PATHS ${PANDA_LIBRARY_DIR} ${PANDA_INSTALL_DIR} ${PANDA_INSTALL_DIR}/lib /usr/local/lib /usr/lib
    )
    LIST(APPEND PANDA_LIBRARIES ${PANDA_LIBRARY_${lib}})
ENDFOREACH(lib)

SET(PANDA_TEST_NAMES panda_testutils)
FOREACH( lib ${PANDA_TEST_NAMES} )
    FIND_LIBRARY(PANDA_LIBRARY_${lib} 
	NAMES ${lib}
	PATHS ${PANDA_LIBRARY_DIR} ${PANDA_LIBRARY_DIR}/test ${PANDA_INSTALL_DIR}/lib /usr/local/lib /usr/lib
    )
    LIST(APPEND PANDA_TEST_LIBRARIES ${PANDA_LIBRARY_${lib}})
ENDFOREACH(lib)

# handle the QUIETLY and REQUIRED arguments and set PANDA_FOUND to TRUE if.
# all listed variables are TRUE
INCLUDE(FindPackageHandleCompat)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(PANDA DEFAULT_MSG PANDA_LIBRARIES PANDA_INCLUDE_DIR)

IF(NOT PANDA_FOUND)
    SET( PANDA_LIBRARIES )
    SET( PANDA_TEST_LIBRARIES )
ENDIF(NOT PANDA_FOUND)

MARK_AS_ADVANCED(PANDA_LIBRARIES PANDA_TEST_LIBRARIES PANDA_INCLUDE_DIR)

