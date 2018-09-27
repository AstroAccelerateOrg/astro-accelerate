# Store the git hash of the current head
# based on blog of David Gobbi (http://www.cognitive-antics.net/?p=816)

if(EXISTS "${PROJECT_SOURCE_DIR}/.git/HEAD")
  file(READ "${PROJECT_SOURCE_DIR}/.git/HEAD" PROJECT_SOURCE_VERSION)
  if("${PROJECT_SOURCE_VERSION}" MATCHES "^ref:")
    string(REGEX REPLACE "^ref: *([^ \n\r]*).*" "\\1" PROJECT_GIT_REF "${PROJECT_SOURCE_VERSION}")
    set(git_version_file "${PROJECT_SOURCE_DIR}/.git/${PROJECT_GIT_REF}")
    if(EXISTS "${git_version_file}")
        file(READ "${git_version_file}" PROJECT_SOURCE_VERSION)
    endif(EXISTS "${git_version_file}")
  endif()
  string(STRIP "${PROJECT_SOURCE_VERSION}" PROJECT_SOURCE_VERSION)
endif()
