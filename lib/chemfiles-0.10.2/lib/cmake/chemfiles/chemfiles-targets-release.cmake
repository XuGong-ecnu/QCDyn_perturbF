#----------------------------------------------------------------
# Generated CMake target import file for configuration "release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "chemfiles" for configuration "release"
set_property(TARGET chemfiles APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(chemfiles PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C;CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libchemfiles.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS chemfiles )
list(APPEND _IMPORT_CHECK_FILES_FOR_chemfiles "${_IMPORT_PREFIX}/lib/libchemfiles.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
