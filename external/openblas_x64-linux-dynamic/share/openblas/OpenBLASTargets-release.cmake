#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "OpenBLAS::OpenBLAS" for configuration "Release"
set_property(TARGET OpenBLAS::OpenBLAS APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(OpenBLAS::OpenBLAS PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libopenblas.so.0.3"
  IMPORTED_SONAME_RELEASE "libopenblas.so.0"
  )

list(APPEND _cmake_import_check_targets OpenBLAS::OpenBLAS )
list(APPEND _cmake_import_check_files_for_OpenBLAS::OpenBLAS "${_IMPORT_PREFIX}/lib/libopenblas.so.0.3" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
