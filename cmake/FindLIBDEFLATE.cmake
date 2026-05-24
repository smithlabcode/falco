# SPDX-License-Identifier: GPL-3.0-or-later; (c) 2025 Andrew D Smith (author)
#[=======================================================================[.rst:
FindLIBDEFLATE
--------------

Find the native libdeflate includes and library.

#]=======================================================================]

# FindLIBDEFLATE.cmake
# Custom CMake module to find libdeflate

# Support preference of static libs by adjusting CMAKE_FIND_LIBRARY_SUFFIXES
# ADS: this is taken from the FindBoost.cmake file
if(LIBDEFLATE_USE_STATIC_LIBS)
  set(_libdeflate_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES
    ${CMAKE_FIND_LIBRARY_SUFFIXES}
  )
  if(WIN32)
    list(INSERT CMAKE_FIND_LIBRARY_SUFFIXES 0 .lib .a)
  else()
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
  endif()
endif()

find_path(LIBDEFLATE_INCLUDE_DIR NAMES libdeflate.h)
find_library(LIBDEFLATE_LIBRARY NAMES deflate libdeflate)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LIBDEFLATE
  REQUIRED_VARS LIBDEFLATE_LIBRARY LIBDEFLATE_INCLUDE_DIR
  VERSION_VAR LIBDEFLATE_VERSION
)

if(LIBDEFLATE_FOUND AND NOT TARGET LIBDEFLATE::LIBDEFLATE)
  add_library(LIBDEFLATE::LIBDEFLATE UNKNOWN IMPORTED)
  set_target_properties(LIBDEFLATE::LIBDEFLATE PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${LIBDEFLATE_INCLUDE_DIR}"
    IMPORTED_LOCATION "${LIBDEFLATE_LIBRARY}"
  )
endif()

# Restore the original find library ordering
# ADS: this is take from the FindBoost.cmake file
if(LIBDEFLATE_USE_STATIC_LIBS)
  set(CMAKE_FIND_LIBRARY_SUFFIXES
    ${_libdeflate_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES}
  )
endif()
