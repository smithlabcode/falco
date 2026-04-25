# SPDX-License-Identifier: MIT; (c) 2025 Andrew D Smith (author)
#[=======================================================================[.rst:
FindHTSLIB
--------

Find the native HTSLib includes and library. Based on the ZLIB module.

#]=======================================================================]

cmake_policy(PUSH)
cmake_policy(SET CMP0159 NEW) # file(STRINGS) with REGEX updates CMAKE_MATCH_<n>

if(HTSLIB_FIND_COMPONENTS AND NOT HTSLIB_FIND_QUIETLY)
  message(AUTHOR_WARNING
    "HTSLib does not provide any COMPONENTS.  Calling\n"
    "  find_package(HTSLIB COMPONENTS ...)\n"
    "will always fail."
  )
endif()

set(_HTSLIB_SEARCHES)

# Search HTSLIB_ROOT first if it is set.
if(HTSLIB_ROOT)
  set(_HTSLIB_SEARCH_ROOT PATHS ${HTSLIB_ROOT} NO_DEFAULT_PATH)
  list(APPEND _HTSLIB_SEARCHES _HTSLIB_SEARCH_ROOT)
endif()

# Normal search.
# Windows stuff
set(_HTSLIB_x86 "(x86)")
set(_HTSLIB_SEARCH_NORMAL
  PATHS "$ENV{ProgramFiles}/htslib"
  "$ENV{ProgramFiles${_HTSLIB_x86}}/htslib")
unset(_HTSLIB_x86)
list(APPEND _HTSLIB_SEARCHES _HTSLIB_SEARCH_NORMAL)

if(HTSLIB_USE_STATIC_LIBS)
  set(HTSLIB_NAMES hts)
  set(HTSLIB_NAMES_DEBUG hts)
else()
  set(HTSLIB_NAMES hts)
  set(HTSLIB_NAMES_DEBUG hts)
endif()

# Try each search configuration.
foreach(search ${_HTSLIB_SEARCHES})
  find_path(HTSLIB_INCLUDE_DIR NAMES htslib ${${search}} PATH_SUFFIXES include)
endforeach()

# Allow HTSLIB_LIBRARY to be set manually, as the location of the htslib library
if(NOT HTSLIB_LIBRARY)
  if(DEFINED CMAKE_FIND_LIBRARY_PREFIXES)
    set(_htslib_ORIG_CMAKE_FIND_LIBRARY_PREFIXES "${CMAKE_FIND_LIBRARY_PREFIXES}")
  else()
    set(_htslib_ORIG_CMAKE_FIND_LIBRARY_PREFIXES)
  endif()
  if(DEFINED CMAKE_FIND_LIBRARY_SUFFIXES)
    set(_htslib_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES "${CMAKE_FIND_LIBRARY_SUFFIXES}")
  else()
    set(_htslib_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES)
  endif()
  # Prefix/suffix of the win32/Makefile.gcc build
  if(WIN32)
    list(APPEND CMAKE_FIND_LIBRARY_PREFIXES "" "lib")
    list(APPEND CMAKE_FIND_LIBRARY_SUFFIXES ".dll.a")
  endif()
  # Support preference of static libs by adjusting CMAKE_FIND_LIBRARY_SUFFIXES
  if(HTSLIB_USE_STATIC_LIBS)
    if(WIN32)
      set(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
    else()
      set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
    endif()
  endif()

  foreach(search ${_HTSLIB_SEARCHES})
    find_library(HTSLIB_LIBRARY_RELEASE NAMES ${HTSLIB_NAMES} NAMES_PER_DIR ${${search}} PATH_SUFFIXES lib)
    find_library(HTSLIB_LIBRARY_DEBUG NAMES ${HTSLIB_NAMES_DEBUG} NAMES_PER_DIR ${${search}} PATH_SUFFIXES lib)
  endforeach()

  # Restore the original find library ordering
  if(DEFINED _htslib_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES)
    set(CMAKE_FIND_LIBRARY_SUFFIXES "${_htslib_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES}")
  else()
    set(CMAKE_FIND_LIBRARY_SUFFIXES)
  endif()
  if(DEFINED _htslib_ORIG_CMAKE_FIND_LIBRARY_PREFIXES)
    set(CMAKE_FIND_LIBRARY_PREFIXES "${_htslib_ORIG_CMAKE_FIND_LIBRARY_PREFIXES}")
  else()
    set(CMAKE_FIND_LIBRARY_PREFIXES)
  endif()

  include(SelectLibraryConfigurations)
  select_library_configurations(HTSLIB)
endif()

unset(HTSLIB_NAMES)
unset(HTSLIB_NAMES_DEBUG)

mark_as_advanced(HTSLIB_INCLUDE_DIR)

if(HTSLIB_INCLUDE_DIR AND EXISTS "${HTSLIB_INCLUDE_DIR}/htslib/hts.h")
  # Example: #define HTS_VERSION 101300
  file(STRINGS "${HTSLIB_INCLUDE_DIR}/htslib/hts.h" HTSLIB_H_LIST REGEX "^#define HTS_VERSION")
  list(GET HTSLIB_H_LIST 0 HTSLIB_H)  # Take the first matching line
  if (HTSLIB_H MATCHES "#define[ \t]+HTS_VERSION[ \t]+\([0-9]+\)")
    set(NUMERIC_VERSION "${CMAKE_MATCH_1}")
    # Extract digits by position in string
    # XYYYZZ => X = major, YYY = minor, ZZ = patch
    string(SUBSTRING "${NUMERIC_VERSION}" 0 1 HTSLIB_VERSION_MAJOR)
    string(SUBSTRING "${NUMERIC_VERSION}" 1 3 HTSLIB_VERSION_MINOR)
    string(SUBSTRING "${NUMERIC_VERSION}" 4 2 HTSLIB_VERSION_PATCH)
  else()
    set(HTSLIB_VERSION_STRING "")
    set(HTSLIB_VERSION_MAJOR "")
    set(HTSLIB_VERSION_MINOR "")
    set(HTSLIB_VERSION_PATCH "")
  endif()

  # Make sure the version numbers don't have leading zeros
  # The minor version is encoded in such a way that it often will
  foreach(part HTSLIB_VERSION_MAJOR HTSLIB_VERSION_MINOR HTSLIB_VERSION_PATCH)
    if(${${part}} MATCHES "^[0]+$")
      set(${part} "0")
    else()
      string(REGEX REPLACE "^0+" "" ${part} "${${part}}")
    endif()
  endforeach()

  # Set canonical variables
  set(HTSLIB_MAJOR_VERSION "${HTSLIB_VERSION_MAJOR}")
  set(HTSLIB_MINOR_VERSION "${HTSLIB_VERSION_MINOR}")
  set(HTSLIB_PATCH_VERSION "${HTSLIB_VERSION_PATCH}")
  # Build the standard version string
  set(HTSLIB_VERSION "${HTSLIB_VERSION_MAJOR}.${HTSLIB_VERSION_MINOR}")
  # Only append patch if it's not "00"
  if(NOT HTSLIB_VERSION_PATCH STREQUAL "00")
    set(HTSLIB_VERSION "${HTSLIB_VERSION}.${HTSLIB_VERSION_PATCH}")
  endif()
endif()

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(
  HTSLIB
  REQUIRED_VARS
  HTSLIB_LIBRARY
  HTSLIB_INCLUDE_DIR
  VERSION_VAR
  HTSLIB_VERSION
  HANDLE_COMPONENTS
)

if(HTSLIB_FOUND)
  set(HTSLIB_INCLUDE_DIRS ${HTSLIB_INCLUDE_DIR})
  if(NOT HTSLIB_LIBRARIES)
    set(HTSLIB_LIBRARIES ${HTSLIB_LIBRARY})
  endif()
  if(NOT TARGET HTSLIB::HTSLIB)
    add_library(HTSLIB::HTSLIB UNKNOWN IMPORTED)
    set_target_properties(HTSLIB::HTSLIB PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${HTSLIB_INCLUDE_DIRS}")
    if(HTSLIB_LIBRARY_RELEASE)
      set_property(TARGET HTSLIB::HTSLIB APPEND PROPERTY
        IMPORTED_CONFIGURATIONS RELEASE)
      set_target_properties(HTSLIB::HTSLIB PROPERTIES
        IMPORTED_LOCATION_RELEASE "${HTSLIB_LIBRARY_RELEASE}")
    endif()
    if(HTSLIB_LIBRARY_DEBUG)
      set_property(TARGET HTSLIB::HTSLIB APPEND PROPERTY
        IMPORTED_CONFIGURATIONS DEBUG)
      set_target_properties(HTSLIB::HTSLIB PROPERTIES
        IMPORTED_LOCATION_DEBUG "${HTSLIB_LIBRARY_DEBUG}")
    endif()
    if(NOT HTSLIB_LIBRARY_RELEASE AND NOT HTSLIB_LIBRARY_DEBUG)
      set_property(TARGET HTSLIB::HTSLIB APPEND PROPERTY
        IMPORTED_LOCATION "${HTSLIB_LIBRARY}")
    endif()
  endif()
endif()

cmake_policy(POP)
