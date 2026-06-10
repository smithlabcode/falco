# SPDX-License-Identifier: MIT; (c) 2026 Andrew D Smith (author)
#[=======================================================================[.rst:
FindISAL
--------

Find the native ISA-L includes and library. Based on the ZLIB module (BSD-3-Clause).

#]=======================================================================]


cmake_policy(PUSH)
cmake_policy(SET CMP0159 NEW) # file(STRINGS) with REGEX updates CMAKE_MATCH_<n>

if(ISAL_FIND_COMPONENTS AND NOT ISAL_FIND_QUIETLY)
  message(AUTHOR_WARNING
    "ISAL does not provide any COMPONENTS.  Calling\n"
    "  find_package(ISAL COMPONENTS ...)\n"
    "will always fail."
    )
endif()

set(_ISAL_SEARCHES)

# Search ISAL_ROOT first if it is set.
if(ISAL_ROOT)
  set(_ISAL_SEARCH_ROOT PATHS ${ISAL_ROOT} NO_DEFAULT_PATH)
  list(APPEND _ISAL_SEARCHES _ISAL_SEARCH_ROOT)
endif()

# Normal search.
set(_ISAL_x86 "(x86)")
set(_ISAL_SEARCH_NORMAL
    PATHS "[HKEY_LOCAL_MACHINE\\SOFTWARE\\GnuWin32\\isa-l;InstallPath]"
          "$ENV{ProgramFiles}/isa-l"
          "$ENV{ProgramFiles${_ISAL_x86}}/isa-l")
unset(_ISAL_x86)
list(APPEND _ISAL_SEARCHES _ISAL_SEARCH_NORMAL)

if(ISAL_USE_STATIC_LIBS)
  set(ISAL_NAMES isal)
  set(ISAL_NAMES_DEBUG isal)
else()
  set(ISAL_NAMES isal)
  set(ISAL_NAMES_DEBUG isal)
endif()

# Try each search configuration.
foreach(search ${_ISAL_SEARCHES})
  find_path(ISAL_INCLUDE_DIR NAMES isa-l.h ${${search}} PATH_SUFFIXES include)
endforeach()

# Allow ISAL_LIBRARY to be set manually, as the location of the isa-l library
if(NOT ISAL_LIBRARY)
  if(DEFINED CMAKE_FIND_LIBRARY_PREFIXES)
    set(_isal_ORIG_CMAKE_FIND_LIBRARY_PREFIXES "${CMAKE_FIND_LIBRARY_PREFIXES}")
  else()
    set(_isal_ORIG_CMAKE_FIND_LIBRARY_PREFIXES)
  endif()
  if(DEFINED CMAKE_FIND_LIBRARY_SUFFIXES)
    set(_isal_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES "${CMAKE_FIND_LIBRARY_SUFFIXES}")
  else()
    set(_isal_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES)
  endif()
  # Prefix/suffix of the win32/Makefile.gcc build
  if(WIN32)
    list(APPEND CMAKE_FIND_LIBRARY_PREFIXES "" "lib")
    list(APPEND CMAKE_FIND_LIBRARY_SUFFIXES ".dll.a")
  endif()
  # Support preference of static libs by adjusting CMAKE_FIND_LIBRARY_SUFFIXES
  if(ISAL_USE_STATIC_LIBS)
    if(WIN32)
      set(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
    else()
      set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
    endif()
  endif()

  foreach(search ${_ISAL_SEARCHES})
    find_library(ISAL_LIBRARY_RELEASE NAMES ${ISAL_NAMES} NAMES_PER_DIR ${${search}} PATH_SUFFIXES lib)
    find_library(ISAL_LIBRARY_DEBUG NAMES ${ISAL_NAMES_DEBUG} NAMES_PER_DIR ${${search}} PATH_SUFFIXES lib)
  endforeach()

  # Restore the original find library ordering
  if(DEFINED _isal_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES)
    set(CMAKE_FIND_LIBRARY_SUFFIXES "${_isal_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES}")
  else()
    set(CMAKE_FIND_LIBRARY_SUFFIXES)
  endif()
  if(DEFINED _isal_ORIG_CMAKE_FIND_LIBRARY_PREFIXES)
    set(CMAKE_FIND_LIBRARY_PREFIXES "${_isal_ORIG_CMAKE_FIND_LIBRARY_PREFIXES}")
  else()
    set(CMAKE_FIND_LIBRARY_PREFIXES)
  endif()

  include(SelectLibraryConfigurations)
  select_library_configurations(ISAL)
endif()

unset(ISAL_NAMES)
unset(ISAL_NAMES_DEBUG)

mark_as_advanced(ISAL_INCLUDE_DIR)

if(ISAL_INCLUDE_DIR AND EXISTS "${ISAL_INCLUDE_DIR}/isa-l.h")
  file(STRINGS "${ISAL_INCLUDE_DIR}/isa-l.h" ISAL_MAJOR_LIST REGEX "^#define ISAL_MAJOR_VERSION")
  list(GET ISAL_MAJOR_LIST 0 ISAL_MAJOR)  # Take the first matching line
  if (ISAL_MAJOR MATCHES "#define[ \t]+ISAL_MAJOR_VERSION[ \t]+\([0-9]+\)")
    set(ISAL_VERSION_MAJOR "${CMAKE_MATCH_1}")
  else()
    set(ISAL_VERSION_MAJOR "0")
  endif()

  file(STRINGS "${ISAL_INCLUDE_DIR}/isa-l.h" ISAL_MINOR_LIST REGEX "^#define ISAL_MINOR_VERSION")
  list(GET ISAL_MINOR_LIST 0 ISAL_MINOR)  # Take the first matching line
  if (ISAL_MINOR MATCHES "#define[ \t]+ISAL_MINOR_VERSION[ \t]+\([0-9]+\)")
    set(ISAL_VERSION_MINOR "${CMAKE_MATCH_1}")
  else()
    set(ISAL_VERSION_MINOR "0")
  endif()

  file(STRINGS "${ISAL_INCLUDE_DIR}/isa-l.h" ISAL_PATCH_LIST REGEX "^#define ISAL_PATCH_VERSION")
  list(GET ISAL_PATCH_LIST 0 ISAL_PATCH)  # Take the first matching line
  if (ISAL_PATCH MATCHES "#define[ \t]+ISAL_PATCH_VERSION[ \t]+\([0-9]+\)")
    set(ISAL_VERSION_PATCH "${CMAKE_MATCH_1}")
  else()
    set(ISAL_VERSION_PATCH "0")
  endif()

  set(ISAL_VERSION "${ISAL_VERSION_MAJOR}.${ISAL_VERSION_MINOR}.${ISAL_VERSION_PATCH}")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ISAL REQUIRED_VARS ISAL_LIBRARY ISAL_INCLUDE_DIR
                                       VERSION_VAR ISAL_VERSION
                                       HANDLE_COMPONENTS)

if(ISAL_FOUND)
    set(ISAL_INCLUDE_DIRS ${ISAL_INCLUDE_DIR})

    if(NOT ISAL_LIBRARIES)
      set(ISAL_LIBRARIES ${ISAL_LIBRARY})
    endif()

    if(NOT TARGET ISAL::ISAL)
      add_library(ISAL::ISAL UNKNOWN IMPORTED)
      set_target_properties(ISAL::ISAL PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${ISAL_INCLUDE_DIRS}")

      if(ISAL_LIBRARY_RELEASE)
        set_property(TARGET ISAL::ISAL APPEND PROPERTY
          IMPORTED_CONFIGURATIONS RELEASE)
        set_target_properties(ISAL::ISAL PROPERTIES
          IMPORTED_LOCATION_RELEASE "${ISAL_LIBRARY_RELEASE}")
      endif()

      if(ISAL_LIBRARY_DEBUG)
        set_property(TARGET ISAL::ISAL APPEND PROPERTY
          IMPORTED_CONFIGURATIONS DEBUG)
        set_target_properties(ISAL::ISAL PROPERTIES
          IMPORTED_LOCATION_DEBUG "${ISAL_LIBRARY_DEBUG}")
      endif()

      if(NOT ISAL_LIBRARY_RELEASE AND NOT ISAL_LIBRARY_DEBUG)
        set_property(TARGET ISAL::ISAL APPEND PROPERTY
          IMPORTED_LOCATION "${ISAL_LIBRARY}")
      endif()
    endif()
endif()

cmake_policy(POP)
