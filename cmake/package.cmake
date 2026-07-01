# MIT License
#
# Copyright (c) 2026 Andrew Smith
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

message(STATUS "Configuring to package for distribution")

set(LIBDEFLATE_USE_STATIC_LIBS on)
set(HTSLIB_USE_STATIC_LIBS on)
set(ZLIB_USE_STATIC_LIBS on)
if(USE_ISAL)
  set(ISAL_USE_STATIC_LIBS on)
endif()
add_link_options(-static-libgcc -static-libstdc++)
if(NOT CMAKE_SYSTEM_NAME STREQUAL "Darwin")  # not macOS
  add_link_options(-s)
endif()

set(CPACK_PACKAGE_VENDOR "Andrew D Smith")
set(CPACK_PACKAGE_CONTACT "Andrew D Smith")
set(CPACK_RESOURCE_FILE_LICENSE ${PROJECT_SOURCE_DIR}/LICENSE)
set(CPACK_GENERATOR "TGZ;STGZ")

if(SYSTEM_NAME_SUFFIX)  # For macOS releases to get the OS version right
  set(CPACK_PACKAGE_FILE_NAME
    ${PROJECT_NAME}-${PROJECT_VERSION}-${SYSTEM_NAME_SUFFIX}
  )
endif()
message(STATUS
  "Package file name (without suffixes): ${CPACK_PACKAGE_FILE_NAME}"
)

if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")  # macOS
  set(CPACK_SYSTEM_NAME "macOS")
endif()

# For custom STGZ generator
set(CPACK_STGZ_HEADER_FILE "cmake/CPack.STGZ_Header.sh.in")

# Only make the Debian installer if building on Linux
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  # Get the system processor architecture
  message(STATUS "System Architecture: ${CMAKE_SYSTEM_PROCESSOR}")

  # Run 'ldd --version' to get glibc version
  execute_process(
    COMMAND sh -c "ldd --version | head -1 | awk '{print \$NF}'"
    OUTPUT_VARIABLE GLIBC_VERSION
    ERROR_VARIABLE LDD_ERROR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  if(LDD_ERROR OR NOT GLIBC_VERSION)
    message(FATAL_ERROR "Failed to detect glibc version")
  else()
    message(STATUS "Detected glibc version: ${GLIBC_VERSION}")
  endif()

  # Check that the DEB can be built
  find_program(DPKG_FOUND dpkg)
  if (NOT DPKG_FOUND)
    message(STATUS "Failed to detect dpkg; skipping DEB packaging")
  else()
    list(APPEND CPACK_GENERATOR "DEB")
    set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Andrew D Smith <andrewds@usc.edu>")
    set(CPACK_DEBIAN_PACKAGE_HOMEPAGE
      "https://github.com/smithlabcode/falco"
    )
    # Set dependency libc6 with the version constraint
    set(CPACK_DEBIAN_PACKAGE_DEPENDS "libc6 (>= ${GLIBC_VERSION})")
    # Check if the architecture is x86_64 and set to amd64 for Debian packaging;
    # this is not the same as for RPMs.
    if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
      set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
    elseif(CMAKE_SYSTEM_PROCESSOR STREQUAL "aarch64")
      set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "arm64")
    else()
      set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE ${CMAKE_SYSTEM_PROCESSOR})
    endif()
    message(STATUS
      "Debian package architecture: ${CPACK_DEBIAN_PACKAGE_ARCHITECTURE}")
  endif()

  # Check that the RPM can be built
  find_program(RPMBUILD_FOUND rpmbuild)
  if (NOT RPMBUILD_FOUND)
    message(STATUS "Failed to detect rpmbuild; skipping RPM packaging")
  else()
    list(APPEND CPACK_GENERATOR "RPM")
    set(CPACK_RPM_PACKAGE_SUMMARY ${PROJECT_DESCRIPTION})
    set(CPACK_RPM_PACKAGE_LICENSE "MIT")
    set(CPACK_RPM_PACKAGE_NAME ${PROJECT_NAME})
    set(CPACK_RPM_FILE_NAME RPM_DEFAULT)
    set(CPACK_RPM_FILE_NAME ${CPACK_PACKAGE_VERSION})
    set(CPACK_RPM_PACKAGE_RELEASE 1)
    set(CPACK_RPM_PACKAGE_URL ${PROJECT_HOMEPAGE})
    set(CPACK_RPM_PACKAGE_DESCRIPTION ${PROJECT_DESCRIPTION})
    # RPM automatically deduces dependencies in binaries
  endif()
endif()

# For Source (only package source from a clean clone)
set(CPACK_SOURCE_IGNORE_FILES /build/ /.git/ /.github/ /.gitignore)
set(CPACK_SOURCE_GENERATOR "TGZ")

include(CPack)
