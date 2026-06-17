# MIT License
#
# Copyright (c) 2026 Andrew Smith
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to
# deal in the Software without restriction, including without limitation the
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
# sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
# IN THE SOFTWARE.

# to find the version of cmake do
# $ cmake --version
cmake_minimum_required(VERSION 3.30)
project(
  falco2
  VERSION 0.1.0
  LANGUAGES CXX)

# Set language version used
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED yes)
set(CMAKE_CXX_EXTENSIONS off)  # prevents std=gnu++23
set(CMAKE_EXPORT_COMPILE_COMMANDS on)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

configure_file(config.h.in config.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_compile_options(
  -Wall
  -Wextra
  -Wpedantic
  -Werror
  -Wfatal-errors
)

option(STATIC_ANALYSIS "Enable static analysis" off)
if(STATIC_ANALYSIS)
  include(cmake/static_analysis.cmake)
endif()

file(READ data/LICENSE LICENSE_TEXT)
configure_file(data/LICENSES.inc.in license.h @ONLY)

option(PACKAGE "Enabling packaging mode (experts only)" off)
if(PACKAGE)
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
endif()

if(ORIGINAL_DUPS)
  add_compile_options(-DORIGINAL_DUPS)
endif()

find_package(LIBDEFLATE REQUIRED)
find_package(HTSLIB REQUIRED)
find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)
if(USE_ISAL)
  find_package(ISAL REQUIRED)
  add_compile_options(-DHAVE_ISAL)
endif()

if(MAKE_HTML)
  find_package(fmt)
  add_compile_options(-DMAKE_HTML)
endif()

if(CMAKE_BUILD_TYPE STREQUAL Release)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT result OUTPUT output)
  if(result)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
  else()
    message(STATUS "IPO is not supported in your build environment")
  endif()
endif()

add_library(tile_processor OBJECT src/tile_processor.cpp)
target_link_libraries(tile_processor PRIVATE HTSLIB::HTSLIB)

add_library(contaminants OBJECT src/contaminants.cpp)

add_library(kmer_counter OBJECT src/kmer_counter.cpp)
target_link_libraries(kmer_counter PRIVATE HTSLIB::HTSLIB)

add_library(duplication_results OBJECT src/duplication_results.cpp)
target_link_libraries(duplication_results PRIVATE HTSLIB::HTSLIB)

add_library(fastq_file OBJECT src/fastq_file.cpp)
target_link_libraries(fastq_file PRIVATE HTSLIB::HTSLIB)
if(USE_ISAL)
  target_link_libraries(fastq_file PUBLIC ISAL::ISAL)
endif()

add_library(bam_file OBJECT src/bam_file.cpp)
target_link_libraries(bam_file PRIVATE HTSLIB::HTSLIB)

add_library(falco_file_format OBJECT src/falco_file_format.cpp)
target_link_libraries(falco_file_format PRIVATE HTSLIB::HTSLIB)

add_library(adapter_matcher OBJECT src/adapter_matcher.cpp)
target_link_libraries(adapter_matcher PRIVATE HTSLIB::HTSLIB)

add_library(format_output OBJECT src/format_output.cpp)
target_link_libraries(format_output PRIVATE HTSLIB::HTSLIB)

add_library(falco_utils OBJECT src/falco_utils.cpp)
target_link_libraries(falco_utils PRIVATE HTSLIB::HTSLIB)

add_executable(falco2 src/falco2.cpp)
target_link_libraries(falco2 PRIVATE
  tile_processor
  contaminants
  adapter_matcher
  duplication_results
  kmer_counter
  fastq_file
  bam_file
  format_output
  falco_file_format
  falco_utils
  HTSLIB::HTSLIB
  LIBDEFLATE::LIBDEFLATE
  ZLIB::ZLIB
  Threads::Threads
)
if(MAKE_HTML)
  target_link_libraries(falco2 PRIVATE fmt::fmt)
endif()

install(TARGETS falco2 RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
