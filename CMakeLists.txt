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
  falco
  VERSION 2.0.0
  LANGUAGES CXX)

# Set language version used
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED yes)
set(CMAKE_CXX_EXTENSIONS off)  # prevents std=gnu++23
set(CMAKE_EXPORT_COMPILE_COMMANDS on)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
list(APPEND CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/cmake")

if(STATIC_ANALYSIS)
  include(cmake/static_analysis.cmake)
endif()

add_library(tile_processor OBJECT src/tile_processor.cpp)
add_library(kmer_counter OBJECT src/kmer_counter.cpp)
add_library(fastq_file OBJECT src/fastq_file.cpp)
add_library(duplication_results OBJECT src/duplication_results.cpp)
add_library(adapter_matcher OBJECT src/adapter_matcher.cpp)

add_executable(falco2 src/falco2.cpp)
target_link_libraries(falco2 PUBLIC
  tile_processor
  adapter_matcher
  duplication_results
  fastq_file
  kmer_counter
)

# find_package(HTSLIB REQUIRED)
# if(HTSLIB_FOUND)
#   target_include_directories(falco2 PUBLIC src/bamxx)
#   target_link_libraries(falco2 PUBLIC HTSLIB::HTSLIB)
# else()
#   message(FATAL_ERROR "HTSLib required but not found")
# endif()
