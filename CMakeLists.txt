# SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

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

include(CheckIncludeFileCXX)
include(CheckFunctionExists)
include(CheckCXXCompilerFlag)
include(GNUInstallDirs)

include(cmake/configure_data.cmake)

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

option(FALCO_TESTS "Enable testing (ctest or make test)" on)
if(FALCO_TESTS)
  include(cmake/tests.cmake)
endif()

file(READ data/LICENSE LICENSE_TEXT)
if(USE_ISAL)
  file(READ data/ISA_license_text.txt ISAL_LICENSE_TEXT)
  string(APPEND LICENSE_TEXT ${ISAL_LICENSE_TEXT})
endif()
configure_file(data/LICENSES.inc.in license.h @ONLY)

option(PACKAGE "Enabling packaging mode (experts only)" off)
if(PACKAGE)
  include(cmake/package.cmake)
endif()

if(ORIGINAL_DUPS)
  add_compile_options(-DORIGINAL_DUPS)
endif()

find_package(LIBDEFLATE REQUIRED)
find_package(HTSLIB 1.18 REQUIRED)
find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)
if(USE_ISAL)
  find_package(ISAL REQUIRED)
  add_compile_options(-DHAVE_ISAL)
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

add_library(adapter_set OBJECT src/adapter_set.cpp)
add_library(falco_utils OBJECT src/falco_utils.cpp)
add_library(base_groups OBJECT src/base_groups.cpp)

add_library(bgzf_block OBJECT src/bgzf_block.cpp)
target_link_libraries(bgzf_block PUBLIC LIBDEFLATE::LIBDEFLATE)

add_library(bgzf_reader OBJECT src/bgzf_reader.cpp)
target_link_libraries(bgzf_reader PUBLIC bgzf_block)

add_library(fastq_file OBJECT src/fastq_file.cpp)
target_link_libraries(fastq_file PUBLIC HTSLIB::HTSLIB)

add_library(fastq_stdin OBJECT src/fastq_stdin.cpp)

add_library(fastq_gz_file OBJECT src/fastq_gz_file.cpp)
target_link_libraries(fastq_gz_file PUBLIC HTSLIB::HTSLIB)
if(USE_ISAL)
  target_link_libraries(fastq_gz_file PUBLIC ISAL::ISAL)
endif()

add_library(fastq_bgzf_file OBJECT src/fastq_bgzf_file.cpp)
target_link_libraries(fastq_bgzf_file PUBLIC
  HTSLIB::HTSLIB
  LIBDEFLATE::LIBDEFLATE
)

add_library(bam_header OBJECT src/bam_header.cpp)

add_library(bam_file OBJECT src/bam_file.cpp)
target_link_libraries(bam_file PUBLIC HTSLIB::HTSLIB)

add_library(samrec OBJECT src/samrec.cpp)
add_library(sam_file OBJECT src/sam_file.cpp)
target_link_libraries(sam_file PUBLIC HTSLIB::HTSLIB)

add_library(falco_file_format OBJECT src/falco_file_format.cpp)
target_link_libraries(falco_file_format PUBLIC HTSLIB::HTSLIB)

add_library(quality_score OBJECT src/quality_score.cpp)
target_link_libraries(quality_score PUBLIC falco_file_format)

add_library(tile_processor OBJECT src/tile_processor.cpp)
target_link_libraries(tile_processor PUBLIC
  falco_file_format
  HTSLIB::HTSLIB
)

add_library(falco_grade OBJECT src/falco_grade.cpp)
target_link_libraries(falco_grade PUBLIC
  falco_utils
  quality_score
)

add_library(kmer_counter OBJECT src/kmer_counter.cpp)
target_link_libraries(kmer_counter PUBLIC
  falco_grade
  falco_utils
)

add_library(contaminants OBJECT src/contaminants.cpp)

add_library(duplication_results OBJECT src/duplication_results.cpp)
target_link_libraries(duplication_results PUBLIC
  falco_grade
  falco_utils
  contaminants
)

add_library(adapter_matcher OBJECT src/adapter_matcher.cpp)
target_link_libraries(adapter_matcher PUBLIC
  falco_utils
  falco_grade
  adapter_set
)

add_library(report OBJECT src/report.cpp)
target_link_libraries(report PUBLIC
  tile_processor
  quality_score
  kmer_counter
  falco_utils
)

add_library(html OBJECT src/html.cpp)
target_link_libraries(html PUBLIC
  tile_processor
  quality_score
  kmer_counter
  falco_utils
)

add_library(results_collector INTERFACE)
target_link_libraries(results_collector INTERFACE fastq_file)

add_library(results_summary OBJECT src/results_summary.cpp)
target_link_libraries(results_summary PUBLIC
  adapter_matcher
  duplication_results
  falco_grade
  falco_utils
  html
  kmer_counter
  quality_score
  report
  tile_processor
  results_collector
)

add_library(falco_analyzer OBJECT src/falco_analyzer.cpp)
# ADS: not sure why the line below is needed on some systems, but dependencies
# in this file are a mess
target_link_libraries(falco_analyzer PUBLIC HTSLIB::HTSLIB)

add_library(get_binary_dir OBJECT src/get_binary_dir.cpp)
add_library(run_mode OBJECT src/run_mode.cpp)
add_library(falco_config OBJECT src/falco_config.cpp)
target_link_libraries(falco_config PUBLIC run_mode)

add_executable(falco src/falco.cpp)
target_link_libraries(falco PUBLIC
  tile_processor
  contaminants
  adapter_matcher
  adapter_set
  duplication_results
  kmer_counter
  fastq_file
  fastq_stdin
  fastq_gz_file
  fastq_bgzf_file
  bam_file
  sam_file
  samrec
  bam_header
  falco_file_format
  falco_utils
  falco_grade
  falco_analyzer
  results_summary
  quality_score
  falco_config
  html
  report
  get_binary_dir
  run_mode
  base_groups
  bgzf_reader
  bgzf_block
  HTSLIB::HTSLIB
  LIBDEFLATE::LIBDEFLATE
  ZLIB::ZLIB
  Threads::Threads
)

include(cmake/install.cmake)
