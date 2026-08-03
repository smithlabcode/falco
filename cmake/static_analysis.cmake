# SPDX-License-Identifier: MIT; (c) 2026 Andrew D Smith

# StaticAnalysis
message(STATUS "Enabling static analysis")
# If no specific static analysis is requested, do them all
if(NOT RUN_CPPCHECK AND NOT RUN_IWYU AND NOT RUN_CPPLINT AND NOT RUN_CLANG_TIDY)
  set(RUN_CPPCHECK on)
  set(RUN_IWYU on)
  set(RUN_CPPLINT on)
  set(RUN_CLANG_TIDY on)
endif()

set(STATIC_ANALYSIS_CHECKS "")
if(RUN_CPPCHECK)
  list(APPEND STATIC_ANALYSIS_CHECKS "cppcheck")
endif()
if(RUN_CPPLINT)
  list(APPEND STATIC_ANALYSIS_CHECKS "cpplint")
endif()
if(RUN_IWYU)
  list(APPEND STATIC_ANALYSIS_CHECKS "iwyu")
endif()
if(RUN_CLANG_TIDY)
  list(APPEND STATIC_ANALYSIS_CHECKS "clang-tidy")
endif()

message(STATUS "Requested static analysis: ${STATIC_ANALYSIS_CHECKS}")

# cpplint: all options are in the config file
if ("cpplint" IN_LIST STATIC_ANALYSIS_CHECKS)
  find_program(FOUND_CPPLINT cpplint)
  if(FOUND_CPPLINT)
    message(STATUS "Enabling cpplint analysis")
    set(CMAKE_CXX_CPPLINT cpplint --quiet)
  else()
    message(STATUS "Could not find cpplint; disabling cpplint")
  endif()
endif()

# include-what-you-use: config is a mappings file
if ("iwyu" IN_LIST STATIC_ANALYSIS_CHECKS)
  find_program(FOUND_IWYU include-what-you-use)
  if(FOUND_IWYU)
    message(STATUS "Enabling include-what-you-use analysis")
    set(CMAKE_CXX_INCLUDE_WHAT_YOU_USE
      include-what-you-use
      -Xiwyu
      --comment_style=long
      -Xiwyu
      --quoted_includes_first
      -Xiwyu
      --mapping_file=${PROJECT_SOURCE_DIR}/iwyu.json
    )
  else()
    message(STATUS "Could not find iwyu; disabling iwyu")
  endif()
endif()

# cppcheck: options on the command line as there is no config file
if ("cppcheck" IN_LIST STATIC_ANALYSIS_CHECKS)
  find_program(FOUND_CPPCHECK cppcheck)
  if(FOUND_CPPCHECK)
    message(STATUS "Enabling cppcheck analysis")
    set(CMAKE_CXX_CPPCHECK
      cppcheck
      --quiet
      --enable=all
      --inline-suppr
      --max-configs=1
      --suppressions-list=${PROJECT_SOURCE_DIR}/.cppcheck_suppress
    )
  else()
    message(STATUS "Could not find cppcheck; disabling cppcheck")
  endif()
endif()

# clang-tidy: need to make sure version is at least 20
if ("clang-tidy" IN_LIST STATIC_ANALYSIS_CHECKS)
  find_program(CLANG_TIDY_EXECUTABLE NAMES clang-tidy)
  # Minimum required version
  set(MIN_CLANG_TIDY_VERSION "20.0.0")
  if(CLANG_TIDY_EXECUTABLE)
    execute_process(
      COMMAND
      bash -c
      "${CLANG_TIDY_EXECUTABLE} --version | grep version | tr -cd '0-9.\n'"
      OUTPUT_VARIABLE CLANG_TIDY_VERSION
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    # Compare the version numbers
    if(CLANG_TIDY_VERSION VERSION_GREATER_EQUAL MIN_CLANG_TIDY_VERSION)
      message(STATUS "Enabling clang-tidy (version: ${CLANG_TIDY_VERSION})")
      set(CMAKE_CXX_CLANG_TIDY
        clang-tidy
        --quiet
        --allow-no-checks
        -p ${PROJECT_BINARY_DIR}
      )
    else()
      message(STATUS "Not enabling clang-tidy (min version not found")
    endif()
  else()
    message(STATUS "Could not find clang-tidy; disabling clang-tidy")
  endif()
endif()
