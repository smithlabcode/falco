# SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

## DATADIR is used in sources and assigned in config.h
set(DATADIR "${CMAKE_INSTALL_DATADIR}")
set(BUILD_DATA_DIR "${PROJECT_BINARY_DIR}/${DATADIR}/falco")
file(MAKE_DIRECTORY ${BUILD_DATA_DIR})
set(FastQC_DATA_DIR "${PROJECT_SOURCE_DIR}/data/Configuration")
set(
  FastQC_DATA_FILES
  limits.txt
  adapter_list.txt
  contaminant_list.txt
  LICENSE
  README.md
)

foreach(FastQC_FILE_NAME ${FastQC_DATA_FILES})
  configure_file(
    ${FastQC_DATA_DIR}/${FastQC_FILE_NAME}
    ${BUILD_DATA_DIR}/${FastQC_FILE_NAME}
    COPYONLY
  )
endforeach()
