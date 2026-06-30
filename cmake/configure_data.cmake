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

## DATADIR is used in sources and assigned in config.h
set(DATADIR "${CMAKE_INSTALL_DATADIR}")
set(BUILD_DATA_DIR "${PROJECT_BINARY_DIR}/${DATADIR}/falco2")
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
