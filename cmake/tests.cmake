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

# Add test scripts and data

include(CTest)

## Ensure test data is available in the build dir
file(CREATE_LINK
  ${PROJECT_SOURCE_DIR}/data/test_data
  ${PROJECT_BINARY_DIR}/test_data
  SYMBOLIC
)

## Ensure test scripts are available in the build dir
file(CREATE_LINK
  ${PROJECT_SOURCE_DIR}/data/test_scripts
  ${PROJECT_BINARY_DIR}/test_scripts
  SYMBOLIC
)

find_program(BASH_PROGRAM bash)

## Add each test
add_test(NAME "BGZF input" COMMAND bash test_scripts/fastq_bgzf.sh)
add_test(NAME "BGZF multiple" COMMAND bash test_scripts/fastq_bgzf_mult.sh)
add_test(NAME "BAM input" COMMAND bash test_scripts/bam.sh)
add_test(NAME "BAM multiple" COMMAND bash test_scripts/bam_mult.sh)
add_test(NAME "GZ input" COMMAND bash test_scripts/fastq_gz.sh)
add_test(NAME "BGZF threads" COMMAND bash test_scripts/fastq_bgzf_threads.sh)
add_test(NAME "Plain FASTQ" COMMAND bash test_scripts/fastq_plain.sh)
add_test(NAME "Tiles" COMMAND bash test_scripts/tiles.sh)
add_test(NAME "Kmers" COMMAND bash test_scripts/kmers.sh)
