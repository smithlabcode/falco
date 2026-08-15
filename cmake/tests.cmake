# SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

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
add_test(NAME "Groups" COMMAND bash test_scripts/groups.sh)
add_test(NAME "SAM input" COMMAND bash test_scripts/sam.sh)
add_test(NAME "Orig dups" COMMAND bash test_scripts/orig_dups.sh)
add_test(NAME "Preseq output" COMMAND bash test_scripts/preseq.sh)
