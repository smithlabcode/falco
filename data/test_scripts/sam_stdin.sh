#!/usr/bin/env bash
# SPDX-License-Identifier: MIT

# Test for SAM input from stdin

prog=./falco
infile1=test_data/bam_1.bam
name=sam_stdin
outdir=sam_stdin_out
if [[ -e "${infile1}" ]]; then
    mkdir -p ${outdir}
    samtools view -h ${infile1} | \
    ${prog} -o ${outdir} --stdin sam ${name}
    x=$(md5sum --ignore-missing -c test_data/md5sum.txt | \
            grep "${outdir}" | \
            grep -c "OK$")
    if [[ "${x}" != "1" ]]; then
        exit 1;
    fi
    rm -r ${outdir}
else
    echo "${infile} not found; skipping remaining tests";
    exit 77;
fi
