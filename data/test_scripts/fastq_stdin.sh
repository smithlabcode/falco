#!/usr/bin/env bash
# SPDX-License-Identifier: MIT

# Test for fastq input from stdin

prog=./falco
infile1=test_data/fastq_bgzip_1.fq.gz
name=fastq_stdin
outdir=fastq_stdin_out
if [[ -e "${infile1}" ]]; then
    mkdir -p ${outdir}
    gunzip -c ${infile1} | \
    ${prog} -o ${outdir} --stdin fq ${name}
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
