#!/usr/bin/env bash
# SPDX-License-Identifier: MIT

# Basic test on BGZF compressed fastq file using more than one thread

prog=./falco
infile=test_data/fastq_bgzip_1.fq.gz
outdir=fastq_bgzf_threads_out
if [[ -e "${infile}" ]]; then
    ${prog} -t 4 -o ${outdir} ${infile}
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
