#!/usr/bin/env bash
# SPDX-License-Identifier: MIT

prog=./falco
infile1=test_data/fastq_bgzip_1.fq.gz
infile2=test_data/fastq_bgzip_2.fq.gz
outdir=fastq_bgzf_mult_out
if [[ -e "${infile1}" && -e "${infile2}" ]]; then
    ${prog} -o ${outdir} ${infile1} ${infile2}
    x=$(md5sum --ignore-missing -c test_data/md5sum.txt | \
            grep "${outdir}" | \
            grep -c "OK$")
    if [[ "${x}" != "2" ]]; then
        exit 1;
    fi
    rm -r ${outdir}
else
    echo "${infile1} or ${infile2} not found; skipping remaining tests";
    exit 77;
fi
