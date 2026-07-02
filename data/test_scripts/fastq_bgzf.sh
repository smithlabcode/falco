#!/usr/bin/env bash
# SPDX-License-Identifier: MIT

# Most basic test on BGZF compressed fastq file

prog=./falco
infile=test_data/fastq_bgzip_1.fq.gz
outdir=fastq_bgzf_out
if [[ -e "${infile}" ]]; then
    ${prog} -o ${outdir} ${infile}
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
