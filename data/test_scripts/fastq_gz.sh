#!/usr/bin/env bash
# SPDX-License-Identifier: MIT

# Most basic test on regular gzip fastq file

prog=./falco
infile1=test_data/fastq_bgzip_1.fq.gz
infile2=fastq_gz_1.fq.gz
outdir=fastq_gz_out
if [[ -e "${infile1}" ]]; then
    mkdir -p ${outdir}
    gunzip -c ${infile1} | gzip > ${outdir}/${infile2}
    ${prog} -o ${outdir} ${outdir}/${infile2}
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
