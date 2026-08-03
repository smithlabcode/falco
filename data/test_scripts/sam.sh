#!/usr/bin/env bash
# SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

# Most basic test on regular SAM file

if [[ ! -e $(type -P samtools) ]]; then
    echo "samtools not found; skipping remaining tests";
    exit 77;
fi

prog=./falco
infile1=test_data/bam_1.bam
infile2=sam_1.sam
outdir=sam_out
if [[ -e "${infile1}" ]]; then
    mkdir -p ${outdir}
    samtools view -h -o ${infile2} ${infile1}
    ${prog} -o ${outdir} ${infile2}
    x=$(md5sum --ignore-missing -c test_data/md5sum.txt | \
            grep "${outdir}" | \
            grep -c "OK$")
    if [[ "${x}" != "1" ]]; then
        exit 1;
    fi
    rm ${infile2}
    rm -r ${outdir}
else
    echo "${infile1} not found; skipping remaining tests";
    exit 77;
fi
