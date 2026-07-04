#!/usr/bin/env bash
# SPDX-License-Identifier: MIT

prog=./falco
infile=test_data/bam_1.bam
outdir=bam_out
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
