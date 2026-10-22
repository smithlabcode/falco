#!/usr/bin/env bash
# SPDX-License-Identifier: MIT

prog=./falco
infile1=test_data/fastq_bgzip_1.fq.gz
infile2=fastq_1_iupac.fq
outdir=fastq_iupac_out
if [[ -e "${infile1}" ]]; then
    mkdir -p ${outdir}
    # ADS: on the first sequence line, modify all C to K
    gunzip -c ${infile1}| sed -e "2s/C/K/g" > ${outdir}/${infile2}
    ${prog} --iupac -o ${outdir} ${outdir}/${infile2}
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
