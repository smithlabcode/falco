#!/usr/bin/env bash
# SPDX-License-Identifier: MIT

# Test k-mer analysis with BAM input and one thread

prog=./falco
infile1=test_data/bam_1.bam
infile2=test_data/fastq_bgzip_1.fq.gz
infile3=test_data/fastq_bgzip_2.fq.gz
outdir=kmers_out
if [[ -e "${infile1}" ||-e "${infile2}" || -e "${infile3}" ]]; then
    ${prog} --kmers -o ${outdir} ${infile1} ${infile2} ${infile3}
    x=$(md5sum --ignore-missing -c test_data/md5sum.txt | \
            grep "${outdir}" | \
            grep -c "OK$")
    if [[ "${x}" != "3" ]]; then
        exit 1;
    fi
    rm -r ${outdir}
else
    echo "input files for kmer tests not found; skipping remaining tests";
    exit 77;
fi
