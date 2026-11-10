#!/usr/bin/env bash
# SPDX-License-Identifier: MIT

# Test that IUPAC symbols are counted as N

# Need samtools to check that this works for sam/bam formats
if [[ ! -e $(type -P samtools) ]]; then
    echo "samtools not found; skipping remaining tests";
    exit 77;
fi

prog=./falco
infile1=test_data/iupac_dups.fq.gz
infile2=iupac_dups_bam.bam
outdir=iupac_dups_out
if [[ -e "${infile1}" ]]; then
    ${prog} --orig-dups -o ${outdir} ${infile1}
    samtools view -h -o ${outdir}/${infile2} ${infile1}
    ${prog} --orig-dups -o ${outdir} ${outdir}/${infile2}
    x=$(md5sum --ignore-missing -c test_data/md5sum.txt | \
            grep "${outdir}" | \
            grep -c "OK$")
    if [[ "${x}" != "2" ]]; then
        exit 1;
    fi
    rm -r ${outdir}
else
    echo "${infile} not found; skipping remaining tests";
    exit 77;
fi
