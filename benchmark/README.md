Benchmarking tools
===================

This directory shows the script and one example output of how benchmarking was
performed. the shell files in this directory can be run given an input
directory, and will run the three compared tools, redirecting the `time` output
to the `outs` directory.

Downloading benchmark fastqs
============================

FASTQ files can be downloaded using the `fastq-dump` program in the [SRA
Toolkit](https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft) (alternatively
available on [conda](https://anaconda.org/bioconda/sra-tools)).

With an active internet connection, the `fastq-dump` command in your `PATH`
variable and at least 300GB of disk space, run the following to download the
fastq files onto the test directory:
```
$ bash download_files.sh
```
QC software download links
==========================
The following tools need to be installed to run the benchmarking:
 * [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc)
 * [fastp](https://github.com/OpenGene/fastp/releases)
 * [HTQC](https://sourceforge.net/projects/htqc)

Command to run benchmarking
===========================
Once files are downloaded and programs are installed and added to your local
`PATH` variable, you can reproduce the benchmarking in the paper by running the
following three commands:
```
$ bash run_all_falco_tests.sh
$ bash run_all_fastp_tests.sh
$ bash run_all_fastqc_tests.sh
$ bash run_all_htqc_tests.sh
```

This will output the real, user and sys runtimes for each tool in each dataset.

List of SRR accessions tested
=============================
 * SRR10124060
 * SRR10143153
 * SRR3897196
 * SRR9624732
 * SRR1853178
 * SRR6387347
 * SRR891268
 * SRR1772703
 * SRR9878537
 * SRR6059706

The human genome nanopore file can be downloaded from the [Human Whole Genome
Sequencing
Project](https://github.com/nanopore-wgs-consortium/NA12878/blob/master/nanopore-human-genome/rel_3_4.md)
(file
[FAB49164](http://s3.amazonaws.com/nanopore-human-wgs/rel3-nanopore-wgs-4045668814-FAB49164.fastq.gz))
and extracted into the tests/fastq directory
