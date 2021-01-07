Benchmarking tools
===================

This directory shows the script and one example output of how benchmarking was
performed. the shell files in this directory can be run given an input
directory, and will run the three compared tools, redirecting the `time` output
to the `outs` directory.

### Downloading benchmark fastqs

FASTQ files can be downloaded using the `fastq-dump` program in the [SRA
Toolkit](https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft) (alternatively
available on [conda](https://anaconda.org/bioconda/sra-tools)).

With an active internet connection, the `fastq-dump` command in your `PATH`
variable and at least 300GB of disk space, run the following to download the
fastq files onto the test directory:
```
$ bash download_files.sh
```
### QC software download links
The following tools need to be installed to run the benchmarking:
 * [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc)
 * [fastp](https://github.com/OpenGene/fastp/releases)
 * [HTQC](https://sourceforge.net/projects/htqc)

### Command to run benchmarking
Once files are downloaded and programs are installed and added to your local
`PATH` variable, you can reproduce the benchmarking in the paper by running the
following three commands:
```
$ ./run_all_falco_tests.sh
$ ./run_all_fastp_tests.sh
$ ./run_all_fastqc_tests.sh
$ ./run_all_htqc_tests.sh
```

This will output the real, user and sys runtimes for each tool in each dataset.

#### List of SRR accessions tested
 The URLS below link to the `.sra` file that can be then converted to
 FASTQ using the `fastq-dump` command. Details about each dataset can
 be found at https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=[srr],
 where [srr] can be replaced by each accession number (e.g.
`?run=SRR1853178`)
 * [SRR10124060](https://sra-download.ncbi.nlm.nih.gov/traces/sra4/SRR/009886/SRR10124060)
 * [SRR10143153](https://sra-download.ncbi.nlm.nih.gov/traces/sra68/SRR/009905/SRR10143153)
 * [SRR3897196](https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-9/SRR3897196/SRR3897196.1)
 * [SRR9624732](https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-15/SRR9624732/SRR9624732.1)
 * [SRR1853178](https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1853178/SRR1853178.1)
 * [SRR6387347](https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR6387347/SRR6387347.1)
 * [SRR891268](https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR891268/SRR891268.1)
 * [SRR1772703](https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-2/SRR1772703/SRR1772703.1)
 * [SRR9878537](https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-15/SRR9878537/SRR9878537.1)
 * [SRR6059706](https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/SRR6059706/SRR6059706.1)

The human genome nanopore file can be downloaded from the [Human Whole
Genome Sequencing
Project](https://github.com/nanopore-wgs-consortium/NA12878/blob/master/nanopore-human-genome/rel_3_4.md)
(file
[FAB49164](http://s3.amazonaws.com/nanopore-human-wgs/rel3-nanopore-wgs-4045668814-FAB49164.fastq.gz))
and extracted into the tests/fastq directory
