# Benchmarking tools
This directory shows the script and one example output of how benchmarking was
performed. the `run_benchmarking.sh` shell file can be run given an input
directory, and will run the three compared tools, redirecting the `time` output
to the `outs` directory.

Downloading benchmark fastqs
=============================
Fastq files can be downloaded using the `fastq-dump` program in the [SRA
Toolkit](https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft) (alternatively
available on [conda](https://anaconda.org/bioconda/sra-tools)).

List of SRR accessions tested
=============================
SRR10124060
SRR10143153
SRR3897196
SRR9624732
SRR1853178
SRR6387347
SRR891268
SRR1772703
SRR9878537
SRR6059706

Command to download all fastq files
===================================
To populate the `tests` directory, just run the command below. Note that you
need at least 300GB of disk space:
```
files="SRR10124060 SRR10143153 SRR3897196 SRR9624732 SRR1853178 SRR6387347
SRR891268 SRR1772703 SRR9878537"
for i in $files
do
  echo "Downloading ${i}..."
  fastq-dump --outdir  --gzip --skip-technical --readids \
               --read-filter pass --dumpbase --split-3 --clip \
               ${i}
done
```

Command to run benchmarking
===========================
Once files are downloaded you can reproduce the benchmarking in the paper by
running `./run_benchmark.sh`.
