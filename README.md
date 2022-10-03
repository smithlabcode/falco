# Falco: FastQC Alternative Code

[![GitHub Downloads](https://img.shields.io/github/downloads/smithlabcode/falco/total.svg?style=social&logo=github&label=Download)](https://github.com/smithlabcode/falco/releases/latest)
[![DOI](https://zenodo.org/badge/214499063.svg)](https://zenodo.org/badge/latestdoi/214499063)
[![Install on conda](https://anaconda.org/bioconda/falco/badges/platforms.svg)](https://anaconda.org/bioconda/falco)
[![Install on conda](https://anaconda.org/bioconda/falco/badges/license.svg)](https://anaconda.org/bioconda/falco)
[![Install on conda](https://img.shields.io/conda/dn/bioconda/falco?color=red&label=conda%20downloads&style=flat-square)](https://anaconda.org/bioconda/falco)

This program is an emulation of the popular
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc)
software to check large sequencing reads for common problems.

Installing falco
================

## Installing through conda
If you use [anaconda](https://anaconda.org) to manage your packages,
and the `conda` binary is in your path, you can install the most
recent release of `falco` by running
```
$ conda install -c bioconda falco
```

`falco` can be found inside the `bin` directory of your anaconda
installer.

## Installing from source (code release)

Compilation from source can be done by downloading a `falco` release
from the [releases](https://github.com/smithlabcode/falco/releases)
section above. Upon downloading the release (in `.tar.gz` or `.zip`
format), and inflating the downloaded file to a directory
(e.g. `falco`), move to the target directory the file was inflated
(e.g. `cd falco`). You should see a `configure` file in it. In this
directory, run

```
$ ./configure CXXFLAGS="-O3 -Wall"
$ make all
$ make install
```
if you wish to install the falco binaries on a specific directory, you can use
the `--prefix` argument when running `./configure`, for instance:

```
$ ./configure CXXFLAGS="-O3 -Wall" --prefix=/path/to/installation/directory
```

The `falco` binary will be found in the `bin` directory inside the
specified prefix.

## Installing from a cloned repository

We strongly recommend using `falco` through stable releases as
described above, as the latest commits might contain undocumented
bugs. For the more advanced users who wish to test the most recent
code, `falco` can be installed by first cloning the repository

```
$ git clone https://github.com/smithlabcode/falco.git
$ cd falco
```

Once inside the generated repsotory directory, run
```
$ make all
$ make install
```

This should create a `bin` directory inside the cloned repository
containing `falco`.

### Required C++ dependencies

[zlib](https://zlib.net) is required to read gzip compressed FASTQ
files. It is usually installed by default in most UNIX computers and
is part of the htslib setup, but it can also be installed with
standard package managers like apt, brew or conda.

On Ubuntu, zlib C++ libraries can be installed with `apt`:
```
$ sudo apt install zlib1g zlib1g-dev
```

### Optional C++ dependencies

[htslib](https://github.com/samtools/htslib) is required to process
bam files. If not provided, bam files will be treated as unrecognized
file formats.

If htslib is installed, falco can be compiled with it by simply replacing the
configure command above with the `--enable-hts` flag:

```
$ ./configure CXXFLAGS="-O3 -Wall" --enable-hts
```

If `falco` was cloned from the repository, run the following commands
to allow BAM file reading:

```
$ make HAVE_HTSLIB=1 all
$ make HAVE_HTSLIB=1 install
```

If successfully compiled, `falco` can be used in BAM files the same way as it is
used with fastq and sam files.

Running falco
=============

Run falco in with the following command, where the `example.fq` file
provided can be replaced with the path to any FASTQ file you want to run
`falco`
```
$ falco example.fq
```

This will generate three files in the same directory as the input fastq file:

* `fastqc_data.txt` is a text file with a summary of the QC metrics

* `fastqc_report.html` is the visual HTML report showing plots of the
   QC metrics summarized in the text summary.

* `summary.txt`: A tab-separated file describing whether the
  pass/warn/fail result for each module. If multiple files are
  provided, only one summary file is generated, with one of the
  columns being the file name associated to each module result.

The full list of arguments and options can be seen by running `falco`
without any arguments, as well as `falco -?` or `falco --help`. This
will print the following list:

```
Usage: falco [OPTIONS] <seqfile1> <seqfile2> ...
Options:
  -h, --help               Print this help file and exit  
  -v, --version            Print the version of the program and exit  
  -o, --outdir             Create all output files in the specified 
                           output directory. FALCO-SPECIFIC: If the 
                           directory does not exists, the program will 
                           create it. If this option is not set then 
                           the output file for each sequence file is 
                           created in the same directory as the 
                           sequence file which was processed.  
      --casava             [IGNORED BY FALCO] Files come from raw 
                           casava output. Files in the same sample 
                           group (differing only by the group number) 
                           will be analysed as a set rather than 
                           individually. Sequences with the filter flag 
                           set in the header will be excluded from the 
                           analysis. Files must have the same names 
                           given to them by casava (including being 
                           gzipped and ending with .gz) otherwise they 
                           won't be grouped together correctly.  
      --nano               [IGNORED BY FALCO] Files come from nanopore 
                           sequences and are in fast5 format. In this 
                           mode you can pass in directories to process 
                           and the program will take in all fast5 files 
                           within those directories and produce a 
                           single output file from the sequences found 
                           in all files.  
      --nofilter           [IGNORED BY FALCO] If running with --casava 
                           then don't remove read flagged by casava as 
                           poor quality when performing the QC 
                           analysis.  
      --extract            [ALWAYS ON IN FALCO] If set then the zipped 
                           output file will be uncompressed in the same 
                           directory after it has been created. By 
                           default this option will be set if fastqc is 
                           run in non-interactive mode.  
  -j, --java               [IGNORED BY FALCO] Provides the full path to 
                           the java binary you want to use to launch 
                           fastqc. If not supplied then java is assumed 
                           to be in your path.  
      --noextract          [IGNORED BY FALCO] Do not uncompress the 
                           output file after creating it. You should 
                           set this option if you do not wish to 
                           uncompress the output when running in 
                           non-interactive mode.  
      --nogroup            Disable grouping of bases for reads >50bp. 
                           All reports will show data for every base in 
                           the read. WARNING: When using this option, 
                           your plots may end up a ridiculous size. You 
                           have been warned!  
      --min_length         [NOT YET IMPLEMENTED IN FALCO] Sets an 
                           artificial lower limit on the length of the 
                           sequence to be shown in the report. As long 
                           as you set this to a value greater or equal 
                           to your longest read length then this will 
                           be the sequence length used to create your 
                           read groups. This can be useful for making 
                           directly comaparable statistics from 
                           datasets with somewhat variable read 
                           lengths.  
  -f, --format             Bypasses the normal sequence file format 
                           detection and forces the program to use the 
                           specified format. Valid formats are bam, sam, 
                           bam_mapped, sam_mapped, fastq, fq, fastq.gz 
                           or fq.gz.  
  -t, --threads            [NOT YET IMPLEMENTED IN FALCO] Specifies the 
                           number of files which can be processed 
                           simultaneously. Each thread will be 
                           allocated 250MB of memory so you shouldn't 
                           run more threads than your available memory 
                           will cope with, and not more than 6 threads 
                           on a 32 bit machine [1] 
  -c, --contaminants       Specifies a non-default file which contains 
                           the list of contaminants to screen 
                           overrepresented sequences against. The file 
                           must contain sets of named contaminants in 
                           the form name[tab]sequence. Lines prefixed 
                           with a hash will be ignored. Default: 
                           Configuration/contaminant_list.txt 
  -a, --adapters           Specifies a non-default file which contains 
                           the list of adapter sequences which will be 
                           explicity searched against the library. The 
                           file must contain sets of named adapters in 
                           the form name[tab]sequence. Lines prefixed 
                           with a hash will be ignored. Default: 
                           Configuration/adapter_list.txt 
  -l, --limits             Specifies a non-default file which contains 
                           a set of criteria which will be used to 
                           determine the warn/error limits for the 
                           various modules. This file can also be used 
                           to selectively remove some modules from the 
                           output all together. The format needs to 
                           mirror the default limits.txt file found in 
                           the Configuration folder. Default: 
                           Configuration/limits.txt 
  -k, --kmers              [IGNORED BY FALCO AND ALWAYS SET TO 7] 
                           Specifies the length of Kmer to look for in 
                           the Kmer content module. Specified Kmer 
                           length must be between 2 and 10. Default 
                           length is 7 if not specified.  
  -q, --quiet              Supress all progress messages on stdout and 
                           only report errors.  
  -d, --dir                [IGNORED: FALCO DOES NOT CREATE TMP FILES] 
                           Selects a directory to be used for temporary 
                           files written when generating report images. 
                           Defaults to system temp directory if not 
                           specified.  
  -s, -subsample           [Falco only] makes falco faster (but 
                           possibly less accurate) by only processing 
                           reads that are multiple of this value (using 
                           0-based indexing to number reads). [1] 
  -b, -bisulfite           [Falco only] reads are whole genome 
                           bisulfite sequencing, and more Ts and fewer 
                           Cs are therefore expected and will be 
                           accounted for in base content.  
  -r, -reverse-complement  [Falco only] The input is a 
                           reverse-complement. All modules will be 
                           tested by swapping A/T and C/G  
      -skip-data           [Falco only] Do not create FastQC data text 
                           file.  
      -skip-report         [Falco only] Do not create FastQC report 
                           HTML file.  
      -skip-summary        [Falco only] Do not create FastQC summary 
                           file  
  -D, -data-filename       [Falco only] Specify filename for FastQC 
                           data output (TXT). If not specified, it will 
                           be called fastq_data.txt in either the input 
                           file's directory or the one specified in the 
                           --output flag. Only available when running 
                           falco with a single input.  
  -R, -report-filename     [Falco only] Specify filename for FastQC 
                           report output (HTML). If not specified, it 
                           will be called fastq_report.html in either 
                           the input file's directory or the one 
                           specified in the --output flag. Only 
                           available when running falco with a single 
                           input.  
  -S, -summary-filename    [Falco only] Specify filename for the short 
                           summary output (TXT). If not specified, it 
                           will be called fastq_report.html in either 
                           the input file's directory or the one 
                           specified in the --output flag. Only 
                           available when running falco with a single 
                           input.  
  -K, -add-call            [Falco only] add the command call call to 
                           FastQC data output and FastQC report HTML 
                           (this may break the parse of fastqc_data.txt 
                           in programs that are very strict about the 
                           FastQC output format).  

Help options:
  -?, -help                print this help message  
      -about               print about message  

PROGRAM: falco
A high throughput sequence QC analysis tool
```

Citing falco
============

If `falco` was helpful for your research, you can cite us as follows:

*de Sena Brandine G and Smith AD. Falco: high-speed FastQC emulation for
quality control of sequencing data. F1000Research 2021, 8:1874
(https://doi.org/10.12688/f1000research.21142.2)*

**Please do not cite this manuscript if you used FastQC directly and not falco!**

Copyright and License Information
=================================

Copyright (C) 2019-2022 Guilherme de Sena Brandine and
                        Andrew D. Smith

Authors: Guilherme de Sena Brandine and Andrew D. Smith

This is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your
option) any later version.

This software is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.
