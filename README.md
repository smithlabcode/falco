# Falco: FastQC Alternative Code
This program is an emulation of the popular
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc) software to
check large sequencing reads for common problems.

Installing falco
================

## Installing through conda
If you use [anaconda](https://anaconda.org) to manage your packages, you can
install the most recent release of `falco` by running:
```
conda install -c bioconda falco
```

## Installing from source (code release)

Compilation from source can be done by downloading a `falco` release from the
[releases](https://github.com/smithlabcode/falco/releases)
section above. Upon downloading, inflating and moving to the source
directory, installation can be done through the following commands:

```
$ ./configure CXXFLAGS="-O3 -Wall"
$ make all
$ make install
```
if you wish to install the falco binaries on a specific directory, you can use the `--prefix` argument when running `./configure`, for instance:

```
$ ./configure CXXFLAGS="-O3 -Wall" --prefix=/path/to/installation/directory
```

## Installing from a cloned repository
We strongly recommend using `falco` through stable releases as described above,
as the latest commits might contain undocumented bugs. However, if users wish
to test the most recent code, the
[GNU autotools](https://www.gnu.org/software/automake) software is required
for compilation. Simply run `autoreconf -i` to generate the `configure` file,
then run the commands listed in the previous section to check for dependencies
and compile the program.

### Required C++ dependencies

[zlib](https://zlib.net) is required to read gzip compressed FASTQ files. It is
usually installed by default in most UNIX computers and is part of the htslib
setup, but it can also be installed with standard package managers like
apt, brew or conda.

On Ubuntu, zlib C++ libraries can be installed with `apt`:
```
$ sudo apt install zlib1g zlib1g-dev
```

### Optional C++ dependencies

[htslib](https://github.com/samtools/htslib) is required to process bam
files. If not provided, bam files will be treated as unrecognized file
formats.

If htslib is installed, falco can be compiled with it by simply replacing the
configure command above with the `--enable-hts` flag:

```
$ ./configure CXXFLAGS="-O3 -Wall" --enable-hts
```

If successfully compiled, `falco` can be used in bam files the same way as it is
used with fastq and sam files.

Running falco
=============

Run falco in a hypothetical fastq file `example.fastq` with the following
command:
```
$ falco example.fastq
```

This will generate three files in the same directory as the input fastq file:
 * ``fastqc_data.txt`` is a text file with a summary of the QC
   metrics
 * ``fastqc_report.html`` is the visual HTML report showing plots of the
   QC metrics summarized in the text summary.
* ``summary.txt``: A tab-separated file describing whether the pass/warn/fail
  result for each module. If multiple files are provided, only one summary file
  is generated, with one of the columns being the file name associated to each
  module result.

the full list of arguments and options can be seen by running `falco` without any arguments, as well as `falco -?` or `falco --help`. This will print the following list:

```
Usage: ./bin/falco [OPTIONS] <seqfile1> <seqfile2> ...

Options:
  -h, --help                print this help file and exit
  -v, --version             print the program version and exit
  -o, --outdir              Create all output files in the specified
                            output directory. If notprovided, files
                            will be created in the same directory as
                            the input file.
  -C, --casava              Files come from raw casava output
                            (currently ignored)
  -n, --nano                Files come from fast5 nanopore sequences
  -F, --nofilter            If running with --casava do not sequences
                            (currently ignored)
  -e, --noextract           If running with --casava do not remove poor
                            quality sequences (currently ignored)
  -g, --nogroup             Disable grouping of bases for reads >50bp
  -f, --format              Force file format
  -t, --threads             Specifies number of threads to process
                            simultaneos files in parallel (currently
                            set for compatibility with fastqc. Not yet
                            supported!)
  -c, --contaminants        Non-default filer with a list of
                            contaminants
  -a, --adapters            Non-default file with a list of adapters
  -l, --limits              Non-default file with limits and warn/fail
                            criteria
  -T, --skip-text           Skip generating text file (Default = false)
  -H, --skip-html           Skip generating HTML file (Default = false)
  -S, --skip-short-summary  Skip short summary(Default = false)
  -q, --quiet               print more run info
  -d, --dir                 directory in which to create temp files
  -A, --advanced-mode       advanced mode: adds more information to the
                            FastQC output depending on non-fastqc user
                            flags
  -B, --bisulfite           reads are whole genome bisulfite
                            sequencing, and more Ts and fewer Cs are
                            therefore expected and will be accounted
                            for in base content (advanced mode)
  -R, --reverse-complement  The input is a reverse-complement. All
                            modules will be tested by swapping A/T and
                            C/G

Help options:
  -?, -help                 print this help message
      -about                print about message

PROGRAM: ./bin/falco
A high throughput sequence QC analysis tool
```
Copyright and License Information
=================================

Copyright (C) 2019
University of Southern California,

Current Authors: Guilherme de Sena Brandine, Andrew D. Smith

This is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your
option) any later version.

This software is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.
