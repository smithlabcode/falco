# Falco: FastQC Alternative Code
This program is an emulation of the popular
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc) software to
check large sequencing reads for common problems.

Installation
============
Installation can be done with the standard autotools commands. By default, we
will look for the zlib (`-lz`) library in your
`LIBRARY_PATH` environment variables,as well as their sources in your
`CPLUS_INCLUDE_PATH` variable. If those are installed, compiling can be done by
running the following commands, where `$` is your shell:
```
$ ./configure
$ make all
$ make install
```

Required C++ dependencies
============
[zlib](https://zlib.net) is required to read gzipped fastq files. It is
usually installed by default in most UNIX computers and is part of the htslib
setup, but it can also be installed with apt or brew. If not available,
fastq.gz files will be considered unrecognized file formats.

On Ubuntu, zlib C++ libraries can be installed with `apt`:
```
$ sudo apt install zlib1g zlib1g-dev
```


Optional C++ dependencies
============
[htslib](https://github.com/samtools/htslib) is required to process bam
files. If not provided, bam files will be treated as unrecognized file
formats.

If htslib is installed, `falco` can be compiled with it by simply replacing the
configure command above with the `--enable-hts` flag:

```
$ ./configure --enable-hts
```
If successfully compiled, `falco` can be used in bam files the same way as it is
used with fastq and sam files.

Running
=======

Run `falco` in a hypothetical fastq file `example.fastq` with the following
command:
```
$ falco example.fastq
```

This will generate three files in the same directory as the input fastq file:
 * ``example.fastq_qc_summary.txt`` is a text file with a summary of the QC
   metrics
 * ``example.fastq_report.html`` is the visual HTML report showing plots of the
   QC metrics summarized in the text summary.
* ``summary.txt``: A tab-separated file describing whether the pass/warn/fail
  result for each module. If multiple files are provided, only one summary file
  is generated, with one of the columns being the file name associated to each
  module result.

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
