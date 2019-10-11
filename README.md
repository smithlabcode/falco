# Falco: FastQC Alternative Code
This program is an emulation of the popular
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc) software to
check large sequencing reads for common problems.

Full installation
============
Installation can be done by first compiling and then installing. By default, we
will look for the zlib (`-lz`) and htslib (`-lhts`) libraries in your
`LIBRARY_PATH` environment variables, as well as their sources in your
`CPLUS_INCLUDE_PATH` variable. If those are installed, compiling can be done by
running:
```
make all
```

Minimal installation
============
If you only wish to run the program on uncompressed files (FASTQ and SAM) you
can disable the zlib and hts dependencies by running:

```
make NO_HTS=1 NO_ZLIB=1
```
You can also disable only one or the other. Most unix machines have the zlib
dependency by default, so it is likely that you will not need the `NO_ZLIB=1`
flag.

In either case (full or minimum), installation can be done by running:
```
make install
```
This will create a **bin** directory with the **falco** executable inside, which
can either be added to your PATH variable or run locally.

Optional C++ dependencies
============
 * [htslib](https://github.com/samtools/htslib) is required to process bam
   files. If not provided, bam files will be treated as unrecognized file
   formats.
 * [zlib](https://zlib.net) is required to read gzipped fastq files. It is
   usually installed by default in most UNIX computers and is part of the htslib
   setup, but it can also be installed with apt or brew. If not available,
   fastq.gz files will be considered unrecognized file formats.

Running
=======

Run an example as follows:
```
falco example.fastq
```

This will generate two files :
 * ``example.fastq_qc_summary.txt`` is a text file with a summary of the QC
   metrics
 * ``example.fast_report.html`` is the visual HTML report showing plots of the
   QC metrics summarized in the text summary.

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
