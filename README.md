[![GitHub Downloads](https://img.shields.io/github/downloads/smithlabcode/falco/total.svg?style=social&logo=github&label=Download)](https://github.com/smithlabcode/falco/releases/latest)
[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.214499063-blue)](https://zenodo.org/badge/latestdoi/214499063)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Install on conda](https://anaconda.org/bioconda/falco/badges/platforms.svg)](https://anaconda.org/bioconda/falco)
[![Install on conda](https://img.shields.io/conda/dn/bioconda/falco?color=red&label=conda%20downloads&style=flat-square)](https://anaconda.org/bioconda/falco)

# Falco

Falco v2.0 has arrived. If you want to see some recent benchmarks you can find
them [here](https://github.com/smithlabcode/falco/discussions).

Falco was conceived as an emulation of the popular
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc) software to
check large sequencing reads for common problems. Falco was rewritten for
version 2.0 in order to facilitate incoroprating new functionality moving
forward.

## Quick start

You will be able to find binaries for Linux and macOS with the releases.

Example:
```
falco -o output input.fq
```
This generates 3 files in the direcotry named `output` (creating it if needed):
* `fastqc_data.txt` is a text file with a summary of the QC metrics.
* `fastqc_report.html` is the visual HTML report showing plots of the QC metrics
  summarized in the text summary.
* `summary.txt`: A tab-separated file describing whether the pass/warn/fail
  result for each analysis done.

## Installing through conda

If you use [anaconda](https://anaconda.org) to manage your packages,
and the `conda` binary is in your path, you can install the most
recent release of `falco` by running
```
conda install -c bioconda falco
```

## Building

- Falco has moved to using c++23 (GCC >= 14.2.0 or LLVM-Clang >= 20.0.0; on
  macOS GCC >= 15).
- We moved from autotools to cmake.
- Dependencies:
  * [HTSLib](https://github.com/samtools/htslib): used for identifying file formats.
  * [Zlib](https://github.com/madler/zlib): used by HTSLib for regular gzip files.
  * [libdeflate](https://github.com/ebiggers/libdeflate): also used by HTSLib for BGZF files.
  * [ISA-L](https://github.com/intel/isa-l) (highly recommended; not required):
    speeds up parsing regular gzip files.
  * Other dependencies are in the source and listed with `falco --licenses`
- The "workflows" have exactly working steps for builds but are likely more than
  you need: `.github/workflows/linux-release.yml` and `.github/workflows/macos-release.yml`

From the root of the repo if you have all dependencies:
```
cmake -B build -DCMAKE_CXX_COMPILER=g++ -DCMAKE_BUILD_TYPE=Release
cmake --build build -j8  # use -j for more cores when building
```
Then `falco` will be in the `build` directory. Please do not build with
`-DCMAKE_BUILD_TYPE=Build` because I wrote the source with the intention of
allowing the compiler do most of the optimizations and without using `Release`
falco will become slow. You can use a different compiler for the
`-DCMAKE_CXX_COMPILER` option, but I've found a compiler often needs to be
specified directly.

### Linux detailed instructions

I'm explaining this via a clean Ubuntu instance in docker:
```
docker pull ubuntu:latest
docker run -it ubuntu:latest bash
```
Inside the docker:
```
export DEBIAN_FRONTEND=noninteractive &&
apt-get update &&
apt-get install -y --no-install-recommends \
    libssl-dev \
    zlib1g-dev \
    libdeflate-dev \
    libisal-dev \
    libhts-dev \
    ca-certificates \
    git \
    g++-15 \
    cmake \
    make \
    samtools && # samtools is for running tests
git clone https://github.com/smithlabcode/falco.git &&
cd falco &&
cmake -B build -DUSE_ISAL=on -DCMAKE_CXX_COMPILER=g++-15 -DCMAKE_BUILD_TYPE=Release &&
cmake --build build -j8 &&  # use -j for more cores when building
ctest --test-dir build
```
If you want instructions that also include building the dependencies from source,
you can find them here: `.github/workflows/linux-release.yml`

### macOS detailed instructions

I don't have the same ability to test with clean OS images for macOS
(suggestions welcome). The best I can do is use the GitHub macOS runners, which
already have some of the dependencies installed. Here is what works:
```
brew install libdeflate isa-l htslib samtools &&  # samtools for testing
git clone https://github.com/smithlabcode/falco.git &&
cd falco &&
cmake -B build -DUSE_ISAL=on -DCMAKE_CXX_COMPILER=g++-15 -DCMAKE_BUILD_TYPE=Release &&
cmake --build build -j8 &&
ctest --test-dir build
```
ZLib is already installed on macOS, HTSLib installs libdeflate as a dependency
and samtools installs both as dependency. To see what's already installed on
GitHub's macOS look
[here](https://github.com/actions/runner-images/blob/main/images/macos/macos-26-Readme.md).

## Changes in Falco v2.0

### Tiles results

I found that the method for tile analysis is a bit unstable, and the tile grade
can be slightly unstable. The only way to notice this is to process reads from
the same input file in different orders. This happens as a side effect of
analyzing reads concurrently with threads. Here is my understanding of how
FastQC works, and how I implemented falco v2.0. Please comment if you see
anything incorrect.

- Tile analysis is done for 1/10 of the reads (though FastQC includes all among
  the first 10k reads).
- Accumulating results: For each counted read, for each position in the read,
  the quality score contributes to that tile's mean for the given position.
- Summarizing tile results: For each read position, the mean over tiles' quality
  scores is taken. Then for each tile, for each read position, the value is
  centered by subtracting the mean (the 'centered' tile values).
- The summary stat for deciding the grade is based on the minimum value among
  all centered values, across all tiles and across all positions.

Based on the assumptions above, the grade uses an extreme value statistic. When
introducing multithreading to falco, the order of reads analyzed changes between
runs, so the 1/10 reads contributing to the tile analysis also changes between
runs. I've noticed that this can lead to differences between runs, and in some
cases this has changed the grade between pass/warn and warn/fail. So it is
possible the grade can differ between runs for the same data.

### Duplication results

I changed how falco evaluates "duplcation". Although the format of the output is
the same, the numbers differ dramatically. The motivation for the change is to
produce more useful output, and to soon build
[preseq](https://github.com/smithlabcode/preseq) into falco.

The original duplication analysis method is still implemented in falco v2.0, and
can be turned on with `--orig-dups`. Here's how the old and new analyses differ.

#### Original method
- 100,000 unique reads are hashed and counted (first 50nt of each read).
- These are taken from the first 100,000 reads in the input.
- After the first 100,000 unique has been reached, reads that match one
  previously hashed are counted.

#### New method
- 1,000,000 reads are hashed and counted (first 50nt of each read).
- These are taken approximately uniformly throughout the input.
- Although there is no randomization, if multiple threads are used the results
  will appear as though they are randomly sampled due to fluctaions in thread
  speed changing which 1M reads are hashed.

## Citing falco

If falco was helpful for your research, you can cite us as follows:

de Sena Brandine G and Smith AD. Falco: high-speed FastQC emulation for quality
control of sequencing data. F1000Research 2021, 8:1874
(https://doi.org/10.12688/f1000research.21142.2)

## Copyright and License Information

```
MIT License

Copyright (c) 2026 Andrew D Smith and Guilherme de Sena Brandine

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
