# falco changelog

## falco 2.0.1 (2026-08-03)

Changes:
* GC content distribution across reads now smoothing more similarly to falco v1
* Minor improvements to HTML, including fixing axis ranges
* Adding `--bisulfite` more from falco v1
* Original duplicates mode now enabled with command line flag `--orig-dups`
* Bugfix: rare issue with parsing SAM quality scores now fixed (BAM was ok)

## falco 2.0.0 (2026-08-03)

Changes:
* Source is 100% rewritten, so individual changes not listed here.
* Output files are the same: `summary.txt`, `fastqc_data.txt`,
  `fastqc_report.html`
* All output is written to a subdirectory of the output directory, named based
  on input file names.
* Command line arguments are 100% different from falco v1.
* File formats accepted: FASTQ (plain, gzip and BGZF) and BAM/SAM
* Falco v2.0.0 is changing to a MIT license
