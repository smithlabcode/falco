# falco changelog

## falco 2.0.4 (2026-11-21)

Changes:
* Improvements to user interface (argument parsing and verbose output to
  console)
* Falco now accepts BAM/SAM files with reads that don't have quality scores
* GC content can now be output with the original FastQC smoothing
* All non-ACGTN are converted and handled as N for robustness to IUPAC symbols
* Fixed bug in decompressing some BGZF files
* Fixed bugs related to grade cutoffs (by @MaximilianStammnitz)

## falco 2.0.3 (2026-09-27)

Changes:
* Bugfix: possible race condition for very small input files (or attempting to
  use massive numbers of threads) that so far has only been seen with very small
  data files of the kind used for unit tests.

## falco 2.0.2 (2026-09-09)

Changes:
* Input for preseq: falco can output a file (`outdir/preseq_hist.txt`) that
  can be directly input to preseq for library complexity analysis.
* Building on macOS: code updates for easier building using Apple Clang on
  macOS.
* Standard input: Falco can now accept piped or redirected input from stdin in
  either SAM or FASTQ format (thanks to @kerrycobb for initiating this).
* HTML report: for reads longer than 1000 bp, plots in the HTML will have
  bases grouped to prevent generating very large, slow and unreadable HTML.

## falco 2.0.1 (2026-08-12)

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
