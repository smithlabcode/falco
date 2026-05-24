# falco2
A revision of falco

# Note about threads

The following results categories can be impacted by thread use:

- Tile quality scores: The reads used to make estimates about tile qualities can
  differ depending on the number of threads used. The trends should not change
  and the magnitude of any differences should be small unless the data is very
  small.

- Kmer counts: Similar to tile quality scores, the identities of reads used for
  kmer statistics will differ depending on the number of threads used. The
  number of reads used for this analysis will not change.

# What's not yet working

- `fastqc_data.txt` slightly differs on tiles (seems to be arbitrary) and in GC
  content, which I'm not sure I want to directly replicate.
- Use of threads means that the "Duplication levels" are very hard to replicate
  exactly, though that might not matter very much.
- Counting duplicate reads truncates at 50nt, where previously anything above
  75nt was truncated to 50, while reads of length [51, 75] were left as their
  original length, which was a bit unusual.
- There is no `summary.txt` output file yet. The same info is in the
  `fastqc_data.txt`.
- There is no HTML formatting of the output.
