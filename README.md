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
