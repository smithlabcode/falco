/* Copyright (C) 2026 Andrew D. Smith
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "bam_stream_reader.hpp"

#include "StreamReader.hpp"

template <class T>
[[nodiscard]] inline T
min8(const T a, const T b) {
  return a < b ? a : b;
}

void
BamReader::process_sequence_base_from_buffer(FastqStats &stats) {
  // count Ns even if asked not to report them
  if (base_from_buffer == 'N') {
    ++stats.n_base_count[read_pos];
    num_bases_after_n = 1;  // start over the current kmer
    return;
  }

  const auto two_bit = actg_to_2bit(base_from_buffer);

  // ATGC bases, increments basic statistic counts
  cur_gc_count += (two_bit & 1);
  ++stats.base_count[(read_pos << Constants::bit_shift_base) | two_bit];

  if (do_sliding_window) {  // Update k-mer sequence
    cur_kmer = (cur_kmer << Constants::bit_shift_base) | two_bit;
    // registers k-mer if seen at least k nucleotides since the last n
    if (do_kmer && do_kmer_read && num_bases_after_n >= Constants::kmer_size) {
      ++stats.kmer_count[(read_pos << Constants::bit_shift_kmer) |
                         (cur_kmer & Constants::kmer_mask)];
      ++stats.pos_kmer_count[read_pos];
    }
    // GS: slow, need to use fsm
    if (do_adapter_optimized && num_bases_after_n == adapter_size) {
      cur_kmer &= adapter_mask;
      for (std::size_t i = 0; i < num_adapters; ++i)
        if (cur_kmer == adapters[i] && !adapters_found[i]) {
          const auto shifted = read_pos << Constants::bit_shift_adapter;
          ++stats.pos_adapter_count[shifted | i];
          adapters_found[i] = true;
        }
    }
    num_bases_after_n += (num_bases_after_n != adapter_size);
  }
}

// puts base either on buffer or leftover
void
BamReader::put_base_in_buffer(const std::size_t pos) {
  base_from_buffer = seq_nt16_str[bam_seqi(cur_char, pos)];
  if (still_in_buffer) {
    buffer[read_pos] = base_from_buffer;
    return;
  }
  if (leftover_ind > std::size(leftover_buffer))
    leftover_buffer.resize(leftover_ind);
  leftover_buffer[leftover_ind] = base_from_buffer;
}

[[nodiscard]] inline std::size_t
get_truncate_point(const std::size_t read_pos) {
  return read_pos <= Constants::unique_reads_max_length
           ? read_pos
           : Constants::unique_reads_truncate;
}

// Specially made to work directly with bam1_t
void
BamReader::read_sequence_line(FastqStats &stats) {
  if (!do_read)
    return;

  const std::size_t seq_len = aln.b->core.l_qseq;
  // MN: TODO: make sure everything works in this scope
  if (do_adapters_slow) {
    std::string seq_line_str(seq_len, '\0');
    for (std::size_t i = 0; i < seq_len; i++)
      seq_line_str[i] = seq_nt16_str[bam_seqi(cur_char, i)];
    for (std::size_t i = 0; i < num_adapters; ++i) {
      const std::size_t adapt_index = seq_line_str.find(adapter_seqs[i]);
      if (adapt_index < stats.SHORT_READ_THRESHOLD) {
        const auto x = (adapt_index + std::size(adapter_seqs[i]) - 1)
                       << Constants::bit_shift_adapter;
        ++stats.pos_adapter_count[x | i];
      }
    }
  }

  const auto trunc_point = get_truncate_point(seq_len);
  sequence_to_hash.clear();
  sequence_to_hash.resize(trunc_point);

  // restart line counters
  read_pos = 0;
  cur_gc_count = 0;
  truncated_gc_count = 0;
  num_bases_after_n = 1;
  still_in_buffer = true;
  next_truncation = 100;
  do_kmer_read = (stats.num_reads == next_kmer_read);
  adapters_found.reset();

  /*********************************************************/
  /********** THIS LOOP MUST BE ALWAYS OPTIMIZED ***********/
  /*********************************************************/
  // In the following loop, cur_char does not change, but rather i changes
  // and we access bases using bam_seqi(cur_char, i) in
  // put_base_in_buffer.
  for (std::size_t i = 0; i < seq_len; i++, ++read_pos) {
    // if we reached the buffer size, stop using it and start using leftover
    if (read_pos == buffer_size) {
      still_in_buffer = false;
      leftover_ind = 0;
    }

    // make sure we have sufficient space left in buffer
    if (!still_in_buffer && leftover_ind == stats.num_extra_bases)
      stats.allocate_new_base(tile_ignore);

    // puts base either on buffer or leftover
    // BamReader::put_base_in_buffer(i);
    base_from_buffer = seq_nt16_str[bam_seqi(cur_char, i)];
    if (i < trunc_point)
      sequence_to_hash[i] = base_from_buffer;

    // statistics updated base by base use buffer
    if (still_in_buffer)
      process_sequence_base_from_buffer(stats);
    else {  // use dynamic allocation
      process_sequence_base_from_leftover(stats);
      ++leftover_ind;  // increment leftover pos if no longer in buffer
    }

    // Truncate GC counts to multiples of 100
    if (do_gc_sequence && read_pos == next_truncation) {
      truncated_gc_count = cur_gc_count;
      truncated_length = read_pos;
      next_truncation += 100;
    }
  }

  // statistics summarized after the read
  postprocess_sequence_line(stats);
}

void
BamReader::read_quality_line(FastqStats &stats) {
  static constexpr auto max_qual_value = 255;
  if (!do_read || aln.b->core.qual == max_qual_value) {
    read_fast_forward_line_eof();
    return;
  }

  // reset quality counts
  read_pos = 0;
  cur_quality = 0;
  still_in_buffer = true;

  const std::size_t seq_len = aln.b->core.l_qseq;
  for (std::size_t i = 0; i < seq_len; ++cur_char, i++) {
    if (read_pos == buffer_size) {
      still_in_buffer = false;
      leftover_ind = 0;
    }

    get_base_from_buffer();

    // MN: Adding quality_zero to emulate the behavior of the original function.
    stats.lowest_char =
      min8(stats.lowest_char, static_cast<char>(*cur_char + 33));

    // Converts quality ascii to zero-based
    quality_value = *cur_char;

    if (still_in_buffer)  // Fast bases from buffer
      process_quality_base_from_buffer(stats);
    else {  // slow bases from dynamic allocation
      process_quality_base_from_leftover(stats);
      ++leftover_ind;
    }

    // sums quality value so we can bin the average at the end
    cur_quality += quality_value;

    // flag to start reading and writing outside of buffer
    ++read_pos;
  }

  // Average quality approximated to the nearest integer. Used to make a
  // histogram in the end of the summary.
  if (read_pos != 0)
    ++stats.quality_count[cur_quality / read_pos];  // avg quality histogram
}

/*******************************************************/
/*************** READ BAM RECORD ***********************/
/*******************************************************/

/* ADS: The table below (due to Masaru Nakajima) converts 2 bases packed in a
 * byte to their reverse complement. The input is therefore a unit8_t
 * representing 2 bases. It is assumed that the input uint8_t value is of form
 * "xx" or "x-", where 'x' a 4-bit number representing one of {A, C, G, T, N}
 * and '-' is 0000.  For example, the ouptut for "AG" is "CT". The format "x-"
 * is used at the end of an odd-length sequence. The output of "A-" is "-T", and
 * the output of "C-" is "-G", and so forth. The calling function must take care
 * when handling odd-length sequences.
 */
static const std::uint8_t byte_revcom_table[] = {
  // clang-format off
    0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,
    8, 136,  72,   0,  40,   0,   0,   0,
   24,   0,   0,   0,   0,   0,   0, 248,
    4, 132,  68,   0,  36,   0,   0,   0,
   20,   0,   0,   0,   0,   0,   0, 244,
    0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,
    2, 130,  66,   0,  34,   0,   0,   0,
   18,   0,   0,   0,   0,   0,   0, 242,
    0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,
    1, 129,  65,   0,  33,   0,   0,   0,
   17,   0,   0,   0,   0,   0,   0, 241,
    0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,
   15, 143,  79,   0,  47,   0,   0,   0,
   31,   0,   0,   0,   0,   0,   0, 255
};
// clang-format on

static inline void
revcom_byte_then_reverse(unsigned char *a, unsigned char *b) {
  unsigned char *p1, *p2;
  for (p1 = a, p2 = b - 1; p2 > p1; ++p1, --p2) {
    *p1 = byte_revcom_table[*p1];
    *p2 = byte_revcom_table[*p2];
    *p1 ^= *p2;
    *p2 ^= *p1;
    *p1 ^= *p2;
  }
  if (p1 == p2)
    *p1 = byte_revcom_table[*p1];
}

static inline void
revcomp_seq_by_byte(bam1_t *aln) {
  const int32_t l_qseq = aln->core.l_qseq;
  auto seq = bam_get_seq(aln);
  const std::int32_t n_bytes = (l_qseq + 1) >> 1;  // ceil(l_qseq/2.0)
  auto seq_end = seq + n_bytes;
  revcom_byte_then_reverse(seq, seq_end);
  if (l_qseq % 2 == 1) {  // for odd-length sequences
    for (int32_t i = 0; i < n_bytes - 1; ++i) {
      // swap 4-bit chunks within consecutive bytes like this:
      // (----aaaa bbbbcccc dddd....) => (aaaabbbb ccccdddd ....)
      seq[i] = (seq[i] << 4) | (seq[i + 1] >> 4);
    }
    seq[n_bytes - 1] <<= 4;  // final byte is just shifted
  }
}

static inline void
reverse_quality_scores(bam1_t *aln) {
  auto s = bam_get_qual(aln);
  if (s[0] == 0xff)
    return;
  std::reverse(s, bam_get_aux(aln));  // aux is just past qual
}

// set sam separator as tab
BamReader::BamReader(FalcoConfig &config, const std::size_t buffer_size) :
  StreamReader(config, buffer_size, '\t', '\n'), tp(config.threads),
  hts(config.filename), hdr(hts) {
  tp.set_io(hts);
}

[[nodiscard]] std::size_t
BamReader::load() {
  // if (!(hts = hts_open(std::data(filename), "r")))
  //   throw std::runtime_error("cannot load bam file : " + filename);
  // if (!(hdr = sam_hdr_read(hts)))
  //   throw std::runtime_error("failed to read header from file: " + filename);
  // if (!(b = bam_init1()))
  //   throw std::runtime_error("failed to read record from file: " + filename);
  // GS TODO: implement this properly
  return 1;
}

// We will check eof on the >> operator
inline bool
BamReader::is_eof() {
  return (cur_char == last - 1);
}

void
BamReader::postprocess_fastq_record(FastqStats &stats) {
  if (do_sequence_hash) {
    const auto itr = stats.sequence_count.find(sequence_to_hash);
    if (itr != std::cend(stats.sequence_count)) {
      itr->second++;
      stats.count_at_limit += continue_storing_sequences;
    }
    else if (continue_storing_sequences) {
      stats.sequence_count.insert({{sequence_to_hash, 1}});
      stats.count_at_limit = stats.num_reads;
      ++stats.num_unique_seen;
      continue_storing_sequences =
        (stats.num_unique_seen < Constants::unique_reads_stop_counting);
    }
  }
  // count tile if applicable
  next_tile_read += (!tile_ignore && do_tile_read) ? num_reads_for_tile : 0;
  next_kmer_read += do_kmer_read ? num_reads_for_kmer : 0;
}

bool
BamReader::read_entry(FastqStats &stats, std::size_t &num_bytes_read) {
  // const auto bam_to_string = [&](const auto &h, const auto &a) {
  //   kstring_t ks = KS_INITIALIZE;
  //   const int ret = sam_format1(h.h, a.b, &ks);
  //   if (ret < 0)
  //     throw std::runtime_error("Can't format record: " +
  //                              std::string(bam_get_qname(a.b)));
  //   std::string s{ks.s};
  //   ks_free(&ks);
  //   return s;
  // };

  if (hts.read(hdr, aln)) {
    auto b = aln.b;
    if (bam_is_rev(b)) {
      revcomp_seq_by_byte(b);
      reverse_quality_scores(b);
      b->core.flag ^= BAM_FREVERSE;
    }
    num_bytes_read = 0;
    do_read = (stats.num_reads == next_read);

    // Read tile line
    cur_char = bam_get_qname(b);
    last = cur_char + b->m_data;
    const std::size_t first_padding_null =
      b->core.l_qname - b->core.l_extranul - 1;
    // Turn "QUERYNAME\0\0\0" into "QUERYNAME\t\0\0" (assuming
    // field_separtor = '\t') to be compatible with read_fast_forward_line().
    cur_char[first_padding_null] = field_separator;
    read_tile_line(stats);

    // Read sequence line
    std::size_t seq_len = b->core.l_qseq;
    cur_char = reinterpret_cast<char *>(bam_get_seq(b));
    BamReader::read_sequence_line(stats);

    // Read quality line
    cur_char = reinterpret_cast<char *>(bam_get_qual(b));

    // Set the first byte after qual to line_separator so that read_quality_line
    // stops at the end of qual
    cur_char[seq_len] = line_separator;

    BamReader::read_quality_line(stats);

    if (do_read)
      postprocess_fastq_record(stats);

    next_read += do_read * read_step;
    ++stats.num_reads;
    return true;
  }
  // If I could not read another line it means it's eof
  return false;
}
