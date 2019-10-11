#include "FastqStats.hpp"
#include <algorithm>
#include <cmath>
using std::string;
using std::vector;
using std::unordered_map;
using std::sort;
using std::min;
using std::max;
using std::ostream;

/******************* IMPLEMENTATION OF FASTQC FUNCTIONS **********************/
// FastQC extrapolation of counts to the full file size
double get_corrected_count(size_t count_at_limit,
                           size_t num_reads,
                           size_t dup_level,
                           size_t num_obs) {
  // See if we can bail out early
  if (count_at_limit == num_reads)
    return num_obs;

  // If there aren't enough sequences left to hide another sequence with this
  // count the we can also skip the calculation
  if (num_reads - num_obs < count_at_limit)
    return num_obs;

  // If not then we need to see what the likelihood is that we had
  // another sequence with this number of observations which we would
  // have missed. We'll start by working out the probability of NOT seeing a
  // sequence with this duplication level within the first count_at_limit
  // sequences of num_obs.  This is easier than calculating
  // the probability of seeing it.
  double p_not_seeing = 1.0;

  // To save doing long calculations which are never going to produce anything
  // meaningful we'll set a limit to our p-value calculation.  This is the
  // probability below which we won't increase our count by 0.01 of an
  // observation.  Once we're below this we stop caring about the corrected
  // value since it's going to be so close to the observed value thatwe can
  // just return that instead.
  double limit_of_caring = 1.0 - (num_obs/(num_obs + 0.01));
  for (size_t i = 0; i < count_at_limit; ++i) {
    p_not_seeing *= static_cast<double>((num_reads-i)-dup_level) /
                         static_cast<double>(num_reads-i);

    if (p_not_seeing < limit_of_caring) {
      p_not_seeing = 0;
      break;
    }
  }

  // Now we can assume that the number we observed can be
  // scaled up by this proportion
  return num_obs/(1 - p_not_seeing);
}

// Function to create a smooth gc model similar to FastQC's
inline void
smooth_gc_model(std::array <double, 101> &smooth_gc_count,
                std::array <size_t, 101> gc_count) {

  // impute based on adjacent values
  for (size_t i = 1; i < 100; ++i) {
    if (gc_count[i] == 0) {
      gc_count[i] = (gc_count[i-1] + gc_count[i+1])/2.0;
    }
  }

  // smooth with a regular sliding window
  for (size_t i = 1; i < 100; ++i) {
    smooth_gc_count[i] = (gc_count[i-1] + gc_count[i] + gc_count[i+1])/3.0;
  }
}

// Function to calculate the deviation of a histogram with 100 bins from a
// theoretical normal distribution with same mean and standard deviation
double
sum_deviation_from_normal(const std::array <double, 101> &smooth_gc_count,
                           std::array <double, 101> &theoretical) {
  double mode = 0.0,  num_reads = 0.0;

  // Mode and total number of reads (note that we "smoothed" the gc content by
  // filling up the 0 reads, so the number of reads is not the same as the
  // number of reads in the whole fastq
  for (size_t i = 0; i < 101; ++i) {
    mode += i*smooth_gc_count[i];
    num_reads += smooth_gc_count[i];
  }
  mode /= num_reads;

  // stadard deviation from the modde
  double stdev = 0.0;
  for (size_t i = 0; i < 101; ++i)
    stdev += (mode - i) * (mode - i) * smooth_gc_count[i];
  stdev = sqrt(stdev / (num_reads - 1));

  // theoretical sampling from a normal distribution with mean = mode and stdev
  // = stdev to the mode from the sampled gc content from the data
  double ans = 0.0, theoretical_sum = 0.0, z;
  theoretical.fill(0);
  for (size_t i = 0; i < 101; ++i) {
    z = i - mode;
    theoretical[i] = exp(- (z*z)/ (2 * stdev *stdev));
    theoretical_sum += theoretical[i];
  }

  // Normalize theoretical so it sums to the total of readsq
  for (size_t i = 0; i < 101; ++i) {
    theoretical[i] = theoretical[i] * num_reads / theoretical_sum;
  }
  for (size_t i = 0; i < 101; ++i)
    ans += fabs(smooth_gc_count[i] - theoretical[i]);

  // Fractional deviation
  return 100.0 * ans / num_reads;
}

/****************************************************************/
/******************** FASTQ STATS FUNCTIONS *********************/
/****************************************************************/
// Default constructor
FastqStats::FastqStats() {
  total_bases = 0;
  num_extra_bases = 0;
  avg_read_length = 0;
  total_gc = 0;
  avg_gc = 0;
  num_reads = 0;
  num_reads_kmer = 0;
  min_read_length = 0;
  max_read_length = 0;
  num_poor = 0;

  num_unique_seen = 0;
  count_at_limit = 0;

  // Initialize IO arrays
  base_count.fill(0);
  n_base_count.fill(0);
  read_length_freq.fill(0);
  quality_count.fill(0);
  gc_count.fill(0);
  position_quality_count.fill(0);

  // Defines k-mer mask, length and allocates vector
  kmer_mask = (1ll << (2*kmer_size)) - 1;
  kmer_count = vector<size_t>(min(kNumBases, kKmerMaxBases)
                                  * (kmer_mask + 1), 0);
}

// When we read new bases, dynamically allocate new space for their statistics
void 
FastqStats::allocate_new_base(const bool ignore_tile) {
  for (size_t i = 0; i < kNumNucleotides; ++i) {
    long_base_count.push_back(0);
  }
  long_n_base_count.push_back(0);

  // space for quality boxplot
  for (size_t i = 0; i < kNumQualityValues; ++i)
    long_position_quality_count.push_back(0);

  long_read_length_freq.push_back(0);

  // space for tile quality in each position.
  if (!ignore_tile) {
    for (auto  &v : tile_position_quality)
      v.second.push_back(0);
    for (auto  &v : tile_position_count)
      v.second.push_back(0);
  }

  // Successfully allocated space for a new base
  ++num_extra_bases;
}

// Calculates all summary statistics and pass warn fails
void
FastqStats::summarize(Config &config) {
  /******************* BASIC STATISTICS **********************/
  pass_basic_statistics = "pass";  // in fastqc, basic statistics is always pass
  // Average read length
  avg_read_length = 0;
  for (size_t i = 0; i < max_read_length; ++i) {
    if (i < kNumBases)
      total_bases += i * read_length_freq[i];
    else
      total_bases += i * long_read_length_freq[i - kNumBases];
  }

  avg_read_length = total_bases / num_reads;

  // counts bases G and C in each base position
  avg_gc = 0;

  // GC %
  avg_gc = 100 * total_gc / static_cast<double>(total_bases);

  // Poor quality reads
  num_poor = 0;
  for (size_t i = 0; i < config.kPoorQualityThreshold; ++i)
    num_poor += quality_count[i];

  // Cumulative read length frequency
  size_t cumulative_sum = 0;
  for (size_t i = 0; i < max_read_length; ++i) {
    if (i < kNumBases) {
      cumulative_sum += read_length_freq[i];
      if (read_length_freq[i] > 0)
        if (min_read_length == 0)
          min_read_length = i + 1;
    }
    else
      cumulative_sum += long_read_length_freq[i - kNumBases];
  }

  for (size_t i = 0; i < max_read_length; ++i) {
    if (i < kNumBases) {
      cumulative_read_length_freq[i] = cumulative_sum;
      cumulative_sum -= read_length_freq[i];
    }
    else {
      long_cumulative_read_length_freq.push_back(cumulative_sum);
      cumulative_sum -= long_read_length_freq[i - kNumBases];
    }
  }

  /******************* PER BASE SEQUENCE QUALITY **********************/
  if (config.do_quality_base) {
    pass_per_base_sequence_quality = "pass";

    // Quality quantiles for base positions
    size_t cur;  // for readability, get the quality x position count
    double ldecile_thresh,
           lquartile_thresh,
           median_thresh,
           uquartile_thresh,
           udecile_thresh;

    size_t cur_ldecile = 0,
           cur_lquartile = 0,
           cur_median = 0,
           cur_uquartile = 0,
           cur_udecile = 0;

    double cur_mean;
    for (size_t i = 0; i < max_read_length; ++i) {
      cur_mean = 0;
      size_t counts = 0;

      // Number of counts I need to see to know in which bin each *ile is
      if (i < kNumBases) {
        ldecile_thresh = 0.1 * cumulative_read_length_freq[i];
        lquartile_thresh = 0.25 * cumulative_read_length_freq[i];
        median_thresh = 0.5 * cumulative_read_length_freq[i];
        uquartile_thresh = 0.75 * cumulative_read_length_freq[i];
        udecile_thresh = 0.9 * cumulative_read_length_freq[i];
      } else {
        ldecile_thresh = 0.1 * long_cumulative_read_length_freq[i - kNumBases];
        lquartile_thresh =
          0.25 * long_cumulative_read_length_freq[i - kNumBases];
        median_thresh = 0.5 * long_cumulative_read_length_freq[i - kNumBases];
        uquartile_thresh =
          0.75 * long_cumulative_read_length_freq[i - kNumBases];
        udecile_thresh = 0.9 * long_cumulative_read_length_freq[i - kNumBases];
      }
      // Iterate through quality values to find quantiles in increasing order
      for (size_t j = 0; j < kNumQualityValues; ++j) {
        if (i < kNumBases)
          cur = position_quality_count[(i << kBitShiftQuality) | j];
        else
          cur = long_position_quality_count[
            ((i - kNumBases) << kBitShiftQuality) | j];

        // Finds in which bin of the histogram reads are
        if (counts < ldecile_thresh && counts + cur >= ldecile_thresh)
          cur_ldecile = j;

        if (counts < lquartile_thresh && counts + cur >= lquartile_thresh)
          cur_lquartile = j;

        if (counts < median_thresh && counts + cur >= median_thresh)
          cur_median = j;

        if (counts < uquartile_thresh && counts + cur >= uquartile_thresh)
          cur_uquartile = j;

        if (counts < udecile_thresh && counts + cur >= udecile_thresh)
          cur_udecile = j;

        cur_mean += cur*j;
        counts += cur;
      }

      // Normalize mean
      if (i < kNumBases)
        cur_mean = cur_mean / cumulative_read_length_freq[i];
      else
        cur_mean = cur_mean / long_cumulative_read_length_freq[i - kNumBases];

      if (i < kNumBases) {
        mean[i] = cur_mean;
        ldecile[i] = cur_ldecile;
        lquartile[i] = cur_lquartile;
        median[i] = cur_median;
        uquartile[i] = cur_uquartile;
        udecile[i] = cur_udecile;
      } else {
        long_mean.push_back(cur_mean);
        long_ldecile.push_back(cur_ldecile);
        long_lquartile.push_back(cur_lquartile);
        long_median.push_back(cur_median);
        long_uquartile.push_back(cur_uquartile);
        long_udecile.push_back(cur_udecile);
      }

      // Pass warn fail criteria
      if (pass_per_base_sequence_quality != "fail") {
        if (cur_lquartile < config.limits["quality_base_lower"]["error"])
          pass_per_base_sequence_quality = "fail";
        else if (cur_lquartile < config.limits["quality_base_lower"]["warn"])
          pass_per_base_sequence_quality = "warn";

        if (cur_median < config.limits["quality_base_median"]["error"])
          pass_per_base_sequence_quality = "fail";
        else if (cur_median < config.limits["quality_base_median"]["warn"])
          pass_per_base_sequence_quality = "warn";
      }
    }
  }
  /******************* PER SEQUENCE QUALITY SCORE **********************/
  if (config.do_quality_sequence) {
    pass_per_sequence_quality_scores = "pass";
    size_t mode_val = 0;
    size_t mode_ind = 0;
    for (size_t i = 0; i < kNumQualityValues; ++i) {
      if (quality_count[i] > mode_val) {
        mode_val = quality_count[i];
        mode_ind = i;
      }
    }

    if (mode_ind < config.limits["quality_sequence"]["warn"])
      pass_per_sequence_quality_scores = "warn";

    else if (mode_ind < config.limits["quality_sequence"]["error"])
      pass_per_sequence_quality_scores = "fail";
  }
  /******************* PER BASE SEQUENCE CONTENT **********************/
  if (config.do_sequence) {
    pass_per_base_sequence_content = "pass";
    double a, t, g, c, n;
    double total;
    double max_diff = 0.0;
    for (size_t i = 0; i < max_read_length; ++i) {
      if (i < kNumBases) {
        a = base_count[(i << kBitShiftNucleotide)];
        c = base_count[(i << kBitShiftNucleotide) | 1];
        t = base_count[(i << kBitShiftNucleotide) | 2];
        g = base_count[(i << kBitShiftNucleotide) | 3];
        n = n_base_count[i];
      } else {
        a = long_base_count[((i - kNumBases) << kBitShiftNucleotide)];
        c = long_base_count[((i - kNumBases) << kBitShiftNucleotide) | 1];
        t = long_base_count[((i - kNumBases) << kBitShiftNucleotide) | 2];
        g = long_base_count[((i - kNumBases) << kBitShiftNucleotide) | 3];
        n = long_n_base_count[i - kNumBases];
      }

      // turns above values to percent
      total = static_cast<double>(a + c + t + g + n);
      a = 100.0*a / total;
      c = 100.0*c / total;
      t = 100.0*t / total;
      g = 100.0*g / total;
      n = 100.0*n / total;
      if (i < kNumBases) {
        g_pct[i] = g;
        a_pct[i] = a;
        t_pct[i] = t;
        c_pct[i] = c;
        n_pct[i] = n;
      } else {
        long_g_pct.push_back(g);
        long_a_pct.push_back(a);
        long_t_pct.push_back(t);
        long_c_pct.push_back(c);
        long_n_pct.push_back(n);
      }

      max_diff = max(max_diff, fabs(a-c));
      max_diff = max(max_diff, fabs(a-t));
      max_diff = max(max_diff, fabs(a-g));
      max_diff = max(max_diff, fabs(c-t));
      max_diff = max(max_diff, fabs(c-g));
      max_diff = max(max_diff, fabs(t-g));

      if (pass_per_base_sequence_content != "fail") {
        if (max_diff > config.limits["sequence"]["error"])
          pass_per_base_sequence_content = "fail";
        else if (max_diff > config.limits["sequence"]["warn"])
          pass_per_base_sequence_content = "warn";
      }
    }
  }

  /******************* PER SEQUENCE GC CONTENT *****************/
  if (config.do_gc_sequence) {
    pass_per_sequence_gc_content = "pass";
    // Here we smooth the gc content using a sliding window:
    smooth_gc_model(smooth_gc_count, gc_count);

    // Calculate pass warn fail statistics
    double gc_deviation = sum_deviation_from_normal(smooth_gc_count,
                                                    theoretical_gc_count);
    if (gc_deviation >= config.limits["gc_sequence"]["error"]) {
      pass_per_sequence_gc_content = "fail";
    }
    else if (gc_deviation >= config.limits["gc_sequence"]["warn"]) {
      pass_per_sequence_gc_content = "warn";
    }
  }
  /******************* PER BASE N CONTENT **********************/
  if (config.do_n_content) {
    pass_per_base_n_content = "pass";
    double cur_n_pct;
    for (size_t i = 0; i < max_read_length; ++i) {
      if (pass_per_base_n_content != "fail") {
        if (i < kNumBases) {
          cur_n_pct = n_pct[i];
        }
        else {
          cur_n_pct = long_n_pct[i - kNumBases];
        }

        if (cur_n_pct > config.limits["n_content"]["error"]) {
          pass_per_base_n_content = "fail";
        }

        else if (cur_n_pct > config.limits["n_content"]["warn"]) {
          pass_per_base_n_content = "warn";
        }
      }
    }
  }
  /************** SEQUENCE LENGTH DISTRIBUTION *****************/
  if (config.do_sequence_length) {
    pass_sequence_length_distribution = "pass";
    size_t freq_of_avg;

    if (avg_read_length < kNumBases) {
      freq_of_avg = read_length_freq[avg_read_length];
    }
    else {
      freq_of_avg = long_read_length_freq[avg_read_length - kNumBases];
    }

    if (config.limits["sequence_length"]["error"] == 1) {
      if (freq_of_avg != num_reads) {
        pass_sequence_length_distribution = "warn";
      }

      if (read_length_freq[0] > 0) {
        pass_sequence_length_distribution = "fail";
      }
    }
  }

  /************** DUPLICATE SEQUENCES **************************/
  if (config.do_duplication) {
    pass_duplicate_sequences = "pass";

    double seq_total = 0.0;
    double seq_dedup = 0.0;

    // Key is frequenccy (r), value is number of times we saw a sequence
    // with that frequency
    unordered_map <size_t, size_t> counts_by_freq;
    for (auto v : sequence_count) {
      if (counts_by_freq.count(v.second) == 0) {
        counts_by_freq[v.second] = 0;
      }
      counts_by_freq[v.second]++;
    }

    // Now we run the fastqc corrected extrapolation
    for (auto v : counts_by_freq) {
      counts_by_freq[v.first] = get_corrected_count(count_at_limit,
                                                    num_reads,
                                                    v.first,
                                                    v.second);
    }

    // Group in blocks similarly to fastqc
    for (auto v : counts_by_freq) {
      size_t dup_slot = v.first - 1;
      if (v.first >= 10000) dup_slot = 15;
      else if (v.first >= 5000) dup_slot = 14;
      else if (v.first >= 1000) dup_slot = 13;
      else if (v.first >= 500) dup_slot = 12;
      else if (v.first >= 100) dup_slot = 11;
      else if (v.first >= 50) dup_slot = 10;
      else if (v.first >= 10) dup_slot = 9;

      percentage_deduplicated[dup_slot] += v.second;
      percentage_total[dup_slot] += v.second * v.first;

      seq_total += v.second * v.first;
      seq_dedup += v.second;
    }

    // "Sequence duplication estimate" in the summary
    total_deduplicated_pct = 100.0 * seq_dedup / seq_total;

    // Convert to percentage
    for (auto &v : percentage_deduplicated)
      v = 100.0 * v / seq_dedup;  // Percentage of unique sequences in bin

     // Convert to percentage
    for (auto &v : percentage_total)
      v = 100.0 * v / seq_total;  // Percentage of sequences in bin

    // pass warn fail criteria : unique reads must be >80%
    // (otherwise warn) or >50% (otherwisefail)
    if (total_deduplicated_pct <= config.limits["duplication"]["error"]) {
      pass_duplicate_sequences = "fail";
    }
    else if (total_deduplicated_pct <= config.limits["duplication"]["warn"]) {
      pass_duplicate_sequences = "warn";
    }
  }
  /************** OVERREPRESENTED SEQUENCES ********************/
  if (config.do_overrepresented) {
    pass_overrepresented_sequences = "pass";

    // Keep only sequences that pass the input cutoff
    for (auto it = sequence_count.begin(); it != sequence_count.end(); ++it) {
      if (it->second > num_reads * config.kOverrepMinFrac) {
        overrep_sequences.push_back(*it);
      }

      // implment pass warn fail for overrep sequences
      if (pass_overrepresented_sequences != "fail") {
        // get percentage that overrep reads represent
        double pct = 100.0 * it->second / num_reads;
        if (pct > config.limits["overrepresented"]["error"]) {
          pass_overrepresented_sequences = "fail";
        }
        else if (pct > config.limits["overrepresented"]["warn"]) {
          pass_overrepresented_sequences = "warn";
        }
      }
    }

    // Sort strings by frequency
    sort(overrep_sequences.begin(), overrep_sequences.end(),
         [](auto &a, auto &b){
            return a.second > b.second;
         });
  }
  /************** ADAPTER CONTENT ******************************/
  if (config.do_adapter) {
    pass_adapter_content = "pass";

    // Cumulative count of adapter kmers by position
    size_t jj;
    for (size_t i = 0; i < min(kNumBases, kKmerMaxBases); ++i) {
      if (cumulative_read_length_freq[i] > 0) {
        // Makes the count of kmers by position cumulative by copying
        // the previous position count
        if (i == 0) {
          kmer_by_base[i] = vector<double> (config.adapters.size(), 0.0);
        }
        else {
          kmer_by_base[i] = vector<double> (kmer_by_base[i-1].begin(),
                                            kmer_by_base[i-1].end());
        }
        jj = 0;

        // Get count for adapter's k-prefix
        for (auto v : config.adapters) {
          size_t kmer_ind = (i << kBitShiftKmer) | v.second;
          kmer_by_base[i][jj] += kmer_count[kmer_ind];
          ++jj;
        }
      }
    }

    for (size_t i = 0; i < min(kNumBases, kKmerMaxBases); ++i) {
      if (cumulative_read_length_freq[i] > 0) {
        jj = 0;
        for (auto v : config.adapters) {
          kmer_by_base[i][jj] = kmer_by_base[i][jj] * 100.0 / num_reads_kmer;

          // Update pass warn fail
          if (pass_adapter_content != "fail") {
            if (kmer_by_base[i][jj] > config.limits["adapter"]["error"])
              pass_adapter_content = "fail";
            else if (kmer_by_base[i][jj] > config.limits["adapter"]["warn"])
              pass_adapter_content = "warn";
          }
          ++jj;
        }
      }
    }
  }

  /************** PER TILE SEQUENCE QUALITY ********************/
  if (config.do_tile) {
    // First thing is to average all tile qualities in all positions and
    // subtract the mean for that base position
    double mean_in_base;
    for (auto &v : tile_position_quality) {
      for (size_t i = 0; i < max_read_length; ++i) {
        if (i < kNumBases) mean_in_base = mean[i];
        else mean_in_base = long_mean[i - kNumBases];
        v.second[i] = v.second[i] / tile_position_count[v.first][i]
                      - mean_in_base;
      }
    }

    // Now check the deviation based on the config file
    pass_per_tile_sequence_quality = "pass";
    for (size_t i = 0; i < max_read_length; ++i) {
      for (auto &v : tile_position_quality) {
        if (pass_per_tile_sequence_quality != "fail") {
          if (v.second[i] <= -config.limits["tile"]["error"]) {
            pass_per_tile_sequence_quality = "fail";
          }
          else if (v.second[i] <= -config.limits["tile"]["warn"]) {
            pass_per_tile_sequence_quality = "warn";
          }
        }
      }
    }
  }
}

/****************** WRITE STATS ***********************/
void
FastqStats::write(ostream &os, const Config &config) {
  // Header
  os << "##FastQC\t0.11.8\n";

  // Basic statistics
  os << ">>Basic Statistics\t" << pass_basic_statistics << "\n";
  os << "#Measure\tValue\n";
  os << "Filename\t" << config.filename_stripped << "\n";
  os << "Total Sequences\t" << num_reads << "\n";
  os << "Sequences flagged as poor quality \t" << num_poor << "\n";

  if (config.do_sequence)
    os << "%GC \t" << avg_gc << "\n";
  os << ">>END_MODULE\n";

  // Per base quality
  if (config.do_quality_base) {
    os << ">>Per base sequence quality\t" <<
           pass_per_base_sequence_quality << "\n";

    os << "#Base\tMean\tMedian\tLower Quartile\tUpper Quartile" <<
          "\t10th Percentile 90th Percentile\n";
    for (size_t i = 0; i < max_read_length; ++i) {
      if (i < kNumBases) {
        // Write distribution to new line
        os << i + 1 << "\t"
           << mean[i] << "\t"
           << median[i] << "\t"
           << lquartile[i] << "\t"
           << uquartile[i] << "\t"
           << ldecile[i] << "\t"
           << udecile[i] << "\n";
      } else {
        os << i + 1 << "\t"
           << long_mean[i - kNumBases] << "\t"
           << long_median[i - kNumBases] << "\t"
           << long_lquartile[i - kNumBases] << "\t"
           << long_uquartile[i - kNumBases] << "\t"
           << long_ldecile[i - kNumBases] << "\t"
           << long_udecile[i - kNumBases] << "\n";
      }
    }
    os << ">>END_MODULE\n";
  }

  // Per sequence quality scores
  if (config.do_quality_sequence) {
    os << ">>Per sequence quality scores\t" <<
          pass_per_sequence_quality_scores << "\n";

    os << "#Quality\tCount\n";

    for (size_t i = 0; i < kNumQualityValues; ++i) {
      if (quality_count[i] > 0)
        os << i << "\t" << quality_count[i] << "\n";
    }
    os << ">>END_MODULE\n";

    // Per base sequence content
    os << ">>Per base sequence content\t" <<
          pass_per_base_sequence_content << "\n";

    os << "#Base\tG\tA\tT\tC\n";

    for (size_t i = 0; i < max_read_length; ++i) {
      if (i < kNumBases) {
        os << i+1 << "\t" <<
              g_pct[i] << "\t" <<
              a_pct[i] << "\t" <<
              t_pct[i] << "\t" <<
              c_pct[i] << "\n";
      } else {
        os << i+1 << "\t" <<
              long_g_pct[i - kNumBases] << "\t" <<
              long_a_pct[i - kNumBases] << "\t" <<
              long_t_pct[i - kNumBases] << "\t" <<
              long_c_pct[i - kNumBases] << "\n";
      }
    }
    os << ">>END_MODULE\n";
  }

  // Per tile sequence quality
  // get tile values and sort
  vector<size_t> tiles_sorted;
  for (auto v : tile_count) {
    tiles_sorted.push_back(v.first);
  }

  sort(tiles_sorted.begin(), tiles_sorted.end());
  if (config.do_tile) {
    os << ">>Per tile sequence quality\t" <<
          pass_per_tile_sequence_quality << "\n";

    // prints tiles sorted by value
    for (size_t i = 0; i < tiles_sorted.size(); ++i) {
      for (size_t j = 0; j < max_read_length; ++j) {
        os << tiles_sorted[i] << "\t" << j + 1 << "\t"
           << tile_position_quality[tiles_sorted[i]][j];
        os << "\n";
      }
    }

    os << ">>END_MODULE\n";
  }

  // Per sequence gc content
  if (config.do_gc_sequence) {
    os << ">>Per sequence gc content\t" << pass_per_sequence_gc_content << "\n";
    os << "#GC Content\tCount\n";
    for (size_t i = 0; i <= 100; ++i) {
      if (gc_count[i] > 0) {
        os << i << "\t" << gc_count[i] << "\n";
      }
    }
    os << ">>END_MODULE\n";

    // Per base N content
    os << ">>Per base N concent\t" << pass_per_base_n_content << "\n";
    os << "#Base\tN-Count\n";

    for (size_t i = 0; i < max_read_length; ++i) {
      if (i < kNumBases)
        os << i+1 << "\t" << n_pct[i] << "\n";
      else
        os << i+1 << "\t" << long_n_pct[i - kNumBases] << "\n";
    }

    os << ">>END_MODULE\n";
  }

  // Sequence length distribution
  if (config.do_sequence_length) {
    os << "Sequence Length Distribution\t" <<
          pass_sequence_length_distribution << "\n";

    os << "Length\tCount\n";
    for (size_t i = 0; i < max_read_length; ++i) {
      if (i < kNumBases) {
        if (read_length_freq[i] > 0) {
          os << i+1 << "\t" << read_length_freq[i] << "\n";
        }
      } else {
        if (long_read_length_freq[i - kNumBases] > 0) {
          os << i + 1  << "\t" << long_read_length_freq[i - kNumBases] << "\n";
        }
      }
    }
    os << ">>END_MODULE\n";
  }

  // Sequence duplication levels
  if (config.do_duplication) {
    os << ">>Sequence Duplication Levels\t" <<
           pass_duplicate_sequences << "\n";

    os << ">>Total Deduplicated Percentage\t" <<
           total_deduplicated_pct << "\n";

    os << "#Duplication Level  Percentage of deduplicated  "
       << "Percentage of total\n";

    for (size_t i = 0; i < 9; ++i) {
      os << i+1 << "\t" << percentage_deduplicated[i] << "\t"
         << percentage_total[i] << "\n";
    }

    os << ">10\t" << percentage_deduplicated[9]
       << "\t" << percentage_total[9] << "\n";
    os << ">50\t" << percentage_deduplicated[10]
       << "\t" << percentage_total[10] << "\n";
    os << ">100\t" << percentage_deduplicated[11]
       << "\t" << percentage_total[11] << "\n";
    os << ">500\t" << percentage_deduplicated[12]
       << "\t" << percentage_total[12] << "\n";
    os << ">1k\t" << percentage_deduplicated[13]
       << "\t" << percentage_total[13] << "\n";
    os << ">5k\t" << percentage_deduplicated[14]
       << "\t" << percentage_total[14] << "\n";
    os << ">10k+\t" << percentage_deduplicated[15]
       << "\t" << percentage_total[15] << "\n";
    os << ">>END_MOUDLE\n";
  }

  // Overrepresented sequences
  if (config.do_overrepresented) {
    os << ">>Overrepresented sequences\t" <<
          pass_overrepresented_sequences << "\n";
    os << "#Sequence\tCount\tPercentage\tPossible Source\n";

    for (auto seq : overrep_sequences) {
        os << seq.first << "\t" << seq.second <<  "\t" <<
          100.0 * seq.second / num_reads << "\t"
          << config.get_matching_contaminant(seq.first) << "\n";
    }
    os << ">>END_MODULE\n";
    os << ">>Adapter Content\t" << pass_adapter_content << "\n";

    os << "#Position\t";
    for (auto v : config.adapters) {
      os << v.first << "\t";
    }

    os << "\n";

    // Number of kmers counted
    size_t jj;
    for (size_t i = 0; i < min(kNumBases, kKmerMaxBases); ++i) {
      if (cumulative_read_length_freq[i] > 0) {
        os << i + 1 << "\t";
        jj = 0;
        for (auto v : config.adapters) {
          os << kmer_by_base[i][jj];
          ++jj;

          if (jj != config.adapters.size()) {
            os << " ";
          }
        }
        os << "\n";
      }
    }
    os << ">>END_MODULE\n";
  }
}

