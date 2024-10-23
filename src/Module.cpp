/* Copyright (C) 2019 Guilherme De Sena Brandine and
 *                    Andrew D. Smith
 * Authors: Guilherme De Sena Brandine, Andrew Smith
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 */
// clang-format off
#include "Module.hpp"
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <sstream>
#include <cstdlib>
#include <cassert>

using std::string;
using std::vector;
using std::array;
using std::unordered_map;
using std::sort;
using std::min;
using std::max;
using std::ostream;
using std::pair;
using std::transform;
using std::setprecision;
using std::runtime_error;
using std::make_pair;
using std::ostringstream;
using std::istringstream;
using std::getline;

template<typename T> using num_lim = std::numeric_limits<T>;

/*****************************************************************************/
/******************* AUX FUNCTIONS *******************************************/
/*****************************************************************************/
// convert string to uppercase, to be used in making the short summaries where
// PASS, WARN and FAIL are uppercased
static string
toupper(const string &s) {
  string out;
  transform(s.begin(), s.end(), std::back_inserter(out),(int (*)(int))toupper);
  return out;
}


/*****************************************************************************/
/******************* BASEGROUPS COPIED FROM FASTQC ***************************/
/*****************************************************************************/

/************* NO BASE GROUP ****************/
void
make_default_base_groups(vector<BaseGroup> &base_groups,
                         const size_t &num_bases) {
  base_groups.clear();
  for (size_t i = 0; i < num_bases; ++i)
    base_groups.push_back(BaseGroup(i,i));
}


/************* EXP BASE GROUP **************/
void
make_exponential_base_groups(vector<BaseGroup> &base_groups,
                 const size_t &num_bases) {
  size_t starting_base = 0,
         end_base,
         interval = 1;

  base_groups.clear();
  for (; starting_base < num_bases;) {
    end_base = starting_base + interval - 1;
    if (end_base >= num_bases)
      end_base = num_bases;

    base_groups.push_back(BaseGroup(starting_base, end_base));
    starting_base += interval;
    if (starting_base == 9 && num_bases > 75)
      interval = 5;
    if (starting_base == 49 && num_bases > 200)
      interval = 10;
    if (starting_base == 99 && num_bases > 300)
      interval = 50;
    if (starting_base == 499 && num_bases > 1000)
      interval = 100;
    if (starting_base == 1000 && num_bases > 2000)
      interval = 500;
  }
}


/************* LINEAR BASE GROUP *************/
// aux function to get linear interval
size_t
get_linear_interval(const size_t num_bases) {
  // The the first 9bp as individual residues since odd stuff
  // can happen there, then we find a grouping value which gives
  // us a total set of groups below 75.  We limit the intervals
  // we try to sensible whole numbers.
  vector<size_t> baseValues = {2,5,10};
  size_t multiplier = 1;
  while (true) {
    for (size_t b = 0; b<baseValues.size(); b++) {
      size_t interval = baseValues[b] * multiplier;
      size_t group_count = 9 + ((num_bases-9)/interval);

      if ((num_bases-9) % interval != 0)
        group_count++;

      if (group_count < 75)
        return interval;
    }
    multiplier *= 10;
    if (multiplier == 10000000)
      throw runtime_error("Couldn't find a sensible interval grouping for len ="
                      + std::to_string(num_bases));
  }
}


void
make_linear_base_groups(vector<BaseGroup> &base_groups,
                        const size_t num_bases) {

  // For lengths below 75bp we just return everything.
  if (num_bases <= 75) {
    make_default_base_groups(base_groups, num_bases);
    return;
  }

  // We need to work out what interval we're going to use.
  const size_t interval = get_linear_interval(num_bases);
  size_t starting_base = 1;

  while (starting_base <= num_bases) {
    size_t end_base = starting_base + interval - 1;

    if (starting_base < 10)
      end_base = starting_base;

    if (starting_base == 10 && interval > 10)
      end_base = interval - 1;

    if (end_base > num_bases)
      end_base = num_bases;

    BaseGroup bg = BaseGroup(starting_base - 1, end_base - 1);
    base_groups.push_back(bg);

    if (starting_base < 10)
      starting_base++;
    else if (starting_base == 10 && interval > 10)
      starting_base = interval;
    else
      starting_base += interval;

  }
}

template<class T>
T& operator << (T &out, const BaseGroup &group) {
  if (group.start == group.end) out << group.start + 1;
  else out << group.start + 1 << "-" << group.end + 1;
  return out;
}

// FastQC extrapolation of counts to the full file size
double get_corrected_count(size_t count_at_limit,
                           size_t num_reads,
                           size_t dup_level,
                           size_t num_obs) {
  // See if we can bail out early (ADS: can we know if num_reads <=
  // count_at_limit always holds?)
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
  return num_obs/std::max(num_lim<double>::min(), 1.0 - p_not_seeing);
}

// Function to calculate the deviation of a histogram with 100 bins from a
// theoretical normal distribution with same mode and standard deviation
double
sum_deviation_from_normal(const array <double, 101> &gc_count,
                          array <double, 101> &theoretical) {
  /******************* BEGIN COPIED FROM FASTQC **********************/
  const size_t num_gc_bins = 101;

  // Sum of all gc counts in all histogram bins
  double total_count = 0.0;

  // We use the mode to calculate the theoretical distribution
  // so that we cope better with skewed distributions.
  size_t first_mode = 0;
  double mode_count = 0.0;

  for (size_t i = 0; i < num_gc_bins; ++i) {
    total_count += gc_count[i];
    if (gc_count[i] > mode_count) {
      mode_count = gc_count[i];
      first_mode = i;
    }
  }

  // The mode might not be a very good measure of the centre
  // of the distribution either due to duplicated vales or
  // several very similar values next to each other.  We therefore
  // average over adjacent points which stay above 95% of the modal
  // value

  double mode = 0;
  size_t mode_duplicates = 0;
  bool fell_off_top = true;

  for (size_t i = first_mode; i < num_gc_bins; ++i) {
    if (gc_count[i] > gc_count[first_mode] - (gc_count[first_mode]/10.0)) {
      mode += i;
      mode_duplicates++;
    }
    else {
      fell_off_top = false;
      break;
    }
  }

  bool fell_off_bottom = true;
  for (int i = first_mode - 1; i >= 0; --i) {
    if (gc_count[i] > gc_count[first_mode]
                          - (gc_count[first_mode]/10.0)) {
      mode += i;
      mode_duplicates++;
    }
    else {
      fell_off_bottom = false;
      break;
    }
  }

  if (fell_off_bottom || fell_off_top) {
    // If the distribution is so skewed that 95% of the mode
    // is off the 0-100% scale then we keep the mode as the
    // centre of the model
    mode = first_mode;
  } else {
    // ADS: check if we need to avoid divide-by-zero here
    mode /= std::max(static_cast<size_t>(1), mode_duplicates);
  }

  // We can now work out a theoretical distribution
  double stdev = 0.0;
  for (size_t i = 0; i < num_gc_bins; ++i) {
    stdev += (i - mode) * (i - mode) * gc_count[i];
  }

  // ADS: check if we need to avoid divide-by-zero here
  stdev = stdev / std::max(num_lim<double>::min(), total_count - 1.0);
  stdev = sqrt(stdev);

  /******************* END COPIED FROM FASTQC **********************/
  // theoretical sampling from a normal distribution with mean = mode and stdev
  // = stdev to the mode from the sampled gc content from the data
  double ans = 0.0, theoretical_sum = 0.0, z;
  theoretical.fill(0);
  // ADS: lonely magic below; what is the 100?
  for (size_t i = 0; i <= 100; ++i) {
    z = i - mode;
    // ADS: check if we need to avoid divide-by-zero here
    theoretical[i] = exp(- (z*z)/ (2.0 * stdev *stdev));
    theoretical_sum += theoretical[i];
  }

  // Normalize theoretical so it sums to the total of readsq
  for (size_t i = 0; i <= 100; ++i) {
    // ADS: check if we need to avoid divide-by-zero here
    theoretical[i] = theoretical[i] * total_count /
      std::max(num_lim<double>::min(), theoretical_sum);
  }

  for (size_t i = 0; i <= 100; ++i) {
    ans += fabs(gc_count[i] - theoretical[i]);
  }
  // Fractional deviation (ADS: check if we need to avoid
  // divide-by-zero here)
  return 100.0 * ans / std::max(num_lim<double>::min(), total_count);
}

/***************************************************************/
/********************* ABSTRACT MODULE *************************/
/***************************************************************/
Module::Module(const string &_module_name) : module_name(_module_name) {
  // make placeholders
  placeholder = module_name;

  // removes spaces
  placeholder.erase(remove_if(placeholder.begin(),
        placeholder.end(), isspace), placeholder.end());

  // lowercases it
  transform(placeholder.begin(),
            placeholder.end(), placeholder.begin(),
          [](unsigned char c){ return std::tolower(c); });

  // makes html placeholders
  placeholder_name = "{{" + placeholder + "name" + "}}";
  placeholder_data = "{{" + placeholder + "data" + "}}";
  placeholder_cs = "{{" + placeholder + "cs" + "}}";
  placeholder_ce = "{{" + placeholder + "ce" + "}}";
  placeholder_grade = "{{pass" + placeholder + "}}";
  grade = "pass";
  summarized = false;
}

Module::~Module() {

}

void
Module::write(ostream &os) {
  if (!summarized)
    throw runtime_error("Attempted to write module before summarizing : "
                        + module_name);
  os << ">>" << module_name << "\t" << grade << "\n";
  write_module(os);
  os << ">>END_MODULE\n";
}

void
Module::write_short_summary(ostream &os, const string &filename) {
  if (!summarized)
    throw runtime_error("Attempted to write module before summarizing : "
                        + module_name);
  os << toupper(grade) << "\t"
     << module_name << "\t"
     << filename << "\n";
}

// Guarantees the summarized flag is only set to true when module
// data has been summarized
void
Module::summarize(FastqStats &stats) {
  summarize_module(stats);
  make_grade();
  html_data = make_html_data();
  summarized = true;
}

/***************************************************************/
/********************* SUMMARIZE FUNCTIONS *********************/
/***************************************************************/

/******************* BASIC STATISTICS **************************/
const string ModuleBasicStatistics::module_name = "Basic Statistics";
ModuleBasicStatistics::
ModuleBasicStatistics(const FalcoConfig &config)
: Module(ModuleBasicStatistics::module_name) {
  is_nanopore = config.nanopore;
  filename_stripped = config.filename_stripped;
}

void
ModuleBasicStatistics::summarize_module(FastqStats &stats) {
  // Total sequences
  total_sequences = stats.num_reads;

  // min and max read length
  min_read_length = stats.min_read_length;
  max_read_length = stats.max_read_length;

  // These seem to always be the same on FastQC
  // File type
  file_type = "Conventional base calls";
  if (is_nanopore) {
    file_encoding = "Oxford Nanopore";
    //stats.encoding_offset = ?
  }
  else {
    // File encoding
    const char lowest_char = stats.lowest_char;

    //copied from FastQC:
    static const size_t SANGER_ENCODING_OFFSET = 33;
    static const char ILLUMINA_1_3_ENCODING_OFFSET = 64;
    if (lowest_char < 33) {
      throw runtime_error("No known encoding with chars < 33. Yours was " +
                          std::to_string(lowest_char) + ")");
    }
    else if (lowest_char < 64) {
      file_encoding = "Sanger / Illumina 1.9";
      stats.encoding_offset = SANGER_ENCODING_OFFSET - Constants::quality_zero;
    }
    else if (lowest_char == ILLUMINA_1_3_ENCODING_OFFSET+1) {
      file_encoding = "Illumina 1.3";
      stats.encoding_offset = ILLUMINA_1_3_ENCODING_OFFSET - Constants::quality_zero;
    }
    else if (lowest_char <= 126) {
      file_encoding = "Illumina 1.5";
      stats.encoding_offset = ILLUMINA_1_3_ENCODING_OFFSET - Constants::quality_zero;
    }

    // found a char but it's not the default one
    else if (lowest_char != stats.lowest_char) {
      throw runtime_error("No known encodings with chars > 126 (Yours was " +
                          std::to_string(lowest_char) + ")");
    }
  }

  // Poor quality reads
  num_poor = 0;

  // Average read length
  avg_read_length = 0;
  size_t total_bases = 0;
  for (size_t i = 0; i < max_read_length; ++i) {
    if (i < FastqStats::SHORT_READ_THRESHOLD)
      total_bases += i * stats.read_length_freq[i];
    else
      total_bases += i * stats.long_read_length_freq[i - FastqStats::SHORT_READ_THRESHOLD];
  }

  avg_read_length =
    total_bases / std::max(static_cast<size_t>(1), total_sequences);

  // counts bases G and C in each base position
  avg_gc = 0;

  // GC %
  // GS: TODO delete gc calculation during stream and do it using the total G
  // counts in all bases
  avg_gc = 100 * stats.total_gc / std::max(1.0, static_cast<double>(total_bases));

}

// It's always a pass
void
ModuleBasicStatistics::make_grade() {
}

void
ModuleBasicStatistics::write_module(ostream &os) {
  os << "#Measure\tValue\n";
  os << "Filename\t" << filename_stripped << "\n";
  os << "File type\t" << file_type << "\n";
  os << "Encoding\t" << file_encoding << "\n";
  os << "Total Sequences\t" << total_sequences << "\n";
  os << "Sequences flagged as poor quality\t" << num_poor << "\n";
  os << "Sequence length\t";
  if (min_read_length == max_read_length) {
    os << min_read_length;
  } else {
    os << min_read_length << "-" << max_read_length;
  }
  os << "\n";
  os << "%GC\t" << static_cast<size_t>(avg_gc) << "\n";
}

string
ModuleBasicStatistics::make_html_data() {
  ostringstream data;
  data << "<table><thead><tr><th>Measure</th><th>Value"
       << "</th></tr></thead><tbody>";
  data << "<tr><td>Filename</td><td>" << filename_stripped
       << "</td></tr>";
  data << "<tr><td>File type</td><td>" << file_type
       << "</td></tr>";
  data << "<tr><td>Encoding</td><td>" << file_encoding
       << "</td></tr>";
  data << "<tr><td>Total Sequences</td><td>" << total_sequences << "</td></tr>";
  data << "<tr><td>Sequences Flagged As Poor Quality</td><td>"
       << num_poor << "</td></tr>";
  data << "<tr><td>Sequence length</td><td>";
  if (min_read_length != max_read_length) {
    data << min_read_length << " - " << max_read_length;
  }
  else {
    data << max_read_length;
  }
  data << "</td></tr>";
  data << "<tr><td>%GC:</td><td>" << avg_gc << "</td></tr>";
  data << "</tbody></table>";
  return data.str();
}

void
ModuleBasicStatistics::read_data_line(const std::string &line) {
  string lhs,rhs;
  istringstream iss(line);

  // get text before and after tab
  getline(iss, lhs, '\t');
  getline(iss, rhs, '\t');

  if (lhs ==  "Filename")
    filename_stripped = rhs;
  else if (lhs ==  "File type")
    file_type = rhs;
  else if (lhs ==  "Encoding")
    file_encoding = rhs;
  else if (lhs ==  "Total Sequences")
    total_sequences = atoi(rhs.c_str());
  else if (lhs ==  "Sequences flagged as poor quality")
    num_poor = atoi(rhs.c_str());

  else if (lhs ==  "Sequence length") {
    // non-constant sequence length
    if (rhs.find("-") != string::npos) {
      istringstream seq_iss (rhs);
      string min_l, max_l;
      getline(seq_iss, min_l, '-');
      getline(seq_iss, max_l, '-');
      min_read_length = atoi(min_l.c_str());
      max_read_length = atoi(max_l.c_str());
    }
  }
  else if (lhs == "%GC")
    avg_gc = atoi(rhs.c_str());
  else {
    throw runtime_error("malformed basic statistic" + lhs);
  }
}

/******************* PER BASE SEQUENCE QUALITY **********************/
const string
ModulePerBaseSequenceQuality::module_name = "Per base sequence quality";

ModulePerBaseSequenceQuality::ModulePerBaseSequenceQuality
(const FalcoConfig &config):
Module(ModulePerBaseSequenceQuality::module_name){
  auto base_lower = config.limits.find("quality_base_lower");
  auto base_median = config.limits.find("quality_base_median");

  base_lower_warn = (base_lower->second).find("warn")->second;
  base_lower_error = (base_lower->second).find("error")->second;
  base_median_warn = (base_median->second).find("warn")->second;
  base_median_error = (base_median->second).find("error")->second;

  do_group = !config.nogroup;
}

void
ModulePerBaseSequenceQuality::summarize_module(FastqStats &stats) {
  // Quality quantiles for base positions
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

  size_t ldecile_group_sum = 0,
         lquartile_group_sum = 0,
         median_group_sum = 0,
         uquartile_group_sum = 0,
         udecile_group_sum = 0;
  double mean_group_sum;

  size_t cur;
  size_t cur_sum;
  size_t counts;
  double cur_mean;
  num_bases = stats.max_read_length;

  // first, do the groups
  if (do_group) {
    make_linear_base_groups(base_groups, num_bases);
  }

  else make_default_base_groups(base_groups, num_bases);
  num_groups = base_groups.size();

  // Reserves space I know I will use
  group_mean = vector<double>(num_groups, 0.0);
  group_ldecile = vector<double>(num_groups, 0.0);
  group_lquartile = vector<double>(num_groups, 0.0);
  group_median = vector<double>(num_groups, 0.0);
  group_uquartile = vector<double>(num_groups, 0.0);
  group_udecile = vector<double>(num_groups, 0.0);

  // temp
  vector <size_t> histogram(128, 0);
  size_t bases_in_group = 0;

  for (size_t group = 0; group < num_groups; ++group) {
    mean_group_sum = 0;
    ldecile_group_sum = 0;
    lquartile_group_sum = 0;
    median_group_sum = 0;
    uquartile_group_sum = 0;
    udecile_group_sum = 0;

    // Find quantiles for each base group
    for (size_t i = base_groups[group].start;
              i  <= base_groups[group].end; ++i) {

      // reset group values
      bases_in_group = 0;
      for (size_t j = 0; j < 128; ++j)
        histogram[j] = 0;

      for (size_t j = 0; j < FastqStats::kNumQualityValues; ++j) {
        // get value
        if (i < FastqStats::SHORT_READ_THRESHOLD) {
          cur = stats.position_quality_count[
            (i << FastqStats::kBitShiftQuality) | j];
        }
        else {
          cur = stats.long_position_quality_count [
            ((i - FastqStats::SHORT_READ_THRESHOLD) << FastqStats::kBitShiftQuality) | j];
        }

        // Add to Phred histogram
        histogram[j] += cur;
      }

      // Number of bases seen in position i
      if (i < FastqStats::SHORT_READ_THRESHOLD) {
        bases_in_group += stats.cumulative_read_length_freq[i];
      } else {
        bases_in_group +=
          stats.long_cumulative_read_length_freq[i - FastqStats::SHORT_READ_THRESHOLD];
      }

      ldecile_thresh = 0.1 * bases_in_group;
      lquartile_thresh = 0.25 * bases_in_group;
      median_thresh = 0.5 * bases_in_group;
      uquartile_thresh = 0.75 * bases_in_group;
      udecile_thresh = 0.9 * bases_in_group;

      // now go again through the counts in each quality value to find the
      // quantiles
      cur_sum = 0;
      counts = 0;

      for (size_t j = 0; j < FastqStats::kNumQualityValues; ++j) {
        // Finds in which bin of the histogram reads are
        cur = histogram[j];
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
        cur_sum += cur * j;
        counts += cur;
      }

      cur_mean = static_cast<double>(cur_sum) /
                 static_cast<double>(bases_in_group);
      const size_t offset = stats.encoding_offset;
      mean_group_sum += cur_mean - offset;
      ldecile_group_sum += cur_ldecile - offset;
      lquartile_group_sum += cur_lquartile - offset;
      median_group_sum += cur_median - offset;
      uquartile_group_sum += cur_uquartile - offset;
      udecile_group_sum += cur_udecile - offset;
    }

    const size_t base_positions = base_groups[group].end - base_groups[group].start + 1;
    assert(base_positions != static_cast<size_t>(0));
    group_mean[group] = mean_group_sum / base_positions;
    group_ldecile[group] = static_cast<double>(ldecile_group_sum) / base_positions;
    group_lquartile[group] = static_cast<double>(lquartile_group_sum) / base_positions;
    group_median[group] = static_cast<double>(median_group_sum) / base_positions;
    group_uquartile[group] = static_cast<double>(uquartile_group_sum) / base_positions;
    group_udecile[group] = static_cast<double>(udecile_group_sum) / base_positions;
  }
}

void
ModulePerBaseSequenceQuality::make_grade() {
  num_warn = 0;
  num_error = 0;
  for (size_t i = 0; i < num_groups; ++i) {
    // there was enough data to make this assessment
    if (group_lquartile[i] > 0) {
      if (grade != "fail") {
        if (group_lquartile[i] < base_lower_error ||
            group_median[i] < base_median_error) {
          num_error++;
        } else if (group_lquartile[i] < base_lower_warn ||
                   group_median[i] < base_median_warn) {
          num_warn++;
        }
      }
    }
  }
  if (num_error > 0)
    grade = "fail";
  else if (num_warn > 0)
    grade = "warn";
}

void
ModulePerBaseSequenceQuality::write_module(ostream &os) {
  os << "#Base\tMean\tMedian\tLower Quartile\tUpper Quartile" <<
        "\t10th Percentile\t90th Percentile\n";

  // GS: TODO make base groups
  for (size_t i = 0; i < num_groups; ++i) {
    os << base_groups[i] << "\t"
       << group_mean[i] << "\t"
       << group_median[i] << "\t"
       << group_lquartile[i] << "\t"
       << group_uquartile[i] << "\t"
       << group_ldecile[i] << "\t"
       << group_udecile[i] << "\n";
  }
}

// Plotly data
string
ModulePerBaseSequenceQuality::make_html_data() {
  ostringstream data;
  for (size_t i = 0; i < num_groups; ++i) {
    data << "{y : [";

      data << group_ldecile[i] << ", "
           << group_lquartile[i] << ", "
           << group_median[i] << ", "
           << group_uquartile[i] << ", "
           << group_udecile[i] << "], ";
    data << "type : 'box', name : ' ";
    if (base_groups[i].start == base_groups[i].end)
      data << base_groups[i].start + 1;
    else
      data << base_groups[i].start + 1 << "-" << base_groups[i].end + 1;
    data << "bp', ";
    data << "marker : {color : '";

    // I will color the boxplot based on whether it passed or failed
    if (group_median[i] < base_median_error ||
        group_lquartile[i] < base_lower_error)
      data << "red";
    else if (group_median[i] < base_median_warn ||
             group_lquartile[i] < base_lower_warn)
      data << "yellow";
    else
      data << "green";
    data << "'}}";
    if (i < num_bases - 1) {
      data << ", ";
    }
  }
  return data.str();
}

void
ModulePerBaseSequenceQuality::read_data_line(const string &line) {

}

/************** PER TILE SEQUENCE QUALITY ********************/
const string
ModulePerTileSequenceQuality::module_name = "Per tile sequence quality";
ModulePerTileSequenceQuality::
ModulePerTileSequenceQuality(const FalcoConfig &config) :
Module(ModulePerTileSequenceQuality::module_name) {
  auto grade_tile = config.limits.find("tile")->second;
  grade_warn = grade_tile.find("warn")->second;
  grade_error = grade_tile.find("error")->second;
}

void
ModulePerTileSequenceQuality::summarize_module(FastqStats &stats) {
  max_read_length = stats.max_read_length;
  tile_position_quality = stats.tile_position_quality;
  // First I calculate the number of counts for each position
  vector<size_t> position_counts(max_read_length, 0);

  for (auto v : stats.tile_position_quality) {
    for (size_t i = 0; i < v.second.size(); ++i) {
      position_counts[i] +=
        stats.tile_position_count.find(v.first)->second[i];
    }
  }

  // Now I calculate the sum of all tile qualities in each position
  vector<double> mean_in_base(max_read_length, 0.0);
  for (auto v : tile_position_quality) {
    for(size_t i = 0; i < v.second.size(); ++i) {
      mean_in_base[i] += v.second[i];
    }
  }

  // Now transform sum into mean
  for (size_t i = 0; i < max_read_length; ++i)
    if (position_counts[i] > 0.0)
      mean_in_base[i] = mean_in_base[i] / position_counts[i];
    else
      mean_in_base[i] = 0.0;

  for (auto &v : tile_position_quality) {
    const size_t lim = v.second.size();
    for (size_t i = 0; i < lim; ++i) {
      // transform sum of all qualities in mean
      const auto itr = stats.tile_position_count.find(v.first);
      if (itr == cend(stats.tile_position_count))
        throw runtime_error("failure ModulePerTileSequenceQuality::summarize_module");
      const size_t count_at_pos = itr->second[i];

      if (count_at_pos > 0)
        v.second[i] = v.second[i] / count_at_pos;

      // subtract the global mean
      v.second[i] -= mean_in_base[i];
    }
  }
  // sorts tiles
  tiles_sorted.clear();
  for (auto v : tile_position_quality)
    tiles_sorted.push_back(v.first);

  sort(tiles_sorted.begin(), tiles_sorted.end());
}

void
ModulePerTileSequenceQuality::make_grade() {
  grade = "pass";
  for (auto &v : tile_position_quality) {
    for (size_t i = 0; i < v.second.size(); ++i) {
      if (grade != "fail") {
        if (v.second[i] <= -grade_error) {
          grade = "fail";
        }
        else if (v.second[i] <= -grade_warn) {
          grade = "warn";
        }
      }
    }
  }
}

void
ModulePerTileSequenceQuality::write_module(ostream &os) {

  // prints tiles sorted by value
  os << "#Tile\tBase\tMean\n";
  for (size_t i = 0; i < tiles_sorted.size(); ++i) {
    for (size_t j = 0; j < max_read_length; ++j) {

      if (tile_position_quality[tiles_sorted[i]].size() >= j) {
        os << tiles_sorted[i] << "\t" << j + 1 << "\t"
           << tile_position_quality[tiles_sorted[i]][j];
        os << "\n";
      }
    }
  }
}

inline double
round_quantile(const double val, const double num_quantiles) {
  // ADS: check if we need to worry about divide by zero here
  return static_cast<int>(val * num_quantiles) / num_quantiles;
}

string
ModulePerTileSequenceQuality::make_html_data() {
  // find quantiles based on data for a standardized color scale
  double min_val = std::numeric_limits<double>::max();
  double max_val = std::numeric_limits<double>::min();

  ostringstream data;
  data << "{x : [";
  for (size_t i = 0; i < max_read_length; ++i) {
    data << i+1;
    if (i < max_read_length - 1)
      data << ",";
  }

  // Y : Tile
  data << "], y: [";
  bool first_seen = false;
  for (size_t i = 0; i < tiles_sorted.size(); ++i) {
    if (!first_seen) first_seen = true;
    else data << ",";
    data << tiles_sorted[i];
  }

  // Z: quality z score
  data << "], z: [";
  first_seen = false;
  for (size_t i = 0; i < tiles_sorted.size(); ++i) {
    if (!first_seen) first_seen = true;
    else data << ", ";

    // start new array with all counts
    data << "[";
    for (size_t j = 0; j < max_read_length; ++j) {
      const double val =  tile_position_quality[tiles_sorted[i]][j];

      data <<  val;
      if (j < max_read_length - 1) data << ",";

      max_val = max(max_val, val);
      min_val = min(min_val, val);
    }
    data << "]";
  }
  data << "]";
  data << ", type : 'heatmap',";

  // fixed color scale
  data << "colorscale: [";

  // We will now discretize the quantiles so plotly understands
  // the color scheme
  static const double num_quantiles = 20.0;
  // ADS: not sure if we need to worry about divide by zero here?
  double mid_point = round_quantile(min_val/(min_val - max_val), num_quantiles);

  // - 10: red
  data << "[0.0, 'rgb(210,65,83)'],";

  // 0: light blue
  data << "[" << mid_point << ", 'rgb(178,236,254)'],";

  // + 10: dark blue
  data << "[1.0, 'rgb(34,57,212)']";
  data << "],";
  data << "showscale : true";
  data << "}";

  return data.str();
}

/******************* PER SEQUENCE QUALITY SCORE **********************/
const string
ModulePerSequenceQualityScores::module_name = "Per sequence quality scores";
ModulePerSequenceQualityScores::
ModulePerSequenceQualityScores(const FalcoConfig &config) :
Module(ModulePerSequenceQualityScores::module_name) {
  mode_val = 0;
  mode_ind = 0;
  offset = 0;

  auto mode_limits = config.limits.find("quality_sequence");
  mode_warn = (mode_limits->second).find("warn")->second;
  mode_error = (mode_limits->second).find("error")->second;
}

void
ModulePerSequenceQualityScores::summarize_module(FastqStats &stats) {
  // Need to copy this to write later
  quality_count = stats.quality_count;
  offset = stats.encoding_offset;

  // get mode for grade
  for (size_t i = 0; i < FastqStats::kNumQualityValues; ++i) {
    if (stats.quality_count[i] > mode_val) {
      mode_val = stats.quality_count[i];
      mode_ind = i - offset;
    }
  }
}

void
ModulePerSequenceQualityScores::make_grade() {
  if (mode_ind < mode_warn) {
    grade = "warn";
  }

  if (mode_ind < mode_error) {
    grade = "fail";
  }
}

void
ModulePerSequenceQualityScores::write_module(ostream &os) {
  os << "#Quality\tCount\n";
  for (size_t i = 0; i < FastqStats::kNumQualityValues; ++i) {
    if (quality_count[i] > 0)
      os << i - offset << "\t" << quality_count[i] << "\n";
  }
}

string
ModulePerSequenceQualityScores::make_html_data() {
  ostringstream data;
  data << "{x : [";
  bool seen_first = false;
  for (size_t i = 0; i < FastqStats::kNumQualityValues; ++i) {
    if (quality_count[i] > 0){
      if (seen_first)
        data << ",";
      else
        seen_first = true;
      data << i - offset;
    }
  }

  // Y values: frequency with which they were seen
  data << "], y : [";
  seen_first = false;
  for (size_t i = 0; i < FastqStats::kNumQualityValues; ++i) {
    if (quality_count[i] > 0){
      if (seen_first)
        data << ",";
      else
        seen_first = true;
      data << quality_count[i];
    }
  }
  data << "], type: 'line', line : {color : 'red'}, "
       << "name : 'Sequence quality distribution'}";

  return data.str();
}

/******************* PER BASE SEQUENCE CONTENT **********************/
const string
ModulePerBaseSequenceContent::module_name = "Per base sequence content";
ModulePerBaseSequenceContent::
ModulePerBaseSequenceContent(const FalcoConfig &config) :
Module(ModulePerBaseSequenceContent::module_name) {
  auto sequence_limits = config.limits.find("sequence")->second;
  sequence_warn = sequence_limits.find("warn")->second;
  sequence_error = sequence_limits.find("error")->second;
  is_bisulfite = config.is_bisulfite;
  is_reverse_complement = config.is_reverse_complement;
  do_group = !config.nogroup;
}

void
ModulePerBaseSequenceContent::summarize_module(FastqStats &stats) {
  double a_group, t_group, g_group, c_group, n_group;
  double a_pos{}, t_pos{}, g_pos{}, c_pos{}, n_pos{};
  double total; //a+c+t+g+n
  max_diff = 0.0;

  num_bases = stats.max_read_length;

  if (do_group) {
    make_linear_base_groups(base_groups, num_bases);
  }

  else make_default_base_groups(base_groups, num_bases);
  num_groups = base_groups.size();

  a_pct = vector<double>(num_groups, 0.0);
  c_pct = vector<double>(num_groups, 0.0);
  g_pct = vector<double>(num_groups, 0.0);
  t_pct = vector<double>(num_groups, 0.0);

  for (size_t group = 0; group < num_groups; ++group) {
    // Find quantiles for each base group
    a_group = c_group = g_group = t_group = n_group = 0.0;
    for (size_t i = base_groups[group].start;
              i  <= base_groups[group].end; ++i) {
      if (i < FastqStats::SHORT_READ_THRESHOLD) {
        a_pos = stats.base_count[(i << FastqStats::kBitShiftNucleotide)];
        c_pos = stats.base_count[(i << FastqStats::kBitShiftNucleotide) | 1];
        t_pos = stats.base_count[(i << FastqStats::kBitShiftNucleotide) | 2];
        g_pos = stats.base_count[(i << FastqStats::kBitShiftNucleotide) | 3];
        n_pos = stats.n_base_count[i];
      } else {
        a_pos = stats.long_base_count[
              ((i - FastqStats::SHORT_READ_THRESHOLD) << FastqStats::kBitShiftNucleotide)
            ];
        c_pos = stats.long_base_count[
              ((i - FastqStats::SHORT_READ_THRESHOLD) << FastqStats::kBitShiftNucleotide) | 1
            ];
        t_pos = stats.long_base_count[
              ((i - FastqStats::SHORT_READ_THRESHOLD) << FastqStats::kBitShiftNucleotide) | 2
            ];
        g_pos = stats.long_base_count[
              ((i - FastqStats::SHORT_READ_THRESHOLD) << FastqStats::kBitShiftNucleotide) | 3
            ];
        n_pos = stats.long_n_base_count[
                i - FastqStats::SHORT_READ_THRESHOLD
            ];
      }
      a_group += a_pos; c_group += c_pos; g_group += g_pos; t_group += t_pos;
      n_group += n_pos;

      const double total_pos =
        static_cast<double>(a_pos + c_pos + g_pos + t_pos + n_pos);
      a_pos = 100.0 * a_pos / std::max(num_lim<double>::min(), total_pos);
      c_pos = 100.0 * c_pos / std::max(num_lim<double>::min(), total_pos);
      g_pos = 100.0 * g_pos / std::max(num_lim<double>::min(), total_pos);
      t_pos = 100.0 * t_pos / std::max(num_lim<double>::min(), total_pos);

      // for WGBS, we only test non-bisulfite treated bases
      if (!is_reverse_complement)
        max_diff = max(max_diff, fabs(a_pos - g_pos));
      else
        max_diff = max(max_diff, fabs(c_pos - t_pos));

      if (!is_bisulfite) {
        max_diff = max(max_diff, fabs(a_pos - c_pos));
        max_diff = max(max_diff, fabs(a_pos - t_pos));
        max_diff = max(max_diff, fabs(c_pos - g_pos));
        max_diff = max(max_diff, fabs(t_pos - g_pos));

        if (!is_reverse_complement)
          max_diff = max(max_diff, fabs(c_pos - t_pos));
        else
          max_diff = max(max_diff, fabs(a_pos - g_pos));
      }
      // WGBS specific base content count
      else {
        max_diff = max(max_diff, fabs((c_pos + t_pos) - (a_pos + g_pos)));
      }
    }

    // turns above values to percent
    total = static_cast<double>(a_group + c_group + t_group + g_group + n_group);
    a_pct[group] = 100.0*a_group / std::max(num_lim<double>::min(), total);
    c_pct[group] = 100.0*c_group / std::max(num_lim<double>::min(), total);
    g_pct[group] = 100.0*g_group / std::max(num_lim<double>::min(), total);
    t_pct[group] = 100.0*t_group / std::max(num_lim<double>::min(), total);
  }
}

void
ModulePerBaseSequenceContent::make_grade() {
  if (max_diff > sequence_error) {
    grade = "fail";
  }
  else if (max_diff > sequence_warn) {
    grade = "warn";
  }
}

void
ModulePerBaseSequenceContent::write_module(ostream &os) {
  os << "#Base\tG\tA\tT\tC";
  if (is_bisulfite)
    os << "\tC+T\tA+G";

  os << "\n";
  for (size_t i = 0; i < num_groups; ++i) {
    os << base_groups[i] << "\t" <<
          g_pct[i] << "\t" <<
          a_pct[i] << "\t" <<
          t_pct[i] << "\t" <<
          c_pct[i];

    if (is_bisulfite) {
      os << "\t" << c_pct[i]+t_pct[i] << "\t"
         << a_pct[i]+g_pct[i];
    }
    os << "\n";
  }
}

string
ModulePerBaseSequenceContent::make_html_data() {
  ostringstream data;
  // ATGC
  for (size_t base = 0; base < 4; ++base) {
    // start line
    data << "{";

    // X values : base position
    data << "x : [";
    for (size_t i = 0; i < num_groups; ++i) {
      if (base_groups[i].start == base_groups[i].end)
        data << base_groups[i].start + 1;
      else
        data << "\"" << base_groups[i].start + 1 << "-"
                     << base_groups[i].end + 1 << "\"";
      if (i < num_groups - 1)
        data << ", ";
    }

    // Y values: frequency with which they were seen
    data << "], y : [";
    for (size_t i = 0; i < num_groups; ++i) {
      if (base == 0) data << a_pct[i];
      if (base == 1) data << c_pct[i];
      if (base == 2) data << t_pct[i];
      if (base == 3) data << g_pct[i];
      if (i < num_groups - 1)
        data << ", ";
    }
    data << "], mode : 'lines', name : '" + size_t_to_seq(base, 1) + "', ";

    // color
    data << "line :{ color : '";
    if (base == 0)
      data << "green";
    else if (base == 1)
      data << "blue";
    else if (base == 2)
      data << "red";
    else
      data << "black";
    data << "'}";
    // end color

    // end line
    data << "}";
    if (base < 4)
      data << ", ";
  }

  // bisulfite dashed lines
  if (is_bisulfite) {
    for (size_t line = 0; line <= 1; ++line) {
      data << "{";
      data << "x : [";
      for (size_t i = 0; i < num_groups; ++i) {
        if (base_groups[i].start == base_groups[i].end)
          data << base_groups[i].start + 1;
        else
          data << "\"" << base_groups[i].start + 1 << "-"
                       << base_groups[i].end + 1 << "\"";
        if (i < num_groups - 1)
          data << ", ";
      }

      // Y values: frequency with which they were seen
      data << "], y : [";
      for (size_t i = 0; i < num_groups; ++i) {
        if (line == 0)
          data << a_pct[i] + g_pct[i];
        else
          data << c_pct[i] + t_pct[i];

        if (i < num_groups - 1)
          data << ", ";
      }
      data << "], mode : 'lines', name : '";
      if (line == 0) {
        data << "A+G";
      } else {
        data << "C+T";
      }
      data <<  "', ";

      // color and dash
      data << "line :{ color : '";
      if (line == 0)
        data << "#CCCCCC";
      else
        data << "#999999";
      data << "', dash : 'dash'}";
      // end color

      // end line
      data << "}";
      if (line == 0)
        data << ", ";

    }
  }

  return data.str();
}

/******************* PER SEQUENCE GC CONTENT *****************/
const string
ModulePerSequenceGCContent::module_name = "Per sequence GC content";
ModulePerSequenceGCContent::
ModulePerSequenceGCContent(const FalcoConfig &config) :
Module(ModulePerSequenceGCContent::module_name) {
  auto gc_vars = config.limits.find("gc_sequence")->second;
  gc_warn = gc_vars.find("warn")->second;
  gc_error = gc_vars.find("error")->second;
}

void
ModulePerSequenceGCContent::summarize_module(FastqStats &stats) {
  gc_count = stats.gc_count;
  gc_deviation = sum_deviation_from_normal(gc_count,
                                           theoretical_gc_count);
}

void
ModulePerSequenceGCContent::make_grade() {
  if (gc_deviation >= gc_error) {
    grade = "fail";
  }
  else if (gc_deviation >= gc_warn) {
    grade = "warn";
  }
}

void
ModulePerSequenceGCContent::write_module(ostream &os) {
  os << "#GC Content\tCount\n";
  for (size_t i = 0; i <= 100; ++i)
    os << i << "\t" << gc_count[i] << "\n";
}

string
ModulePerSequenceGCContent::make_html_data() {
  ostringstream data;
  // Actual count
  data << "{x : [";
  for (size_t i = 0; i < 101; ++i) {
    data << i + 1;
    if (i < 101)
      data << ", ";
  }

  // Y values: frequency with which they were seen
  data << "], y : [";
  for (size_t i = 0; i < 101; ++i) {
    data << gc_count[i];
    if (i < 101)
      data << ", ";
  }
  data << "], type: 'line', line : {color : 'red'},name : 'GC distribution'}";

  // Theoretical count
  data << ", {x : [";
  for (size_t i = 0; i < 101; ++i) {
    data << i + 1;
    if (i < 101)
      data << ", ";
  }

  // Y values: frequency with which they were seen
  data << "], y : [";
  for (size_t i = 0; i < 101; ++i) {
    data << theoretical_gc_count[i];
    if (i < 101)
      data << ", ";
  }
  data << "], type: 'line', line : {color : 'blue'},"
       << "name : 'Theoretical distribution'}";

  return data.str();
}

/******************* PER BASE N CONTENT **********************/
const string
ModulePerBaseNContent::module_name = "Per base N content";
ModulePerBaseNContent::
ModulePerBaseNContent(const FalcoConfig &config) :
Module(ModulePerBaseNContent::module_name) {
  auto grade_n = config.limits.find("n_content")->second;
  grade_n_warn = grade_n.find("warn")->second;
  grade_n_error = grade_n.find("error")->second;
  do_group = !config.nogroup;
}

void
ModulePerBaseNContent::summarize_module(FastqStats &stats) {
  num_bases = stats.max_read_length;

  // first, do the groups
  if (do_group) {
    make_linear_base_groups(base_groups, num_bases);
  }

  else make_default_base_groups(base_groups, num_bases);
  num_groups = base_groups.size();

  n_pct = vector<double>(num_groups, 0.0);
  size_t this_n_cnt = 0;
  size_t this_n_total = 0;
  double this_n_pct = 0.0;

  max_n_pct = 0.0;
  for (size_t group = 0; group < num_groups; ++group) {
    size_t group_n_cnt = 0, group_n_total = 0;
    for (size_t i = base_groups[group].start;
                i <= base_groups[group].end; ++i) {
      this_n_cnt = (i < FastqStats::SHORT_READ_THRESHOLD) ?
        (stats.n_base_count[i]) : (stats.long_n_base_count[i - FastqStats::SHORT_READ_THRESHOLD]);

      this_n_total = (i < FastqStats::SHORT_READ_THRESHOLD) ? (stats.cumulative_read_length_freq[i]) :
                     (stats.long_cumulative_read_length_freq[i - FastqStats::SHORT_READ_THRESHOLD]);
      this_n_pct = this_n_cnt / std::max(num_lim<double>::min(),
                                         static_cast<double>(this_n_total));
      max_n_pct = max(max_n_pct, this_n_pct);
      group_n_cnt += this_n_cnt;
      group_n_total += this_n_total;
    }
    n_pct[group] = 100.0*group_n_cnt / std::max(num_lim<double>::min(),
                                                static_cast<double>(group_n_total));
  }
}

void
ModulePerBaseNContent::make_grade() {
  if(grade != "fail") {
    if (max_n_pct > grade_n_error) {
      grade = "fail";
    }

    else if (max_n_pct > grade_n_warn) {
      grade = "warn";
    }
  }
}

void
ModulePerBaseNContent::write_module(ostream &os) {
  os << "#Base\tN-Count\n";
  for (size_t i = 0; i < num_groups; ++i) {
    os << base_groups[i] << "\t" << n_pct[i] << "\n";
  }
}

string
ModulePerBaseNContent::make_html_data() {
  ostringstream data;
  // base position
  data << "{x : [";
  for (size_t i = 0; i < num_groups; ++i) {
    if (base_groups[i].start == base_groups[i].end)
      data << "\"" << base_groups[i].start + 1 << "\"";
    else
      data << "\"" << base_groups[i].start + 1 << "-" <<
                      base_groups[i].end + 1 << "\"";
    if (i < num_groups - 1)
      data << ", ";
  }

  // Y values: frequency with which they were seen
  data << "], y : [";
  for (size_t i = 0; i < num_groups; ++i) {
    data << n_pct[i];

    if (i < num_groups - 1)
      data << ", ";
  }
  data << "], type: 'line', line : {color : 'red'}, "
       << "name : 'Fraction of N reads per base'}";

  return data.str();
}

/************** SEQUENCE LENGTH DISTRIBUTION *****************/
const string
ModuleSequenceLengthDistribution::module_name = "Sequence Length Distribution";
ModuleSequenceLengthDistribution::
ModuleSequenceLengthDistribution(const FalcoConfig &config) :
Module(ModuleSequenceLengthDistribution::module_name) {
  auto length_grade = config.limits.find("sequence_length")->second;
  do_grade_error = (length_grade.find("error")->second != 0);
  do_grade_warn = (length_grade.find("warn")->second != 0);
}

void
ModuleSequenceLengthDistribution::summarize_module(FastqStats &stats) {
  max_read_length = stats.max_read_length;

  empty_reads = stats.empty_reads;
  is_all_same_length = true;
  // store the read lengths
  sequence_lengths = vector<size_t>(max_read_length, 0);

  size_t num_nonzero = 0;
  for (size_t i = 0; i < max_read_length; ++i) {
    if (i < FastqStats::SHORT_READ_THRESHOLD) {
      sequence_lengths[i] = stats.read_length_freq[i];
    } else {
      sequence_lengths[i] = stats.long_read_length_freq[
                              i - FastqStats::SHORT_READ_THRESHOLD
                            ];
    }

    if (sequence_lengths[i] > 0) {
      num_nonzero++;
      if (num_nonzero > 1)
        is_all_same_length = false;
    }
  }
}

void
ModuleSequenceLengthDistribution::make_grade() {
  if (do_grade_warn) {
    if (!is_all_same_length) {
      grade = "warn";
    }
  }
  if (do_grade_error) {
    if (empty_reads > 0) {
      grade = "fail";
    }
  }
}

void
ModuleSequenceLengthDistribution::write_module(ostream &os) {
  os << "#Length\tCount\n";
  if (empty_reads > 0)
    os << "0\t" << empty_reads << ".0\n";
  for (size_t i = 0; i < max_read_length; ++i) {
    if (sequence_lengths[i] > 0) {
      os << i+1 << "\t" << sequence_lengths[i] << ".0\n";
    }
  }
}

string
ModuleSequenceLengthDistribution::make_html_data() {
  ostringstream data;
  // X values : avg quality phred scores
  data << "{x : [";
  bool first_seen = false;

  if (empty_reads > 0) {
    first_seen = true;
    data << "\"0 bp\"";
  }
  for (size_t i = 0; i < max_read_length; ++i) {
    if (sequence_lengths[i] > 0) {
      if (first_seen)
        data << ",";
      data << "\"" << i+1 << " bp\"";
      first_seen = true;
    }
  }

  // Y values: frequency with which they were seen
  data << "], y : [";
  first_seen = false;

  if (empty_reads > 0) {
    first_seen = true;
    data << empty_reads;
  }
  for (size_t i = 0; i < max_read_length; ++i) {
    if (sequence_lengths[i] > 0) {
      if (first_seen)
        data << ",";
      data << sequence_lengths[i];
      first_seen = true;
    }
  }

  // Put the sequence value in the text
  data << "], text : [";
  first_seen = false;
  for (size_t i = 0; i < max_read_length; ++i) {
    if (sequence_lengths[i] > 0) {
      if (first_seen)
        data << ",";
      data << i+1;
      first_seen = true;
    }
  }

  data << "], type: 'bar', marker : {color : 'rgba(55,128,191,1.0)',"
       << "line : {width : 2}}, "
       << "name : 'Sequence length distribution'}";

  return data.str();
}

/************** DUPLICATE SEQUENCES **************************/
const string
ModuleSequenceDuplicationLevels::module_name = "Sequence Duplication Levels";
ModuleSequenceDuplicationLevels::
ModuleSequenceDuplicationLevels(const FalcoConfig &config) :
Module(ModuleSequenceDuplicationLevels::module_name) {
  percentage_deduplicated.fill(0);
  percentage_total.fill(0);
  auto grade_dup = config.limits.find("duplication")->second;
  grade_dup_warn = grade_dup.find("warn")->second;
  grade_dup_error = grade_dup.find("error")->second;
}

void
ModuleSequenceDuplicationLevels::summarize_module(FastqStats &stats) {
    seq_total = 0.0;
    seq_dedup = 0.0;

    // Key is frequenccy (r), value is number of times we saw a sequence
    // with that frequency (Nr)
    for (auto v : stats.sequence_count) {
      if (counts_by_freq.count(v.second) == 0) {
        counts_by_freq[v.second] = 0;
      }
      counts_by_freq[v.second]++;
    }

    // Now we change it to the FastQC corrected extrapolation
    for (auto v : counts_by_freq) {
      counts_by_freq[v.first] =
      get_corrected_count(stats.count_at_limit, stats.num_reads,
                          v.first, v.second);
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
    total_deduplicated_pct = 100.0 * seq_dedup / std::max(1.0, seq_total);

    // Convert to percentage
    for (auto &v : percentage_deduplicated)
      v = 100.0 * v / std::max(1.0, seq_dedup);  // Percentage of unique sequences in bin

     // Convert to percentage
    for (auto &v : percentage_total)
      v = 100.0 * v / std::max(1.0, seq_total);  // Percentage of sequences in bin
}

void
ModuleSequenceDuplicationLevels::make_grade() {
  // pass warn fail criteria : unique reads must be >80%
  // (otherwise warn) or >50% (otherwisefail)
  if (total_deduplicated_pct <= grade_dup_error) {
    grade = "fail";
  }
  else if (total_deduplicated_pct <= grade_dup_warn) {
    grade = "warn";
  }
}

void
ModuleSequenceDuplicationLevels::write_module(ostream &os) {
  os << "#Total Deduplicated Percentage\t" <<
         total_deduplicated_pct << "\n";

  os << "#Duplication Level\tPercentage of deduplicated\t"
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
}

string
ModuleSequenceDuplicationLevels::make_html_data() {
  ostringstream data;
  // non-deduplicated
  data << "{x : [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]";

  // total percentage in each bin
  data << ", y : [";
  for (size_t i = 0; i < 16; ++i) {
    data << percentage_total[i];

    if (i < 15)
      data << ", ";
  }
  data << "], type: 'line', line : {color : 'blue'}, "
       << "name : 'total sequences'}";

  // deduplicated
  data << ", {x : [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]";


  // total percentage in deduplicated
  data << ", y : [";
  for (size_t i = 0; i < 16; ++i) {
    data << percentage_deduplicated[i];

    if (i < 15)
      data << ", ";
  }
  data << "], type: 'line', line : {color : 'red'}, "
       << "name : 'deduplicated sequences'}";

  return data.str();
}


/************** OVERREPRESENTED SEQUENCES ********************/
const string
ModuleOverrepresentedSequences::module_name = "Overrepresented sequences";
ModuleOverrepresentedSequences::
ModuleOverrepresentedSequences(const FalcoConfig &config) :
Module(ModuleOverrepresentedSequences::module_name) {
  contaminants = config.contaminants;
  auto grade_overrep = config.limits.find("overrepresented")->second;
  grade_warn = grade_overrep.find("warn")->second;
  grade_error = grade_overrep.find("error")->second;
}

// gets the largest suffix of left which is a prefix of right
// returns 0 if none exist
size_t
get_overlap(const string &left, const string &right) {
  const auto left_st(begin(left));
  const auto left_end(end(left));

  const auto right_st(begin(right));
  const auto right_end(end(right));

  size_t best = 0;
  const auto sz = left.size();
  for (size_t start = 0; start < sz; ++start) {
    size_t cur = 0;
    auto left_it(left_st + start);
    auto right_it(right_st);
    for (; (left_it != left_end) && (right_it != right_end) && (*left_it == *right_it);
        ++left_it, ++right_it, ++cur);

    // does not accept if overlap does not cover a suffix or the
    // entirety of right
    cur = (left_it == left_end || right_it == right_end) ? (cur) : 0;
    best = max(best, cur);
  }
  return best;
}

string
ModuleOverrepresentedSequences::get_matching_contaminant (const string &seq) {
  size_t best = 0;
  string ret;
  for (auto v : contaminants) {
    const size_t cand = max(get_overlap(v.second, seq), get_overlap(seq, v.second));
    if (cand > best) {
      best = cand;
      ret = v.first;
    }
  }

  // If any sequence is a match, return the best one
  if (best >= min(ret.size(), seq.size())/2)
    return ret;
  return "No Hit";
}

void
ModuleOverrepresentedSequences::summarize_module(FastqStats &stats) {
  // Keep only sequences that pass the input cutoff
  num_reads = stats.num_reads;
  for (auto it = stats.sequence_count.begin();
            it != stats.sequence_count.end(); ++it) {
    if (it->second > num_reads * min_fraction_to_overrepresented) {
      overrep_sequences.push_back(*it);
    }
  }

  // Sort strings by frequency
  sort(begin(overrep_sequences), end(overrep_sequences),
       [](const pair<string, size_t> &a, const pair<string, size_t> &b){
         return a.second > b.second;
       });
}

void
ModuleOverrepresentedSequences::make_grade() {
  for (auto seq : overrep_sequences) {
    // implment pass warn fail for overrep sequences
    if (grade != "fail") {
      // get percentage that overrep reads represent
      double pct = 100.0 * seq.second / std::max(static_cast<size_t>(1), num_reads);
      if (pct > grade_error) {
        grade = "fail";
      }
      else if (pct > grade_warn) {
        grade = "warn";
      }
    }
  }
}

void
ModuleOverrepresentedSequences::write_module(ostream &os) {
  if (overrep_sequences.size() > 0)
    os << "#Sequence\tCount\tPercentage\tPossible Source\n";
  for (auto seq : overrep_sequences) {
      os << seq.first << "\t" << seq.second <<  "\t" <<
        100.0 * seq.second / std::max(static_cast<size_t>(1), num_reads) << "\t"
        << get_matching_contaminant(seq.first) << "\n";
  }
}

string
ModuleOverrepresentedSequences::make_html_data() {
  if (overrep_sequences.size() == 0) {
    return "No overrepresented sequences";
  }
  ostringstream data;
  // Header
  data << "<table><thead><tr>";
  data << "<th>Sequence</th>";
  data << "<th>Count</th>";
  data << "<th>Percentage</th>";
  data << "<th>Possible Source</th>";
  data << "</tr></thead><tbody>";

  // All overrep sequences
  for (auto v : overrep_sequences) {
    data << "<tr><td>" << v.first << "</td>";
    data << "<td>" << v.second << "</td>";
    data << "<td>" << 100.0 * v.second / std::max(static_cast<size_t>(1), num_reads) << "</td>";
    data << "<td>" << get_matching_contaminant(v.first)
         << "</td>";
    data << "</tr>";
  }
  data << "</tbody></table>";

  return data.str();
}


/************** ADAPTER CONTENT ***********/
const string
ModuleAdapterContent::module_name = "Adapter Content";
ModuleAdapterContent::
ModuleAdapterContent(const FalcoConfig &config) :
Module(ModuleAdapterContent::module_name) {
  // data parsed from config
  adapter_names = config.adapter_names;
  adapter_seqs = config.adapter_seqs;
  adapter_hashes = config.adapter_hashes;
  shortest_adapter_size = config.shortest_adapter_size;

  // check if they are all the same size
  if (adapter_names.size() != adapter_seqs.size())
    throw runtime_error("Adapter name and sequence vector sizes differ");

  if (adapter_names.size() != adapter_hashes.size())
    throw runtime_error("Adapter name and hash vector sizes differ");

  num_adapters = adapter_names.size();

  // maximum adapter % before pass/warn/fail
  auto grade_adapter = config.limits.find("adapter")->second;
  grade_warn = grade_adapter.find("warn")->second;
  grade_error = grade_adapter.find("error")->second;

  adapter_size = config.adapter_size;
}

void
ModuleAdapterContent::summarize_module(FastqStats &stats) {
  num_bases = max(min(stats.max_read_length, FastqStats::SHORT_READ_THRESHOLD),
                  ((shortest_adapter_size >= 1) ? (shortest_adapter_size - 1) : 0));

  if (num_bases + 1 >= shortest_adapter_size) { // otherwise this module does not make sense
    adapter_pos_pct.clear();
    for (size_t i = 0; i < num_adapters; ++i)
      adapter_pos_pct.push_back(
          vector<double>(num_bases - shortest_adapter_size + 1, 0.0)
      );

    size_t cnt = 0;
    for (size_t i = 0; i < adapter_pos_pct.size(); ++i) {
      for (size_t j = 0; j < adapter_pos_pct[0].size(); ++j) {
        cnt = 0;
        // check if position even makes sense
        if (j + adapter_seqs[i].size() < num_bases + 1 && j + adapter_seqs[i].size() >= 1)
          cnt = stats.pos_adapter_count[
            ((j + adapter_seqs[i].size() - 1) << stats.kBitShiftAdapter) | i
          ];

        if (j == 0) adapter_pos_pct[i][j] = cnt;
        else adapter_pos_pct[i][j] = adapter_pos_pct[i][j-1] + cnt;
      }
    }

    // now convert the counts we got before to percentage
    for (size_t i = 0; i < adapter_pos_pct.size(); ++i) {
      for (size_t j = 0; j < adapter_pos_pct[0].size(); ++j) {
        adapter_pos_pct[i][j] *= 100.0;
        adapter_pos_pct[i][j] /= std::max(num_lim<double>::min(),
                                          static_cast<double>(stats.num_reads));
      }
    }
  }
}

void
ModuleAdapterContent::make_grade() {
  for (size_t i = 0; i < adapter_pos_pct.size(); ++i) {
    for (size_t j = 0; j < adapter_pos_pct[0].size(); ++j) {
      if (grade != "fail") {
        if (adapter_pos_pct[i][j] > grade_error) {
          grade = "fail";
        } else if (adapter_pos_pct[i][j] > grade_warn) {
          grade = "warn";
        }
      }
    }
  }
}

void
ModuleAdapterContent::write_module(ostream &os) {
  // ADS: number of positions with data calculated
  const size_t n_pos_calc =
    (adapter_pos_pct.empty() ? 0 : adapter_pos_pct[0].size());

  os << "#Position";

  // adapter names
  for (size_t i = 0; i < num_adapters; ++i)
    os << "\t" << adapter_names[i];
  os << "\n";

  // matrix of percentages
  for (size_t i = 0; i < n_pos_calc; ++i) {
    os << i + 1;
    for (size_t j = 0; j < num_adapters; ++j)
      os << "\t" << adapter_pos_pct[j][i];
    os << "\n";
  }
  // ADS: now be sure to print the full read length
  for (size_t i = n_pos_calc; i < num_bases; ++i) {
    os << i + 1;
    for (size_t j = 0; j < num_adapters; ++j)
      // ADS: since this is cumulative, pad with final entry; but only
      // if that exists...
      if (!adapter_pos_pct[j].empty())
        os << "\t" << adapter_pos_pct[j].back();
    os << "\n";
  }
}

string
ModuleAdapterContent::make_html_data() {
  bool seen_first = false;
  ostringstream data;

  // ADS: number of positions with data calculated
  const size_t n_pos_calc =
    (adapter_pos_pct.empty() ? 0 : adapter_pos_pct[0].size());

  for (size_t i = 0; i < num_adapters; ++i) {
    if (!seen_first) {
      seen_first = true;
    }
    else {
      data << ",";
    }
    data << "{";

    // X values : read position
    data << "x : [";
    for (size_t j = 0; j < n_pos_calc; ++j) {
      data << j+1;
      if (j + 1 < num_bases) data << ",";
    }
    for (size_t j = n_pos_calc; j < num_bases; ++j) {
      data << j+1;
      if (j + 1 < num_bases) data << ",";
    }
    data << "]";

    // Y values : cumulative adapter frequency
    data << ", y : [";
    for (size_t j = 0; j < n_pos_calc; ++j) {
      data << adapter_pos_pct[i][j];
      if (j + 1 < num_bases)
        data << ",";
    }
    // ADS: Only write the final entry of the adaptor position
    // percentage if that exists. I don't know how we could get here
    // if it doesn't, but it's happening.
    if (!adapter_pos_pct[i].empty()) {
      for (size_t j = n_pos_calc; j < num_bases; ++j) {
        data << adapter_pos_pct[i].back();
        if (j + 1 < num_bases)
          data << ",";
      }
    }

    data << "]";
    data << ", type : 'line', ";
    data << "name : \"" << adapter_names[i] << "\"}";
  }
  return data.str();
}

/************** KMER CONTENT ******************************/
const string
ModuleKmerContent::module_name = "Kmer Content";

// ADS: Defining const static integer class variables here for
// correctness. Optimizer has been ignoring the issue. Hopefully it
// still will when turned on, and allow non-optimized code for
// debugging to compile.
const size_t ModuleKmerContent::MIN_OBS_EXP_TO_REPORT;
const size_t ModuleKmerContent::MAX_KMERS_TO_REPORT;
const size_t ModuleKmerContent::MAX_KMERS_TO_PLOT;

ModuleKmerContent::
ModuleKmerContent(const FalcoConfig &config) :
Module(ModuleKmerContent::module_name) {
  auto grade_kmer = config.limits.find("kmer")->second;
  grade_warn = grade_kmer.find("warn")->second;
  grade_error = grade_kmer.find("error")->second;
}

void
ModuleKmerContent::summarize_module(FastqStats &stats) {
  kmer_size = Constants::kmer_size;

  // 4^kmer size
  num_kmers = (1 << (2 * kmer_size));
  if (stats.max_read_length < FastqStats::SHORT_READ_THRESHOLD)
    num_kmer_bases = stats.max_read_length;
  else
    num_kmer_bases = FastqStats::SHORT_READ_THRESHOLD;

  // copies counts of all kmers per base position from stats
  pos_kmer_count = stats.pos_kmer_count;

  // Allocates space for all statistics
  obs_exp_max = vector<double>(num_kmers, 0.0);
  where_obs_exp_is_max = vector<size_t>(num_kmers, 0);
  total_kmer_counts = vector<size_t>(num_kmers, 0);

  // temp variables
  size_t observed_count;
  double expected_count;
  double obs_exp_ratio;
  num_seen_kmers = 0;

  // Here we get the total count of all kmers and the number of observed kmers
  for (size_t kmer = 0; kmer < num_kmers; ++kmer) {
    for (size_t i = kmer_size - 1; i < num_kmer_bases; ++i) {
      observed_count = stats.kmer_count[
                          (i << Constants::bit_shift_kmer) | kmer
                       ];
      total_kmer_counts[kmer] += observed_count;
    }
    if (total_kmer_counts[kmer] > 0) ++num_seen_kmers;
  }

  double dividend = static_cast<double>(num_seen_kmers);
  for (size_t kmer = 0; kmer < num_kmers; ++kmer) {
    for (size_t i = kmer_size - 1; i < num_kmer_bases; ++i) {
      observed_count =
        stats.kmer_count[(i << Constants::bit_shift_kmer) | kmer];

      expected_count = pos_kmer_count[i] / std::max(num_lim<double>::min(), dividend);
      // ADS: below, denom can't be zero if not above?
      obs_exp_ratio = (expected_count > 0) ? (observed_count / expected_count) : 0;

      if (i == 0 || obs_exp_ratio > obs_exp_max[kmer]) {
        obs_exp_max[kmer] = obs_exp_ratio;
        where_obs_exp_is_max[kmer] = i + 1 - kmer_size;
      }
    }

    if (obs_exp_max[kmer] > MIN_OBS_EXP_TO_REPORT) {
      kmers_to_report.push_back(make_pair(kmer, obs_exp_max[kmer]));
    }
  }

  sort (begin(kmers_to_report), end(kmers_to_report),
        [](const pair<size_t, double> &a, const pair<size_t, double> &b) {
          return a.second > b.second;
        });
}

void
ModuleKmerContent::make_grade() {
  const size_t lim = min(kmers_to_report.size(), MAX_KMERS_TO_REPORT);
  grade = "pass";

  // the worst kmer is at the top
  if (lim > 0) {
    const size_t kmer = kmers_to_report[0].first;
    const double obs_exp = obs_exp_max[kmer];
    if (obs_exp >= grade_error)
      grade = "fail";
    else if (obs_exp >= grade_warn)
      grade = "warn";
  }
}

void
ModuleKmerContent::write_module(ostream &os) {
  os << "#Sequence" << '\t'
     << "Count" << '\t'
     << "PValue" << '\t'
     << "Obs/Exp Max" << '\t'
     << "Max Obs/Exp Position" << '\n';
  const size_t lim = min(kmers_to_report.size(), MAX_KMERS_TO_REPORT);
  for (size_t i = 0; i < lim; ++i) {
    const size_t kmer = kmers_to_report[i].first;
    os << size_t_to_seq(kmer, kmer_size) << "\t"
       << total_kmer_counts[kmer] << "\t"
       << "0.0" << "\t"
       << obs_exp_max[kmer] << "\t"
       << where_obs_exp_is_max[kmer] << "\n";
  }
}

string
ModuleKmerContent::make_html_data() {
  bool seen_first = false;
  ostringstream data;
  const size_t lim = min(kmers_to_report.size(), MAX_KMERS_TO_PLOT);

  // get xlim to plot: whatever the largest position with some
  // reported k-mer is
  size_t xlim = 0;
  for (size_t i = 0; i < lim; ++i)
    xlim = max(xlim, where_obs_exp_is_max[kmers_to_report[i].first]);
  xlim += kmer_size;

  for (size_t i = 0; i < lim; ++i) {
    const size_t kmer = kmers_to_report[i].first;
    const double log_obs_exp = log(kmers_to_report[i].second)/log(2.0);
    if (!seen_first)
      seen_first = true;
    else
      data << ",";
    data << "{";

    // X values : read position
    data << "x : [";
    for (size_t j = 0; j < xlim; ++j) {
      data << j+1;
      if (j < xlim - 1) data << ",";
    }
    data << "]";

    // Y values : A peak wherever the k-mer is seen the most
    data << ", y : [";
    for (size_t j = 0; j < xlim; ++j) {
      data << ((j == (where_obs_exp_is_max[kmer])) ? (log_obs_exp) : 0);
      if (j < xlim - 1)
        data << ",";
    }

    data << "]";
    data << ", type : 'line', ";
    data << "name : \"" << size_t_to_seq(kmer, Constants::kmer_size) << "\"}";
  }
  return data.str();
}
// clang-format on
