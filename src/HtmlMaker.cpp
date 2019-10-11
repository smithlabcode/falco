#include "HtmlMaker.hpp"
#include <algorithm>
#include <sstream>
#include <chrono>
#include <fstream>
#include <algorithm>

using std::ostringstream;
using std::string;
using std::vector;
using std::sort;
using std::ifstream;
using std::runtime_error;
using std::chrono::system_clock;
using std::min;

HtmlMaker::HtmlMaker(string filepath) {
  sourcecode = "";
  ifstream in(filepath);
  if (!in) {
    throw runtime_error("HTML layout not found: " + filepath);
  }

  // pass the whole source code template to a string
  string line;
  while (getline(in, line))
    sourcecode += line + "\n";
}

void
HtmlMaker::replace_placeholder_with_data(const string &placeholder,
                                           const string &data) {
  auto pos = sourcecode.find(placeholder);

  // Placeholder not found
  if (pos == string::npos) {
    throw runtime_error("placeholder not found: " + placeholder);
  }

  // at least one placeholder found
  while (pos != string::npos) {
    sourcecode.replace(pos, placeholder.size(), data);
    pos = sourcecode.find(placeholder, pos + 1);
  }
}

// Comments out html pieces if analyses were skipped
void
HtmlMaker::comment_if_skipped(string ph_begin,
                                string ph_end,
                                bool done) {
  // put html comments if analysis was skipped
  if (!done) {
    replace_placeholder_with_data(ph_begin, "<!--");
    replace_placeholder_with_data(ph_end, "-->");
  }

  // otherwise delete placeholder
  else {
    replace_placeholder_with_data(ph_begin, "");
    replace_placeholder_with_data(ph_end, "");
  }
}

void
HtmlMaker::comment_out(const Config &config) {
  comment_if_skipped("{{skipbasesequence_s}}",
                     "{{skipbasesequence_e}}",
                     config.do_quality_base);
  comment_if_skipped("{{skiptile_s}}",
                     "{{skiptile_e}}",
                     config.do_tile);

  comment_if_skipped("{{skipsequencequal_s}}",
                     "{{skipsequencequal_e}}",
                     config.do_quality_sequence);

  comment_if_skipped("{{skipbasecontent_s}}",
                     "{{skipbasecontent_e}}",
                     config.do_sequence);

  comment_if_skipped("{{skipsequencegc_s}}",
                     "{{skipsequencegc_e}}",
                     config.do_gc_sequence);

  comment_if_skipped("{{skipbasencontent_s}}",
                     "{{skipbasencontent_e}}",
                     config.do_n_content);

  comment_if_skipped("{{skipseqlength_s}}",
                     "{{skipseqlength_e}}",
                     config.do_n_content);

  comment_if_skipped("{{skipseqdup_s}}",
                     "{{skipseqdup_e}}",
                     config.do_duplication);

  comment_if_skipped("{{skipoverrep_s}}",
                     "{{skipoverrep_e}}",
                     config.do_overrepresented);

  comment_if_skipped("{{skipadapter_s}}",
                     "{{skipadapter_e}}",
                     config.do_adapter);
}

void
HtmlMaker::put_file_details(const Config &config) {
  // Put filename in filename placeholder
  replace_placeholder_with_data("{{filename}}",
                                config.filename_stripped);

  // Put date on date placeholder
  auto tmp = system_clock::to_time_t(system_clock::now());
  string time_fmt = string(ctime(&tmp));
  replace_placeholder_with_data("{{date}}", time_fmt);
}

void
HtmlMaker::put_pass_warn_fail(const FastqStats &stats) {
  replace_placeholder_with_data("{{passbasic}}",
                                stats.pass_basic_statistics);
  replace_placeholder_with_data("{{passbasesequence}}",
                                stats.pass_per_base_sequence_quality);
  replace_placeholder_with_data("{{passpertile}}",
                                stats.pass_per_tile_sequence_quality);
  replace_placeholder_with_data("{{passpersequencequal}}",
                                stats.pass_per_sequence_quality_scores);
  replace_placeholder_with_data("{{passperbasecontent}}",
                                stats.pass_per_base_sequence_content);
  replace_placeholder_with_data("{{passpersequencegc}}",
                                stats.pass_per_sequence_gc_content);
  replace_placeholder_with_data("{{passperbasencontent}}",
                                stats.pass_per_base_n_content);
  replace_placeholder_with_data("{{passseqlength}}",
                                stats.pass_sequence_length_distribution);
  replace_placeholder_with_data("{{passseqdup}}",
                                stats.pass_duplicate_sequences);
  replace_placeholder_with_data("{{passoverrep}}",
                                stats.pass_overrepresented_sequences);
  replace_placeholder_with_data("{{passadapter}}",
                                stats.pass_adapter_content);
}
void
HtmlMaker::make_basic_statistics(const FastqStats &stats,
                                   Config &config) {
  string placeholder = "{{BASICSTATSDATA}}";
  ostringstream data;
  data << "<table><thead><tr><th>Measure</th><th>Value"
       << "</th></tr></thead><tbody>";
  data << "<tr><td>Filename</td><td>" << config.filename_stripped
       << "</td></tr>";
  data << "<tr><td>Total Sequences</td><td>" << stats.num_reads << "</td></tr>";
  data << "<tr><td>Sequences Flagged As Poor Quality</td><td>"
       << stats.num_poor << "</td></tr>";
  data << "<tr><td>Sequence length</td><td>";
  if (stats.min_read_length != stats.max_read_length) {
    data << stats.min_read_length << " - " << stats.max_read_length;
  }
  else {
    data << stats.max_read_length;
  }
  data << "</td></tr>";

  if (config.do_sequence) {
    data << "<tr><td>%GC:</td><td>" << stats.avg_gc << "</td></tr>";
  }
  data << "</tbody></table>";

  replace_placeholder_with_data(placeholder, data.str());
}

void
HtmlMaker::make_position_quality_data(const FastqStats &stats,
                                        const Config &config) {
  ostringstream data;
  const string placeholder = "{{SEQBASEQUALITYDATA}}";
  if (config.do_quality_base) {
    size_t cur_median;
    for (size_t i = 0; i < stats.max_read_length; ++i) {
      data << "{y : [";

     if (i < stats.kNumBases) {
       cur_median = stats.median[i];
       data << stats.ldecile[i] << ", "
            << stats.lquartile[i] << ", "
            << stats.median[i] << ", "
            << stats.uquartile[i] << ", "
            << stats.udecile[i] << "], ";
     } else {
       cur_median = stats.long_median[i - stats.kNumBases];
       data << stats.long_ldecile[i - stats.kNumBases] << ", "
            << stats.long_lquartile[i - stats.kNumBases] << ", "
            << stats.long_median[i - stats.kNumBases] << ", "
            << stats.long_uquartile[i - stats.kNumBases] << ", "
            << stats.long_udecile[i - stats.kNumBases] << "], ";
     }
     data << "type : 'box', name : ' " << i << "', ";
     data << "marker : {color : '";
     if (cur_median > 30)
       data << "green";
     else if (cur_median > 20)
       data << "yellow";
     else
       data << "red";
     data << "'}}";
     if (i < stats.max_read_length - 1)
       data << ", ";
    }
  }
  replace_placeholder_with_data(placeholder, data.str());
}

void
HtmlMaker::make_tile_quality_data(FastqStats &stats,
                                    const Config &config) {
  ostringstream data;
  const string placeholder = "{{TILEQUALITYDATA}}";

  if (config.do_tile) {
    // first thing is to get the sorted tile values again
    vector<size_t> tiles_sorted;
    for (auto v : stats.tile_count)
      tiles_sorted.push_back(v.first);
    sort(tiles_sorted.begin(), tiles_sorted.end());

    // X: position
    data << "{x : [";
    for (size_t i = 0; i < stats.max_read_length; ++i) {
      data << i+1;
      if (i < stats.max_read_length - 1)
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
      for (size_t j = 0; j < stats.max_read_length; ++j) {
        data << stats.tile_position_quality[tiles_sorted[i]][j];
        if (j < stats.max_read_length - 1) data << ",";
      }
      data << "]";
    }
    data << "]";
    data << ", type : 'heatmap' }";
  }
  replace_placeholder_with_data(placeholder, data.str());
}

void
HtmlMaker::make_sequence_quality_data(const FastqStats &stats,
                                        const Config &config) {
  ostringstream data;
  const string placeholder = "{{SEQQUALITYDATA}}";

  if (config.do_quality_sequence) {
    // X values : avg quality phred scores
    data << "{x : [";
    for (size_t i = 0; i < 41; ++i) {
      data << i + stats.kBaseQuality;
      if (i < 40)
        data << ", ";
    }

    // Y values: frequency with which they were seen
    data << "], y : [";
    for (size_t i = 0; i < 41; ++i) {
      data << stats.quality_count[i];
      if (i < 40)
        data << ", ";
    }
    data << "], type: 'line', line : {color : 'red'}, "
         << "name : 'Sequence quality distribution'}";
  }
  replace_placeholder_with_data(placeholder, data.str());
}

void 
HtmlMaker::make_base_sequence_content_data(const FastqStats &stats,
                                             const Config &config) {
  ostringstream data;
  const string placeholder = "{{BASESEQCONTENTDATA}}";

  if (config.do_sequence) {
    // ATGC
    for (size_t base = 0; base < stats.kNumNucleotides; ++base) {
      // start line
      data << "{";

      // X values : base position
      data << "x : [";
      for (size_t i = 0; i < stats.max_read_length; ++i) {
        data << i+1;
        if (i < stats.max_read_length - 1)
          data << ", ";
      }

      // Y values: frequency with which they were seen
      data << "], y : [";
      for (size_t i = 0; i < stats.max_read_length; ++i) {
        if (base == 0) {
          if (i < stats.kNumBases) data << stats.a_pct[i];
          else data << stats.long_a_pct[i - stats.kNumBases];
        }
        if (base == 1) {
          if (i < stats.kNumBases) data << stats.c_pct[i];
          else data << stats.long_c_pct[i - stats.kNumBases];
        }
        if (base == 2) {
          if (i < stats.kNumBases) data << stats.t_pct[i];
          else data << stats.long_t_pct[i - stats.kNumBases];
        }
        if (base == 3) {
          if (i < stats.kNumBases) data << stats.g_pct[i];
          else data << stats.long_g_pct[i - stats.kNumBases];
        }
        if (i < stats.max_read_length - 1)
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
      if (base < stats.kNumNucleotides - 1)
        data << ", ";
    }
  }
  replace_placeholder_with_data(placeholder, data.str());
}

void
HtmlMaker::make_sequence_gc_content_data(const FastqStats &stats,
                                           const Config &config) {
  ostringstream data;
  const string placeholder = "{{SEQGCCONTENTDATA}}";
  if (config.do_gc_sequence) {
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
      data << stats.smooth_gc_count[i];
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
      data << stats.theoretical_gc_count[i];
      if (i < 101)
        data << ", ";
    }
    data << "], type: 'line', line : {color : 'blue'},"
         << "name : 'Theoretical distribution'}";
  }
  replace_placeholder_with_data(placeholder, data.str());
}

void
HtmlMaker::make_base_n_content_data(const FastqStats &stats,
                                      const Config &config) {
  ostringstream data;
  const string placeholder = "{{BASENCONTENTDATA}}";

  if (config.do_n_content) {
    // base position
    data << "{x : [";
    for (size_t i = 0; i < stats.max_read_length; ++i) {
      data << i + 1;
      if (i < stats.max_read_length - 1)
        data << ", ";
    }

    // Y values: frequency with which they were seen
    data << "], y : [";
    for (size_t i = 0; i < stats.max_read_length; ++i) {
      if (i < stats.kNumBases)
        data << stats.n_pct[i];
      else
        data << stats.long_n_pct[i - stats.kNumBases];

      if (i < stats.max_read_length - 1)
        data << ", ";
    }
    data << "], type: 'line', line : {color : 'red'}, "
         << "name : 'Fraction of N reads per base'}";
  }
  replace_placeholder_with_data(placeholder, data.str());
}

void
HtmlMaker::make_sequence_length_data(const FastqStats &stats,
                                       const Config &config) {
  ostringstream data;
  const string placeholder = "{{SEQLENDATA}}";

  if (config.do_sequence_length) {
    // X values : avg quality phred scores
    data << "{x : [";
    bool first_seen = false;
    size_t val;
    for (size_t i = 0; i < stats.max_read_length; ++i) {
      val = 0;
      if (i < stats.kNumBases) {
        if (stats.read_length_freq[i] > 0) {
          val = stats.read_length_freq[i];
        }
      } else {
        if (stats.long_read_length_freq[i - stats.kNumBases] > 0) {
          val = stats.long_read_length_freq[i - stats.kNumBases];
        }
      }

      if (val > 0) {
        if (first_seen)
          data << ",";
        data << "\"" << i+1 << " bp\"";
        first_seen = true;
      }
    }

    // Y values: frequency with which they were seen
    data << "], y : [";
    first_seen = false;
    for (size_t i = 0; i < stats.max_read_length; ++i) {
      val = 0;
      if (i < stats.kNumBases) {
        if (stats.read_length_freq[i] > 0) {
          val = stats.read_length_freq[i];
        }
      } else {
        if (stats.long_read_length_freq[i - stats.kNumBases] > 0) {
          val = stats.long_read_length_freq[i - stats.kNumBases];
        }
      }

      if (val > 0) {
        if (first_seen)
          data << ",";
        data << val;
        first_seen = true;
      }
    }

    // Put the sequence value in the text
    data << "], text : [";
    first_seen = false;
    for (size_t i = 0; i < stats.max_read_length; ++i) {
      val = 0;
      if (i < stats.kNumBases) {
        if (stats.read_length_freq[i] > 0) {
          val = stats.read_length_freq[i];
        }
      } else {
        if (stats.long_read_length_freq[i - stats.kNumBases] > 0) {
          val = stats.long_read_length_freq[i - stats.kNumBases];
        }
      }

      if (val > 0) {
        if (first_seen) data << ",";
        data << i+1;
        first_seen = true;
      }
    }


    data << "], type: 'bar', marker : {color : 'rgba(55,128,191,1.0)',"
         << "line : {width : 2}}, "
         << "name : 'Sequence length distribution'}";
  }
  replace_placeholder_with_data(placeholder, data.str());
}


void
HtmlMaker::make_sequence_duplication_data(const FastqStats &stats,
                                            const Config &config) {
  ostringstream data;
  const string placeholder = "{{SEQDUPDATA}}";

  if (config.do_duplication) {
    // non-deduplicated
    data << "{x : [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]";

    // total percentage in each bin
    data << ", y : [";
    for (size_t i = 0; i < 16; ++i) {
      data << stats.percentage_total[i];

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
      data << stats.percentage_deduplicated[i];

      if (i < 15)
        data << ", ";
    }
    data << "], type: 'line', line : {color : 'red'}, "
         << "name : 'total sequences'}";
  }
  replace_placeholder_with_data(placeholder, data.str());
}


void
HtmlMaker::make_overrepresented_sequences_data(const FastqStats &stats,
                                                 const Config &config) {
  string placeholder = "{{OVERREPSEQDATA}}";
  ostringstream data;

  if (config.do_overrepresented) {
    // Header
    data << "<table><thead><tr>";
    data << "<th>Sequence</th>";
    data << "<th>Count</th>";
    data << "<th>Percentage</th>";
    data << "<th>Possible Source</th>";
    data << "</tr></thead><tbody>";

    // All overrep sequences
    for (auto v : stats.overrep_sequences) {
      data << "<tr><td>" << v.first << "</td>";
      data << "<td>" << v.second << "</td>";
      data << "<td>" << 100.0 * v.second / stats.num_reads << "</td>";
      data << "<td>" << config.get_matching_contaminant(v.first)
          << "</td>";

      data << "</tr>";
    }
    data << "</tbody></table>";
  }
  replace_placeholder_with_data(placeholder, data.str());
}

void
HtmlMaker::make_adapter_content_data(FastqStats &stats,
                                       Config &config) {
  string placeholder = "{{ADAPTERDATA}}";
  ostringstream data;

  if (config.do_adapter) {
    // Number of bases to make adapter content
    size_t num_bases =  min(stats.kNumBases, stats.kKmerMaxBases);
    bool seen_first = false;

    size_t jj = 0;
    for (auto v : config.adapters) {
      if (!seen_first) {
        seen_first = true;
      }
      else {
        data << ",";
      }
      data << "{";

      // X values : read position
      data << "x : [";
      for (size_t i = 0; i < num_bases; ++i) {
        if (stats.cumulative_read_length_freq[i] > 0) {
          data << i+1;
          if (i < num_bases - 1)
            data << ",";
        }
      }
      data << "]";

      // Y values : cumulative adapter frequency
      data << ", y : [";
      for (size_t i = 0; i < num_bases; ++i) {
        if (stats.cumulative_read_length_freq[i] > 0) {
          data << stats.kmer_by_base[i][jj];
          if (i < num_bases - 1)
            data << ",";
        }
      }

      data << "]";
      data << ", type : 'line', ";
      data << "name : '" << v.first << "'}";
      ++jj;
    }
  }
  replace_placeholder_with_data(placeholder, data.str());
}

void
HtmlMaker::write(FastqStats &stats, Config &config) {
  // Filename and date
  put_file_details(config);

  // Put html comments around things to skip
  comment_out(config);

  // Put colors in pass, warn and fail
  put_pass_warn_fail(stats);

  // Put data on tables and plotly javascripts
  make_basic_statistics(stats, config);
  make_position_quality_data(stats, config);
  make_tile_quality_data(stats, config);
  make_sequence_quality_data(stats, config);
  make_base_sequence_content_data(stats, config);
  make_sequence_gc_content_data(stats, config);
  make_base_n_content_data(stats, config);
  make_sequence_length_data(stats, config);
  make_sequence_duplication_data(stats, config);
  make_overrepresented_sequences_data(stats, config);
  make_adapter_content_data(stats, config);
}


