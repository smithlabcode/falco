#ifndef _HTMLMAKER_HPP
#define _HTMLMAKER_HPP
#include <string>
#include <fstream>
#include "Config.hpp"
#include "FastqStats.hpp"

/*******************************************************/
/*************** HTML MAKER ****************************/
/*******************************************************/
class HtmlMaker {
 private:
  void replace_placeholder_with_data(const std::string &placeholder,
                                     const std::string &data);

  // Replace placeholders with html comments if analysis was skipped
  void comment_if_skipped(std::string ph_begin, std::string ph_end, bool skip);

  // First thing to do: file details
  void put_file_details(const Config &config);

  // Comment out analyses that were skipped
  void comment_out(const Config &config);

  // Functions to replace pass warn fail
  void put_pass_warn_fail(const FastqStats &stats);

  // Function to replace template placeholders with data
  void make_basic_statistics(const FastqStats &stats,
                             Config &config);

  void make_position_quality_data(const FastqStats &stats,
                                  const Config &config);

  void make_tile_quality_data(FastqStats &stats,
                              const Config &config);

  void make_sequence_quality_data(const FastqStats &stats,
                                  const Config &config);
  void make_base_sequence_content_data(const FastqStats &stats,
                                       const Config &config);

  void make_sequence_gc_content_data(const FastqStats &stats,
                                     const Config &config);

  void make_base_n_content_data(const FastqStats &stats,
                                const Config &config);

  void make_sequence_length_data(const FastqStats &stats,
                                 const Config &config);

  void make_sequence_duplication_data(const FastqStats &stats,
                                      const Config &config);

  void make_overrepresented_sequences_data(const FastqStats &stats,
                                           const Config &config);

  void make_adapter_content_data(FastqStats &stats,
                                 Config &config);

 public:
  std::string sourcecode;
  explicit HtmlMaker(std::string filepath);
  void write(FastqStats &stats, Config &config);
};
#endif
