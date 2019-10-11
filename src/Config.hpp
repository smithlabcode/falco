#ifndef _CONFIG_HPP
#define _CONFIG_HPP
#include <string>
#include <vector>
#include <unordered_map>
#include <utility>
#include "aux.hpp"

/*************************************************************
 ******************** CUSTOM CONFIGURATION *******************
 *************************************************************/

// config from options, constants, magic numbers, etc
struct Config {
 private:
  std::string strip_path (std::string full_path) const;
 public:
  /************************************************************
   *************** MY UNIVERSAL CONSTANTS *********************
   ************************************************************/
  // threshold for a sequence to be considered  poor quality
  size_t kPoorQualityThreshold;

  /************ OVERREPRESENTATION ESTIMTES **********/
  // fraction of the number of slow reads a sequence needs to be seen to be
  // considered a candiate for overrepresentation
  double kOverrepMinFrac;

  /************************************************************
   *************** FASTQC OPTION PARSER************************
   ************************************************************/
  bool casava;  // files from raw casava output
  bool nanopore;  // fast5 format
  bool nofilter;  // if running with --casava flag
  bool extract;  // if set the zipped file will be uncompressed
  bool nogroup;  // disable grouping of bases for reads >50bp
  bool compressed;  // whether or not to inflate file
  bool quiet;
  size_t min_length;  // lower limit in sequence length to be shown in report
  size_t threads;  // number of threads to read multiple files in parallel
  size_t kmer_size;  // kmer size
  std::string format;  // force file format
  std::string contaminants_file;  // custom contaminants file
  std::string adapters_file;  // adapters file
  std::string limits_file;  // file with limits and options and custom analyses
  std::string html_file;  // file with limits and options and custom analyses
  std::string tmpdir;  // dir for temp files when generating report images

  // config on how to handle reads
  bool do_duplication,
       do_kmer,
       do_n_content,
       do_overrepresented,
       do_quality_base,
       do_sequence,
       do_gc_sequence,
       do_quality_sequence,
       do_tile,
       do_adapter,
       do_sequence_length;

  /************************************************************
   *************** FASTQC LIMITS *******************************
   ************************************************************/
  // These will become const bools in the stream reader
  std::unordered_map <std::string,
                      std::unordered_map <std::string, double> > limits;
  static const std::vector <std::string> values_to_check;

  /*************** CONTAMINANTS *****************/
  std::vector <std::pair <std::string, std::string> > 
               contaminants;  // first = name, scond = seq

  std::vector <std::pair <std::string, size_t> > 
               adapters;  // kmer of the adapter prefix

  /*************** DEFINE FILE TYPE ************/

  // IO
  bool is_sam, is_bam, is_fastq, is_fastq_gz;
  std::string filename;
  std::string filename_stripped;

  /*********** FUNCTIONS TO READ FILES *************/
  Config();  // set magic defaults
  void define_file_format();
  void read_limits();  // populate limits hash map
  void read_adapters();
  void read_contaminants();

  void setup();
  std::string get_matching_contaminant(std::string seq) const;
};
#endif
