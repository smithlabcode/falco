/* MIT License
 *
 * Copyright (c) 2026 Andrew D Smith
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "tile_processor.hpp"

#include "falco_file_format.hpp"
#include "falco_grade.hpp"
#include "falco_utils.hpp"
#include "quality_score.hpp"

#include "boost/boost_unordered.hpp"

#include <htslib/bgzf.h>
#include <htslib/sam.h>

#define FMT_HEADER_ONLY
#include "fmt/format.h"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <format>
#include <iterator>
#include <map>
#include <memory>
#include <numeric>
#include <ranges>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

[[nodiscard]] auto
get_name_fastq(const std::string &filename) -> std::string {
  std::unique_ptr<BGZF, int (*)(BGZF *)> in(bgzf_open(std::data(filename), "r"),
                                            &bgzf_close);
  if (!in)
    throw std::runtime_error("failed to open gz file: " + filename);
  kstring_t str = KS_INITIALIZE;
  const auto r = bgzf_getline(in.get(), '\n', &str);
  if (r < 0)
    throw std::runtime_error("failed to read line from: " + filename);
  std::string line(str.s, str.l);
  ks_free(&str);
  return line;
}

[[nodiscard]] auto
get_name_bam(const std::string &filename) -> std::string {
  std::unique_ptr<samFile, int (*)(samFile *)> in(
    hts_open(std::data(filename), "r"), &hts_close);
  if (!in)
    throw std::runtime_error("failed to open BAM/SAM file: " + filename);
  std::unique_ptr<sam_hdr_t, void (*)(sam_hdr_t *)> h(sam_hdr_read(in.get()),
                                                      &sam_hdr_destroy);
  std::unique_ptr<bam1_t, void (*)(bam1_t *)> b(bam_init1(), &bam_destroy1);
  const auto r = sam_read1(in.get(), h.get(), b.get());  // -1 on EOF
  if (r < -1)
    throw std::runtime_error("failed reading bam record");
  return bam_get_qname(b);
}

[[nodiscard]] auto
tile_processor::get_centered() -> tile_processor::tiles_centered_t {
  std::vector<double> means(max_read_len);
  std::vector<double> n_tiles_for_size(max_read_len);
  for (const auto &tile_quals : quals | std::views::values) {
    for (const auto [i, q] : std::views::enumerate(tile_quals))
      means[i] += as_frac(q.first, q.second);
    if (!tile_quals.empty())
      ++n_tiles_for_size[std::size(tile_quals) - 1ul];
  }
  std::partial_sum(std::crbegin(n_tiles_for_size), std::crend(n_tiles_for_size),
                   std::rbegin(n_tiles_for_size));

  std::ranges::transform(
    means, n_tiles_for_size, std::begin(means),
    [&](const auto m, const auto l) { return as_frac(m, l); });

  // using map to get sorted order by tile id
  const auto get_cent = [](const auto &a, const auto mean) {
    const auto [qual_sum, n_obs] = a;
    return as_frac(qual_sum, n_obs) - mean;
  };

  tile_processor::tiles_centered_t centered;
  for (const auto &[id, vals] : quals) {
    const auto cent_vals = std::views::zip_transform(get_cent, vals, means);
    centered.emplace(id, cent_vals | std::ranges::to<std::vector>());
  }
  return centered;
};

auto
tile_processor::adjust_fastq_qual_encoding(const falco::encoding enc) -> void {
  const auto qual_offset = get_quality_score_offset(enc);
  for (auto &tile_quals : quals | std::views::values)
    for (auto &q : tile_quals)
      q.first -= q.second * qual_offset;
}

auto
tile_processor::trim() -> void {
  for (auto &tile_quals : quals | std::views::values) {
    auto first_trailing_zero = 0L;
    for (const auto [i, q] : std::views::enumerate(tile_quals))
      first_trailing_zero = q.second > 0 ? i + 1 : first_trailing_zero;
    tile_quals.resize(first_trailing_zero);
  }
}

auto
tile_processor::finalize(const file_info &info) -> void {
  trim();
  if (!is_mapped_reads(info.format))
    adjust_fastq_qual_encoding(info.encoding);
}

auto
tile_processor::apply_groups(const run_mode &mode) -> void {
  if (do_groups(mode)) {
    const auto groups = get_default_base_groups(max_read_len, do_groups(mode));
    for (auto &quals_for_tile : std::views::values(quals))
      apply_base_groups(groups, quals_for_tile, [](auto &a, const auto &b) {
        a.first += b.first;
        a.second += b.second;
      });
  }
}

auto
tile_processor::operator+=(const tile_processor &rhs)
  -> const tile_processor & {
  const auto pair_plus = [](const auto &a, const auto &b) {
    return std::pair{a.first + b.first, a.second + b.second};
  };
  for (const auto &[i, q] : rhs.quals) {
    if (quals.contains(i)) {
      quals[i].resize(std::max(std::size(quals[i]), std::size(q)));
      std::ranges::transform(quals[i], q, std::begin(quals[i]), pair_plus);
    }
    else
      quals.emplace(i, q);
  }
  return *this;
}

auto
get_tile_info(const std::string &filename) -> std::uint32_t {
  // colon cutoffs taken from FastQC
  static constexpr auto colon_cutoff_1 = 6;
  static constexpr auto colon_cutoff_1_val = 4;
  static constexpr auto colon_cutoff_2 = 4;
  static constexpr auto colon_cutoff_2_val = 2;

  std::unique_ptr<htsFile, int (*)(htsFile *)> fp(
    hts_open(std::data(filename), "r"), &hts_close);
  if (!fp)
    throw std::runtime_error("failed to open file: " + filename);

  const auto hts_fmt = hts_get_format(fp.get());
  if (hts_fmt->format != fastq_format && hts_fmt->format != bam &&
      hts_fmt->format != sam)
    throw std::runtime_error(std::format("unknown file format: {}",
                                         std::to_underlying(hts_fmt->format)));

  const auto line = (hts_fmt->format == bam || hts_fmt->format == sam)
                      ? get_name_bam(filename)
                      : get_name_fastq(filename);

  const auto colons_found = std::ranges::count(line, ':');
  return colons_found >= colon_cutoff_1
           ? colon_cutoff_1_val
           : (colons_found >= colon_cutoff_2 ? colon_cutoff_2_val : 0);
}

[[nodiscard]] auto
get_grade_tile(const tile_processor::tiles_centered_t &centered)
  -> std::string {
  static constexpr auto label = "tile";
  assert(!centered.empty());
  const auto neg_min_cent_qual =
    -std::ranges::min(centered | std::views::values | std::views::join);
  return grader_set::get_grade(label, neg_min_cent_qual);
}
