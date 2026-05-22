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
#include "falco_utils.hpp"

#include <htslib/bgzf.h>
#include <htslib/sam.h>

#include <algorithm>
#include <cstdint>
#include <format>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

std::uint32_t tile_processor::preceding_colons = 0;

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

auto
tile_processor::set_preceding_colons(const std::string &filename)
  -> std::uint32_t {
  // colon cutoffs taken from FastQC
  static constexpr auto colon_cutoff_1 = 6;
  static constexpr auto colon_cutoff_1_val = 4;
  static constexpr auto colon_cutoff_2 = 4;
  static constexpr auto colon_cutoff_2_val = 2;

  std::unique_ptr<htsFile, int (*)(htsFile *)> fp(
    hts_open(std::data(filename), "r"), &hts_close);
  if (!fp)
    throw std::runtime_error("failed to open file: " + filename);

  const auto fmt = hts_get_format(fp.get());
  if (fmt->format != fastq_format && fmt->format != bam && fmt->format != sam)
    throw std::runtime_error(
      std::format("unknown file format: {}", std::to_underlying(fmt->format)));

  const auto line = (fmt->format == bam || fmt->format == sam)
                      ? get_name_bam(filename)
                      : get_name_fastq(filename);

  const auto colons_found = std::ranges::count(line, ':');
  // clang-format off
  preceding_colons =
    colons_found >= colon_cutoff_1 ? colon_cutoff_1_val :
    (colons_found >= colon_cutoff_2 ? colon_cutoff_2_val : 0);
  // clang-format on
  return preceding_colons;
}

[[nodiscard]] auto
tile_processor::string(const std::uint32_t len) const -> std::string {
  static constexpr auto max_precision{std::numeric_limits<double>::digits10};
  static constexpr auto start_tag = ">>Per tile sequence quality\t{}\n";
  static constexpr auto header = "#Tile\t"
                                 "Base\t"
                                 "Mean\n";
  auto r = std::format(start_tag, "pass");
  r += header;

  auto means = std::vector<double>(max_read_len, 0.0);
  for (auto i = 0u; i < max_read_len; ++i)
    for (const auto &q : quals | std::views::values)
      means[i] += as_frac(q[i].first, q[i].second);

  std::ranges::transform(means, std::begin(means),
                         [&](const auto x) { return x / std::size(quals); });

  // also sorting tiles
  std::map<std::uint32_t, std::vector<double>> centered;
  const auto cent = [](const auto &a, const auto b) {
    return as_frac(a.first, a.second) - b;
  };
  for (const auto &[id, vals] : quals)
    centered.emplace(id, std::views::zip_transform(cent, vals, means) |
                           std::ranges::to<std::vector>());
  for (const auto &[i, q] : centered | std::views::take(len))
    for (auto j = 0u; j < std::size(q); ++j)
      r += std::format("{}\t{}\t{:.{}f}\n", i, j + 1, q[j], max_precision);
  return r + end_module_tag;
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
