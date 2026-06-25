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

#define FMT_HEADER_ONLY
#include "fmt/format.h"
#include "fmt/ranges.h"

#include <htslib/bgzf.h>
#include <htslib/sam.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <format>
#include <iterator>
#include <limits>
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
get_centered(const auto max_read_len, const auto &quals)
  -> std::map<std::uint32_t, std::vector<double>> {
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
  const auto cent = [](const auto &a, const auto mean) {
    const auto [qual_sum, n_obs] = a;
    return as_frac(qual_sum, n_obs) - mean;
  };
  std::map<std::uint32_t, std::vector<double>> centered;
  for (const auto &[id, vals] : quals)
    centered.emplace(id, std::views::zip_transform(cent, vals, means) |
                           std::ranges::to<std::vector>());
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

[[nodiscard]] static inline auto
get_max_size(const auto &x) {
  assert(!x.empty());
  const auto sz = [](const auto &y) { return std::size(y); };
  return std::ranges::max(std::views::transform(x | std::views::values, sz));
}

auto
tile_processor::finalize(const run_mode &mode, const file_info &info) -> void {
  assert(!quals.empty() && centered.empty());
  trim();
  if (!is_mapped_reads(info.format))
    adjust_fastq_qual_encoding(info.encoding);
  if (do_groups(mode)) {
    const auto groups = get_default_base_groups(max_read_len, do_groups(mode));
    for (auto &quals_for_tile : std::views::values(quals))
      apply_base_groups(groups, quals_for_tile, [](auto &a, const auto &b) {
        a.first += b.first;
        a.second += b.second;
      });
  }
  // ADS: max_read_len not useful after grouping
  centered = get_centered(get_max_size(quals), quals);
  quals.clear();  // ADS: force a crash if this is used after centering
}

[[nodiscard]] auto
tile_processor::get_grade() const -> std::string {
  static constexpr auto label = "tile";
  assert(!centered.empty());
  const auto neg_min_cent_qual =
    -std::ranges::min(centered | std::views::values | std::views::join);
  return grader_set::get_grade(label, neg_min_cent_qual);
}

[[nodiscard]] auto
tile_processor::get_report(const std::vector<base_group_t> &groups,
                           const std::string &grade) const -> std::string {
  static constexpr auto max_precision{std::numeric_limits<double>::digits10};
  static constexpr auto start_tag = ">>Per tile sequence quality\t{}\n";
  static constexpr auto header = "#Tile\t"
                                 "Base\t"
                                 "Mean\n";
  if (quals.empty() && centered.empty())
    return {};  // ADS: in case this fun is called for files w/o tile info
  auto r = std::format(start_tag, grade);
  r += header;
  for (const auto &[tile_id, q] : centered)
    for (auto j = 0u; j < std::size(q); ++j)
      r += std::format("{}\t{}\t{:.{}f}\n", tile_id, make_group_tag(groups[j]),
                       q[j], max_precision);
  return r + end_module_tag;
}

[[nodiscard]] auto
tile_processor::get_html(const std::vector<base_group_t> &groups,
                         const file_grades &grades) const -> std::string {
  static constexpr auto label = "tile";
  static constexpr auto n_quants = 20.0;
  // ADS: ??? (-10: red, 0: light blue, +10: dark blue)
  static constexpr auto tiles_plot_fmt =
    R"""([{{
x: [{}],
y: [{}],
z: [{}],
type: "heatmap",
colorscale: [
[0.0, "rgb(210,65,83)"],
[{}, "rgb(178,236,254)"],
[1.0, "rgb(34,57,212)"]
],
showscale: true,
}}],
{{
margin: {{t: 0}},
showlegend: false,
xaxis: {{title: "Base position"}},
yaxis: {{title: "tile", type: "category"}},
}}
)""";
  static constexpr auto plot_fmt = R"(<div id="{}"></div>
<script>Plotly.newPlot("{}",
{}
);</script>
)";
  assert(centered.empty() || get_max_size(centered) == std::size(groups));
  if (centered.empty())
    return {};  // ADS: in case we are here but tile analysis was not done
  const auto tag = [&](const auto &g) { return make_group_tag_quoted(g); };

  const auto [minval, maxval] =
    std::ranges::minmax(centered | std::views::values | std::views::join);
  // discretize quantiles so plotly understands color scheme (???)
  const auto midpoint = minval / (minval - maxval);
  const auto mid_discrete = std::round(n_quants * midpoint) / n_quants;

  // z: one formatted vector of quality z scores per tile
  const auto format1 = [](const auto &c) {
    // ADS: note the format and precision (3) here
    return fmt::format("[{:.3f}]", fmt::join(c, ","));  // inner lists
  };
  const auto z = std::views::transform(centered | std::views::values, format1);

  const auto grade = grades.grade(label);
  const auto title = grades.get_title(label);
  return fmt::format(
    html_module_fmt, grade, label, title, grade,
    fmt::format(
      plot_fmt, "tiles_plot", "tiles_plot",
      fmt::format(tiles_plot_fmt,
                  fmt::join(groups | std::views::transform(tag), ","),  // x
                  fmt::join(centered | std::views::keys, ","),          // y
                  fmt::join(z, ","),                                    // z
                  mid_discrete)));
}

auto
tile_processor::operator+=(const tile_processor &rhs)
  -> const tile_processor & {
  assert(centered.empty());
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
