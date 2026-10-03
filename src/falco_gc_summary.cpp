// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#include "falco_gc_summary.hpp"

#include "falco_utils.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <compare>
#include <cstdint>
#include <ctime>  // for std::localtime
#include <format>
#include <iomanip>  // for std::put_time
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

[[nodiscard]] auto
get_theoretical_distribution(const std::vector<double> &gc,
                             const std::uint64_t total_count)
  -> std::vector<double> {
  static constexpr auto mode_width = 0.90;  // ADS: where does this come from?
  // calculate deviation of a hist from a "normal" with same mode and sd
  const auto n_bins = std::size(gc);
  const auto gc_beg = std::cbegin(gc);
  if (total_count <= 1)  // we will be dividing by (total_count - 1)
    return std::vector<double>(std::size(gc));

  // get mode
  const auto mode_itr = std::ranges::max_element(gc);
  const std::uint64_t mode_pos = std::distance(gc_beg, mode_itr);
  const auto mode_val = static_cast<double>(*mode_itr);

  // ADS: in case mode is not sharp average nearby values (not clear on why)
  const auto gt_cut = [&](const double x) { return x < mode_width * mode_val; };
  const auto right_itr = std::find_if(mode_itr, std::cend(gc), gt_cut);
  const auto left_itr =
    std::find_if(std::reverse_iterator(mode_itr), std::crend(gc), gt_cut);
  const auto n = std::distance(std::reverse_iterator(right_itr), left_itr);
  const auto mode = static_cast<double>(mode_pos) + as_frac(n - 1, 2.0);

  // theoretical distribution
  const auto cntr_sq = [m = mode](const auto v) { return (v - m) * (v - m); };
  const auto sd_term = [&](const auto &x) {
    return cntr_sq(std::get<0>(x)) * std::get<1>(x);
  };
  const auto id_gc =
    falco::views::enumerate(gc) | std::views::transform(sd_term);
  const auto sd = std::sqrt(as_frac(
    std::reduce(std::cbegin(id_gc), std::cend(id_gc)), total_count - 1));
  const auto to_normal = [&](const auto val) {
    return std::exp(-cntr_sq(val) / (2.0 * sd * sd));
  };
  auto normed = std::views::iota(0U, n_bins) |
                std::views::transform(to_normal) |
                std::ranges::to<std::vector>();
  const auto denom = std::reduce(std::cbegin(normed), std::cend(normed));
  std::ranges::transform(normed, std::begin(normed),
                         [&](const auto x) { return x * total_count / denom; });
  return normed;
}

[[nodiscard]] auto
sum_deviation_from_normal(const std::vector<double> &gc) -> double {
  const auto gc_beg = std::cbegin(gc);
  const auto gc_end = std::cend(gc);
  const auto total_count = std::reduce(gc_beg, gc_end);
  if (total_count <= 1)  // we will be dividing by (total_count - 1)
    return 0.0;
  const auto theor =
    get_theoretical_distribution(gc, static_cast<std::uint64_t>(total_count));
  const auto diff = [](const auto a, const auto b) { return std::fabs(a - b); };
  const auto r = std::transform_reduce(gc_beg, gc_end, std::cbegin(theor), 0.0,
                                       std::plus{}, diff);
  return pct(as_frac(r, total_count));
}

[[nodiscard]] auto
smooth_gc_content(const std::vector<double> &data,
                  const std::int64_t window_size) -> std::vector<double> {
  const auto get_mean = [&](std::ranges::viewable_range auto &&r) {
    return as_frac(std::reduce(std::cbegin(r), std::cend(r)), std::size(r));
  };
  assert(std::cmp_less(window_size, std::size(data)));
  std::vector<double> smoothed;
  for (auto w = 1; w < (window_size + 1) / 2; ++w)
    smoothed.push_back(get_mean(
      std::ranges::subrange(std::cbegin(data), std::cbegin(data) + w)));
#if __cpp_lib_ranges_slide
  for (const auto &window : data | std::views::slide(window_size))
    smoothed.push_back(get_mean(window));
#else
  const auto lim = std::cend(data) - window_size;
  for (auto itr = std::cbegin(data); itr != lim; ++itr)
    smoothed.push_back(get_mean(std::ranges::subrange(itr, itr + window_size)));
#endif
  const auto d_end = std::cend(data);
  for (auto w = (window_size + 1) / 2; w > 1; --w)
    smoothed.push_back(get_mean(std::ranges::subrange(d_end - w + 1, d_end)));
  return smoothed;
}

[[nodiscard]] static auto
get_last_bin(const std::int64_t histogram_size, const double next_percent)
  -> std::int64_t {
  const auto one_past_last = std::min(
    histogram_size, static_cast<std::int64_t>(std::ceil(next_percent)));
  return one_past_last - 1;
}

static auto
fill_histogram(
  const falco::gc_content_t &gc,
  const std::int32_t gc_idx,
  const double mass,
  std::vector<double> &hist  // cppcheck-suppress constParameterReference
  ) -> void {
  // ADS: below, not sure best way to do this for all edge cases
  assert(!hist.empty());
  const auto curr_percent = gc_idx * mass;
  const auto next_percent = (gc_idx + 1) * mass;
  const auto first_bin = static_cast<std::int64_t>(std::floor(curr_percent));
  const auto last_bin = get_last_bin(std::ssize(hist), next_percent);
  assert(std::cmp_less(last_bin, std::size(hist)));
  const auto has_splits = first_bin != last_bin;
  const auto mass_left =
    has_splits ? static_cast<double>(first_bin) + 1.0 - curr_percent : mass;
  const auto mass_right =
    has_splits ? next_percent - static_cast<double>(last_bin) : mass;
  const auto hist_begin = std::begin(hist) + first_bin;
  const auto hist_end = std::cbegin(hist) + last_bin;
  const auto gc_val = as_frac(gc[gc_idx], mass);
  for (auto hist_itr = hist_begin; hist_itr <= hist_end; ++hist_itr)
    *hist_itr += gc_val * ((hist_itr == hist_begin) ? mass_left
                           : (hist_itr == hist_end) ? mass_right
                                                    : 1.0);
}

[[nodiscard]] auto
combine_gc_content_for_lengths(const std::vector<falco::gc_content_t> &gcs)
  -> std::vector<double> {
  static constexpr auto histogram_size = 101;
  std::vector<double> hist(histogram_size);
  // ADS: iterate over each gc content vector (per length or length group)
  for (const auto &gc : gcs) {
    // ignore any length groups that have no reads
    if (std::ranges::none_of(gc, [](const auto x) { return x > 0.0; }))
      continue;
    // 'mass' is the fraction of each bin in the universal gc-content
    // histogram that is contributed by each full bin in the current
    // length-specific gc-content histogram
    const auto mass = as_frac(histogram_size, std::size(gc));
    for (auto gc_idx = 0; gc_idx < std::ssize(gc); ++gc_idx)
      fill_histogram(gc, gc_idx, mass, hist);
  }
  return hist;
}
