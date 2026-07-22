// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#include "nlohmann/json.hpp"

#include <cstdint>
#include <string>

struct bam_header {
  using iterator = char *;
  using const_iterator = const char *;

  static constexpr auto magic = "BAM\1";
  static constexpr auto magic_size = 4;
  std::uint8_t magic_bytes_remaining{4};
  std::uint8_t l_text_bytes_remaining{4};
  std::uint32_t l_text{};
  std::uint8_t n_ref_bytes_remaining{4};
  std::uint32_t n_ref{};
  std::uint8_t l_name_bytes_remaining{4};
  std::uint32_t name_bytes_remaining{};
  std::uint8_t l_ref_bytes_remaining{4};

  [[nodiscard]] auto
  to_string() const -> std::string {
    static constexpr auto n_indent = 4;
    return nlohmann::json(*this).dump(n_indent);
  }

  [[nodiscard]] auto
  update(iterator itr, const iterator end) -> iterator;

  [[nodiscard]] auto
  ref_incomplete() const -> bool {
    return l_name_bytes_remaining || name_bytes_remaining ||
           l_ref_bytes_remaining;
  }

  [[nodiscard]] auto
  ref_complete() const -> bool {
    return !ref_incomplete();
  }

  operator bool() const { return ref_incomplete() || n_ref; }

  auto
  decrement_ref() -> void {
    --n_ref;
    if (n_ref) {
      l_name_bytes_remaining = 4;
      name_bytes_remaining = 0;
      l_ref_bytes_remaining = 4;
    }
  }

  NLOHMANN_DEFINE_TYPE_INTRUSIVE(bam_header, magic_bytes_remaining,
                                 l_text_bytes_remaining, l_text,
                                 n_ref_bytes_remaining, n_ref,
                                 l_name_bytes_remaining, name_bytes_remaining,
                                 l_ref_bytes_remaining);
};
