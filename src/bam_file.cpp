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

#include "bam_file.hpp"
#include "falco_utils.hpp"

#include <htslib/bgzf.h>   // for BGZF
#include <htslib/hfile.h>  // for htell

#include <cstdint>
#include <filesystem>
#include <memory>
#include <string>
#include <system_error>

[[nodiscard]] auto
estimate_n_reads_bam(const std::string &filename) -> std::uint64_t {
  static constexpr auto max_n_reads = 128 * 1024;
  std::unique_ptr<htsFile, int (*)(htsFile *)> f(
    hts_open(std::data(filename), "r"), &hts_close);
  if (!f)
    throw std::system_error(std::make_error_code(std::errc(errno)),
                            "failed to open file: " + filename);
  std::unique_ptr<sam_hdr_t, void (*)(sam_hdr_t *)> h(sam_hdr_read(f.get()),
                                                      &sam_hdr_destroy);
  if (!h)
    throw std::system_error(std::make_error_code(std::errc(errno)),
                            "failed to read header: " + filename);

  const auto format = hts_get_format(f.get());
  if (!format)
    throw std::runtime_error("failed to identify file format: " + filename);

  // NOLINTNEXTLINE (cppcoreguidelines-pro-type-union-access)
  const auto &fp = format->format == bam ? f->fp.bgzf->fp : f->fp.hfile;
  const auto pos_after_header = htell(fp);
  bam1_t rec{};
  std::uint64_t n_reads{};
  int r{};
  while (n_reads++ < max_n_reads &&
         (r = sam_read1(f.get(), h.get(), &rec)) >= 0)
    ;
  if (r < -1)  // error
    throw std::system_error(std::make_error_code(std::errc(errno)),
                            "error reading bam record from: " + filename);
  const auto pos_after_reads = htell(fp);
  const auto n_compressed_bytes = pos_after_reads - pos_after_header;
  const auto filesize = std::filesystem::file_size(filename);
  return static_cast<std::uint64_t>(
    as_frac(n_reads * (filesize - pos_after_header), n_compressed_bytes));
}
