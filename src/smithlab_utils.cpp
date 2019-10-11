/*
 *    Part of SMITHLAB software
 *
 *    Copyright (C) 2011 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "smithlab_utils.hpp"
#include <cstring>
#include <string>
#include <cmath>

using std::vector;
using std::string;

char have_smithlab_cpp() {return 1;}

/* Performs Hochberg step-up p-value adjustment.
 */
void smithlab::correct_pvals(const size_t n_tests, vector<double> &pvals ) {
  assert(!pvals.empty());

  vector<std::pair<double, vector<size_t> > > idx;
  for (size_t i = 0; i < pvals.size(); ++i)
    idx.push_back(std::make_pair(pvals[i], vector<size_t>(1, i)));
  std::sort(idx.begin(), idx.end());

  size_t j = 0;
  for (size_t i = 1; i < idx.size(); ++i)
    if (idx[j].first == idx[i].first)
      idx[j].second.push_back(idx[i].second.front());
    else std::swap(idx[++j], idx[i]);
  idx.erase(idx.begin() + j + 1, idx.end());

  double running_count = idx.front().second.size();
  idx.front().first *= (n_tests/running_count);
  for (size_t i = 1; i < idx.size(); ++i) {
    running_count += idx[i].second.size();
    idx[i].first = std::max(idx[i-1].first, idx[i].first*(n_tests/running_count));
  }

  for (size_t i = 0; i < idx.size(); ++i) {
    for (size_t j = 0; j < idx[i].second.size(); ++j)
      pvals[idx[i].second[j]] = idx[i].first;
  }

}

/* Benjamini-Hochberg step-up FDR cutoff
 */
double
smithlab::get_fdr_cutoff(const size_t n_tests, vector<double> &pvals,
                         const double alpha) {
  if (alpha <= 0) return std::numeric_limits<double>::max();
  else if (alpha > 1) return std::numeric_limits<double>::min();

  std::sort(pvals.begin(), pvals.end());

  const size_t n_pvals = pvals.size();
  size_t i = 0;
  while (i < n_pvals - 1 && pvals[i+1] <
         alpha*static_cast<double>(i+1)/n_tests)
    ++i;
  return pvals[i];
}

char
complement(int i) {
  static const int b2c_size = 20;
  static const char b2c[] = {
    //A,  b,  C,  d,  e,  f,  g,  h,  i,  j,  k,  l,  m,  n,  o,  p,  q,  r,  s,  T
    'T','N','G','N','N','N','C','N','N','N','N','N','N','N','N','N','N','N','N','A'
  };
  static const char b2cl[] = {
    //A,  b,  C,  d,  e,  f,  g,  h,  i,  j,  k,  l,  m,  n,  o,  p,  q,  r,  s,  T
    't','n','g','n','n','n','c','n','n','n','n','n','n','n','n','n','n','n','n','a'
  };
  if (i - 'A' >= 0 && i - 'A' < b2c_size)
    return b2c[i - 'A'];
  else if (i - 'a' >= 0 && i - 'a' < b2c_size)
    return b2cl[i - 'a'];
  else return 'N';
}


std::vector<std::string>
smithlab::split(std::string s, const char *delim, bool get_empty_fields) {
  std::vector <std::string> parts;
  size_t i = 0, j = 0, dlen = strlen(delim);
  while (i < s.length()) {
    bool prev_in_delim = false;
    if (i > 0)
      for (size_t k = 0; k < dlen; ++k)
        if (s[i-1] == delim[k])
          prev_in_delim = true;
    bool curr_in_delim = false;
    for (size_t k = 0; k < dlen; ++k)
      if (s[i] == delim[k])
        curr_in_delim = true;
    if (curr_in_delim && (get_empty_fields || !prev_in_delim)) {
      parts.push_back(s.substr(j, i - j));
      j = i + 1;
    }
    i++;
  }
  if (i > j || (get_empty_fields && i == j))
    parts.push_back(s.substr(j, i - j));
  return parts;
}

std::string
smithlab::strip(const std::string& s) {
  const size_t len = s.length();
  size_t i = 0;
  while (i < len && isspace(s[i])) {
    i++;
  }
  size_t j = len;
  do {
    j--;
  } while (j >= i && isspace(s[j]));
  j++;
  if (i == 0 && j == len)
    return s;
  else return s.substr(i, j - i);
}

void
smithlab::split_whitespace(const string &s, vector<string> &v) {
  v.clear();
  std::istringstream iss(s);
  copy(std::istream_iterator<string>(iss),
       std::istream_iterator<string>(),
       std::back_inserter(v));
}

std::vector<std::string>
smithlab::split_whitespace_quoted(std::string to_split) {
  static const char *non_word_chars = " \t\"'";
  to_split = smithlab::strip(to_split);

  std::vector<std::string> words;
  size_t start_pos = 0, end_pos = 0;
  const size_t length_of_to_split = to_split.length();
  while (start_pos < length_of_to_split) {
    /** find next position that is not a word character */
    end_pos = to_split.find_first_of(non_word_chars, end_pos);
    if (end_pos > to_split.length()) { /** If we hit the end: done */
      words.push_back(to_split.substr(start_pos, end_pos - start_pos));
      break;
    }
    /** unescaped, unquoted white space: definitely a word delimiter */
    if (to_split[end_pos] == ' ' || to_split[end_pos] == '\t') {
      words.push_back(to_split.substr(start_pos, end_pos - start_pos));
      end_pos = to_split.find_first_not_of(" \t", end_pos);
      start_pos = end_pos;
    }
    /** preserve whatever is being escaped; will become part of the
        current word */
    else if (to_split[end_pos] == '\\')
      end_pos = to_split.find_first_not_of(non_word_chars, end_pos + 2);
    else {
      const std::string current_delim = "\\" + to_split.substr(end_pos, 1);
      do { // slurp doubly- or singly-quoted string
        end_pos = to_split.find_first_of(current_delim, end_pos + 1);
        if (end_pos == std::string::npos) {
          end_pos = length_of_to_split;
          break;
        }
        if (to_split[end_pos] == '\\')
          ++end_pos;
        else break;
      } while (true);
      ++end_pos;
    }
    if (end_pos >= length_of_to_split) {
      words.push_back(to_split.substr(start_pos,
                                      end_pos - start_pos));
      break;
    }
  }
  return words;
}

double
smithlab::log_sum_log_vec(const std::vector<double> &vals, size_t limit) {
  const std::vector<double>::const_iterator x =
    max_element(vals.begin(), vals.begin() + limit);
  const double max_val = *x;
  const size_t max_idx = x - vals.begin();
  double sum = 1.0;
  for (size_t i = 0; i < limit; ++i)
    if (i != max_idx)
      sum += std::exp(vals[i] - max_val);
  return max_val + std::log(sum);
}

/****
 * @summary: remove empty (or only whitespace) strings from a vector of string
 */
std::vector<std::string>
smithlab::squash(const std::vector<std::string>& v) {
  std::vector<std::string> res;
  for (size_t i = 0; i < v.size(); i++) {
    std::string t = v[i];
    strip(t);
    if (t != "") res.push_back(t);
  }
  return res;
}

void
ProgressBar::report(std::ostream &out, const size_t i) {
  prev = std::round((100.0*i)/total);
  const size_t x =
    std::min(static_cast<size_t>(bar_width*(prev/100.0)), bar_width);
  fill_n(begin(bar), x, '=');
  out << left_tag << mid_tag << "|" << bar << "|"
      << std::setw(3) << prev << right_tag;
  if (i >= total)
    out << '\n';
}
