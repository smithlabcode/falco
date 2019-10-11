#ifndef _AUX_HPP
#define _AUX_HPP
#include <string>
#include <chrono>
#include <ctime>
#include <cmath>
#include <iostream>
#include <algorithm>

/*************************************************************
 ******************** AUX FUNCTIONS **************************
 *************************************************************/
// converts 64 bit integer to a sequence std::string by reading 2 bits at a time and
// converting back to ACTG
static inline std::string
size_t_to_seq(size_t v, const size_t seq_length) {
  std::string ans;
  for (size_t i = 0; i < seq_length; ++i) {
    switch (v & 3) {
      case 0: ans.push_back('A'); break;
      case 1: ans.push_back('C'); break;
      case 2: ans.push_back('T'); break;
      case 3: ans.push_back('G'); break;
    }   
    v >>= 2;
  }

  std::reverse(ans.begin(), ans.end());
  return ans;
}

// Converts A,T,G,C to 2-bit values
static inline size_t
actg_to_2bit(const char &c) {
  return ((c >> 1) & 3);
}

#endif
