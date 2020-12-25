// This file is part of PGM-index <https://github.com/gvinciguerra/PGM-index>.
// Copyright (c) 2018 Giorgio Vinciguerra.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#pragma once

#include "piecewise_linear_model.hpp"
#include <algorithm>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

namespace pgm {

#define PGM_SUB_EPS(x, epsilon) ((x) <= (epsilon) ? 0 : ((x) - (epsilon)))
#define PGM_ADD_EPS(x, epsilon, size) ((x) + (epsilon) + 2 >= (size) ? (size) : (x) + (epsilon) + 2)

/**
 * A struct that stores the result of a query to a @ref PGMIndex, that is, a range [@ref lo, @ref hi)
 * centered around an approximate position @ref pos of the sought key.
 */
struct ApproxPos {
  size_t pos; ///< The approximate position of the key.
  size_t lo;  ///< The lower bound of the range.
  size_t hi;  ///< The upper bound of the range.
};

/**
 * A space-efficient index that enables fast search operations on a sorted sequence of @c n numbers.
 *
 * A search returns a struct @ref ApproxPos containing an approximate position of the sought key in the sequence and
 * the bounds of a range of size 2*Epsilon+1 where the sought key is guaranteed to be found if present.
 * If the key is not present, the range is guaranteed to contain a key that is not less than (i.e. greater or equal to)
 * the sought key, or @c n if no such key is found.
 * In the case of repeated keys, the index finds the position of the first occurrence of a key.
 *
 * The @p Epsilon template parameter should be set according to the desired space-time trade-off. A smaller value
 * makes the estimation more precise and the range smaller but at the cost of increased space usage.
 *
 * Internally the index uses a succinct piecewise linear mapping from keys to their position in the sorted order.
 * This mapping is represented as a sequence of linear models (segments) which, if @p EpsilonRecursive is not zero, are
 * themselves recursively indexed by other piecewise linear mappings.
 *
 * @tparam K the type of the indexed keys
 * @tparam Epsilon controls the size of the returned search range
 * @tparam EpsilonRecursive controls the size of the search range in the internal structure
 */
template <typename K, size_t Epsilon = 64, size_t EpsilonRecursive = 4> class PGMIndex {
  static_assert(Epsilon > 0);
  static_assert(EpsilonRecursive > 0);
  struct Segment;

  size_t n;                           ///< The number of elements this index was built on.
  K first_key;                        ///< The smallest element.
  std::vector<Segment> segments;      ///< The segments composing the index, stored contiguously.
  std::vector<size_t> levels_offsets; ///< The starting position of each level in segments[].
  const std::vector<K> &data;         ///< The data vector this index was built on.

  template <typename RandomIt> void build(RandomIt first, RandomIt last, size_t epsilon, size_t epsilon_recursive) {
    if (n == 0)
      return;

    first_key = *first;
    levels_offsets.push_back(0);
    segments.reserve(n / (epsilon * epsilon));

    // The back_check function pushes a sentinel segment with a large key at the end of each level. This avoids costly
    // bounds checking / branches in segment_for_key()
    auto back_check = [this, last](size_t n_segments, size_t prev_level_size) {
      if (segments.back().slope == 0) {
        // Technicality to ensure that keys > *(last-1) are approximated to a position == prev_level_size
        segments.emplace_back(*std::prev(last) + 1, 0, prev_level_size);
        ++n_segments;
      }

      // Add the sentinel segment with intercept = prev_level_size
      segments.emplace_back(std::numeric_limits<K>::max(), 0, prev_level_size);
      return n_segments;
    };

    // Build first level
    auto in_fun = [this, first](auto i) {
      auto x = first[i];
      if (i > 0 && i + 1u < n && x == first[i - 1] && x != first[i + 1] && x + 1 != first[i + 1])
        return std::pair<K, size_t>(x + 1, i); // technicality to predict the correct position when keys are repeated
      return std::pair<K, size_t>(x, i);
    };
    auto out_fun = [this](auto cs) { segments.emplace_back(cs); };
    auto prev_level_size = n;
    auto n_segments = make_segmentation_par(prev_level_size, epsilon, in_fun, out_fun);
    prev_level_size = back_check(n_segments, prev_level_size);
    levels_offsets.push_back(levels_offsets.back() + prev_level_size + 1);

    // Build upper levels
    while (prev_level_size > 1) {
      auto offset = levels_offsets[levels_offsets.size() - 2];
      auto in_fun_rec = [this, offset](auto i) { return std::pair<K, size_t>(segments[offset + i].key, i); };
      n_segments = make_segmentation(prev_level_size, epsilon_recursive, in_fun_rec, out_fun);
      prev_level_size = back_check(n_segments, prev_level_size);
      levels_offsets.push_back(levels_offsets.back() + prev_level_size + 1);
    }

    levels_offsets.pop_back();
    std::reverse(levels_offsets.begin(), levels_offsets.end());
  }

  /** Returns the segment responsible for a given key, that is, the rightmost segment having key <= the sought key. */
  auto segment_for_key(const K &key) const {
    auto it = segments.begin() + levels_offsets.front();

    for (auto l = 1; l < height(); ++l) {
      auto level_begin = segments.begin() + levels_offsets[l];
      auto pos = std::min<size_t>((*it)(key), std::next(it)->intercept);
      auto lo = level_begin + PGM_SUB_EPS(pos, EpsilonRecursive + 1);
      for (; std::next(lo)->key <= key; ++lo) // no need for bounds checking: each level has a final sentinel segment
        continue;
      it = lo;
    }

    return it;
  }

public:
  /** Constructs the index on the given sorted keys. */
  PGMIndex(const std::vector<K> &data) : n(data.size()), first_key(), segments(), levels_offsets(), data(data) {
    build(data.begin(), data.end(), Epsilon, EpsilonRecursive);
  }

  /** Returns the approximate position and the range where @p key can be found. */
  ApproxPos search(const K &key) const {
    auto k = std::max(first_key, key);
    auto it = segment_for_key(k);
    auto pos = std::min<size_t>((*it)(k), std::next(it)->intercept);
    auto lo = PGM_SUB_EPS(pos, Epsilon);
    auto hi = PGM_ADD_EPS(pos, Epsilon, n);
    return {pos, lo, hi};
  }

  K nextGEQ(const K &key) const { // This function is required for the challenge
    auto range = search(key);
    auto lo = data.begin() + range.lo;
    auto hi = data.begin() + range.hi;
    return *std::lower_bound(lo, hi, key);
  }

  /*  Returns the number of levels of the index. */
  size_t height() const { return levels_offsets.size(); }

  /* Returns the size of the index in bytes. */
  size_t size_in_bytes() const { // This function is required for the challenge
    return sizeof(this) + segments.size() * sizeof(Segment);
  }
};

#pragma pack(push, 1)

template <typename K, size_t Epsilon, size_t EpsilonRecursive> struct PGMIndex<K, Epsilon, EpsilonRecursive>::Segment {
  K key;             ///< The first key that the segment indexes.
  double slope;      ///< The slope of the segment.
  int32_t intercept; ///< The intercept of the segment.

  Segment() = default;

  Segment(K key, double slope, double intercept) : key(key), slope(slope), intercept(intercept){};

  explicit Segment(const typename OptimalPiecewiseLinearModel<K, size_t>::CanonicalSegment &cs)
      : key(cs.get_first_x()) {
    auto [cs_slope, cs_intercept] = cs.get_floating_point_segment(key);
    slope = cs_slope;
    intercept = std::round(cs_intercept);
  }

  friend inline bool operator<(const K &k, const Segment &s) { return k < s.key; }

  /** Returns the approximate position of the specified key. */
  inline size_t operator()(const K &k) const {
    auto pos = int64_t(slope * (k - key)) + intercept;
    return pos > 0 ? size_t(pos) : 0ull;
  }
};

#pragma pack(pop)
}