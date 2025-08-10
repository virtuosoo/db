//===----------------------------------------------------------------------===//
//
//                         BusTub
//
// hyperloglog.cpp
//
// Identification: src/primer/hyperloglog.cpp
//
// Copyright (c) 2015-2025, Carnegie Mellon University Database Group
//
//===----------------------------------------------------------------------===//

#include "primer/hyperloglog.h"
#define MAX_NBITS 16

namespace bustub {

/** @brief Parameterized constructor. */
template <typename KeyType>
HyperLogLog<KeyType>::HyperLogLog(int16_t n_bits) : cardinality_(0) {
  if (n_bits < 0) {
    n_bits = 0;
  }
  n_bits_ = n_bits;
  bucket_count_ = 1 << n_bits;
  bucket_values_.resize(bucket_count_, 0);
}

/**
 * @brief Function that computes binary.
 *
 * @param[in] hash
 * @returns binary of a given hash
 */
template <typename KeyType>
auto HyperLogLog<KeyType>::ComputeBinary(const hash_t &hash) const -> std::bitset<BITSET_CAPACITY> {
  std::bitset<BITSET_CAPACITY> bset;
  for (size_t i = 0; i < BITSET_CAPACITY; i++) {
    if ((hash >> i) & 1) {
      bset[i] = 1;
    }
  }
  return bset;
}

/**
 * @brief Function that computes leading zeros.
 *
 * @param[in] bset - binary values of a given bitset
 * @returns leading zeros of given binary set
 */
template <typename KeyType>
auto HyperLogLog<KeyType>::PositionOfLeftmostOne(const std::bitset<BITSET_CAPACITY> &bset) const -> uint64_t {
  uint64_t cnt = 0;
  for (int i = BITSET_CAPACITY - 1 - n_bits_; i >= 0; i--) {
    if (bset[i] == 0) {
      cnt++;
    } else {
      break;
    }
  }
  return cnt + 1;
}

/**
 * @brief Adds a value into the HyperLogLog.
 *
 * @param[in] val - value that's added into hyperloglog
 */
template <typename KeyType>
auto HyperLogLog<KeyType>::AddElem(KeyType val) -> void {
  std::scoped_lock slk(m);

  hash_t hash = CalculateHash(val);
  std::bitset<BITSET_CAPACITY> bset = ComputeBinary(hash);
  uint64_t leftmost1 = PositionOfLeftmostOne(bset);
  
  size_t index = 0;
  for (int i = BITSET_CAPACITY - 1; i >= BITSET_CAPACITY - n_bits_; i--) {
    index += (1 << (n_bits_ - (BITSET_CAPACITY - i))) * bset[i];
  }
  bucket_values_[index] = std::max(bucket_values_[index], leftmost1);
}

/**
 * @brief Function that computes cardinality.
 */
template <typename KeyType>
auto HyperLogLog<KeyType>::ComputeCardinality() -> void {
  std::scoped_lock slk(m);

  double sum = 0.0;
  for (size_t i = 0; i < bucket_count_; i++) {
    sum += std::pow(2.0, -1.0 * bucket_values_[i]);
  }

  cardinality_ = std::floor(CONSTANT * bucket_count_ * bucket_count_ / sum);
}

template class HyperLogLog<int64_t>;
template class HyperLogLog<std::string>;

}  // namespace bustub
