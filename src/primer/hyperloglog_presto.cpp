//===----------------------------------------------------------------------===//
//
//                         BusTub
//
// hyperloglog_presto.cpp
//
// Identification: src/primer/hyperloglog_presto.cpp
//
// Copyright (c) 2015-2025, Carnegie Mellon University Database Group
//
//===----------------------------------------------------------------------===//

#include "primer/hyperloglog_presto.h"

namespace bustub {

/** @brief Parameterized constructor. */
template <typename KeyType>
HyperLogLogPresto<KeyType>::HyperLogLogPresto(int16_t n_leading_bits) : cardinality_(0), baseline_(0) {
  if (n_leading_bits < 0) {
    n_leading_bits = 0;
  }
  n_leading_bits_ = n_leading_bits;
  bucket_count_ = 1 << n_leading_bits_;
  dense_bucket_.resize(1 << n_leading_bits, std::bitset<DENSE_BUCKET_SIZE>());
}

/**
 * @brief Function that computes binary.
 *
 * @param[in] hash
 * @returns binary of a given hash
 */
template <typename KeyType>
auto HyperLogLogPresto<KeyType>::ComputeBinary(const hash_t &hash) const -> std::bitset<BITSET_CAPACITY> {
  std::bitset<BITSET_CAPACITY> bset;
  for (size_t i = 0; i < BITSET_CAPACITY; i++) {
    if ((hash >> i) & 1) {
      bset[i] = 1;
    }
  }
  return bset;
}

/** @brief Function that finds the position of leftmost one.
 *
 * @param[in] bset
 * @returns position of leftmost one
 */
template <typename KeyType>
auto HyperLogLogPresto<KeyType>::NumberOfRightmostZero(const std::bitset<BITSET_CAPACITY> &bset) const -> uint64_t {
  u_int64_t res = 0;
  for (int i = 0; i < BITSET_CAPACITY - n_leading_bits_; i++) {
    if (bset[i] == 0) {
      res++;
    } else {
      break;
    }
  }
  return res;
}

template <typename KeyType>
auto HyperLogLogPresto<KeyType>::GetDeltaOfIndex(uint16_t idx) const -> uint16_t {
  if (idx >= bucket_count_) {
    return 0;  // 索引越界保护
  }

  uint16_t delta = 0;
  auto denseBucket = dense_bucket_[idx];
  for (int i = 0; i < DENSE_BUCKET_SIZE; ++i) {
    delta += (1 << i) * denseBucket[i];
  }
  if (overflow_bucket_.find(idx) != overflow_bucket_.end()) {
    auto overflow = overflow_bucket_.at(idx);
    for (int i = 0; i < OVERFLOW_BUCKET_SIZE; ++i) {
      delta += (1 << (i + DENSE_BUCKET_SIZE)) * overflow[i];
    }
  }
  return delta;
}

template <typename KeyType>
auto HyperLogLogPresto<KeyType>::SetDeltaOfIndex(uint16_t idx, uint16_t delta) -> void {
  dense_bucket_[idx].reset();
  for (int i = 0; i < DENSE_BUCKET_SIZE; ++i) {
    if ((delta >> i) & 1) {
      dense_bucket_[idx][i] = 1;
    }
  }
  if (delta > MAX_DENSE_NUM) {
    std::bitset<OVERFLOW_BUCKET_SIZE> overflow;
    for (int i = 0; i < OVERFLOW_BUCKET_SIZE; ++i) {
      if ((delta >> (i + DENSE_BUCKET_SIZE)) & 1) {
        overflow[i] = 1;
      }
    }
    overflow_bucket_[idx] = overflow;
  } else {
    if (overflow_bucket_.find(idx) != overflow_bucket_.end()) {
      overflow_bucket_.erase(idx);
    }
  }
}

template <typename KeyType>
auto HyperLogLogPresto<KeyType>::RecalculateBaselineAndDelta() -> void {
  uint16_t new_baseline = UINT16_MAX;
  for (uint16_t i = 0; i < bucket_count_; ++i) {
    new_baseline = std::min(new_baseline, static_cast<uint16_t>(GetDeltaOfIndex(i) + baseline_));
  }
  if (new_baseline > baseline_) {
    for (uint16_t i = 0; i < bucket_count_; ++i) {
      SetDeltaOfIndex(i, GetDeltaOfIndex(i) + baseline_ - new_baseline);
    }
    baseline_ = new_baseline;
  }
}

/** @brief Element is added for HLL calculation. */
template <typename KeyType>
auto HyperLogLogPresto<KeyType>::AddElem(KeyType val) -> void {
  std::scoped_lock slk(m);

  hash_t hash = CalculateHash(val);
  std::bitset<BITSET_CAPACITY> bset = ComputeBinary(hash);
  uint64_t zeros = NumberOfRightmostZero(bset);

  size_t index = 0;
  for (int i = BITSET_CAPACITY - 1; i >= BITSET_CAPACITY - n_leading_bits_; i--) {
    index += (1 << (n_leading_bits_ - (BITSET_CAPACITY - i))) * bset[i];
  }

  std::cout << "add elem, index: " << index << " hash: " << hash << " zeros: " << zeros << '\n';

  if (zeros > GetDeltaOfIndex(index) + baseline_) {
    SetDeltaOfIndex(index, zeros - baseline_);
    // RecalculateBaselineAndDelta();
  }
}

/** @brief Function to compute cardinality. */
template <typename T>
auto HyperLogLogPresto<T>::ComputeCardinality() -> void {
  std::scoped_lock slk(m);

  double sum = 0.0;
  for (uint16_t i = 0; i < bucket_count_; i++) {
    sum += std::pow(2.0, -1.0 * (GetDeltaOfIndex(i) + baseline_));
  }
  cardinality_ = std::floor(CONSTANT * bucket_count_ * bucket_count_ / sum);
}

template class HyperLogLogPresto<int64_t>;
template class HyperLogLogPresto<std::string>;
}  // namespace bustub
