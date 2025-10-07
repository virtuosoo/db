//===----------------------------------------------------------------------===//
//
//                         BusTub
//
// lru_k_replacer.h
//
// Identification: src/include/buffer/lru_k_replacer.h
//
// Copyright (c) 2015-2025, Carnegie Mellon University Database Group
//
//===----------------------------------------------------------------------===//

#pragma once

#include <limits>
#include <list>
#include <mutex>  // NOLINT
#include <optional>
#include <set>
#include <unordered_map>
#include <vector>

#include "common/config.h"
#include "common/macros.h"

namespace bustub {

enum class AccessType { Unknown = 0, Lookup, Scan, Index };

class LRUKNode {
 public:
  LRUKNode(frame_id_t fid, size_t k, size_t time_stamp) : fid_(fid), last_access_timestamp_(time_stamp), k_(k) {}

  auto operator<(const LRUKNode &other) const -> bool {
    auto this_key = GetSortingKey();
    auto other_key = other.GetSortingKey();
    if (this_key != other_key) {
      return this_key < other_key;
    }

    return last_access_timestamp_ < other.last_access_timestamp_;
  }
  /** History of last seen K timestamps of this page. Least recent timestamp stored in front. */
  // Remove maybe_unused if you start using them. Feel free to change the member variables as you want.

  frame_id_t fid_;
  size_t last_access_timestamp_;
  size_t k_;
  bool is_evictable_{false};
  size_t access_count_{1};

  auto GetSortingKey() const -> int64_t {
    if (access_count_ < k_) {
      return 0;
    } 
    return last_access_timestamp_;
  }
};

/**
 * LRUKReplacer implements the LRU-k replacement policy.
 *
 * The LRU-k algorithm evicts a frame whose backward k-distance is maximum
 * of all frames. Backward k-distance is computed as the difference in time between
 * current timestamp and the timestamp of kth previous access.
 *
 * A frame with less than k historical references is given
 * +inf as its backward k-distance. When multiple frames have +inf backward k-distance,
 * classical LRU algorithm is used to choose victim.
 */
class LRUKReplacer {
 public:
  explicit LRUKReplacer(size_t num_frames, size_t k);

  DISALLOW_COPY_AND_MOVE(LRUKReplacer);

  /**
   * TODO(P1): Add implementation
   *
   * @brief Destroys the LRUReplacer.
   */
  ~LRUKReplacer() = default;

  auto Evict() -> std::optional<frame_id_t>;

  void RecordAccess(frame_id_t frame_id, AccessType access_type = AccessType::Unknown);

  void SetEvictable(frame_id_t frame_id, bool set_evictable);

  void Remove(frame_id_t frame_id);

  auto Size() -> size_t;

 private:
  // TODO(student): implement me! You can replace these member variables as you like.
  // Remove maybe_unused if you start using them.
  std::unordered_map<frame_id_t, LRUKNode> node_store_;
  std::set<LRUKNode> evictable_nodes_;
  size_t current_timestamp_{0};
  [[maybe_unused]] size_t replacer_size_;
  size_t num_frames_;
  size_t k_;
  std::mutex latch_;
};

}  // namespace bustub
