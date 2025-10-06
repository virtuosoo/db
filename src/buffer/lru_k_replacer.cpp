//===----------------------------------------------------------------------===//
//
//                         BusTub
//
// lru_k_replacer.cpp
//
// Identification: src/buffer/lru_k_replacer.cpp
//
// Copyright (c) 2015-2025, Carnegie Mellon University Database Group
//
//===----------------------------------------------------------------------===//

#include "buffer/lru_k_replacer.h"
#include "common/exception.h"

namespace bustub {

/**
 *
 * TODO(P1): Add implementation
 *
 * @brief a new LRUKReplacer.
 * @param num_frames the maximum number of frames the LRUReplacer will be required to store
 */
LRUKReplacer::LRUKReplacer(size_t num_frames, size_t k) : num_frames_(num_frames), k_(k) {}

/**
 * TODO(P1): Add implementation
 *
 * @brief Find the frame with largest backward k-distance and evict that frame. Only frames
 * that are marked as 'evictable' are candidates for eviction.
 *
 * A frame with less than k historical references is given +inf as its backward k-distance.
 * If multiple frames have inf backward k-distance, then evict frame whose oldest timestamp
 * is furthest in the past.
 *
 * Successful eviction of a frame should decrement the size of replacer and remove the frame's
 * access history.
 *
 * @return the frame ID if a frame is successfully evicted, or `std::nullopt` if no frames can be evicted.
 */
auto LRUKReplacer::Evict() -> std::optional<frame_id_t> { 
    std::scoped_lock lock(latch_);

    if (evictable_nodes_.empty()) {
        return std::nullopt; 
    }

    auto it = evictable_nodes_.begin();
    frame_id_t frame_id = it->fid_;
    evictable_nodes_.erase(it);
    node_store_.erase(frame_id);
    return frame_id;
}

/**
 * TODO(P1): Add implementation
 *
 * @brief Record the event that the given frame id is accessed at current timestamp.
 * Create a new entry for access history if frame id has not been seen before.
 *
 * If frame id is invalid (ie. larger than replacer_size_), throw an exception. You can
 * also use BUSTUB_ASSERT to abort the process if frame id is invalid.
 *
 * @param frame_id id of frame that received a new access.
 * @param access_type type of access that was received. This parameter is only needed for
 * leaderboard tests.
 */
void LRUKReplacer::RecordAccess(frame_id_t frame_id, [[maybe_unused]] AccessType access_type) {
    std::scoped_lock lock(latch_);
    BUSTUB_ASSERT(static_cast<size_t>(frame_id) <= num_frames_, "frame_id must <= nums_frames");
    
    current_timestamp_++;
    auto it = node_store_.find(frame_id);

    if (it == node_store_.end()) {
        node_store_.emplace(frame_id, LRUKNode(frame_id, k_, current_timestamp_));
    } else {
        LRUKNode& node = it->second;

        // 如果节点是可淘汰的，先从 set 中移除旧状态的它
        if (node.is_evictable_) {
            evictable_nodes_.erase(node);
        }

        // 直接在 map 的节点上进行修改
        node.access_count_++;
        node.last_access_timestamp_ = current_timestamp_;

        // 如果节点是可淘汰的，将更新后的它插入 set
        if (node.is_evictable_) {
            evictable_nodes_.insert(node);
        }
    }
}

/**
 * TODO(P1): Add implementation
 *
 * @brief Toggle whether a frame is evictable or non-evictable. This function also
 * controls replacer's size. Note that size is equal to number of evictable entries.
 *
 * If a frame was previously evictable and is to be set to non-evictable, then size should
 * decrement. If a frame was previously non-evictable and is to be set to evictable,
 * then size should increment.
 *
 * If frame id is invalid, throw an exception or abort the process.
 *
 * For other scenarios, this function should terminate without modifying anything.
 *
 * @param frame_id id of frame whose 'evictable' status will be modified
 * @param set_evictable whether the given frame is evictable or not
 */
void LRUKReplacer::SetEvictable(frame_id_t frame_id, bool set_evictable) {
    std::scoped_lock lock(latch_);
    BUSTUB_ASSERT(static_cast<size_t>(frame_id) <= num_frames_, "frame_id must <= nums_frames");

    auto it = node_store_.find(frame_id);
    if (it == node_store_.end()) {
        return;
    }

    LRUKNode &node = it->second;
    if (node.is_evictable_ == set_evictable) {
        return;
    }

    node.is_evictable_ = set_evictable;
    if (set_evictable) {
        evictable_nodes_.insert(node);
    } else {
        evictable_nodes_.erase(node);
    }
}

/**
 * TODO(P1): Add implementation
 *
 * @brief Remove an evictable frame from replacer, along with its access history.
 * This function should also decrement replacer's size if removal is successful.
 *
 * Note that this is different from evicting a frame, which always remove the frame
 * with largest backward k-distance. This function removes specified frame id,
 * no matter what its backward k-distance is.
 *
 * If Remove is called on a non-evictable frame, throw an exception or abort the
 * process.
 *
 * If specified frame is not found, directly return from this function.
 *
 * @param frame_id id of frame to be removed
 */
void LRUKReplacer::Remove(frame_id_t frame_id) {
    std::scoped_lock lock(latch_);
    BUSTUB_ASSERT(static_cast<size_t>(frame_id) <= num_frames_, "frame_id must <= nums_frames");

    auto it = node_store_.find(frame_id);
    if (it == node_store_.end()) {
        return;
    }
    
    LRUKNode& node = it->second;
    BUSTUB_ASSERT(!node.is_evictable_, "frame must be evictable");
    
    node_store_.erase(frame_id);
    evictable_nodes_.erase(node);
}

/**
 * TODO(P1): Add implementation
 *
 * @brief Return replacer's size, which tracks the number of evictable frames.
 *
 * @return size_t
 */
auto LRUKReplacer::Size() -> size_t { return evictable_nodes_.size(); }

}  // namespace bustub
