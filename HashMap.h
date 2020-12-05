//
// Created by Oleg Antokhiv on 26.11.2020.
//

#ifndef HASHMAP__HASHMAP_H_
#define HASHMAP__HASHMAP_H_

#include <functional>
#include <vector>
#include <cmath>
#include <limits>
template<class Key, class T, class Hash = std::hash<Key>,
    class KeyEqual = std::equal_to<void>,
    class Allocator = std::allocator<std::pair<Key, T>>>
class HashMap {
 public:
  using key_type = Key;
  using mapper_type = T;
  using value_type = std::pair<Key, T>;
  using size_type = std::size_t;
  using hasher = Hash;
  using key_equal = KeyEqual;
  using allocator_type = Allocator;
  using buckets = std::vector<value_type, allocator_type>;

  template<typename IterVal>
  struct hm_iterator {
    using difference_type = std::ptrdiff_t;
    using value_type = IterVal;
    using pointer = value_type *;
    using reference = value_type &;
    using iterator_category = std::forward_iterator_tag;

    bool operator==(const hm_iterator &other) const {
      return other.hm_ == hm_ && other.idx_ == idx_;
    }
    bool operator!=(const hm_iterator &other) const {
      return !(other == *this);
    }

    hm_iterator &operator++() {
      ++idx_;
      advance_past_empty();
      return *this;
    }

    reference operator*() const { return hm_->buckets_[idx_]; }
    pointer operator->() const { return &hm_->buckets_[idx_]; }

   private:
    explicit hm_iterator(HashMap *hm) : hm_(hm) { advance_past_empty(); }
    explicit hm_iterator(HashMap *hm, size_type idx) : hm_(hm), idx_(idx) {}
    template<typename OtherHashMap, typename OtherIterVal>
    hm_iterator(const hm_iterator<OtherIterVal> &other)
        : hm_(other.hm_), idx_(other.idx_) {}

    void advance_past_empty() {
      while (idx_ < hm_->buckets_.size() &&
          key_equal()(hm_->buckets_[idx_].first, hm_->empty_key_)) {
        ++idx_;
      }
    }

    HashMap *hm_ = nullptr;
    typename HashMap::size_type idx_ = 0;

  };
  using iterator = hm_iterator<value_type>;
  using const_iterator = hm_iterator<const value_type>;
 public:
  HashMap(size_type bucket_count, key_type empty_key,
          const allocator_type &alloc = allocator_type())
      : empty_key_(empty_key), buckets_(alloc) {
    buckets_.reserve(bucket_count);
    capacity = bucket_count;
    buckets_.resize(bucket_count, std::make_pair(empty_key_, T()));
    h1 = bucket_count;
    h2 = get_left_prime_number(bucket_count);
  }

  HashMap(const HashMap &other, size_type bucket_count)
      : HashMap(bucket_count, other.empty_key_, other.get_allocator()) {
    for (auto it = other.begin(); it != other.end(); ++it) {
      insert(*it);
    }
  }

  HashMap(const std::vector<Key, T> &x, key_type empty_key,
          const allocator_type &alloc = allocator_type()) :
      empty_key_(empty_key), buckets_(alloc), capacity(x.capacity()) {

    for (const auto &i : x) {
      insert(i);
    }

  }

  allocator_type get_allocator() const noexcept {
    return buckets_.get_allocator();
  }

  iterator begin() noexcept {
    return iterator(this);
  }
  const_iterator cbegin() const noexcept {
    return const_iterator(this);
  }

  iterator end() {
    return iterator(this, buckets_.size());
  }

  const_iterator cend() {
    return const_iterator(this, buckets_.size());
  }

  bool empty() const noexcept {
    return size_ == 0;
  };

  void clear() noexcept {
    for (auto &x : buckets_) {
      if (x.first != empty_key_) {
        x.first = empty_key_;
      }
    }
    size_ = 0;
  }

  HashMap<Key, T> &operator=(const HashMap<Key, T> &other) {
    empty_key_ = other.empty_key_;
    size_ = other.size_;
    capacity = other.capacity;
    h1 = other.h1;
    h2 = other.h2;
    buckets_ = other.buckets_;
    return *this;
  }

  iterator insert(std::pair<Key, mapper_type> &x){

  }
 private:
  size_t hash_to_idx(const key_type &key, int i = 0) const noexcept(noexcept(hasher()(key))) {
    size_t hash = hasher()(key);
    return (hash % h1 + i * hash % h2) % size_;
  }

  size_t get_left_prime_number(size_t x) const noexcept {
    for (size_t i = x - 1; i > 1; i--) {
      bool prime_flag = true;

      for (size_t j = 2; j < ceil(sqrt(i)) && prime_flag; j++) {
        if (i % j == 0) {
          prime_flag = false;
        }
      }

      if (prime_flag) {
        return i;
      }
    }
    return 1;
  }
  void rehash() {
    HashMap<Key, T> temp(buckets_.size() * 2, empty_key_);
    for (auto &i: buckets_) {
      if (i.first != empty_key_) {
        temp.insert(i);
      }
    }
    *this = temp;

  }
  iterator emplace_impl(const key_type &key, const T &value) {
    if (load_factor > max_load_factor){
      rehash();
    }
    size_t i = 0;
    for (size_t idx = hash_to_idx(key);; idx = hash_to_idx(key, ++i)) {
      if (key_equal()(buckets_[idx].first, empty_key_) ||
          key_equal()(buckets_[idx].first, key)) {
        buckets_[idx].second = mapper_type(value);
        buckets_[idx].first = key;
        size_++;
        return iterator(this, idx);
      }
      if (i > 0 && idx == hash_to_idx()) {
        rehash();
      }
    }
  }
 private:
  key_type empty_key_;
  buckets buckets_;
  size_t size_ = 0;
  size_t h1;
  size_t h2;
  const double max_load_factor = 0.75;
  double load_factor = 0;
  size_t capacity = 0;
};

#endif //HASHMAP__HASHMAP_H_
