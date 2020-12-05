//
// Created by Oleg Antokhiv on 26.11.2020.
//

#ifndef HASHMAP__HASH_MAP_H_
#define HASHMAP__HASH_MAP_H_

#include <functional>
#include <vector>
#include <cmath>
#include <limits>
template<class Key, class T, class Hash = std::hash<Key>,
    class KeyEqual = std::equal_to<void>,
    class Allocator = std::allocator<std::pair<Key, T>>>
class hash_map {
 public:
  using key_type = Key;
  using mapped_type = T;
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
    explicit hm_iterator(hash_map *hm) : hm_(hm) { advance_past_empty(); }
    explicit hm_iterator(hash_map *hm, size_type idx) : hm_(hm), idx_(idx) {}
    template<typename OtherHashMap, typename OtherIterVal>
    hm_iterator(const hm_iterator<OtherIterVal> &other)
        : hm_(other.hm_), idx_(other.idx_) {}

    void advance_past_empty() {
      while (idx_ < hm_->buckets_.size() &&
          hm_->check[idx_]) {
        ++idx_;
      }
    }

    hash_map *hm_ = nullptr;
    typename hash_map::size_type idx_ = 0;
    friend hash_map;
  };
  using iterator = hm_iterator<value_type>;
  using const_iterator = hm_iterator<const value_type>;

 public:
  hash_map() = default;
  hash_map(size_type bucket_count,
           const allocator_type &alloc = allocator_type())
      : buckets_(alloc) {
    buckets_.reserve(bucket_count);
    capacity = bucket_count;
    buckets_.resize(bucket_count, std::make_pair(key_type(), T()));
    h1 = bucket_count;
    h2 = bucket_count;
  }

  hash_map(const allocator_type &alloc)
      : buckets_(alloc) {
    size_t bucket_count = 0;
    buckets_.reserve(bucket_count);
    capacity = bucket_count;
    buckets_.resize(bucket_count, std::make_pair(key_type(), T()));
    h1 = bucket_count;
    h2 = bucket_count;
  }

  hash_map(const hash_map &other, size_type bucket_count)
      : hash_map(bucket_count, other.get_allocator()) {
    for (auto it = other.begin(); it != other.end(); ++it) {
      insert(*it);
    }
  }

  hash_map(const hash_map &other, const allocator_type &alloc)
      : hash_map(other.capacity, alloc) {
    for (auto it = other.begin(); it != other.end(); ++it) {
      insert(*it);
    }
  }


  hash_map(const std::vector<std::pair<Key, T>> &x,
           const allocator_type &alloc = allocator_type()) :
      buckets_(alloc), capacity(x.capacity(), load_factor_(x.size() / x.capacity())) {

    for (const auto &i : x) {
      insert(i);
    }

  }

  hash_map(std::initializer_list<std::pair<const key_type, mapped_type>> &list,
           const size_t &bucket_count = 0) : hash_map(bucket_count) {
    for (const auto &i : list) {
      insert(i);
    }
  }
  template<class Iterator>
  hash_map(Iterator first, Iterator last, size_t bucket_count = 0)
      : hash_map(bucket_count) {
    auto i = first;
    while (first != last) {
      insert(*i);
      i++;
    }
  }

  allocator_type get_allocator() const noexcept {
    return buckets_.get_allocator();
  }

  iterator begin() const {
    return iterator(this);
  }
  const_iterator cbegin() const noexcept {
    return const_iterator(this);
  }

  iterator end() const {
    return iterator(this, buckets_.size());
  }

  const_iterator cend() {
    return const_iterator(this, buckets_.size());
  }

  bool empty() const noexcept {
    return size_ == 0;
  };

  hash_map<Key, T> &operator=(const hash_map<Key, T> &other) {
    size_ = other.size_;
    capacity = other.capacity;
    load_factor_ = other.load_factor_;
    h1 = other.h1;
    h2 = other.h2;
    buckets_ = other.buckets_;
    check = other.check;
    return *this;
  }

  iterator insert(const std::pair<Key, mapped_type> &x) {

    return emplace_impl(x.first, x.second);

  }

  size_t erase(const key_type &key) {
    auto it = find_impl(key);
    if (it != end()) {
      erase_impl(it);
      return 1;
    }
    return -1;
  }

  iterator at(const key_type &key) {
    auto it = find_impl(key);
    return at_impl(it);
  }

  mapped_type operator[](const key_type &key) {
    return emplace_impl(key, mapped_type())->second;
  }

  size_t size() {
    return size_;
  }

  size_t bucket_count() {
    return capacity;
  }

  double_t load_factor() {
    return load_factor_;
  }

  void max_load_factor(double_t x) {
    max_load_factor_ = x;
  }
 private:
  size_t hash_to_idx(const key_type &key, int i = 0) const noexcept(noexcept(hasher()(key))) {
    size_t hash = hasher()(key);
    return (hash % h1 + i * (1 + (hash % (capacity - 1)))) % capacity;
  }

  void rehash() {
    std::vector<value_type, allocator_type> temp(capacity *
    2);
    std::vector<bool> temp_check(capacity *
    2, true);
    capacity = temp.capacity();
    load_factor_ = size_ / capacity;
    h1 = capacity;
    h2 = capacity;
    buckets_.swap(temp);
    check.swap(temp_check);
    for (size_t i = 0; i < temp.size(); i++) {
      if (!temp_check[i]) {
        insert(temp[i]);
      }
    }
  }

  iterator emplace_impl(const key_type &key, const T &value) {
    if (load_factor_ > max_load_factor_ || capacity == 0) {
      if (capacity == 0) {
        capacity += 1;
      }
      rehash();
    }
    size_t i = 0;
    for (size_t idx = hash_to_idx(key);; idx = hash_to_idx(key, ++i)) {
      if (check[idx] ||
          key_equal()(buckets_[idx].first, key)) {
        buckets_[idx].second = mapped_type(value);
        buckets_[idx].first = key;
        check[idx] = false;
        size_++;
        load_factor_ = double_t(size_) / double_t(capacity);
        return iterator(this, idx);
      }
      if (i > 0 && idx == hash_to_idx(key)) {
        rehash();
        i = 0;
      }
    }
  }

  iterator find_impl(const key_type &key) {
    size_t i = 0;
    for (size_t ind = hash_to_idx(key);; ind = hash_to_idx(key, ++i)) {
      if (key_equal()(buckets_[ind].first, key)) {
        return iterator(this, ind);
      }
      if (i > 0 && ind == hash_to_idx(key)) {
        return end();
      }
    }
  }

  void erase_impl(const iterator &it) {
    check[it.idx_] = true;
    buckets_[it.idx] = {key_type(), mapped_type()};
    size_--;
    load_factor_ = size_ / capacity;
  }

  iterator at_impl(const iterator &it) {
    if (it != end()) {
      return it;
    }
    throw std::out_of_range("hash_map::at");
  }
 private:
  buckets buckets_;
  std::vector<bool> check;
  size_t size_ = 0;
  size_t h1;
  size_t h2;
  double max_load_factor_ = 0.75;
  double load_factor_ = 0.0;
  size_t capacity = 0;
};

#endif //HASHMAP__HASH_MAP_H_
