#pragma once

#include <cmath>
#include <functional>
#include <memory>
#include <utility>
#include <vector>
#include <type_traits>
#include <algorithm>
#include <exception>
#include <climits>

namespace fefu
{
enum NodeState {
  EMPTY,
  CONTAINS,
  DELETED,
  END
};

template <typename ValueType>
struct IterNode {
  IterNode() : ptr(nullptr), state(EMPTY) {}
  ValueType* ptr;
  NodeState state;
};

template<typename T>
class allocator {
 public:
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;
  using pointer = T*;
  using const_pointer = const T*;
  using reference = typename std::add_lvalue_reference<T>::type;
  using const_reference = typename std::add_lvalue_reference<const T>::type;
  using value_type = T;

  allocator() noexcept {}

  allocator(const allocator& other) noexcept : debug(other.debug) {}

  template <class U>
  allocator(const allocator<U>& other) noexcept : debug(other.debug) {}

  ~allocator() {}

  pointer allocate(size_type n) {
    pointer ptr = static_cast<pointer>(::operator new(n * sizeof(value_type)));
    return ptr;
  }

  void deallocate(pointer p, size_type n) noexcept {
    ::operator delete(static_cast<void*>(p), n * sizeof(value_type));
  }

  int debug = 0;
};


template<typename ValueType>
class hash_map_iterator {
 public:
  using iterator_category = std::forward_iterator_tag;
  using value_type = ValueType;
  using difference_type = std::ptrdiff_t;
  using reference = ValueType&;
  using pointer = ValueType*;

  hash_map_iterator() noexcept : node(nullptr) {}
  hash_map_iterator(const hash_map_iterator& other) noexcept : node(other.node) {}

  reference operator*() const {
    if (node == nullptr || node->ptr == nullptr)
      throw std::out_of_range("Iterator is out of range");
    return *(node->ptr);
  }
  pointer operator->() const {
    return node->ptr;
  }


  // prefix ++
  hash_map_iterator& operator++() {
    if (node->ptr == nullptr)
      throw std::out_of_range("Iterator is out of range");
    ++node;
    while (node->ptr != nullptr && node->state != CONTAINS)
      ++node;
    return *this;
  }
  // postfix ++
  hash_map_iterator operator++(int) {
    hash_map_iterator tmp(*this);
    operator++();
    return tmp;
  }

  friend bool operator==(const hash_map_iterator<ValueType>& lhs, const hash_map_iterator<ValueType>& rhs) {
    return (lhs.node == rhs.node);
  }
  friend bool operator!=(const hash_map_iterator<ValueType>& lhs, const hash_map_iterator<ValueType>& rhs) {
    return !(lhs == rhs);
  }

  template<typename A, typename B, typename C, typename D, typename E>
  friend class hash_map;

  template<typename R>
  friend class hash_map_const_iterator;

 private:
  hash_map_iterator(IterNode<ValueType>* iterNode) : node(iterNode) {
    while (node->state != CONTAINS && node->ptr != nullptr)
      ++node;
  }

  IterNode<ValueType>* node;
};

template<typename ValueType>
class hash_map_const_iterator {
 public:
  using iterator_category = std::forward_iterator_tag;
  using value_type = ValueType;
  using difference_type = std::ptrdiff_t;
  using reference = const ValueType&;
  using pointer = const ValueType*;

  hash_map_const_iterator() noexcept : node(nullptr) {}
  hash_map_const_iterator(const hash_map_const_iterator& other) noexcept : node(other.node) {}
  hash_map_const_iterator(const hash_map_iterator<ValueType>& other) noexcept : node(other.node) {}

  reference operator*() const {
    if (node == nullptr || node->ptr == nullptr)
      throw std::out_of_range("Iterator is out of range");
    return *(node->ptr);
  }
  pointer operator->() const {
    return node->ptr;
  }

  // prefix ++
  hash_map_const_iterator& operator++() {
    if (node->ptr == nullptr)
      throw std::out_of_range("Iterator is out of range");
    ++node;
    while (node->ptr != nullptr && node->state != CONTAINS)
      ++node;
    return *this;
  }
  // postfix ++
  hash_map_const_iterator operator++(int) {
    hash_map_const_iterator tmp(*this);
    operator++();
    return tmp;
  }

  friend bool operator==(const hash_map_const_iterator<ValueType>& lhs,
                         const hash_map_const_iterator<ValueType>& rhs) {
    return (lhs.node == rhs.node);
  }
  friend bool operator!=(const hash_map_const_iterator<ValueType>& lhs,
                         const hash_map_const_iterator<ValueType>& rhs) {
    return !(lhs == rhs);
  }

  template<typename A, typename B, typename C, typename D, typename E>
  friend class hash_map;

 private:
  hash_map_const_iterator(IterNode<ValueType>* iterNode) : node(iterNode) {
    while (node->state != CONTAINS && node->ptr != nullptr)
      node++;
  }

  IterNode<ValueType>* node;
};

template<typename K, typename T,
    typename Hash = std::hash<K>,
    typename Pred = std::equal_to<K>,
    typename Alloc = allocator<std::pair<const K, T>>>
class hash_map
{
 public:
  using key_type = K;
  using mapped_type = T;
  using hasher = Hash;
  using key_equal = Pred;
  using allocator_type = Alloc;
  using value_type = std::pair<const key_type, mapped_type>;
  using reference = value_type&;
  using const_reference = const value_type&;
  using iterator = hash_map_iterator<value_type>;
  using const_iterator = hash_map_const_iterator<value_type>;
  using size_type = std::size_t;

  /// Default constructor.
  hash_map() : mData(nullptr), mNodes(1) {}

  /**
   *  @brief  Default constructor creates no elements.
   *  @param n  Minimal initial number of buckets.
   */
  explicit hash_map(size_type n) : mData(mAlloc.allocate(getPowerOfTwo(n))) {
    n = getPowerOfTwo(n);
    mNodes = std::vector<IterNode<value_type>>(n + 1);
    for (size_type i = 0; i < n; i++)
      mNodes[i].ptr = mData +i;
  }

  /**
   *  @brief  Builds an %hash_map from a range.
   *  @param  first  An input iterator.
   *  @param  last  An input iterator.
   *  @param  n  Minimal initial number of buckets.
   *
   *  Create an %hash_map consisting of copies of the elements from
   *  [first,last).  This is linear in N (where N is
   *  distance(first,last)).
   */
  template<typename InputIterator>
  hash_map(InputIterator first, InputIterator last,
           size_type n = 0) : hash_map(n) {
    insert(first, last);
  }

  /// Copy constructor.
  hash_map(const hash_map& other) : mAlloc(other.mAlloc), mHash(other.mHash), mKeyEqual(other.mKeyEqual),
                                    mCount(other.mCount), mDeleted(other.mDeleted), maxLoadFactor(other.maxLoadFactor),
                                    mNodes(other.mNodes) {
    mData = mAlloc.allocate(other.mNodes.size());

    for (size_t i = 0; i < mNodes.size() - 1; i++) {
      if (other.mNodes[i].state == CONTAINS) {
        new(mData + i) value_type(other.mData[i]);
      }
    }
  }

  /// Move constructor.
  hash_map(hash_map&& rvalue) : mAlloc(std::move(rvalue.mAlloc)), mKeyEqual(std::move(rvalue.mKeyEqual)), mHash(std::move(rvalue.mHash)),
                                mCount(std::move(rvalue.mCount)), mDeleted(std::move(rvalue.mDeleted)), maxLoadFactor(std::move(rvalue.maxLoadFactor)),
                                mNodes(std::move(rvalue.mNodes)), mData(std::move(rvalue.mData)) {
    rvalue.mData = nullptr;
  }

  /**
   *  @brief Creates an %hash_map with no elements.
   *  @param a An allocator object.
   */
  explicit hash_map(const allocator_type& a) : mAlloc(a), mNodes(1) {
    mData = nullptr;
  }

  /**
  * @brief Copy constructor with allocator argument.
  * @param  uset  Input %hash_map to copy.
  * @param  a  An allocator object.
  */
  hash_map(const hash_map& umap,
           const allocator_type& a) : mAlloc(a), mCount(umap.mCount), mDeleted(umap.mDeleted), mNodes(umap.mNodes.size()) {

    mData = mAlloc.allocate(umap.mNodes.size());

    for (size_t i = 0; i < mNodes.size() - 1; i++) {
      if (umap.mNodes[i].state == CONTAINS) {
        mNodes[i].state = CONTAINS;
        new(mData + i) value_type(umap.mData[i]);
      }
    }
  }

  /**
  *  @brief  Move constructor with allocator argument.
  *  @param  uset Input %hash_map to move.
  *  @param  a    An allocator object.
  */
  hash_map(hash_map&& umap,
           const allocator_type& a) : mAlloc(a), mHash(std::move(umap.mHash)), mKeyEqual(std::move(umap.mKeyEqual)),
                                      mCount(std::move(umap.mCount)), mDeleted(std::move(umap.mDeleted)), maxLoadFactor(std::move(umap.maxLoadFactor)),
                                      mNodes(std::move(umap.mNodes)) {
    mData = umap.mData;
    umap.mData = nullptr;


  }

  /**
   *  @brief  Builds an %hash_map from an initializer_list.
   *  @param  l  An initializer_list.
   *  @param n  Minimal initial number of buckets.
   *
   *  Create an %hash_map consisting of copies of the elements in the
   *  list. This is linear in N (where N is @a l.size()).
   */
  hash_map(std::initializer_list<value_type> l,
           size_type n = 0) : hash_map(l.begin(), l.end(), n) {}

  ~hash_map() {
    if (mNodes.size() > 0) {
      for (size_type i = 0; i < mNodes.size() - 1; i++) {
        if (mNodes[i].state == CONTAINS)
          mData[i].~value_type();
      }
    }
    mAlloc.deallocate(mData, mNodes.size());
  }

  /// Copy assignment operator.
  hash_map& operator=(const hash_map&src) {
    hash_map(src).swap(*this);
    return *this;
  }

  /// Move assignment operator.
  hash_map& operator=(hash_map&&src) {
    mAlloc.deallocate(mData, mNodes.size());
    mData = nullptr;
    swap(src);
    return *this;
  }

  /**
   *  @brief  %hash_map list assignment operator.
   *  @param  l  An initializer_list.
   *
   *  This function fills an %hash_map with copies of the elements in
   *  the initializer list @a l.
   *
   *  Note that the assignment completely changes the %hash_map and
   *  that the resulting %hash_map's size is the same as the number
   *  of elements assigned.
   */
  hash_map& operator=(std::initializer_list<value_type> l) {
    hash_map(l).swap(*this);
    return *this;
  }

  ///  Returns the allocator object used by the %hash_map.
  allocator_type get_allocator() const noexcept {
    return mAlloc;
  }

  // size and capacity:

  ///  Returns true if the %hash_map is empty.
  bool empty() const noexcept {
    return (mCount == 0);
  }

  ///  Returns the size of the %hash_map.
  size_type size() const noexcept {
    return mCount;
  }

  ///  Returns the maximum size of the %hash_map.
  size_type max_size() const noexcept {
    return std::numeric_limits<size_type>::max();
  }

  // iterators.

  iterator begin() noexcept {
    //return iterator((decltype(iterator::node))&mNodes[0]);
    return iterator(mNodes.data());
  }

  //@{
  /**
   *  Returns a read-only (constant) iterator that points to the first
   *  element in the %hash_map.
   */
  const_iterator begin() const noexcept {
    return const_iterator((decltype(iterator::node))&mNodes[0]);
  }

  const_iterator cbegin() const noexcept {
    return const_iterator((decltype(iterator::node))&mNodes[0]);
  }

  /**
   *  Returns a read/write iterator that points one past the last element in
   *  the %hash_map.
   */
  iterator end() noexcept {
    return iterator((decltype(iterator::node))&mNodes[mNodes.size() - 1]);
  }

  //@{
  /**
   *  Returns a read-only (constant) iterator that points one past the last
   *  element in the %hash_map.
   */
  const_iterator end() const noexcept {
    return const_iterator((decltype(const_iterator::node))&mNodes[mNodes.size() - 1]);
  }

  const_iterator cend() const noexcept {
    return end();
  }
  //@}

  // modifiers.

  /**
   *  @brief Attempts to build and insert a std::pair into the
   *  %hash_map.
   *
   *  @param args  Arguments used to generate a new pair instance (see
   *	        std::piecewise_contruct for passing arguments to each
  *	        part of the pair constructor).
  *
  *  @return  A pair, of which the first element is an iterator that points
  *           to the possibly inserted pair, and the second is a bool that
  *           is true if the pair was actually inserted.
  *
  *  This function attempts to build and insert a (key, value) %pair into
  *  the %hash_map.
  *  An %hash_map relies on unique keys and thus a %pair is only
  *  inserted if its first element (the key) is not already present in the
  *  %hash_map.
  *
  *  Insertion requires amortized constant time.
  */
  template<typename... _Args>
  std::pair<iterator, bool> emplace(_Args&&... args) {
    return insert(value_type(std::forward<_Args&&>(args)...));
  }

  /**
   *  @brief Attempts to build and insert a std::pair into the
   *  %hash_map.
   *
   *  @param k    Key to use for finding a possibly existing pair in
   *                the hash_map.
   *  @param args  Arguments used to generate the .second for a
   *                new pair instance.
   *
   *  @return  A pair, of which the first element is an iterator that points
   *           to the possibly inserted pair, and the second is a bool that
   *           is true if the pair was actually inserted.
   *
   *  This function attempts to build and insert a (key, value) %pair into
   *  the %hash_map.
   *  An %hash_map relies on unique keys and thus a %pair is only
   *  inserted if its first element (the key) is not already present in the
   *  %hash_map.
   *  If a %pair is not inserted, this function has no effect.
   *
   *  Insertion requires amortized constant time.
   */
  template <typename... _Args>
  std::pair<iterator, bool> try_emplace(const key_type& k, _Args&&... args) {
    return tryEmplace(k, std::forward<_Args>(args)...);
  }

  // move-capable overload
  template <typename... _Args>
  std::pair<iterator, bool> try_emplace(key_type&& k, _Args&&... args) {
    return tryEmplace(std::forward<key_type>(k), std::forward<_Args>(args)...);
  }

  //@{
  /**
   *  @brief Attempts to insert a std::pair into the %hash_map.
   *  @param x Pair to be inserted (see std::make_pair for easy
   *	     creation of pairs).
  *
  *  @return  A pair, of which the first element is an iterator that
  *           points to the possibly inserted pair, and the second is
  *           a bool that is true if the pair was actually inserted.
  *
  *  This function attempts to insert a (key, value) %pair into the
  *  %hash_map. An %hash_map relies on unique keys and thus a
  *  %pair is only inserted if its first element (the key) is not already
  *  present in the %hash_map.
  *
  *  Insertion requires amortized constant time.
  */
  std::pair<iterator, bool> insert(const value_type& x) {
    return Insert(x);
  }

  std::pair<iterator, bool> insert(value_type&& x) {
    return Insert(std::move(x));
  }

  //@}

  /**
   *  @brief A template function that attempts to insert a range of
   *  elements.
   *  @param  first  Iterator pointing to the start of the range to be
   *                   inserted.
   *  @param  last  Iterator pointing to the end of the range.
   *
   *  Complexity similar to that of the range constructor.
   */
  template<typename _InputIterator>
  void insert(_InputIterator first, _InputIterator last) {
    for (auto it = first; it != last; it++) {
      insert(*it);
    }
  }

  /**
   *  @brief Attempts to insert a list of elements into the %hash_map.
   *  @param  l  A std::initializer_list<value_type> of elements
   *               to be inserted.
   *
   *  Complexity similar to that of the range constructor.
   */
  void insert(std::initializer_list<value_type> l) {
    insert(l.begin(), l.end());
  }


  /**
   *  @brief Attempts to insert a std::pair into the %hash_map.
   *  @param k    Key to use for finding a possibly existing pair in
   *                the map.
   *  @param obj  Argument used to generate the .second for a pair
   *                instance.
   *
   *  @return  A pair, of which the first element is an iterator that
   *           points to the possibly inserted pair, and the second is
   *           a bool that is true if the pair was actually inserted.
   *
   *  This function attempts to insert a (key, value) %pair into the
   *  %hash_map. An %hash_map relies on unique keys and thus a
   *  %pair is only inserted if its first element (the key) is not already
   *  present in the %hash_map.
   *  If the %pair was already in the %hash_map, the .second of
   *  the %pair is assigned from obj.
   *
   *  Insertion requires amortized constant time.
   */
  template <typename _Obj>
  std::pair<iterator, bool> insert_or_assign(const key_type& k, _Obj&& obj) {
    return insertAssign(k, std::move(obj));
  }

  // move-capable overload
  template <typename _Obj>
  std::pair<iterator, bool> insert_or_assign(key_type&& k, _Obj&& obj) {
    return insertAssign(std::move(k), std::move(obj));
  }

  //@{
  /**
   *  @brief Erases an element from an %hash_map.
   *  @param  position  An iterator pointing to the element to be erased.
   *  @return An iterator pointing to the element immediately following
   *          @a position prior to the element being erased. If no such
   *          element exists, end() is returned.
   *
   *  This function erases an element, pointed to by the given iterator,
   *  from an %hash_map.
   *  Note that this function only erases the element, and that if the
   *  element is itself a pointer, the pointed-to memory is not touched in
   *  any way.  Managing the pointer is the user's responsibility.
   */
  iterator erase(const_iterator position) {
    if (position == this->cend())
      throw std::out_of_range("Cant erase end iterator");
    (*position.node).ptr->~value_type();
    (*position.node).state = DELETED;
    position++;
    mDeleted++;
    mCount--;

    return iterator(position.node);
  }

  // LWG 2059.
  iterator erase(iterator position) {
    return erase(const_iterator(position));
  }
  //@}

  /**
   *  @brief Erases elements according to the provided key.
   *  @param  x  Key of element to be erased.
   *  @return  The number of elements erased.
   *
   *  This function erases all the elements located by the given key from
   *  an %hash_map. For an %hash_map the result of this function
   *  can only be 0 (not present) or 1 (present).
   *  Note that this function only erases the element, and that if the
   *  element is itself a pointer, the pointed-to memory is not touched in
   *  any way.  Managing the pointer is the user's responsibility.
   */
  size_type erase(const key_type& x) {
    auto it = find(x);
    if (it == cend())
      return 0;
    erase(it);
    return 1;
  }

  /**
   *  @brief Erases a [first,last) range of elements from an
   *  %hash_map.
   *  @param  first  Iterator pointing to the start of the range to be
   *                  erased.
   *  @param last  Iterator pointing to the end of the range to
   *                be erased.
   *  @return The iterator @a last.
   *
   *  This function erases a sequence of elements from an %hash_map.
   *  Note that this function only erases the elements, and that if
   *  the element is itself a pointer, the pointed-to memory is not touched
   *  in any way.  Managing the pointer is the user's responsibility.
   */
  iterator erase(const_iterator first, const_iterator last) {
    for (auto it = first; it != last; ++it) {
      erase(it);
    }
    return iterator(last.node);
  }

  template <typename Pred_>
  void erase_if(Pred_ pred) {
    for (auto it = begin(); it != end(); ++it) {
      if (pred(*it)) {
        erase(it);
      }
    }
  }

  /**
   *  Erases all elements in an %hash_map.
   *  Note that this function only erases the elements, and that if the
   *  elements themselves are pointers, the pointed-to memory is not touched
   *  in any way.  Managing the pointer is the user's responsibility.
   */
  void clear() noexcept {
    erase(this->begin(), this->end());
  }

  /**
   *  @brief  Swaps data with another %hash_map.
   *  @param  x  An %hash_map of the same element and allocator
   *  types.
   *
   *  This exchanges the elements between two %hash_map in constant
   *  time.
   *  Note that the global std::swap() function is specialized such that
   *  std::swap(m1,m2) will feed to this function.
   */
  void swap(hash_map& x) {
    using std::swap;

    if constexpr (std::allocator_traits<Alloc>::propagate_on_container_swap::value)
      swap(this->mAlloc, x.mAalloc);

    swap(this->mData, x.mData);
    swap(this->mNodes, x.mNodes);
    swap(this->mCount, x.mCount);
    swap(this->mDeleted, x.mDeleted);
    swap(this->maxLoadFactor, x.maxLoadFactor);
    swap(this->mKeyEqual, x.mKeyEqual);
    swap(this->mHash, x.mHash);
  }

  template<typename _H2, typename _P2>
  void merge(hash_map<K, T, _H2, _P2, Alloc>& source) {
    Merge(source);
  }

  template<typename _H2, typename _P2>
  void merge(hash_map<K, T, _H2, _P2, Alloc>&& source) {
    Merge(std::move(source));
  }

  // observers.

  ///  Returns the hash functor object with which the %hash_map was
  ///  constructed.
  Hash hash_function() const {
    return mHash;
  }

  ///  Returns the key comparison object with which the %hash_map was
  ///  constructed.
  Pred key_eq() const {
    return mKeyEqual;
  }



  //@{
  /**
   *  @brief Tries to locate an element in an %hash_map.
   *  @param  x  Key to be located.
   *  @return  Iterator pointing to sought-after element, or end() if not
   *           found.
   *
   *  This function takes a key and tries to locate the element with which
   *  the key matches.  If successful the function returns an iterator
   *  pointing to the sought after element.  If unsuccessful it returns the
   *  past-the-end ( @c end() ) iterator.
   */
  iterator find(const key_type& x) {
    if (mCount == 0)
      return end();
    size_type indx = search(x);
    if (mNodes[indx].state == CONTAINS) {
      return iterator((decltype(iterator::node))&mNodes[indx]);
    }
    return end();
  }

  const_iterator find(const key_type& x) const {
    size_type indx = search(x);
    if (mNodes[indx].state == CONTAINS) {
      return iterator((decltype(iterator::node))&mNodes[indx]);
    }
    return end();
  }
  //@}

  /**
   *  @brief  Finds the number of elements.
   *  @param  x  Key to count.
   *  @return  Number of elements with specified key.
   *
   *  This function only makes sense for %unordered_multimap; for
   *  %hash_map the result will either be 0 (not present) or 1
   *  (present).
   */
  size_type count(const key_type& x) const {
    return contains(x);
  }

  /**
   *  @brief  Finds whether an element with the given key exists.
   *  @param  x  Key of elements to be located.
   *  @return  True if there is any element with the specified key.
   */
  bool contains(const key_type& x) const {
    size_type indx = search(x);
    return (mNodes[indx].state == CONTAINS);
  }

  //@{
  /**
   *  @brief  Subscript ( @c [] ) access to %hash_map data.
   *  @param  k  The key for which data should be retrieved.
   *  @return  A reference to the data of the (key,data) %pair.
   *
   *  Allows for easy lookup with the subscript ( @c [] )operator.  Returns
   *  data associated with the key specified in subscript.  If the key does
   *  not exist, a pair with that key is created using default values, which
   *  is then returned.
   *
   *  Lookup requires constant time.
   */
  mapped_type& operator[](const key_type& k) {
    return Operator(k);
  }

  mapped_type& operator[](key_type&& k) {
    return Operator(std::move(k));
  }
  //@}

  //@{
  /**
   *  @brief  Access to %hash_map data.
   *  @param  k  The key for which data should be retrieved.
   *  @return  A reference to the data whose key is equal to @a k, if
   *           such a data is present in the %hash_map.
   *  @throw  std::out_of_range  If no such data is present.
   */
  mapped_type& at(const key_type& k) {
    size_type indx = search(k);
    if (mNodes[indx].state != CONTAINS) {
      throw std::out_of_range("This key is not presented in map");
    }

    return mData[indx].second;
  }

  const mapped_type& at(const key_type& k) const {
    size_type indx = search(k);
    if (mNodes[indx].state != CONTAINS) {
      throw std::out_of_range("This key is not presented in map");
    }

    return mData[indx].second;
  }
  //@}

  // bucket interface.

  /// Returns the number of buckets of the %hash_map.
  size_type bucket_count() const noexcept {
    return (mNodes.size() - 1);
  }

  /**
  * @brief  Returns the bucket index of a given element.
  * @param  _K  A key instance.
  * @return  The key bucket index.
  */
  size_type bucket(const key_type& _K) const {
    size_type indx = search(_K);
    if (mNodes[indx].state != CONTAINS) {
      throw std::out_of_range("This key is not presented in map");
    }
    return indx;
  }

  // hash policy.

  /// Returns the average number of elements per bucket.
  float load_factor() const noexcept {
    return (mCount * 1.0f + mDeleted) / (mNodes.size()-1);
  }

  /// Returns a positive number that the %hash_map tries to keep the
  /// load factor less than or equal to.
  float max_load_factor() const noexcept {
    return maxLoadFactor;
  }

  /**
   *  @brief  Change the %hash_map maximum load factor.
   *  @param  z The new maximum load factor.
   */
  void max_load_factor(float z) {
    if (z <= 0.0 || z >= 1)
      throw std::invalid_argument("Load factor must be positive and less than 1");
    maxLoadFactor = z;
    checkForRehash();
  }

  /**
   *  @brief  May rehash the %hash_map.
   *  @param  n The new number of buckets.
   *
   *  Rehash will occur only if the new number of buckets respect the
   *  %hash_map maximum load factor.
   */
  void rehash(size_type n) {

    n = getPowerOfTwo(n);
    mDeleted = 0;
    hash_map newHashMap(n);
    for (size_type i = 0; i < mNodes.size() - 1; i++) {
      if (mNodes[i].state == CONTAINS) {
        mCount++;
        newHashMap.emplace(std::move(mData[i].first), std::move(mData[i].second));
      }
    }
    swap(newHashMap);

  }

  /**
   *  @brief  Prepare the %hash_map for a specified number of
   *          elements.
   *  @param  n Number of elements required.
   *
   *  Same as rehash(ceil(n / max_load_factor())).
   */
  void reserve(size_type n) {
    rehash(std::ceil(n / maxLoadFactor));
  }

  bool operator==(const hash_map& other) const {
    if (this->size() != other.size())
      return false;
    return std::is_permutation(begin(),end(),other.begin(),other.end());
  }

 private:
  Alloc mAlloc;
  hasher mHash;
  key_equal mKeyEqual;
  size_type mCount = 0;
  size_type mDeleted = 0;
  float maxLoadFactor = 0.75f;
  std::vector<IterNode<value_type>> mNodes;
  value_type* mData;
  const size_t increase = 3;

  size_type search(const key_type& k, bool forFind = true) const {
    size_type indx = mHash(k) % bucket_count();
    auto tr = 1;
    while (!mKeyEqual(mData[indx].first, k) && (mNodes[indx].state == CONTAINS) || mNodes[indx].state == DELETED) {
      indx = (indx + tr*tr*tr) % bucket_count();
      tr++;
    }
    return indx;
  }

  bool checkForRehash() {
    if (mNodes.size() < 2) {
      rehash(2);
      return true;
    }
    if (load_factor() >= maxLoadFactor) {
      rehash(mNodes.size() * increase);
      return true;
    }
    return false;
  }

  inline size_type getPowerOfTwo(size_type n) {
    n--;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    n++;
    return n;
  }

  template <typename _T>
  std::pair<iterator, bool> Insert(_T&& el) {
    checkForRehash();
    size_type indx = search(el.first);
    if (mNodes[indx].state == EMPTY || mNodes[indx].state == DELETED) {
      new (mData + indx) value_type(std::forward<_T>(el));
      mNodes[indx].state = CONTAINS;
      mCount++;
      return make_pair(iterator(&mNodes[indx]), true);
    }

    return make_pair(iterator(&mNodes[indx]), false);
  }

  template <typename _T>
  mapped_type& Operator(_T&& k) {
    checkForRehash();
    size_type indx = search(k);
    if (mNodes[indx].state == EMPTY) {
      indx = search(k, false);
      new(mData + indx) value_type(std::forward<_T>(k), mapped_type());
      mNodes[indx].state = CONTAINS;
      mCount++;
    }

    return mData[indx].second;
  }

  template<typename _T>
  void Merge(_T&& source) {
    for (auto it = source.begin(); it != source.end(); it++) {
      if (!contains(it->first)) {
        insert(std::forward<value_type>(*it));
        source.erase(it);
      }
    }
  }

  template<typename... _Args, typename _T>
  std::pair<iterator, bool> tryEmplace(_T&& k, _Args&&... args) {
    checkForRehash();
    size_type indx = search(k);
    if (mNodes[indx].state == EMPTY) {
      new (mData + indx) value_type(std::piecewise_construct,
                                    std::forward_as_tuple(std::forward<_T>(k)),
                                    std::forward_as_tuple(std::forward<_Args>(args)...));
      mNodes[indx].state = CONTAINS;
      mCount++;
      return make_pair(iterator(&mNodes[indx]), true);
    }

    return make_pair(iterator(&mNodes[indx]), false);
  }

  template <typename _Obj, typename _T>
  std::pair<iterator, bool> insertAssign(_T&& k, _Obj&& obj) {
    checkForRehash();

    size_type indx = search(k);
    if (mNodes[indx].state == EMPTY) {
      indx = search(k, false);
      new(mData + indx) value_type(std::forward<_T>(k), mapped_type(std::move(obj)));
      mNodes[indx].state = CONTAINS;
      mCount++;
      return std::make_pair(iterator(&mNodes[indx]), true);
    }
    mData[indx].second = std::move(obj);

    return std::make_pair(iterator(&mNodes[indx]), false);
  }

};

} // namespace fefu