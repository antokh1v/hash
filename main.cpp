#include <vector>
#include <string>
#include <list>
#include <forward_list>
#include <set>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include <stack>
#include <queue>
#include <deque>

#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <tuple>
#include <type_traits>
#include <functional>
#include <utility>
#include <atomic>
#include <thread>
#include <future>
#include <chrono>
#include <iterator>
#include <memory>
#include <climits>

#define CATCH_CONFIG_MAIN
#include "catch.hpp"



// *************************
// *************************
// *************************
// *************************
// *************************

#include "HashMap.h"
using namespace fefu;

// *************************
// *************************
// *************************
// *************************
// *************************


using namespace std;


template<typename T>
class custom_allocator
{
 public:
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;
  using pointer = T*;
  using const_pointer = const T*;
  using reference = typename std::add_lvalue_reference<T>::type;
  using const_reference = typename std::add_lvalue_reference<const T>::type;
  using value_type = T;

 public:
  int x;

  custom_allocator() : x(rand()) { }

  custom_allocator(int a) : x(a) { }

  custom_allocator(const custom_allocator& other) noexcept : x(other.x) { }

  template <class U>
  custom_allocator(const custom_allocator<U>& other) noexcept : x(other.x) { }

  ~custom_allocator() = default;

  pointer allocate(size_type n)
  {
    return static_cast<pointer>(::operator new(n * sizeof(value_type)));
  }

  void deallocate(pointer p, size_type) noexcept
  {
    ::operator delete(p);
  }

  bool operator==(const custom_allocator&) const
  {
    return true;
  }

  bool operator!=(const custom_allocator&) const
  {
    return false;
  }
};

template<auto N, auto... Ns>
struct is_npack_contain : std::bool_constant<((N == Ns) || ...)> { };

template<auto N, auto... Ns>
constexpr bool is_npack_contain_v = is_npack_contain<N, Ns...>::value;

enum ff
{
  dc, cc, mc, ca, ma
};

template<ff... Forbidden>
struct declarator
{
  template<typename Dummy = void>
  constexpr void check_dc()
  {
    if constexpr (is_npack_contain_v<dc, Forbidden...>) static_assert (!std::is_same_v<Dummy, Dummy>, "Default constructor is deleted");
  }

  template<typename Dummy = void>
  constexpr void check_cc()
  {
    if constexpr (is_npack_contain_v<cc, Forbidden...>) static_assert (!std::is_same_v<Dummy, Dummy>, "Copy constructor is deleted");
  }

  template<typename Dummy = void>
  constexpr void check_mc()
  {
    if constexpr (is_npack_contain_v<mc, Forbidden...>) static_assert (!std::is_same_v<Dummy, Dummy>, "Move constructor is deleted");
  }

  template<typename Dummy = void>
  constexpr void check_ca()
  {
    if constexpr (is_npack_contain_v<ca, Forbidden...>) static_assert (!std::is_same_v<Dummy, Dummy>, "Copy assignment is deleted");
  }

  template<typename Dummy = void>
  constexpr void check_ma()
  {
    if constexpr (is_npack_contain_v<ma, Forbidden...>) static_assert (!std::is_same_v<Dummy, Dummy>, "Move assignment is deleted");
  }

  int x;
  declarator(int x) : x(x) { }
  declarator() : x(0) { check_dc(); }
  declarator(const declarator& other) : x(other.x) { check_cc(); }
  declarator(declarator&& other) : x(std::move(other.x)) { check_mc(); }
  declarator& operator=(const declarator&) { check_ca(); return *this; }
  declarator& operator=(declarator&&) { check_ma(); return *this; }
};

using atype = custom_allocator<std::pair<const int, int>>;
using hmint = hash_map<int, int>;
using hminta = hash_map<int, int, std::hash<int>, std::equal_to<int>, atype>;
const double eps = 10e-7;

TEST_CASE("_ctor()")
{
srand(time(nullptr));
hmint m;
REQUIRE(m.size() == 0);
}

TEST_CASE("_ctor(n)")
{
hmint m(2048);
REQUIRE(m.bucket_count() == 2048);
}

TEST_CASE("_ctor(a)")
{
atype a(4);
hminta m(a);
REQUIRE(m.get_allocator().x == 4);
}

TEST_CASE("_ctor(l, n)")
{
std::initializer_list<pair<const int, int>> l{{1, 2}, {3, 4}};
hmint m(l, 2048);
REQUIRE(m[1] == 2);
REQUIRE(m[3] == 4);
REQUIRE(m.bucket_count() == 2048);
}

TEST_CASE("_ctor(b, e, n")
{
std::initializer_list<pair<const int, int>> l{{1, 2}, {3, 4}};
hmint m(l.begin(), l.end(), 2048);
REQUIRE(m[1] == 2);
REQUIRE(m[3] == 4);
REQUIRE(m.bucket_count() == 2048);
}

TEST_CASE("_ctor(cref)")
{
std::initializer_list<pair<const int, int>> l{{1, 2}, {3, 4}};
hminta m1(l.begin(), l.end(), 2048);
hminta m2(m1);
REQUIRE(m1[1] == m2[1]);
REQUIRE(m1[3] == m2[3]);
REQUIRE(m1.size() == m2.size());
REQUIRE(m1.get_allocator().x == m2.get_allocator().x);
}

TEST_CASE("_ctor(cref, a)")
{
std::initializer_list<pair<const int, int>> l{{1, 2}, {3, 4}};
hminta m1(l.begin(), l.end(), 2048);
atype a;
hminta m2(m1, a);
REQUIRE(m1[1] == m2[1]);
REQUIRE(m1[3] == m2[3]);
REQUIRE(m1.size() == m2.size());
REQUIRE(m2.get_allocator().x == a.x);
}

TEST_CASE("_ctor(rref)")
{
std::initializer_list<pair<const int, int>> l{{1, 2}, {3, 4}};
hminta m1(l.begin(), l.end(), 2048);
hminta m2(std::move(m1));
REQUIRE(m2[1] == 2);
REQUIRE(m2[3] == 4);
REQUIRE(m1.get_allocator().x == m2.get_allocator().x);
}

TEST_CASE("_ctor(rref, a)")
{
std::initializer_list<pair<const int, int>> l{{1, 2}, {3, 4}};
hminta m1(l.begin(), l.end(), 2048);
atype a;
hminta m2(std::move(m1), a);
REQUIRE(m2[1] == 2);
REQUIRE(m2[3] == 4);
REQUIRE(m2.get_allocator().x == a.x);
}

TEST_CASE("_=(cref)")
{
std::initializer_list<pair<const int, int>> l{{1, 2}, {3, 4}};
hminta m1(l.begin(), l.end(), 2048);
hminta m2;
m2 = m1;
REQUIRE(m1[1] == m2[1]);
REQUIRE(m1[3] == m2[3]);
REQUIRE(m1.size() == m2.size());
REQUIRE(m1.get_allocator().x != m2.get_allocator().x);
}

TEST_CASE("_=(rref)")
{
std::initializer_list<pair<const int, int>> l{{1, 2}, {3, 4}};
hminta m1(l.begin(), l.end(), 2048);
hminta m2;
m2 = std::move(m1);
REQUIRE(m2[1] == 2);
REQUIRE(m2[3] == 4);
REQUIRE(m2.size() == 2);
REQUIRE(m1.get_allocator().x != m2.get_allocator().x);
}

TEST_CASE("_=(l)")
{
std::initializer_list<pair<const int, int>> l{{1, 2}, {3, 4}};
hmint m;
m = l;
REQUIRE(m[1] == 2);
REQUIRE(m[3] == 4);
REQUIRE(m.size() == 2);
}

TEST_CASE("_empty")
{
std::initializer_list<pair<const int, int>> l{{1, 2}, {3, 4}};
hmint m(l);
REQUIRE(!m.empty());
m = hmint();
REQUIRE(m.empty());
}

TEST_CASE("_begin_end")
{
std::initializer_list<pair<const int, int>> l{{1, 2}, {3, 4}};
hmint m(l);
unordered_map<int, int> sm(m.begin(), m.end());
REQUIRE(sm.size() == 2);
REQUIRE(sm[1] == 2);
REQUIRE(sm[3] == 4);
}

TEST_CASE("_contains_count")
{
std::initializer_list<pair<const int, int>> l{{1, 2}, {3, 4}};
hmint m(l);
REQUIRE(m.contains(1));
REQUIRE(m.count(3) == 1);
}

TEST_CASE("_emplace")
{
hmint m;
m.emplace(1, 2);
REQUIRE(m.contains(1));
}

TEST_CASE("_insert_copy_deleted")
{
hash_map<int, declarator<cc, ca>> m;
pair<const int, declarator<cc, ca>> p(piecewise_construct, std::tuple<int>(1), std::tuple<>());
m.insert(std::move(p));
}

TEST_CASE("_try_emplace")
{
hash_map<int, declarator<dc, cc, ca>> m;
m.try_emplace(1, 4);
REQUIRE(m.at(1).x == 4);
}

TEST_CASE("_swap")
{
hash_map<int, declarator<dc, cc, ca>> m1;
hash_map<int, declarator<dc, cc, ca>> m2;
m1.try_emplace(1, 4);
m2.try_emplace(5, 2);
m1.swap(m2);
REQUIRE(m1.at(5).x == 2);
REQUIRE(m2.at(1).x == 4);
}

TEST_CASE("_merge_ref")
{
hash_map<int, declarator<dc>> m1;
m1.try_emplace(1, 4);
hash_map<int, declarator<dc>> m2;
m2.try_emplace(5, 2);
m1.merge(m2);
REQUIRE(m1.at(5).x == 2);
REQUIRE(m2.size() == 0);
}

TEST_CASE("_merge_rref")
{
hash_map<int, declarator<dc, cc, ca>> m1;
m1.try_emplace(1, 4);
m1.try_emplace(1,6);
hash_map<int, declarator<dc, cc, ca>> m2;
m2.try_emplace(5, 2);
m1.merge(std::move(m2));
REQUIRE(m1.at(5).x == 2);
}

TEST_CASE("_swap_compile_time_check")
{
hash_map<int, declarator<dc, cc, ca, mc, ma>> m1;
hash_map<int, declarator<dc, cc, ca, mc, ma>> m2;
m1.swap(m2);
}

TEST_CASE("_insert_or_assign_compile_time_check")
{
hash_map<int, declarator<dc, cc, ca>> m;
m.insert_or_assign(1, 2);
}

TEST_CASE("_load_factor")
{
hmint m;
for (int i = 0; i < 100; ++i)
m.insert({i, 1});
for (int i = 0; i < 100; ++i)
m.erase(i);
REQUIRE(m.load_factor() > 0.0f);
}

template<typename TimePoint>
class time_meter
{
 public:
  using time_point = TimePoint;

 private:
  std::function<time_point()> m_fnow;
  std::function<double(time_point, time_point)> m_fsec;
  time_point m_begin;
  time_point m_stop;
  bool m_stopped;

 public:
  template<typename FNow, typename FSec>
  time_meter(FNow&& now, FSec&& sec) : m_fnow(std::forward<FNow>(now)), m_fsec(std::forward<FSec>(sec)), m_begin(m_fnow()), m_stopped(false) { }

  double seconds() const
  {
    if (m_stopped)
      return m_fsec(m_begin, m_stop);
    return m_fsec(m_begin, m_fnow());
  }

  void restart()
  {
    m_stopped = false;
    m_begin = m_fnow();
  }

  void stop()
  {
    if (m_stopped)
      return;
    m_stop = m_fnow();
    m_stopped = true;
  }

  void start()
  {
    if (!m_stopped)
      return;
    m_stopped = false;
    m_begin += m_fnow() - m_stop;
  }
};

auto create_tm()
{
  using tm_type = time_meter<std::chrono::high_resolution_clock::time_point>;

  static const auto get_sec = [](tm_type::time_point p1, tm_type::time_point p2)
  {
    return static_cast<double>((p2 - p1).count()) / std::chrono::high_resolution_clock::period::den;
  };

  return tm_type(std::chrono::high_resolution_clock::now, get_sec);
}

struct custom_hash
{
  vector<long long unsigned int> rnd;

  custom_hash()
  {
    for (int i = 0; i < 1000000; ++i)
      rnd.push_back(rand());
  }

  long long unsigned int operator()(int i) const
  {
    return rnd[i];
  }
};

TEST_CASE("_stress")
{
auto tm = create_tm();
hash_map<int, int, custom_hash> m(1000000);
m.max_load_factor(0.1f);
for (int i = 0; i < 1000000; ++i)
{
m.insert({i, i * 3});
}
for (int i = 100; i < 999999; ++i)
{
m.erase(i);
}
for (int i = 0; i < 1000000; ++i)
{
m.insert({i, i * 3});
}
for (int i = 0; i < 1000000; ++i)
{
m.insert({i, i * 3});
}
for (int i = 100; i < 999999; ++i)
{
m.erase(i);
}
for (int i = 100; i < 999999; ++i)
{
m.erase(i);
}
for (int i = 0; i < 100; ++i)
REQUIRE(m[i] == i * 3);
REQUIRE(m[999999] == 2999997);
REQUIRE(m.size() == 101);
for (int i = 0; i < 1000; ++i)
m.insert({rand() % 1000000, rand()});
cout << "**********\n" << "STRESS TIME = " << tm.seconds() << "\n**********" << endl;
}

TEST_CASE("Allocator") {
vector<int, fefu::allocator<int>> tmp;
for (size_t i = 0; i < 10; i++)
tmp.push_back(i);

REQUIRE(tmp.size() == 10);
for (size_t i = 0; i < 10; i++)
REQUIRE(tmp[i] == i);

for (size_t i = 10; i > 0; i--) {
tmp.pop_back();
REQUIRE(tmp.size() == i - 1);
}
}

TEST_CASE("operator[]") {
fefu::hash_map<int, string> hmap(10);
hmap[2] = "abacaba";
CHECK(hmap[2] == "abacaba");
int t = 5;
hmap[t] = "abasab";
CHECK(hmap[t] == "abasab");
hmap[5] = "abc";
CHECK(hmap[t] == "abc");
CHECK(hmap[0] == "");
}

TEST_CASE("at()") {
fefu::hash_map<int, string> hmap(10);
hmap[2] = "abacaba";
CHECK(hmap.at(2) == "abacaba");
hmap[3] = "ab";
CHECK(hmap.at(3) == "ab");
hmap[2] = "abc";
CHECK(hmap.at(2) == "abc");

const fefu::hash_map<int, string> hmap2(hmap);
CHECK(hmap2.at(3) == "ab");
}

TEST_CASE("bucket_count()") {
fefu::hash_map<int, string> hmap(10);
CHECK(hmap.bucket_count() == 16);
}

TEST_CASE("bucket()") {
fefu::hash_map<int, string> hmap(10);
hmap[4] = "abc";
size_t t = hash<int>{}(4) % 16;
CHECK(t == hmap.bucket(4));
}

TEST_CASE("rehash()") {
fefu::hash_map<int, string> hmap(10);
hmap[4] = "abc";
hmap[-1] = "a";
hmap[2] = "bc";
hmap.rehash(100);
REQUIRE(hmap.bucket_count() == 128);
CHECK(hmap[4] == "abc");
CHECK(hmap[-1] == "a");
CHECK(hmap[2] == "bc");
}

TEST_CASE("reserve()") {
fefu::hash_map<int, string> hmap(6);
hmap[1] = "test";
hmap.reserve(12);
CHECK(hmap.bucket_count() == 16);
CHECK(hmap[1] == "test");
hmap[-1] = "test2";
hmap.reserve(4);
CHECK(hmap.bucket_count() == 8);
CHECK(hmap[1] == "test");
CHECK(hmap[-1] == "test2");
}

TEST_CASE("load_factor") {
fefu::hash_map<int, string> hmap(10);
CHECK(hmap.load_factor() == 0.0);
hmap[4] = "ab";
CHECK(abs(hmap.load_factor() - 0.0625) < eps);
CHECK(abs(hmap.max_load_factor() - 0.4) > eps);
hmap.max_load_factor(0.6);
CHECK(abs(hmap.max_load_factor() - 0.6) < eps);
}

TEST_CASE("auto rehash") {
fefu::hash_map<int, string> hmap(20);
for (int i = 0; i < 40; i++) {
hmap[i] = "aba";
}

CHECK(hmap.bucket_count() >= 40);
CHECK(hmap.size() == 40);
}

TEST_CASE("contains()") {
fefu::hash_map<int, string> hmap(10);
CHECK(!hmap.contains(10));
hmap[10] = "ab";
CHECK(hmap.contains(10));
hmap[-1] = "bc";
CHECK(hmap.contains(-1));
}

TEST_CASE("operator==") {
fefu::hash_map<int, string> hmap1(10);
fefu::hash_map<int, string> hmap2(20);
fefu::hash_map<int, string> hmap3(20);



CHECK(hmap1 == hmap2);
CHECK(hmap2 == hmap3);

hmap1[4] = "ab";
hmap1[1] = "bc";

hmap2[1] = "ab";
hmap2[4] = "bc";
CHECK(!(hmap1 == hmap2));
CHECK(!(hmap1 == hmap3));

hmap3[1] = "ab";
hmap3[4] = "bc";
CHECK(hmap3 == hmap2);
CHECK(!(hmap3 == hmap1));
}

TEST_CASE("count()") {
fefu::hash_map<int, string> hmap(10);
CHECK(hmap.count(2) == 0);
hmap[2] = "ab";
hmap[1] = "d";
CHECK(hmap.count(1) == 1);
CHECK(hmap.count(2) == 1);
}

TEST_CASE("hash_function()") {
fefu::hash_map<int, string> hmap(10);
auto h = hmap.hash_function();
auto defaultHash = hash<int>{};
CHECK(h(4) == defaultHash(4));
CHECK(h(-1) == defaultHash(-1));
CHECK(h(INT_MAX) == defaultHash(INT_MAX));

fefu::hash_map<string, string> hmap2(10);
auto h2 = hmap2.hash_function();
auto defaultHash2 = hash<string>{};
CHECK(h2("") == defaultHash2(""));
CHECK(h2("abacaba") == defaultHash2("abacaba"));
}

TEST_CASE("key_eq") {
fefu::hash_map<int, string> hmap(10);
auto k = hmap.key_eq();
CHECK(!k(1, 2));
CHECK(k(0, 0));
CHECK(k(5, 5));

fefu::hash_map<string, string> hmap2(10);
auto k2 = hmap2.key_eq();
CHECK(!k2("aba", "aca"));
CHECK(k2("", ""));
CHECK(k2("test", "test"));
}

TEST_CASE("InputIterator constructor") {
vector<pair<int, string>> data = { pair<int, string>(1, "aba"),
                                   pair<int, string>(2, "caba"),
                                   pair<int, string>(1, "caba"),
                                   pair<int, string>(2, "aba"),
                                   pair<int, string>(1, "aba"),
                                   pair<int, string>(3, "test") };

fefu::hash_map<int, string> hmap(data.begin(), data.end(), 3);

REQUIRE(hmap.size() == data.size() - 3);
CHECK(hmap.at(1) == "aba");
CHECK(hmap.at(2) == "caba");
CHECK(hmap.at(3) == "test");
}


TEST_CASE("allocator constructor, get_allocator()") {
fefu::allocator<pair<const int, string>> t;
t.debug = 2;
fefu::hash_map<int, string> hmap(t);
auto alloc = hmap.get_allocator();

CHECK(typeid(t).name() == typeid(alloc).name());
CHECK(t.debug == 2);
}

TEST_CASE("Move constructor") {
fefu::allocator<pair<const int, string>> t;
t.debug = 3;

fefu::hash_map<int, string> hmap1(10);
hmap1[4] = "abc";
fefu::hash_map<int, string> hmap2(std::move(hmap1));

CHECK(hmap2.bucket_count() == 16);
CHECK(hmap2[4] == "abc");

fefu::hash_map<int, string> hmap3(10);
hmap3[4] = "abc";
fefu::hash_map<int, string> hmap4(std::move(hmap3), t);
CHECK(hmap4.bucket_count() == 16);
CHECK(hmap4.get_allocator().debug == 3);
CHECK(hmap4[4] == "abc");
}

TEST_CASE("Copy constructor") {
fefu::allocator<pair<const int, string>> t;
t.debug = 5;

fefu::hash_map<int, string> hmap1(10);
hmap1[4] = "abc";
fefu::hash_map<int, string> hmap2(hmap1);

CHECK(hmap2.bucket_count() == 16);
CHECK(hmap2[4] == "abc");

fefu::hash_map<int, string> hmap3(hmap1, t);
CHECK(hmap3.bucket_count() == 16);
CHECK(hmap3[4] == "abc");
CHECK(t.debug == 5);
}

TEST_CASE("Init list constructor") {
fefu::hash_map<int, string> hmap = { pair<int, string>(1, "aba"),
                                     pair<int, string>(2, "caba"),
                                     pair<int, string>(1, "caba"),
                                     pair<int, string>(2, "aba"),
                                     pair<int, string>(1, "aba"),
                                     pair<int, string>(3, "test") };

REQUIRE(hmap.size() == 3);
REQUIRE(hmap.contains(1));
REQUIRE(hmap.contains(2));
REQUIRE(hmap.contains(3));
CHECK(hmap[1] == "aba");
CHECK(hmap[2] == "caba");
CHECK(hmap[3] == "test");
}

TEST_CASE("Assignment operators") {
fefu::hash_map<int, string> hmap1(10);
fefu::hash_map<int, string> hmap2(20);
hmap2 = hmap1;
CHECK(hmap1 == hmap2);

fefu::hash_map<int, string> hmap3(30);
hmap3 = std::move(hmap1);
CHECK(hmap3 == hmap2);

hmap2 = { pair<int, string>(1, "aba"),
          pair<int, string>(2, "caba"),
          pair<int, string>(1, "caba"),
          pair<int, string>(2, "aba"),
          pair<int, string>(1, "aba"),
          pair<int, string>(3, "test") };
REQUIRE(hmap2.size() == 3);
REQUIRE(hmap2.contains(1));
REQUIRE(hmap2.contains(2));
REQUIRE(hmap2.contains(3));
CHECK(hmap2[1] == "aba");
CHECK(hmap2[2] == "caba");
CHECK(hmap2[3] == "test");
}

TEST_CASE("Size") {
fefu::hash_map<int, int> hmap(10);

CHECK(hmap.empty());

hmap[1] = 2;
hmap[2] = 3;
hmap[3] = 5;

CHECK(!hmap.empty());
CHECK(hmap.size() == 3);
CHECK(hmap.max_size() == std::numeric_limits<std::size_t>::max());
}


TEST_CASE("Non-const iterators") {
fefu::hash_map<int, string> hmap(20);

CHECK(hmap.begin() == hmap.end());

hmap[1] = "a";
hmap[-1] = "b";
hmap[3] = "c";
hmap[6] = "d";
fefu::hash_map_iterator<std::pair<const int, string>> it = hmap.begin();
auto itBegin = *it;

fefu::hash_map_const_iterator<std::pair<const int, string>> constIt(it);
CHECK(*it == *constIt);

CHECK(hmap[it->first] == it->second);

fefu::hash_map_iterator<std::pair<const int, string>> tmp = it;
auto tmp2 = it++;
CHECK(tmp2 == tmp);
CHECK(tmp != it);

CHECK(hmap[(*it).first] == (*it).second);
++it;
CHECK(hmap[it->first] == it->second);
auto tmp0 = hmap[it->first];
tmp = it++;
CHECK(hmap[tmp->first] == tmp0);
CHECK(++it == hmap.end());
CHECK_THROWS(++it);
}

TEST_CASE("Const iterators") {
const fefu::hash_map<int, string> emptyHmap;
CHECK(emptyHmap.begin() == emptyHmap.end());
CHECK(emptyHmap.cbegin() == emptyHmap.cend());

fefu::hash_map<int, string> hmap0(20);
hmap0[1] = "a";
hmap0[-1] = "b";
hmap0[3] = "c";
hmap0[6] = "d";
const fefu::hash_map<int, string> hmap(hmap0);
auto it = hmap.begin();
CHECK(it == hmap.cbegin());
CHECK(hmap.at(it->first) == it->second);

auto tmp = it;
auto tmp2 = it++;
CHECK(tmp2 == tmp);
CHECK(tmp != it);

CHECK(hmap.at((*it).first) == (*it).second);
++it;
CHECK(hmap.at(it->first) == it->second);
auto tmp0 = hmap.at(it->first);
auto tmp3 = it++;
CHECK(hmap.at(tmp3->first) == tmp0);
CHECK(++it == hmap.end());
CHECK_THROWS(++it);
}

TEST_CASE("erase") {
fefu::hash_map<int, string> hmap(10);
hmap[1] = "a";
hmap[2] = "b";
auto val = hmap.begin()->first;
hmap.erase(hmap.cbegin());
CHECK(hmap.size() == 1);
CHECK(!hmap.contains(val));

hmap[4] = "d";
size_t count = hmap.erase(4);
CHECK(count == 1);
CHECK(!hmap.contains(4));

count = hmap.erase(5);
CHECK(count == 0);

for (int i = 0; i < 5; i++) {
hmap[i] = "test";
}

fefu::hash_map<int, string> hmapCopy(hmap);

auto it = hmap.erase(hmap.begin(), hmap.end());
REQUIRE(hmap.size() == 0);
CHECK(it == hmap.end());

hmapCopy.clear();
CHECK(hmap.size() == 0);
}

TEST_CASE("find") {
fefu::hash_map<int, string> hmap;
CHECK(hmap.find(2) == hmap.end());
hmap[1] = "a";
hmap[2] = "b";
hmap[-1] = "c";

CHECK(hmap.find(3) == hmap.end());
CHECK(hmap.find(2)->second == "b");

const fefu::hash_map<int, string> constHmap(hmap);
CHECK(constHmap.find(3) == constHmap.end());
CHECK(constHmap.find(1)->second == "a");
}

TEST_CASE("insert") {
fefu::hash_map<int, string> hmap;
auto it = hmap.insert(make_pair(0, "abaca"));
CHECK(hmap.contains(0));
CHECK(hmap.at(0) == "abaca");
CHECK(it.first != hmap.end());
CHECK(it.second);

it = hmap.insert(make_pair(0, "cabada"));
CHECK(hmap.at(0) == "abaca");
CHECK(!it.second);

pair<const int, string> constPair(1, "test");
it = hmap.insert(constPair);
CHECK(hmap.contains(1));
CHECK(hmap.at(1) == "test");
CHECK(it.first != hmap.end());
CHECK(it.second);

pair<const int, string> constPair2(1, "null");
it = hmap.insert(constPair2);
CHECK(hmap.at(1) == "test");
CHECK(!it.second);
}

TEST_CASE("insert range") {
fefu::hash_map<int, string> hmap;
hmap.insert(make_pair(0, "abaca"));
hmap.insert(make_pair(1, "test"));
fefu::hash_map<int, string> hmap2(hmap);

vector<pair<int, string>> inputRange = { { 0, "test0" }, { 1, "test1" }, { 2, "test2"}, { 3, "test3" } };
hmap.insert(inputRange.begin(), inputRange.end());
CHECK(hmap.size() == 4);
CHECK(hmap.at(0) == "abaca");
CHECK(hmap.at(1) == "test");
CHECK(hmap.at(2) == "test2");
CHECK(hmap.at(3) == "test3");

hmap2.insert({ { 0, "test0" }, { 1, "test1" }, { 2, "test2"}, { 3, "test3" } });
CHECK(hmap2.size() == 4);
CHECK(hmap2.at(0) == "abaca");
CHECK(hmap2.at(1) == "test");
CHECK(hmap2.at(2) == "test2");
CHECK(hmap2.at(3) == "test3");
}

TEST_CASE("insert_or_assign") {
fefu::hash_map<int, string> hmap;
hmap.insert(make_pair(0, "abaca"));
hmap.insert(make_pair(1, "test"));

auto res = hmap.insert_or_assign(0, "testb");
CHECK(hmap.at(0) == "testb");
CHECK(!res.second);

res = hmap.insert_or_assign(2, "testc");
CHECK(hmap.at(2) == "testc");
CHECK(res.second);

const int a = 2;
res = hmap.insert_or_assign(a, "constTest");
CHECK(hmap.at(2) == "constTest");
CHECK(!res.second);
}

TEST_CASE("insert into deleted") {
fefu::hash_map<int, string> hmap = { pair<int, string>(1, "aba"),
                                     pair<int, string>(2, "caba"),
                                     pair<int, string>(3, "caba"),
                                     pair<int, string>(4, "aba"),
                                     pair<int, string>(5, "aba"),
                                     pair<int, string>(6, "test") };
hmap.erase_if([](const pair<int, string> tmp) { return true; });
for (int i = 0; i < 100; i++) {
hmap[i] = "d";
}
CHECK(hmap.size() >= 100);

hmap.erase_if([](const pair<int, string> tmp) { return true; });
for (int i = 0; i < 100; i++) {
hmap.insert(make_pair(i, "d"));
}
CHECK(hmap.size() >= 100);

hmap.erase_if([](const pair<int, string> tmp) { return true; });
for (int i = 0; i < 100; i++) {
hmap.emplace(make_pair(i, "d"));
}
CHECK(hmap.size() >= 100);

hmap.erase_if([](const pair<int, string> tmp) { return true; });
for (int i = 0; i < 100; i++) {
hmap.try_emplace(std::move(i), "d");
}
CHECK(hmap.size() >= 100);

hmap.erase_if([](const pair<int, string> tmp) { return true; });
for (int i = 0; i < 100; i++) {
hmap.insert_or_assign(i, "d");
}
CHECK(hmap.size() >= 100);

hmap.erase_if([](const pair<int, string> tmp) { return true; });
for (int i = 0; i < 100; i++) {
hmap.insert_or_assign(std::move(i), "d");
}
CHECK(hmap.size() >= 100);

hmap.erase_if([](const pair<int, string> tmp) { return true; });
for (int i = 0; i < 100; i++) {
hmap[std::move(i)] = "d";
}
CHECK(hmap.size() >= 100);
}

TEST_CASE("exceptions") {
fefu::hash_map<int, string> hmap(10);
const fefu::hash_map<int, string> hmap2(hmap);
CHECK_THROWS(hmap.at(2));
CHECK_THROWS(hmap2.at(2));
CHECK_THROWS(hmap.bucket(2));
CHECK_THROWS(hmap.max_load_factor(1.2));
CHECK_THROWS(hmap.max_load_factor(-0.4));

fefu::hash_map_iterator<pair<int, string>> iter;
fefu::hash_map_const_iterator<pair<int, string>> constIter;
CHECK_THROWS(*iter);
CHECK_THROWS(*constIter);

CHECK_THROWS(hmap.erase(hmap.cend()));
}
TEST_CASE("beg = end"){
fefu::hash_map<int,int> hmap;
CHECK(hmap.begin() == hmap.end());
}

