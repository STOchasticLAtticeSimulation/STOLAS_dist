#ifndef INCLUDED_vec_op_hpp_
#define INCLUDED_vec_op_hpp_

#include <vector>


template <class T1, class T2>
std::vector<T1> operator+(const std::vector<T1> &v1, const std::vector<T2> &v2) {
  std::vector<T1> ans = v1;
  for (size_t i = 0, size = ans.size(); i < size; ++i)
    ans[i] += v2[i];
  return ans;
}

template <class T1, class T2>
std::vector<T1> operator-(const std::vector<T1> &v1, const std::vector<T2> &v2) {
  std::vector<T1> ans = v1;
  for (size_t i = 0, size = ans.size(); i < size; ++i)
    ans[i] -= v2[i];
  return ans;
}

template <class T1, class T2>
std::vector<T1>& operator+=(std::vector<T1> &v1, const std::vector<T2> &v2) {
  for (size_t i = 0, size = v1.size(); i < size; ++i)
    v1[i] += v2[i];
  return v1;
}

template <class T1, class T2>
std::vector<T1>& operator-=(std::vector<T1> &v1, const std::vector<T2> &v2) {
  for (size_t i = 0, size = v1.size(); i < size; ++i)
    v1[i] -= v2[i];
  return v1;
}

template <class Tv, class Tc>
std::vector<Tv> operator*(const std::vector<Tv> &v, const Tc &c) {
  std::vector<Tv> ans = v;
  for (Tv &e : ans)
    e *= c;
  return ans;
}

template <class Tc, class Tv>
std::vector<Tv> operator*(const Tc &c, const std::vector<Tv> &v) {
  std::vector<Tv> ans = v;
  for (Tv &e : ans)
    e *= c;
  return ans;
}

template <class Tv, class Tc>
std::vector<Tv>& operator*=(std::vector<Tv> &v, const Tc &c) {
  for (Tv &e : v)
    e *= c;
  return v;
}

template <class Tv, class Tc>
std::vector<Tv> operator/(const std::vector<Tv> & v, const Tc &c) {
  std::vector<Tv> ans = v;
  for (Tv &e : ans)
    e /= c;
  return ans;
}

template <class Tv, class Tc>
std::vector<Tv>& operator/=(std::vector<Tv> &v, const Tc &c) {
  for (Tv &e : v)
    e /= c;
  return v;
}

namespace vec_op {
  void init(double &x) {
    x = 0;
  }
  
  template <class T>
  void init(std::vector<T> &v) {
    for (T &e : v) {
      init(e);
    }
  }
}

#endif
