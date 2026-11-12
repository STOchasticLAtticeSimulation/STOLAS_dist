#ifndef INCLUDED_array_op_hpp_
#define INCLUDED_array_op_hpp_

#include <array>


template <class T1, class T2, std::size_t Nsize>
std::array<T1,Nsize> operator+(const std::array<T1,Nsize> &v1, const std::array<T2,Nsize> &v2) {
  std::array<T1,Nsize> ans = v1;
  for (size_t i = 0; i < Nsize; ++i)
    ans[i] += v2[i];
  return ans;
}

template <class T1, class T2, std::size_t Nsize>
std::array<T1,Nsize> operator-(const std::array<T1,Nsize> &v1, const std::array<T2,Nsize> &v2) {
  std::array<T1,Nsize> ans = v1;
  for (size_t i = 0; i < Nsize; ++i)
    ans[i] -= v2[i];
  return ans;
}

template <class T1, class T2, std::size_t Nsize>
std::array<T1,Nsize>& operator+=(std::array<T1,Nsize> &v1, const std::array<T2,Nsize> &v2) {
  for (size_t i = 0; i < Nsize; ++i)
    v1[i] += v2[i];
  return v1;
}

template <class T1, class T2, std::size_t Nsize>
std::array<T1,Nsize>& operator-=(std::array<T1,Nsize> &v1, const std::array<T2,Nsize> &v2) {
  for (size_t i = 0; i < Nsize; ++i)
    v1[i] -= v2[i];
  return v1;
}

template <class Tv, class Tc, std::size_t Nsize>
std::array<Tv,Nsize> operator*(const std::array<Tv,Nsize> &v, const Tc &c) {
  std::array<Tv,Nsize> ans = v;
  for (Tv &e : ans)
    e *= c;
  return ans;
}

template <class Tc, class Tv, std::size_t Nsize>
std::array<Tv,Nsize> operator*(const Tc &c, const std::array<Tv,Nsize> &v) {
  std::array<Tv,Nsize> ans = v;
  for (Tv &e : ans)
    e *= c;
  return ans;
}

template <class Tv, class Tc, std::size_t Nsize>
std::array<Tv,Nsize>& operator*=(std::array<Tv,Nsize> &v, const Tc &c) {
  for (Tv &e : v)
    e *= c;
  return v;
}

template <class Tv, class Tc, std::size_t Nsize>
std::array<Tv,Nsize> operator/(const std::array<Tv,Nsize> &v, const Tc &c) {
  std::array<Tv,Nsize> ans = v;
  for (Tv &e : ans)
    e /= c;
  return ans;
}

template <class Tv, class Tc, std::size_t Nsize>
std::array<Tv,Nsize>& operator/=(std::array<Tv,Nsize> &v, const Tc &c) {
  for (Tv &e : v)
    e /= c;
  return v;
}

namespace array_op {

  inline void init(double &x) {
    x = 0;
  }

  template <class T, std::size_t Nsize>
  void init(std::array<T,Nsize> &v) {
    for (T &e : v) {
      init(e);
    }
  }

}

#endif