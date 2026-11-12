#ifndef INCLUDED_noise_shell_hpp_
#define INCLUDED_noise_shell_hpp_

#include <cstdint>

// judge real point
inline bool realpoint(int nx, int ny, int nz, int Num) {
  return (nx==0||nx==Num/2) && (ny==0||ny==Num/2) && (nz==0||nz==Num/2);
}

// judge independent complex point
inline bool complexpoint(int nx, int ny, int nz, int Num) {
  int nxt = (nx<=Num/2 ? nx : nx-Num);
  int nyt = (ny<=Num/2 ? ny : ny-Num);
  int nzt = (nz<=Num/2 ? nz : nz-Num);

  return (1<=nxt && nxt!=Num/2 && nyt!=Num/2 && nzt!=Num/2) ||
    (nxt==Num/2 && nyt!=Num/2 && 1<=nzt && nzt!=Num/2) || (nxt!=Num/2 && 1<=nyt && nyt!=Num/2 && nzt==Num/2) || (1<=nxt && nxt!=Num/2 && nyt==Num/2 && nzt!=Num/2) ||
    (nxt==0 && nyt!=Num/2 && 1<=nzt && nzt!=Num/2) ||
    (nxt==Num/2 && nyt==Num/2 && 1<=nzt && nzt!=Num/2) || (nxt==Num/2 && 1<=nyt && nyt!=Num/2 && nzt==Num/2) || (1<=nxt && nxt!=Num/2 && nyt==Num/2 && nzt!=Num/2) ||
    (nxt==0 && 1<=nyt && nyt!=Num/2 && nzt==0) ||
    (nxt==Num/2 && 1<=nyt && nyt!=Num/2 && nzt==0) || (1<=nxt && nxt!=Num/2 && nyt==0 && nzt==Num/2) || (nxt==0 && nyt==Num/2 && 1<=nzt && nzt!=Num/2);
}

inline uint64_t splitmix64_next(uint64_t &state) {
  uint64_t z = (state += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}

inline uint64_t hash_mix(uint64_t x) {
  uint64_t state = x;
  return splitmix64_next(state);
}

inline void point_normals(uint64_t seed, uint64_t level, uint64_t step, uint64_t field, uint64_t idx, double &z0, double &z1) {
  uint64_t s = hash_mix(seed);
  s = hash_mix(s ^ level);
  s = hash_mix(s ^ step);
  s = hash_mix(s ^ field);
  s = hash_mix(s ^ idx);

  uint64_t r1 = splitmix64_next(s);
  uint64_t r2 = splitmix64_next(s);
  double u1 = ((r1 >> 11) + 1) * (1.0/9007199254740992.0);
  double u2 = (r2 >> 11) * (1.0/9007199254740992.0);

  double radius = sqrt(-2.0 * log(u1));
  double theta = 2.0 * M_PI * u2;
  z0 = radius * cos(theta);
  z1 = radius * sin(theta);
}

#endif
