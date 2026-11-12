#ifndef INCLUDED_noise_bias_hpp_
#define INCLUDED_noise_bias_hpp_

#include <vector>
#include "noise_shell.hpp"

// planner flag of the FFTW plan (FFTW_ESTIMATE gives reproducible plans, used by the regression tests)
#ifndef STOLAS_FFTW_PLANNER
#define STOLAS_FFTW_PLANNER FFTW_MEASURE
#endif

inline fftw_complex *inhalf;
inline double *out;
inline fftw_plan plan;

// ---- shell table ---------------------------------------------------------------------------
// The noise (and the bias) lives only on the thin shell |k| in [nsigma-dn/2, nsigma+dn/2], i.e. on
// a tiny part of the lattice. The lattice points are therefore sorted once by
//   m = nx^2 + ny^2 + nz^2   (indices folded to -N/2..N/2)
// so that the shell of every step is a contiguous range of m, i.e. of shellIdx, and the noise
// generation touches only the shell points instead of scanning the whole N^3 lattice.
//   shellIdx[shellOffset[m] .. shellOffset[m+1])  half-cube points (kz<=N/2, the part used by the
//                                                   c2r transform) with that m; stored as the
//                                                   full-cube index i*N*N+j*N+k
//   shellCountFull[m]                              number of points of the FULL cube with that m
constexpr int shellMmax = 3*(NLnoise/2)*(NLnoise/2);
inline std::vector<int> shellOffset;
inline std::vector<int> shellIdx;
inline std::vector<long> shellCountFull;

inline void init_shell_table() {
  shellOffset.assign(shellMmax+2, 0);
  shellCountFull.assign(shellMmax+1, 0);

  for (int nx = 0; nx < NLnoise; nx++) {
    int nxt = (nx<=NLnoise/2 ? nx : nx-NLnoise);
    for (int ny = 0; ny < NLnoise; ny++) {
      int nyt = (ny<=NLnoise/2 ? ny : ny-NLnoise);
      for (int nz = 0; nz < NLnoise; nz++) {
        int nzt = (nz<=NLnoise/2 ? nz : nz-NLnoise);
        int m = nxt*nxt + nyt*nyt + nzt*nzt;
        shellCountFull[m]++;
        if (nz <= NLnoise/2) shellOffset[m+1]++;
      }
    }
  }
  for (int m = 0; m <= shellMmax; m++) shellOffset[m+1] += shellOffset[m];

  shellIdx.resize(shellOffset[shellMmax+1]);
  std::vector<int> pos(shellOffset.begin(), shellOffset.end()-1);
  for (int nx = 0; nx < NLnoise; nx++) {
    int nxt = (nx<=NLnoise/2 ? nx : nx-NLnoise);
    for (int ny = 0; ny < NLnoise; ny++) {
      int nyt = (ny<=NLnoise/2 ? ny : ny-NLnoise);
      for (int nz = 0; nz <= NLnoise/2; nz++) {
        int m = nxt*nxt + nyt*nyt + nz*nz;
        shellIdx[pos[m]++] = nx*NLnoise*NLnoise + ny*NLnoise + nz;
      }
    }
  }
}

// range of m with |sqrt(m) - nsigma| <= dn/2 (the same condition as the original full scan);
// returns the number of lattice points of the FULL cube in the shell, mlo/mhi = -1 if empty
inline long shell_range(double nsigma, int &mlo, int &mhi) {
  mlo = -1; mhi = -1;
  long count = 0;
  for (int m = 0; m <= shellMmax; m++) {
    if (std::abs(std::sqrt((double)m) - nsigma) <= dn/2.) {
      if (mlo < 0) mlo = m;
      mhi = m;
      count += shellCountFull[m];
    }
  }
  return count;
}

// c2r transforms overwrite their input, so the half-cube input is cleared before every fill
inline void clear_inhalf() {
  double *p = reinterpret_cast<double*>(inhalf);
  const long n = 2L*NLnoiseHalfAll;
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (long i = 0; i < n; i++) p[i] = 0.0;
}

template <typename Array>
inline void zero_array(Array &a) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < NLnoiseAll; i++) a[i] = 0.0;
}

inline void init_fftw_global() {
  static bool is_initialized = false;
  if (!is_initialized) {
    init_shell_table();
    fftw_init_threads();
    #ifdef _OPENMP
      fftw_plan_with_nthreads(omp_get_max_threads());
    #endif

    inhalf = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NLnoiseHalfAll);
    out = fftw_alloc_real(NLnoiseAll);

    if(splan && FFTWwisdom){
      plan = fftw_plan_dft_c2r_3d(NLnoise, NLnoise, NLnoise, inhalf, out, FFTW_PATIENT);
      fftw_export_wisdom_to_filename((sdatadir+"/wisdom"+std::to_string(NLnoise)+"_c2r.dat").c_str());
      std::cout << "Make the FFTW plan." << std::endl;
      FFTwisdomFirst=true;
    }
    else if (FFTWwisdom){
      fftw_import_wisdom_from_filename((sdatadir+"/wisdom"+std::to_string(NLnoise)+"_c2r.dat").c_str());
      plan = fftw_plan_dft_c2r_3d(NLnoise, NLnoise, NLnoise, inhalf, out, FFTW_WISDOM_ONLY);
    }
    else{
      std::cout << "Creating the FFTW plan (may take a while for large lattices)..." << std::flush;
      plan = fftw_plan_dft_c2r_3d(NLnoise, NLnoise, NLnoise, inhalf, out, STOLAS_FFTW_PLANNER);
      std::cout << " done." << std::endl;
    }

    is_initialized = true;
  }
}

void dwlist_gen(double N, int seed, size_t step, int Nfield, int InterpolatingNo) {
  double nsigma = sigma*exp(N);
  int mlo, mhi;
  long count = shell_range(nsigma, mlo, mhi);

  if (count==0) {
    zero_array(dwlist[Nfield]);
    return;
  }

  clear_inhalf();
  const double sqrt_count = sqrt((double)count);
  constexpr int NN = NLnoise*NLnoise;
  const long pbeg = shellOffset[mlo], pend = shellOffset[mhi+1];

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (long p = pbeg; p < pend; p++) {
    const int idx = shellIdx[p];
    const int i = idx / NN, j = (idx / NLnoise) % NLnoise, k = idx % NLnoise;
    double z0, z1, re, im;

    if (realpoint(i, j, k, NLnoise)) {
      point_normals((uint64_t)seed, (uint64_t)InterpolatingNo, (uint64_t)step, (uint64_t)Nfield, (uint64_t)idx, z0, z1);
      re = z0;
      im = 0.0;
    } else if (complexpoint(i, j, k, NLnoise)) {
      point_normals((uint64_t)seed, (uint64_t)InterpolatingNo, (uint64_t)step, (uint64_t)Nfield, (uint64_t)idx, z0, z1);
      re = z0 * inv_sqrt2;
      im = z1 * inv_sqrt2;
    } else {
      // the other half of a conjugate pair: complex conjugate of the independent partner point (-k)
      const int ip = (i == 0 ? 0 : NLnoise - i);
      const int jp = (j == 0 ? 0 : NLnoise - j);
      const int kp = (k == 0 ? 0 : NLnoise - k);
      const int pidx = ip * NN + jp * NLnoise + kp;
      if (realpoint(ip, jp, kp, NLnoise)) {
        point_normals((uint64_t)seed, (uint64_t)InterpolatingNo, (uint64_t)step, (uint64_t)Nfield, (uint64_t)pidx, z0, z1);
        re = z0;
        im = -0.0;
      } else if (complexpoint(ip, jp, kp, NLnoise)) {
        point_normals((uint64_t)seed, (uint64_t)InterpolatingNo, (uint64_t)step, (uint64_t)Nfield, (uint64_t)pidx, z0, z1);
        re = z0 * inv_sqrt2;
        im = -(z1 * inv_sqrt2);
      } else {
        re = 0.0;
        im = -0.0;
      }
    }

    const int idxhalf = i * NLnoise * NLnoiseHalf + j * NLnoiseHalf + k;
    inhalf[idxhalf][0] = re / sqrt_count;
    inhalf[idxhalf][1] = im / sqrt_count;
  }

  // the transform writes directly into dwlist (alignas(64), see STOLAS.hpp)
  fftw_execute_dft_c2r(plan, inhalf, dwlist[Nfield].data());
}


void biaslist1D(double N) {
  double nsigma = sigma*exp(N);
  int mlo, mhi;
  long count = shell_range(nsigma, mlo, mhi);

  if (count==0) {
    zero_array(biaslist[0]);
    return;
  }

  clear_inhalf();
  const double val = 1.0/(double)count;
  constexpr int NN = NLnoise*NLnoise;
  const long pbeg = shellOffset[mlo], pend = shellOffset[mhi+1];

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (long p = pbeg; p < pend; p++) {
    const int idx = shellIdx[p];
    const int i = idx / NN, j = (idx / NLnoise) % NLnoise, k = idx % NLnoise;
    const int idxhalf = i * NLnoise * NLnoiseHalf + j * NLnoiseHalf + k;
    inhalf[idxhalf][0] = val;
    inhalf[idxhalf][1] = 0.0;
  }

  fftw_execute_dft_c2r(plan, inhalf, biaslist[0].data());
}

#endif
