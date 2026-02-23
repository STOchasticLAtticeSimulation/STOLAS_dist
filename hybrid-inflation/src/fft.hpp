#ifndef INCLUDED_fft_hpp_
#define INCLUDED_fft_hpp_

#define _USR_MATH_DEFINES
#include <cmath>
#include <complex>
#include <vector>
#include <mutex>

#include <fftw3.h>

std::mutex fftw_plan_mutex;

std::vector<std::vector<std::vector<std::complex<double>>>> fft_fftw(const std::vector<std::vector<std::vector<std::complex<double>>>>& bk);

std::vector<std::vector<std::vector<std::complex<double>>>> fft_fftw(const std::vector<std::vector<std::vector<std::complex<double>>>>& bk) {
  fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NLnoise * NLnoise * NLnoise);
  fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NLnoise * NLnoise * NLnoise);

  // bk (3D std::vector) -> in (1D fftw_complex) に変換
  int idx = 0;
  for (int i = 0; i < NLnoise; ++i) {
    for (int j = 0; j < NLnoise; ++j) {
      for (int k = 0; k < NLnoise; ++k) {
        in[idx][0] = bk[i][j][k].real();
        in[idx][1] = bk[i][j][k].imag();
        idx++;
      }
    }
  }

  // fftw_plan plan = fftw_plan_dft_3d(NLnoise, NLnoise, NLnoise, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_plan plan;
    {
      std::lock_guard<std::mutex> lock(fftw_plan_mutex);
      plan = fftw_plan_dft_3d(NLnoise, NLnoise, NLnoise, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    }

  fftw_execute(plan);

  // out (1D fftw_complex) -> biaslattice (3D std::vector) に変換
  std::vector<std::vector<std::vector<std::complex<double>>>> biaslattice(NLnoise, std::vector<std::vector<std::complex<double>>>(NLnoise, std::vector<std::complex<double>>(NLnoise)));

  idx = 0;
  for (int i = 0; i < NLnoise; ++i) {
    for (int j = 0; j < NLnoise; ++j) {
      for (int k = 0; k < NLnoise; ++k) {
        biaslattice[i][j][k] = std::complex<double>(out[idx][0], out[idx][1]);
        idx++;
      }
    }
  }

  fftw_destroy_plan(plan);
  fftw_free(in);
  fftw_free(out);

  return biaslattice;
}


#endif
