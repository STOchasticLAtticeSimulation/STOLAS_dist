#ifndef INCLUDED_fft_hpp_
#define INCLUDED_fft_hpp_

#define _USR_MATH_DEFINES
#include <cmath>
#include <complex>
#include <fftw3.h>

std::array<std::complex<double>,NLnoiseHalfAll> bkspectrum{};

void fft_1D_real(const std::array<double,NLnoiseAll>& bk) {
  double* in = fftw_alloc_real(NLnoiseAll);
  fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NLnoiseHalfAll);

  for (int i = 0; i < NLnoiseAll; ++i) {
    in[i] = bk[i];
  }

  fftw_plan plan = fftw_plan_dft_r2c_3d(NLnoise, NLnoise, NLnoise, in, out, FFTW_ESTIMATE);
  fftw_execute(plan);

  for (int i = 0; i < NLnoiseHalfAll; ++i) {
    bkspectrum[i] = out[i][0] + II*out[i][1];
  }

  fftw_destroy_plan(plan);
  fftw_free(in);
  fftw_free(out);
}


#endif
