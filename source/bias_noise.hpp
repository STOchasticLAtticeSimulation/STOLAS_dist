#ifndef INCLUDED_bias_noise_hpp_
#define INCLUDED_bias_noise_hpp_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <sys/time.h>

#include "parameters.hpp"
#include "fft.hpp"
#include "vec_op.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif


const std::complex<double> II(0,1);

// random distribution
std::random_device seed;
std::mt19937 engine(seed());
std::normal_distribution<> dist(0., 1.);

// useful macro
#define LOOP for(int i = 0; i < NL; i++) for(int j = 0; j < NL; j++) for(int k = 0; k < NL; k++)


bool innsigma(int nx, int ny, int nz, int Num, double nsigma, double dn); // judge if point is in nsigma sphere shell
bool realpoint(int nx, int ny, int nz, int Num); // judge real point
bool complexpoint(int nx, int ny, int nz, int Num); // judge independent complex point


// functions
bool innsigma(int nx, int ny, int nz, int Num, double nsigma, double dn) {
  int nxt, nyt, nzt; // shifted index
  if (nx<=Num/2) {
    nxt = nx;
  } else {
    nxt = nx-Num;
  }

  if (ny<=Num/2) {
    nyt = ny;
  } else {
    nyt = ny-Num;
  }

  if (nz<=Num/2) {
    nzt = nz;
  } else {
    nzt = nz-Num;
  }

  return std::abs(sqrt(nxt*nxt + nyt*nyt + nzt*nzt) - nsigma) <= dn/2.;
}

bool realpoint(int nx, int ny, int nz, int Num) {
  return (nx==0||nx==Num/2) && (ny==0||ny==Num/2) && (nz==0||nz==Num/2);
}

bool complexpoint(int nx, int ny, int nz, int Num) {
  int nxt, nyt, nzt; // shifted index
  if (nx<=Num/2) {
    nxt = nx;
  } else {
    nxt = nx-Num;
  }

  if (ny<=Num/2) {
    nyt = ny;
  } else {
    nyt = ny-Num;
  }

  if (nz<=Num/2) {
    nzt = nz;
  } else {
    nzt = nz-Num;
  }

  return (1<=nxt && nxt!=Num/2 && nyt!=Num/2 && nzt!=Num/2) ||
    (nxt==Num/2 && nyt!=Num/2 && 1<=nzt && nzt!=Num/2) || (nxt!=Num/2 && 1<=nyt && nyt!=Num/2 && nzt==Num/2) || (1<=nxt && nxt!=Num/2 && nyt==Num/2 && nzt!=Num/2) ||
    (nxt==0 && nyt!=Num/2 && 1<=nzt && nzt!=Num/2) ||
    (nxt==Num/2 && nyt==Num/2 && 1<=nzt && nzt!=Num/2) || (nxt==Num/2 && 1<=nyt && nyt!=Num/2 && nzt==Num/2) || (1<=nxt && nxt!=Num/2 && nyt==Num/2 && nzt==Num/2) ||
    (nxt==0 && 1<=nyt && nyt!=Num/1 && nzt==0) ||
    (nxt==Num/2 && 1<=nyt && nyt!=Num/2 && nzt==0) || (1<=nxt && nxt!=Num/2 && nyt==0 && nzt==Num/2) || (nxt==0 && nyt==Num/2 && 1<=nzt && nzt!=Num/2);
}

#endif