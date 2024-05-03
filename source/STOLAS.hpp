#ifndef INCLUDED_STOLAS_
#define INCLUDED_STOLAS_

#define _USR_MATH_DEFINES
#include <cmath>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <functional>
#include <complex>
#include <sys/time.h>

#include "parameters.hpp"

#define euler_gamma 0.57721566490153286061

#ifdef _OPENMP
#include <omp.h>
#endif

class STOLAS
{
protected:
  const std::string Nfileprefix = sdatadir + "/Nmap_";
  const std::string Hfileprefix = sdatadir + "/H_";
  const std::string pifileprefix = sdatadir + "/pi_";
  const std::string powfileprefix = sdatadir + "/power_";
  const std::string cmpfileprefix = sdatadir + "/compaction_";
  const std::string prbfileprefix = sdatadir + "/probabilities";
  const std::string powsfileprefix = sdatadir + "/powers";
  bool noisefilefail, biasfilefail;

  std::string model;
  int NL, noisefileNo;
  double dN, bias, Nbias, dNbias;
  std::ifstream noisefile, biasfile;
  std::ofstream Nfile, Hfile, pifile, powfile, cmpfile, prbfile, powsfile;
  std::vector<double> phii;
  std::vector<std::vector<double>> noisedata, biasdata, Hdata, pidata;
  std::vector<double> Ndata;
  std::vector<std::vector<std::vector<std::complex<double>>>> Nmap3D;

public:
  STOLAS(){}
  STOLAS(std::string Model, double DN, std::string sourcedir, int NoisefileNo, std::vector<double> Phii, double Bias, double NBias, double DNbias);

  bool checknoisefile();
  bool checkbiasfile();
  bool noisebiassize();
  bool Nfilefail();
  
  double VV(double phi);
  double Vp(double phi);
  double calPphi(std::vector<double> &phi);
 
  void dNmap();
  void animation();
  void spectrum();
  void compaction();

  double ep(double phi, double pi);
  double hubble(double phi, double pi);

  std::vector<double> dphidN(double N, std::vector<double> phi);

  void RK4(double &t, std::vector<double> &x, double dt);
  void RK4Mbias(double &N, std::vector<double> &phi, double dN, double dw, double Bias);
};

#endif
