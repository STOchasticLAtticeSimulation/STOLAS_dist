#ifndef INCLUDED_STOLAS_
#define INCLUDED_STOLAS_

#define _USR_MATH_DEFINES
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <functional>
#include <random>
#include <sys/time.h>
#include <complex>
#include <boost/numeric/odeint.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

constexpr double euler_gamma = 0.57721566490153286061;
constexpr double LOG2 = 0.69314718055994530942;
const std::complex<double> II(0, 1);
std::normal_distribution<> dist(0., 1.);


#include "parameters.hpp"

// Importance sampling switches: weight and compaction need the bias.
constexpr bool doWeight = sbias && sweight; // output ln W
constexpr bool doCompaction = sbias && scompaction; // compute/output the compaction function
constexpr bool recordWeight = doWeight || doCompaction; // compaction uses ln W, so record it
constexpr bool doPower = spower && !sbias; // the spectrum of biased samples is not meaningful, so only without bias
const double biasamp = sbias ? bias : 0.; // amplitude of the bias actually used
const double sqrt_dN = std::sqrt(dN);
constexpr int NLnoiseAll = NLnoise*NLnoise*NLnoise;
constexpr int NLnoiseHalf = NLnoise/2 + 1; // halved (fastest) dimension for r2c/c2r FFTW transforms
constexpr int NLnoiseHalfAll = NLnoise*NLnoise*NLnoiseHalf;
const int totalstep = ceil(log((NLnoise/2-1)/sigma)/dN); // Total number of time step with noise
const int firststep = ceil(log(nsigmareset/sigma)/dN);
const int itpstep = ceil((log(nsigmareset/sigma)-log(nsigmareset/sigma/2.))/dN);
constexpr double imax_double = LOG2*(NLpower-1) / dlogn;
constexpr int imax = int(imax_double) + (imax_double > int(imax_double));
const double inv_sqrt2 = 1./sqrt(2.);
const double sig2 = sigma*sigma;
const double sig3 = sigma*sigma*sigma;


#include "model.hpp"
#include "src/array_op.hpp"
#include "src/fft.hpp"
#include "src/util.hpp"
#include "src/progress.hpp"

const double Hinit = hubble({PHI_INIT,DPHI_INIT});
const double LL = 2.*M_PI/Hinit; // Box size L
const double dx = LL/NLnoise; // Spacing of each lattice

const std::string Nfileprefix = sdatadir + "/" + model + "/Nmap_";
const std::string fieldfileprefix = sdatadir + "/" + model + "/field_";
const std::string trajectoryfileprefix = sdatadir + "/" + model + "/trajectory_";
const std::string animationfileprefix = sdatadir + "/" + model + "/animation/animation_";
const std::string powfileprefix = sdatadir + "/" + model + "/power_";
const std::string cmpfileprefix = sdatadir + "/" + model + "/compaction_";
const std::string prbfileprefix = sdatadir + "/" + model + "/probabilities";
const std::string logwfileprefix = sdatadir + "/" + model + "/logw_";
const std::string mu2fileprefix = sdatadir + "/" + model + "/mu2_";
const std::string k3fileprefix = sdatadir + "/" + model + "/k3_";

bool Nfilefail, superH = false;
bool FFTwisdomFirst = false;

// int noisefiledirNo, noisefileNo;
std::ofstream Nfile, fieldfile, fieldfileA, trajectoryfile, powfile, cmpfile, prbfile, logwfile, Noisefile, mu2file, k3file;
std::array<double,NLnoiseAll> Ndata{};
std::array<double,NLnoiseAll> Nnoise{}; // use for EoN noise
std::array<double,NLnoiseAll> Ntotal{}; // use for averaging

std::array<state_type,NLnoiseAll> phievol{};
std::array<state_type,NLnoiseAll> PhidataAv{}; // use for averaging
std::array<state_type,NLnoiseAll> Phidata{}; // use for zoom


alignas(64) std::array<std::array<double,NLnoiseAll>,NFIELDS+1> biaslist{}; // FFT output: aligned like the FFTW plan arrays
alignas(64) std::array<std::array<double,NLnoiseAll>,NFIELDS+1+1> dwlist{}; // FFT output: aligned like the FFTW plan arrays
std::array<double,NLnoiseAll> laplacian{};
std::array<double,NLnoiseAll> mutwo{};
std::array<double,NLnoiseAll> kthree{};

// output
std::array<double,imax> disc_power{};
std::array<double,100*NLnoise> weightlist{};
std::array<double,100*NLnoise> weightbool{};
std::array<std::array<double,NLnoise/2>,2> zetar{};
std::array<double,NLnoise/2> dzetar{};

#if MODEL==1
  std::array<double,NLnoiseAll> N0list{};
  std::array<bool,NLnoiseAll> brokenlist{};
#elif MODEL==2
  std::array<double,NLnoiseAll> N1list{};
  std::array<double,NLnoiseAll> N2list{};
  std::array<bool,NLnoiseAll> broken1list{};
  std::array<bool,NLnoiseAll> broken2list{};
  const double dNsub = dN / divdN;
  const int ndiv  = (int)divdN;
#endif

#include "src/noise_bias.hpp"
#ifdef STOLAS_REF_NOISE
// regression test build (make REFNOISE=1): use the original noise generation
#include "tests/noise_bias_ref.hpp"
#define dwlist_gen dwlist_gen_ref
#define biaslist1D biaslist1D_ref
#endif
#include "src/zoom.hpp"
#include "src/laplacian.hpp"
#include "src/output.hpp"
#if MODEL==2
#include "src/simd_rk4.hpp"
#endif


// -- functions -----------------------
void initialize(){
  phii[0] = PHI_INIT;
  phii[1] = DPHI_INIT;

  #if MODEL==3
  NFLOOP{
    phii[2*nf] = PSI_INIT;
    phii[2*nf+1] = DPSI_INIT;
  }
  #endif

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  LOOP{
    phievol[NLnoise*NLnoise*i + NLnoise*j + k] = phii;
  }
}


void evolution(int seed, int starttime, int endtime, int InterpolatingNo) {
  double N = (starttime + InterpolatingNo*itpstep)*dN;
  int animationcount = (starttime==firststep ? firststep%aninum : 0); // for animation
  int animationstep = (starttime==firststep ? firststep/((double)aninum) : 0);

  Progress prog("Evolution " + std::to_string(starttime) + "->" + std::to_string(endtime), endtime - starttime);

  for (size_t n=starttime; n<endtime; n++){
    if constexpr (sbias) biaslist1D(n*dN);

    #if MODEL==3
      NFLOOP dwlist_gen(n*dN,seed,n,nf-1,InterpolatingNo);
    #else
      dwlist_gen(n*dN,seed,n,0,InterpolatingNo); // for phi
      dwlist_gen(n*dN,seed,n,1,InterpolatingNo); // for pi
    #endif

    if(snoisemap){
      if(n==totalstep-500) {
        Noisefile.open(sdatadir + "/" + model + "/noisedata/map_" + std::to_string(n) + std::string(".bin"), std::ios::binary);
        Noisefile << std::setprecision(10);
        Noisefile.write(reinterpret_cast<const char*>(&dwlist[0]), sizeof(double) * NLnoiseAll);
        Noisefile.close();
      }
    }

#if MODEL==2
    simd_rk4_drift_batch(phievol, driftedPhievol, ndiv, dNsub);
#endif

    double GaussianFactor = 1./dNbias/sqrt(2.*M_PI) * exp(-(N-Nbias)*(N-Nbias)/2./dNbias/dNbias);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int i=0; i<NLnoiseAll; i++){
      state_type phi = phievol[i];

      #if MODEL==1
        double calPphival = calPphi(N,phi,N0list[i],brokenlist[i]);
        double calPpival = calPpi(N,phi,N0list[i],brokenlist[i]);
        double crosscor = RecalPphipi(N,phi,N0list[i],brokenlist[i]);
        double phiamp = sqrt(calPphival);
      #elif MODEL==2
        double calPphival = calPphi(N,phi,N1list[i],N2list[i],broken1list[i],broken2list[i]);
        double calPpival = calPpi(N,phi,N1list[i],N2list[i],broken1list[i],broken2list[i]);
        double crosscor = RecalPphipi(N,phi,N1list[i],N2list[i],broken1list[i],broken2list[i]);
        double phiamp = sqrt(calPphival);
      #else
        double phiamp = sqrt(calPphi(phi));
      #endif

      double dw = dwlist[0][i];
      if(i==0 && recordWeight){
        weightlist[n] = dw; // save weight data
        weightbool[n] = true;
      }

      double Bias = biaslist[0][i];

      #if MODEL==2
        double phi_old = phi[0]; // for reflective boundary
        phi[0] = driftedPhievol[i][0];
        phi[1] = driftedPhievol[i][1];
      #else
        boost::numeric::odeint::runge_kutta4<state_type> stepper_noise;
        stepper_noise.do_step(dphidN, phi, N, dN);
      #endif

      #if MODEL==2 || MODEL==1
        double dwpi = dwlist[1][i];

        double b2 = crosscor*crosscor;
        double combi1 = sqrt(4.*b2 + pw2(calPphival-calPpival));

        double sqrtlam1 = sqrt(0.5*(calPphival+calPpival + combi1));
        double eig2 = 0.5*(calPphival+calPpival - combi1);
        if(eig2 < 1e-14*(calPphival+calPpival)) eig2 = 0.;

        double sqrtlam2 = sqrt(eig2);

        double denomplus = sqrt(4.*b2 + pw2(calPphival-calPpival + combi1));
        double denomminus = sqrt(4.*b2 + pw2(calPphival-calPpival - combi1));

        double vplus1 = (calPphival-calPpival + combi1)/denomplus;
        double vplus2 = 2.*crosscor/denomplus;
        double vminus1 = (calPphival-calPpival - combi1)/denomminus;
        double vminus2 = 2.*crosscor/denomminus;

        double biaseddw = dw + biasamp * Bias * GaussianFactor * sqrt_dN;

        phi[0] += (sqrtlam1*vplus1*biaseddw + sqrtlam2*vminus1*dwpi) * sqrt_dN;
        phi[1] += (sqrtlam1*vplus2*biaseddw + sqrtlam2*vminus2*dwpi) * sqrt_dN;

      #elif MODEL==3
        NFLOOP{
          phi[2*nf] += phiamp * sqrt_dN * dwlist[nf-1][i];
        }
      #else
        phi[0] += phiamp * dw * sqrt_dN;
        phi[0] += phiamp * biasamp * Bias * GaussianFactor * dN;
      #endif

      #if MODEL==2
        // reflective boundary
        if (phi_old <= phi1 && phi[0] > phi1) {
          phi[0] = 2.0 * phi1 - phi[0];
          phi[1] = -phi[1];
        }
        if (phi_old <= phi2 && phi[0] > phi2) {
          phi[0] = 2.0 * phi2 - phi[0];
          phi[1] = -phi[1];
        }
      #endif

      phievol[i] = phi;

      #if MODEL==1
        if (!brokenlist[i] && phi[0] < 0) {
          brokenlist[i] = true;
          N0list[i] = N;
        }
      #elif MODEL==2
        if (!broken1list[i] && phi[0] < phi1) {
          broken1list[i] = true;
          N1list[i] = N;
        }
        if (!broken2list[i] && phi[0] < phi2) {
          broken2list[i] = true;
          N2list[i] = N;
        }
      #endif

      if(strajectory && i==0){
        double Nd = N;
        save_trajectory(phi, Nd);
      }
    }

    animationcount++;
    if(sanimation && animationcount>aninum-1){
      animation(phievol, seed, animationstep, InterpolatingNo);
      animationstep++;
      animationcount=0;
    }

    N += dN;
    prog.set((long)n - starttime + 1);
  }
  prog.finish();
}


void evolutionNoise(int seed, int averagetime) {
  Progress prog("Noise evolution to EoN", NLnoiseAll);

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=0; i<NLnoiseAll; i++) {
    double N = totalstep*dN;
    std::seed_seq seq{seed,averagetime,i};
    std::mt19937 engine_av(seq);
    std::normal_distribution<> dist_av(0., 1.);

    state_type phi = phievol[i];

    // stepper
    boost::numeric::odeint::runge_kutta4<state_type> stepper_noise;

    #if MODEL==1
      // local copies: brokenlist/N0list must stay unchanged for the next average
      bool broken = brokenlist[i];
      double N0 = N0list[i];
    #endif

    while (EoN(phi)>0) {
      #if MODEL==3
        double psiamp = sqrt(calPpsi(phi));
        stepper_noise.do_step(dphidN, phi, N, dN);
        N += dN;
        NFLOOP{
          double dw = dist_av(engine_av);
          phi[2*nf] += psiamp * dw * sqrt_dN;
        }
      #elif MODEL==2
        double calPphival = calPphi(N,phi,N1list[i],N2list[i],broken1list[i],broken2list[i]);
        double calPpival = calPpi(N,phi,N1list[i],N2list[i],broken1list[i],broken2list[i]);
        double crosscor = RecalPphipi(N,phi,N1list[i],N2list[i],broken1list[i],broken2list[i]);
        double phiamp = sqrt(calPphival);
        
        double Nstep = N;
        for (int dn=0;dn<(int)divdN;dn++) {
          stepper_noise.do_step(dphidN, phi, Nstep, dN/divdN);
          Nstep += dN/divdN;
        }
        N += dN;
        
        double dw = dist_av(engine_av);
        double dwpi = dist_av(engine_av);

        double b2 = crosscor*crosscor;
        double combi1 = sqrt(4.*b2 + pw2(calPphival-calPpival));

        double sqrtlam1 = sqrt(0.5*(calPphival+calPpival + combi1));
        double eig2 = 0.5*(calPphival+calPpival - combi1);
        if(eig2 < 1e-14*(calPphival+calPpival)) eig2 = 0.;

        double sqrtlam2 = sqrt(eig2);

        double denomplus = sqrt(4.*b2 + pw2(calPphival-calPpival + combi1));
        double denomminus = sqrt(4.*b2 + pw2(calPphival-calPpival - combi1));

        double vplus1 = (calPphival-calPpival + combi1)/denomplus;
        double vplus2 = 2.*crosscor/denomplus;
        double vminus1 = (calPphival-calPpival - combi1)/denomminus;
        double vminus2 = 2.*crosscor/denomminus;

        phi[0] += (sqrtlam1*vplus1*dw + sqrtlam2*vminus1*dwpi) * sqrt_dN;
        phi[1] += (sqrtlam1*vplus2*dw + sqrtlam2*vminus2*dwpi) * sqrt_dN;
      #elif MODEL==1
        double calPphival = calPphi(N,phi,N0,broken);
        double calPpival = calPpi(N,phi,N0,broken);
        double crosscor = RecalPphipi(N,phi,N0,broken);

        double Nold = N;
        stepper_noise.do_step(dphidN, phi, N, dN);
        N += dN;

        double dw = dist_av(engine_av);
        double dwpi = dist_av(engine_av);

        double b2 = crosscor*crosscor;
        double combi1 = sqrt(4.*b2 + pw2(calPphival-calPpival));

        double sqrtlam1 = sqrt(0.5*(calPphival+calPpival + combi1));
        double eig2 = 0.5*(calPphival+calPpival - combi1);
        if(eig2 < 1e-14*(calPphival+calPpival)) eig2 = 0.;

        double sqrtlam2 = sqrt(eig2);

        double denomplus = sqrt(4.*b2 + pw2(calPphival-calPpival + combi1));
        double denomminus = sqrt(4.*b2 + pw2(calPphival-calPpival - combi1));

        double vplus1 = (calPphival-calPpival + combi1)/denomplus;
        double vplus2 = 2.*crosscor/denomplus;
        double vminus1 = (calPphival-calPpival - combi1)/denomminus;
        double vminus2 = 2.*crosscor/denomminus;

        phi[0] += (sqrtlam1*vplus1*dw + sqrtlam2*vminus1*dwpi) * sqrt_dN;
        phi[1] += (sqrtlam1*vplus2*dw + sqrtlam2*vminus2*dwpi) * sqrt_dN;

        if (!broken && phi[0] < 0) {
          broken = true;
          N0 = Nold;
        }
      #else
        double phiamp = sqrt(calPphi(phi));
        
        stepper_noise.do_step(dphidN, phi, N, dN);
        N += dN;
        
        double dw = dist(engine_av);
        phi[0] += phiamp * dw * sqrt_dN;
      #endif

      if(strajectory && i==0 && superH) {
        save_trajectory(phi, N+dN*totalstep);
      }
    }

    phievol[i] = phi;
    Nnoise[i] = N - dN*totalstep;
    prog.tick();
  }
  prog.finish();
}


void dNmap(int InterpolatingNo) {
  Progress prog("dN map", NLnoiseAll);

#ifdef _OPENMP
#pragma omp parallel for schedule(guided)
#endif
  for (int i=0; i<NLnoiseAll; i++) {
    int numstep = 0;
    int LatticePoint = i;
    double N = dN*(totalstep + InterpolatingNo*itpstep) + Nnoise[LatticePoint];
    state_type phi = phievol[LatticePoint];

    // Find zero crossing time
    typedef boost::numeric::odeint::runge_kutta_dopri5<state_type>base_stepper_type;
    auto stepper = make_dense_output(1.0e-15, 1.0e-15, base_stepper_type());

    stepper.initialize(phi, N, dN);
    state_type phil{};
    state_type phir{};
    state_type phim{};

    int countstep = 0;
    int caouttrj=0;
    while (true){
      stepper.do_step(dphidN);
      phi = stepper.current_state();
      N = stepper.current_time();

      caouttrj++;
      if(strajectory && i==0 && caouttrj>1000){
        double Nd = N;
        save_trajectory(phi, Nd);
        caouttrj=0;
      }
      
      if (EoI(phi)<0){
        double precphi = 1.e+2;
        double Nl = N;
        double Nr = N - stepper.current_time_step();
        double Nmid = 0;

        while (precphi>Nprec){
          stepper.calc_state(Nl, phil);
          stepper.calc_state(Nr, phir);
          Nmid = (Nl*EoI(phir) - Nr*EoI(phil)) / (EoI(phir)-EoI(phil));
          stepper.calc_state(Nmid, phim);
          if (EoI(phim)>0){
            Nr = Nmid;
          }
          else {
            Nl = Nmid;
          }
          precphi = std::fabs(EoI(phim));
        }

        N = Nmid;
        break;
      }
      countstep++;
      if (countstep>1e4) {
        std::cout << i << " N = " << N << std::endl;
        break;
      }

      if(strajectory && i==0 && superH) {
        save_trajectory(phi, N);
      }
    }

    Ndata[LatticePoint] = N;
    phievol[LatticePoint] = phi;
    prog.tick();
  }
  prog.finish();
}


bool checkNfilefail() {
  return Nfilefail;
}


#endif
