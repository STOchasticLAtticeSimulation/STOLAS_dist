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

#include <boost/numeric/odeint.hpp>

#define euler_gamma 0.57721566490153286061

#include "parameters.hpp"
#include "model.hpp"
#include "src/vec_op.hpp"
#include "src/fft.hpp"
#include "src/util.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif


#include <random>
std::mt19937 engine(1); // Fixed seed
// std::random_device seed; // random seed
// std::mt19937 engine(seed());
std::normal_distribution<> dist(0., 1.);


const std::string Nfileprefix = sdatadir + "/" + model + "/Nmap_";
const std::string fieldfileprefix = sdatadir + "/" + model + "/field_";
const std::string trajectoryfileprefix = sdatadir + "/" + model + "/trajectory_";
const std::string animationfileprefix = sdatadir + "/" + model + "/animation/animation_";
const std::string powfileprefix = sdatadir + "/" + model + "/power_";
const std::string powsfileprefix = sdatadir + "/" + model + "/powers";
bool noisefilefail, Nfilefail, sanisuperH = false;

int Nn, noisefiledirNo, noisefileNo;
std::ofstream Nfile, fieldfile, fieldfileA, trajectoryfile, powfile, powsfile;
std::vector<std::vector<double>> noisedataTr;
std::vector<std::vector<std::vector<double>>> noisedata(NFIELDS);
std::vector<double> Ndata;
std::vector<std::vector<std::vector<std::complex<double>>>> Nmap3D;

state_type phii(2*NFIELDS+2,0.);
std::vector<std::vector<double>> Phidata(NLnoise*NLnoise*NLnoise, std::vector<double>(2*NFIELDS+2,0.));
std::vector<std::vector<double>> phievol(NLnoise*NLnoise*NLnoise, std::vector<double>(2*NFIELDS+2,0.));
std::vector<std::vector<std::vector<double>>> phianimation;

std::vector<double> Nnoise(NLnoise*NLnoise*NLnoise, 0);
std::vector<double> Naverage(NLnoise*NLnoise*NLnoise, 0);
std::vector<double> Ntotal(NLnoise*NLnoise*NLnoise, 0);

#include "src/zoom.hpp"
#include "src/output.hpp"


// -- functions -----------------------
void initialize(){
  Ndata = std::vector<double>(NLnoise*NLnoise*NLnoise,0);
  phii[0] = PHI_INIT;
  phii[1] = DPHI_INIT;
  NFLOOP{
    phii[2*nf] = PSI_INIT;
    phii[2*nf+1] = DPSI_INIT;
  }

  LOOPLONG{
    phievol[NLnoise*NLnoise*i + NLnoise*j + k] = phii;
  }
}


void STOLAS(int NoisefileDirNo, int NoisefileNo, int InterpolatingNo) {
  std::ifstream noisefile;
  
  NFLOOP{
    noisefile.open(noisefiledir + std::to_string(NoisefileDirNo+nf-1) + noisefilenamediv + std::to_string(NoisefileNo) + "_" + std::to_string(InterpolatingNo) + std::string(".bin"), std::ios::binary);
    noisefilefail = noisefile.fail();
    
    if (!noisefile.fail()) {
      while (true) {
        std::vector<double> vv(NL);
        noisefile.read(reinterpret_cast<char*>(vv.data()), NL * sizeof(double));

        if (noisefile.gcount() != NL * sizeof(double)) {
          break;
        }
        noisedataTr.push_back(vv);
      }
      noisedata[nf-1] = transpose(noisedataTr);
      noisedataTr.clear();
      noisefile.close();
    }
  }

}


void dNmap(int NoiseNo, int InterpolatingNo) {
  
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=0; i<NL; i++) {
    int numstep = 0;
    double N = dN*(totalstep + InterpolatingNo*itpstep) + Nnoise[i + NL*NoiseNo];
    state_type phi = phievol[i + NL*NoiseNo];

    // Find zero crossing time
    typedef boost::numeric::odeint::runge_kutta_dopri5<state_type>base_stepper_type;
    auto stepper = make_dense_output(1.0e-15, 1.0e-15, base_stepper_type());

    stepper.initialize(phi, N, dN);
    state_type phil(2*NFIELDS+2);
    state_type phir(2*NFIELDS+2);
    state_type phim(2*NFIELDS+2);

    while (true){
      stepper.do_step(dphidN);
      phi = stepper.current_state();
      N = stepper.current_time();
      
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

      if(strajectory && i==0 && NoiseNo==0 && sanisuperH) {
        save_trajectory(phi, N);
      }
    }

    Ndata[i + NL*NoiseNo] = N;
    phievol[i + NL*NoiseNo] = phi;
  }
}


void evolution(int NoiseNo) {
  
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=0; i<NL; i++) {
    double N = 0;
    int numstep = 0;
    int animationcount = 0; // for animation
    state_type phi = phievol[i + NL*NoiseNo];

    // stepper
    boost::numeric::odeint::runge_kutta4<state_type> stepper_noise;

    for (size_t n=0; n<noisedata[0][0].size(); n++) {
      double psiamp = sqrt(calPpsi(phi));
      stepper_noise.do_step(dphidN, phi, N, dN);
      N += dN;
      double dw = 0.;
      NFLOOP{
        dw = noisedata[nf-1][i][n];// dist(engine);// 
        phi[2*nf] += psiamp * dw * sqrt(dN);
      }

      if(strajectory && i==0 && NoiseNo==0){
        double Nd = N;
        if(sanisuperH) Nd += dN*firststep;
        save_trajectory(phi, Nd);
      }

      animationcount++;
      if(sanimation && animationcount>aninum-1){
        int ef = int(n/aninum);
        if(sanisuperH) ef += int(firststep/aninum);
        phianimation[ef][i + NL*NoiseNo] = phi;
        animationcount=0;
      }
    }

    phievol[i + NL*NoiseNo] = phi;
  }
}


void evolutionNoise(int NoiseNo) {
  
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=0; i<NL; i++) {
    double N = 0;
    double sqrt_dN = std::sqrt(dN);
    state_type phi = phievol[i + NL*NoiseNo];

    // stepper
    boost::numeric::odeint::runge_kutta4<state_type> stepper_noise;

    while (EoN(phi)>0) {
      double psiamp = sqrt(calPpsi(phi));
      stepper_noise.do_step(dphidN, phi, N, dN);
      N += dN;
      NFLOOP{
        double dw = dist(engine);
        phi[2*nf] += psiamp * dw * sqrt_dN;
      }

      if(strajectory && i==0 && NoiseNo==0 && sanisuperH) {
        save_trajectory(phi, N+dN*totalstep);
      }
    }

    phievol[i + NL*NoiseNo] = phi;
    Nnoise[i + NL*NoiseNo] = N;
  }
}


bool checknoisefile() {
  return noisefilefail;
}


bool checkNfilefail() {
  return Nfilefail;
}


#endif
