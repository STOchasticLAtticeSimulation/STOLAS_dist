#include "../source/STOLAS.hpp"

//--------- User may change ---------
// Potential
double STOLAS::VV(double phi) {
  return mm*mm*phi*phi/2.;
}

// Derivative
double STOLAS::Vp(double phi) {
  return mm*mm*phi;
}

// Power spectrum of phi
double STOLAS::calPphi(std::vector<double> &phi) {
  double HH = hubble(phi[0],phi[1]);
  double eps = Vp(phi[0])*Vp(phi[0])/VV(phi[0])/VV(phi[0])/2.;
  double Hi = hubble(phii[0],phii[1]);
  double NoiseNLO = (1. + eps*(6.-4.*euler_gamma-8.*log(2.))) * pow(sigma*Hi/2./HH, -4.*eps);

  return pow(HH/2./M_PI, 2) * NoiseNLO;
}
//-----------------------------------


int main(int argc, char* argv[])
{
  if (argc!=2) {
    std::cout << "Specify the noise file number correctly." << std::endl;
    return -1;
  }
  
  // ---------- start timer ----------
  struct timeval Nv;
  struct timezone Nz;
  double before, after;
  
  gettimeofday(&Nv, &Nz);
  before = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  // --------------------------------------

  int noisefileNo = atoi(argv[1]);
  
  STOLAS stolas(model,dN,sourcedir,noisefileNo,phii,bias,Nbias,dNbias);

  if (!stolas.checknoisefile()) {
    std::cout << "The noise file couldn't be opened." << std::endl;
    return -1;
  }

  if (!stolas.checkbiasfile()) {
    std::cout << "The bias file couldn't be opened." << std::endl;
    return -1;
  }

  if (!stolas.noisebiassize()) {
    std::cout << "The box sizes of the noise and the bias are inconsistent." << std::endl;
    return -1;
  }

  if (stolas.Nfilefail()) {
    std::cout << "The export file couldn't be opened. 'mkdir data'" << std::endl;
    return -1;
  }

  if(szeta) stolas.dNmap();
  if(spower) stolas.spectrum();
  if(sanimation) stolas.animation();
  if(scompaction) stolas.compaction();

  // ---------- stop timer ----------
  gettimeofday(&Nv, &Nz);
  after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  std::cout << after - before << " sec." << std::endl;
  std::cout << std::endl;
  // -------------------------------------
}
