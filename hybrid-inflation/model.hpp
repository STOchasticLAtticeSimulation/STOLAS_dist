#ifndef INCLUDED_model_hpp_
#define INCLUDED_model_hpp_

#define NFIELDS 1 // Number of waterfall
#define NFLOOP for (int nf = 1; nf < NFIELDS+1; nf++)

extern inline double pw2(double x) {return x*x;}
extern inline double pw3(double x) {return x*x*x;}
extern inline double pw4(double x) {return x*x*x*x;}

// type
typedef std::vector<double> state_type;

// Model parameters
const std::string model = "test"; // Name of the model
const bool CUBIC = false;
const double Pzeta = 2.1e-9;
const double Mm = 1.e-2;
const double pc = sqrt(2.)*Mm;
const double Pi = 10.;//11.;//
const double mu1 = pw2(Pi/Mm)/pc; //7.07e+7;//
const double mu2 = 10.7;//6.8;//
const double mu3 = 0.0375;
const double lam = pow(Pzeta*12.*pw2(M_PI/mu1),0.25);

// Initial value of fields
const double PHI_INIT = pc + 3./mu1;
const double DPHI_INIT = 0.;
const double PSI_INIT = 0.;//sqrt(pw4(lam)*Pi/48./sqrt(2.*M_PI*M_PI*M_PI));
const double DPSI_INIT = 0.;

const double EoIvalue = -2.0;
const double EoNvalue = -1.2;

// Potential
double VV(const state_type &phi) {
  double PHI = phi[0];
  double phiterm = (PHI-pc)/mu1 - pw2(PHI-pc)/pw2(mu2);
  if(CUBIC) phiterm += pw3(PHI-pc)/pw3(mu3);

  double wfterm = 0;
  NFLOOP{
    wfterm += pw2(phi[2*nf]);
  }
  return pw4(lam)*(phiterm + pw2(1-wfterm/pw2(Mm)) + 2.*pw2(PHI/Mm/pc)*wfterm);
}

// Derivative
double Vphi(const state_type &phi) {
  double PHI = phi[0];
  double phiterm = 1./mu1 - 2.*(PHI-pc)/pw2(mu2);
  if(CUBIC) phiterm += 3.*pw2(PHI-pc)/pw3(mu3);

  double wfterm = 0;
  NFLOOP{
    wfterm += pw2(phi[2*nf]);
  }
  return pw4(lam)*(phiterm + 4.*PHI*pw2(1./Mm/pc)*wfterm);
}

double Vpsi(const state_type &phi, int NField) {
  double phiterm = pw2(Mm)*(pw2(phi[0]/pc)-1);
  double wfterm = 0;
  NFLOOP{
    wfterm += pw2(phi[2*nf]);
  }
	return 4.*pw4(lam/Mm) * phi[2*NField] * (phiterm + wfterm);
}

double hubble(const state_type &phi) {
  double wfkterm = 0;
  NFLOOP{
    wfkterm += pw2(phi[2*nf+1])/2.;
  }
  return sqrt((phi[1]*phi[1]/2. + wfkterm + VV(phi))/3.);
}

double ep(const state_type &phi) {
  double HH = hubble(phi);
  return phi[1]*phi[1]/2./HH/HH;
}

// pow spectrum of fields
double calPphi(const state_type &phi) {
  double HH = hubble(phi);
  return pw2(HH/2./M_PI);
}

double calPpsi(const state_type &phi) {
  double HH = hubble(phi);
  return pw2(HH/2./M_PI);
}

// The condition at the end of inflation
double EoI(const state_type &phi) {
  double phiterm = pw2(Mm)*(pw2(phi[0]/pc)-1);
  double wfterm = 0;
  NFLOOP{
    wfterm += pw2(phi[2*nf]);
  }
  return -EoIvalue+4.*pw4(lam/Mm) * (phiterm + 3.*wfterm)/VV(phi);
}

double EoN(const state_type &phi) {
  double phiterm = pw2(Mm)*(pw2(phi[0]/pc)-1);
  double wfterm = 0;
  NFLOOP{
    wfterm += pw2(phi[2*nf]);
  }
  return -EoNvalue+4.*pw4(lam/Mm) * (phiterm + 3.*wfterm)/VV(phi);
}

// EoM
void dphidN(const state_type &phi, state_type &dxdt, const double t) {
  double p1 = phi[1]; // pi
  double HH = hubble(phi);

  dxdt[0] = p1/HH;
  dxdt[1] = -3.*p1 - Vphi(phi)/HH;
  NFLOOP{
    double p2 = phi[2*nf+1]; // pi
    dxdt[2*nf] = p2/HH;
    dxdt[2*nf+1] = -3.*p2 - Vpsi(phi,nf)/HH;
  }
  
}

#endif