#ifndef INCLUDED_model_hpp_
#define INCLUDED_model_hpp_

extern inline double pw2(double x) {return x*x;}
extern inline double pw3(double x) {return x*x*x;}
extern inline double pw4(double x) {return x*x*x*x;}


// The model can be chosen here (default) or from the command line: `make MODEL=2`
#ifndef MODEL
#define MODEL 0 // 0 -> chaotic, 1 -> Starobinsky, 2 -> USR, 3 -> hybrid
#endif

// Number of waterfall fields (hybrid only), or `make MODEL=3 NFIELDS=2`
#ifndef NFIELDS
#if MODEL==3
#define NFIELDS 1
#else
#define NFIELDS 0
#endif
#endif

#define NFLOOP for (int nf = 1; nf < NFIELDS+1; nf++)
constexpr int NumFields = 2*NFIELDS+2;

// type
typedef std::array<double,NumFields> state_type;
state_type phii{};

// input from the environment variable SIGMA, e.g. `SIGMA=100. ./STOLAS 0`
// (only MODEL==2 uses it as the amplitude calPzeta = SIGMA*1e-5; for the other
//  models it just labels the output file names and defaults to 0)
inline double read_env_calPzeta0() {
  const char* s = std::getenv("SIGMA");
  if (s == nullptr) {
#if MODEL==2
    std::cerr << "Environment variable SIGMA is not set. Usage: SIGMA=<value> ./STOLAS <noise seed>" << std::endl;
    std::exit(1);
#else
    return 0.;
#endif
  }
  return std::stod(s);
}
double calPzeta0 = read_env_calPzeta0();

#if MODEL==0

// Model parameters
const std::string model = "chaotic"; // Name of the model
const double mm = 0.0211; // Mass of the inflaton
const double PHI_INIT = 11.;
const double DPHI_INIT = -sqrt(2./3)*mm;

// Potential
double VV(const state_type &phi) {
  double PHI = phi[0];
  return mm*mm*phi[0]*phi[0]/2.;
}

// Derivative
double Vphi(const state_type &phi) {
  return mm*mm*phi[0];
}

double hubble(const state_type &phi) {
  return sqrt((phi[1]*phi[1]/2. + VV(phi))/3.);
}

double ep(const state_type &phi) {
  double HH = hubble(phi);
  return phi[1]*phi[1]/2./HH/HH;
}

// Power spectrum of phi
double calPphi(const state_type &phi) {
  double HH = hubble(phi);
  double eps = Vphi(phi)*Vphi(phi)/VV(phi)/VV(phi)/2.;
  double Hi = hubble(phii);
  double NoiseNLO = (1. + eps*(6.-4.*euler_gamma-8.*log(2.))) * pow(sigma*Hi/2./HH, -4.*eps);

  return pow(HH/2./M_PI, 2) * NoiseNLO;
}

// The condition at the end of inflation
inline double EoI(const state_type &phi) {
  return 1-ep(phi);
}

inline double EoN(const state_type &phi) {
  return 0.6-ep(phi);
}

#elif MODEL==1

// Model parameters
const std::string model = "Starobinsky"; // Name of the model
const double H0 = 1e-5; // Hubble parameter of broken point
const double calPRIR = 8.5e-10; // Amplitude of curvature perturbation
const double Lambda = 1080; // Ratio between Ap to Am
const double Ap = sqrt(9./4/M_PI/M_PI*H0*H0*H0*H0*H0*H0/calPRIR); // Gradient of the potential at first stage
const double Am = Ap/Lambda; // Gradient of the potential at second stage
const double V0 = 3*H0*H0; // Amplitude of flat potential
const double phif = -0.0187; // The inflaton value at the end of inflation
const double PHI_INIT = 0.0835;//0.0193;//
const double DPHI_INIT = -5.45e-7;

// Potential
double VV(const state_type &phi) {
  if (phi[0] > 0) {
    return V0 + Ap*phi[0];
  } else {
    return V0 + Am*phi[0];
  }
}

// Derivative
double Vphi(const state_type &phi){ 
  if (phi[0] > 0) {
    return Ap;
  } else {
    return Am;
  }
}

double hubble(const state_type &phi) {
  return sqrt((phi[1]*phi[1]/2. + VV(phi))/3.);
}


// Power spectrum of phi
double calPphi(double &N, const state_type &phi, double N0, bool broken) {
  double alpha = exp(N-N0);
  double alp2 = alpha*alpha;
  double alp3 = alpha*alpha*alpha;

  if (!broken) {
    return pw2(hubble(phi)/M_PI)*(1 + sig2)/4.;
  }
  else {
    return (pw2(hubble(phi)/(M_PI*alp3*Lambda*sig3))*(3*(Lambda - 1.)*cos(2*(-1 + alpha)*sigma)*
        (3*(1 + (-1 + 4*alpha)*sig2 - (-4 + alpha)*alp3*sig2*sig2 + alp2*alp2*sig3*sig3) + 
          Lambda*(-3 + (3 - 12*alpha)*sig2 + (-16 + 7*alpha)*alp3*sig2*sig2 + (-7 + 4*alpha)*alp2*alp2*sig3*sig3)) + 
       (1 + sig2)*(pw2(Lambda)*(9 + 18*alp2*sig2 + 9*alp2*alp2*sig2*sig2 + 2*alp3*alp3*sig3*sig3) - 
          18*Lambda*(1 + alp2*sig2)*(1 + alp2*sig2) + 9*(1 + alp2*sig2)*(1 + alp2*sig2)) - 
       6*(Lambda - 1)*sigma*(-3*(-1 + alpha)*(1 + alpha*sig2 + alp2*sig2 + alp3*sig2*sig2) + 
          Lambda*(-3 - 3*alpha*(-1 + sig2) - 4*alp3*(-1 + sig2)*sig2 + 7*alp2*alp2*sig2*sig2 + 
             alp2*alp3*(-1 + sig2)*sig2*sig2))*sin(2*(-1. + alpha)*sigma)))/8.;
  }
}

// The power spectrum of pi
double calPpi(double &N, const state_type &phi, double N0, bool broken) {
  double alpha = exp(N-N0);
  double alp2 = alpha*alpha;
  double alp3 = alpha*alpha*alpha;

  if (!broken) {
    return pw4(hubble(phi)*sigma)/pw2(2.*M_PI);
  }
  else {
    return (pw4(hubble(phi))/pw2(M_PI*Lambda*sigma*alp3)*(-18*Lambda + 9*pw2(Lambda) + 9 - 36*Lambda*alp2*sig2 + 
       18*alp2*pw2(Lambda)*sig2 + 18*alp2*sig2 - 18*Lambda*alp2*alp2*sig2*sig2 + 9*alp2*alp2*pw2(Lambda)*sig2*sig2 + 
       9*alp2*alp2*sig2*sig2 + 3*(Lambda - 1)*cos(2*(-1 + alpha)*sigma)*(-3*Lambda + 3 + 7*Lambda*alp2*alp2*sig2*sig2 - 3*alp2*alp2*sig2*sig2) + 
       2*alp3*alp3*pw2(Lambda)*sig3*sig3 + 6*alpha*(Lambda - 1)*sigma*
        (3*(1 + alp2*sig2) + Lambda*(-3 - 4*alp2*sig2 + alp2*alp2*sig2*sig2))*sin(2*(-1 + alpha)*sigma)))/8.;
  }
}

// Cross correlation of phi and pi
double RecalPphipi(double &N, const state_type &phi, double N0, bool broken) {
  double alpha = exp(N-N0);
  double alp2 = alpha*alpha;
  double alp3 = alpha*alpha*alpha;
  double prefactor = ((pw3(hubble(phi)*sigma))/(2.*M_PI*M_PI));
  
  if (!broken) {
    return -prefactor*(1./(2.*sigma));
  }
  else {
    return -prefactor*0.25*(1/pw2(alp3*Lambda*sig3)/sigma*(9*pw2(Lambda - 1) + 18*alp2*pw2(Lambda - 1)*sig2 + 9*alp2*alp2*pw2(Lambda - 1)*sig2*sig2 + 
       2*alp3*alp3*pw2(Lambda)*pw2(sig3) + 3*(Lambda - 1)*(cos(2*(-1 + alpha)*sigma)*
           (-3*Lambda + 3 + 6*alpha*(-Lambda + 1)*sig2 + ((-8 + 7*alpha)*Lambda - 3*(-2 + alpha))*alp3*sig2*sig2 + 2*Lambda*alp3*alp2*pw2(sig3)) + 
          sigma*(-3*(-1 + 2*alpha)*(Lambda - 1) + 2*(-4*Lambda + 3)*alp3*sig2 + ((-7 + 2*alpha)*Lambda + 3)*alp2*alp2*sig2*sig2)*
           sin(2*(-1 + alpha)*sigma))));
  }
}

// // Power spectrum of phi
// double calPphi(double &N, const state_type &phi, double N0, bool broken) {
//   if (!broken) {
//     return pow(hubble(phi)/2./M_PI,2);
//   } else {
//     double alpha = exp(N-N0);
//     double alp2 = alpha*alpha;
//     double sig2 = sigma*sigma;
//     double alp3 = alpha*alpha*alpha;
//     double sig3 = sigma*sigma*sigma;

//     return pow(hubble(phi)/2./M_PI,2) *
//       ((1 + sig2)*(9*pow(1 + alp2*sig2,2) - 
// 				18*Lambda*pow(1 + alp2*sig2,2) + 
// 				pow(Lambda,2)*(9 + 18*alp2*sig2 + 9*alp2*alp2*sig2*sig2 + 
// 					       2*alp3*alp3*sig3*sig3)) + 
// 	    3*(-3*(1 + (-1 + 4*alpha)*sig2 - (-4 + alpha)*alp3*sig2*sig2 + 
// 		   alp2*alp2*sig3*sig3) + Lambda*
// 	       (6 + 6*(-1 + 4*alpha)*sig2 + 2*(14 - 5*alpha)*alp3*sig2*sig2 + 
// 		2*(5 - 2*alpha)*alp2*alp2*sig3*sig3) + 
// 	       pow(Lambda,2)*(-3 + (3 - 12*alpha)*sig2 + alp3*(-16 + 7*alpha)*sig2*sig2 + 
// 			      alp2*alp2*(-7 + 4*alpha)*sig3*sig3))*cos(2*(-1 + alpha)*sigma) + 
// 	    6*sigma*(-3*pow(-1 + Lambda,2) + alp2*alp2*(3 - 10*Lambda + 7*pow(Lambda,2))*sig2*sig2 - 
// 		     3*alpha*pow(-1 + Lambda,2)*(-1 + sig2) - 
// 		     alp3*(3 - 7*Lambda + 4*pow(Lambda,2))*sig2*(-1 + sig2) + 
// 		     alp2*alp3*(-1 + Lambda)*Lambda*sig2*sig2*(-1 + sig2))*sin(2*sigma - 2*alpha*sigma))/
//       (2.*alp3*alp3*pow(Lambda,2)*sig3*sig3);
//   }
// }

// // The power spectrum of pi
// double calPpi(double &N, const state_type &phi, double N0, bool broken) {
//   return 0;
// }

// // Cross correlation of phi and pi
// double RecalPphipi(double &N, const state_type &phi, double N0, bool broken) {
//   return 1;
// }

// The condition at the end of inflation
inline double EoI(const state_type &phi) {
  return phi[0] - phif;
}

inline double EoN(const state_type &phi) {
  return phi[0] - phif;
}


#elif MODEL==2

// Model parameters
const std::string model = "USR"; // Name of the model
const double H0 = 1e-5; // Hubble parameter of broken point
const double calPRIR = 8.5e-10; // Amplitude of curvature perturbation
const double B1 = sqrt(9./4/M_PI/M_PI*H0*H0*H0*H0*H0*H0/calPRIR); // Gradient of the potential at first stage
const double B2 = B1/5000;// 0; // Gradient of the potential at second stage
const double B3 = B1; // Gradient of the potential at third stage
const double V0 = 3.*H0*H0; // Amplitude of flat potential
const double phif = -0.3; // The inflaton value at the end of inflation
const double phiN = phif+0.05; // The inflaton value at the end of inflation
const double PHI_INIT = 0.0826;
const double DPHI_INIT = -5.45e-7;
const double calPzeta = calPzeta0 * 1.e-5;//22.8e-5;//50.e-4;//
const double USRrange = H0/(6.*M_PI*sqrt(calPzeta)) - B1/pw2(3.*H0);//-0.0181435;//-0.0180356;
double divdN = 10.;
double Nprotect = 20.;

const double phi1 = 0.;
const double phi2 = phi1 + USRrange;


// Potential
double VV(const state_type &phi) {
  if (phi[0] > phi1) {
    return V0 + B1*(phi[0] - phi1) - B2*phi2;
  }
  else if (phi[0] > phi2 && phi[0] <= phi1) {
    return V0 + B2*(phi[0] - phi1 - phi2);
  }
  else {
    return V0 + B3*(phi[0] - phi2) - B2*phi1;
  }
}

// Derivative
double Vphi(const state_type &phi) {
  if (phi[0] > phi1) {
    return B1;
  }
  else if (phi[0] > phi2 && phi[0] <= phi1) {
    return B2;
  }
  else {
    return B3;
  }
}

double hubble(const state_type &phi) {
  return sqrt((phi[1]*phi[1]/2. + VV(phi))/3.);
}


// Power spectrum of phi
double calPphi(double &N, const state_type &phi, double N1, double N2, bool broken1, bool broken2) {
  double alpha = exp(N-N1);
  double beta = exp(N-N2);
  double alp2 = alpha*alpha;
  double bet2 = beta*beta;
  double alp3 = alpha*alpha*alpha;
  double bet3 = beta*beta*beta;

  if (N - N1 > Nprotect) {
    return pw2(hubble(phi)/M_PI)*(1. + sig2)/4.;
  }

  if (!broken1&&!broken2) {
    return pw2(hubble(phi)/M_PI)*(1 + sig2)/4.;
  }
  else if(broken1&&!broken2) {
    return (pw2(hubble(phi)/(M_PI*alp3*B1*sig3))*(3*(B1 - B2)*cos(2*(-1 + alpha)*sigma)*
        (3*B2*(1 + (-1 + 4*alpha)*sig2 - (-4 + alpha)*alp3*sig2*sig2 + alp2*alp2*sig3*sig3) + 
          B1*(-3 + (3 - 12*alpha)*sig2 + (-16 + 7*alpha)*alp3*sig2*sig2 + (-7 + 4*alpha)*alp2*alp2*sig3*sig3)) + 
       (1 + sig2)*(pw2(B1)*(9 + 18*alp2*sig2 + 9*alp2*alp2*sig2*sig2 + 2*alp3*alp3*sig3*sig3) - 
          18*B1*B2*(1 + alp2*sig2)*(1 + alp2*sig2) + 9*(B2 + B2*alp2*sig2)*(B2 + B2*alp2*sig2)) - 
       6*(B1 - B2)*sigma*(-3*(-1 + alpha)*B2*(1 + alpha*sig2 + alp2*sig2 + alp3*sig2*sig2) + 
          B1*(-3 - 3*alpha*(-1 + sig2) - 4*alp3*(-1 + sig2)*sig2 + 7*alp2*alp2*sig2*sig2 + 
             alp2*alp3*(-1 + sig2)*sig2*sig2))*sin(2*(-1. + alpha)*sigma)))/8.;
  }
  else {
    return (1./(M_PI*alp3*B1*bet3*sig3*sig3*M_PI*alp3*B1*bet3*sig3*sig3)*
     std::norm(hubble(phi)*(-((1. - II*sigma)*exp(2.*II*sigma)*
             ((3*(B2 + B2*alp2*sig2) + B1*(-3 - 3*alp2*sig2 - 2.*II*alp3*sig3))*
                (-3*B3*alp3*(1 + bet2*sig2) + 2.*II*B1*bet3*bet3*sig3 - 2.*II*B2*bet3*bet3*sig3 + 
                  B2*alp3*(3 + 3*bet2*sig2 + 2.*II*bet3*sig3)) + 
               9*(B1 - B2)*(B2 - B3)*alp3*exp(2.*II*(alpha - beta)*sigma)*(II + alpha*sigma)*(II + alpha*sigma)*(-II + beta*sigma)*(-II + beta*sigma))) + 
          3.*II*(1. + II*sigma)*exp(2.*II*beta*sigma)*
           ((B1 - B2)*exp(2.*II*(alpha - beta)*sigma)*(-3.*II*B3*alp3*(1 + bet2*sig2) + 2*B1*bet3*bet3*sig3 - 
                2*B2*bet3*bet3*sig3 + B2*alp3*(3.*II + 3.*II*bet2*sig2 + 2*bet3*sig3))*
              (II + alpha*sigma)*(II + alpha*sigma) + (B2 - B3)*alp3*
              (3.*II*(-B1 + B2) + 3.*II*(-B1 + B2)*alp2*sig2 + 2*B1*alp3*sig3)*(II + beta*sigma)*(II + beta*sigma)))*
        1./(B2*alp3 + B1*bet3 - B2*bet3)))/64.;
  }
}

// The power spectrum of pi
double calPpi(double &N, const state_type &phi, double N1, double N2, bool broken1, bool broken2) {
  double alpha = exp(N-N1);
  double beta = exp(N-N2);
  double alp2 = alpha*alpha;
  double bet2 = beta*beta;
  double alp3 = alpha*alpha*alpha;
  double bet3 = beta*beta*beta;

  if (N - N1 > Nprotect) {
    return pw4(hubble(phi)*sigma)/pw2(2.*M_PI);
  }

  if (!broken1&&!broken2) {
    return pw4(hubble(phi)*sigma)/pw2(2.*M_PI);
  }
  else if(broken1&&!broken2) {
    return (pw4(hubble(phi))/pw2(M_PI*B1*sigma*alp3)*(-18*B1*B2 + 9*pw2(B1) + 9*pw2(B2) - 36*B1*B2*alp2*sig2 + 
       18*alp2*pw2(B1)*sig2 + 18*alp2*pw2(B2)*sig2 - 18*B1*B2*alp2*alp2*sig2*sig2 + 9*alp2*alp2*pw2(B1)*sig2*sig2 + 
       9*alp2*alp2*pw2(B2)*sig2*sig2 + 3*(B1 - B2)*cos(2*(-1 + alpha)*sigma)*(-3*B1 + 3*B2 + 7*B1*alp2*alp2*sig2*sig2 - 3*B2*alp2*alp2*sig2*sig2) + 
       2*alp3*alp3*pw2(B1)*sig3*sig3 + 6*alpha*(B1 - B2)*sigma*
        (3*(B2 + B2*alp2*sig2) + B1*(-3 - 4*alp2*sig2 + alp2*alp2*sig2*sig2))*sin(2*(-1 + alpha)*sigma)))/8.;
  }
  else {
    return (pw4(hubble(phi))/pw2(M_PI*B1*alp3*bet3*sig2*sig2)*
     std::norm((exp(2.*II*sigma)*(-3*(B2 + B2*alp2*sig2) + B1*(3. + 3*alp2*sig2 + 2.*II*alp3*sig3))*
           (3.*II*B3*alp3*(1 + bet2*sig2) + 2*B1*bet3*bet3*sig3 - 2*B2*bet3*bet3*sig3 + 
             B2*alp3*(-3.*II - 3.*II*bet2*sig2 + 2*bet3*sig3)) + 
          3*(B1 - B2)*exp(2.*II*alpha*sigma)*(-3.*II*B3*alp3*(1 + bet2*sig2) + 2*B1*bet3*bet3*sig3 - 
             2*B2*bet3*bet3*sig3 + B2*alp3*(3.*II + 3.*II*bet2*sig2 + 2*bet3*sig3))*
           (II + alpha*sigma)*(II + alpha*sigma) + 9.*II*(B1 - B2)*(B2 - B3)*alp3*exp(2.*II*(1 + alpha - beta)*sigma)*
           (II + alpha*sigma)*(II + alpha*sigma)*(-II + beta*sigma)*(-II + beta*sigma) + 
          3*(B2 - B3)*alp3*exp(2.*II*beta*sigma)*(3.*II*(-B1 + B2) + 3.*II*(-B1 + B2)*alp2*sig2 + 
             2*B1*alp3*sig3)*(II + beta*sigma)*(II + beta*sigma))/(B2*(alp3 - bet3) + B1*bet3)))/64.;
  }
}

// Cross correlation of phi and pi
double RecalPphipi(double &N, const state_type &phi, double N1, double N2, bool broken1, bool broken2) {
  double alpha = exp(N-N1);
  double beta = exp(N-N2);
  double alp2 = alpha*alpha;
  double bet2 = beta*beta;
  double alp3 = alpha*alpha*alpha;
  double bet3 = beta*beta*beta;
  double prefactor = ((pw3(hubble(phi)*sigma))/(2.*M_PI*M_PI));

  if (N - N1 > Nprotect) {
    return -prefactor*(1./(2.*sigma));
  }
  
  if (!broken1&&!broken2) {
    return -prefactor*(1./(2.*sigma));
  }
  else if(broken1&&!broken2) {
    return -prefactor*0.25*(1/pw2(alp3*B1*sig3)/sigma*(9*pw2(B1 - B2) + 18*alp2*pw2(B1 - B2)*sig2 + 9*alp2*alp2*pw2(B1 - B2)*sig2*sig2 + 
       2*alp3*alp3*pw2(B1)*pw2(sig3) + 3*(B1 - B2)*(cos(2*(-1 + alpha)*sigma)*
           (-3*B1 + 3*B2 + 6*alpha*(-B1 + B2)*sig2 + ((-8 + 7*alpha)*B1 - 3*(-2 + alpha)*B2)*alp3*sig2*sig2 + 2*B1*alp3*alp2*pw2(sig3)) + 
          sigma*(-3*(-1 + 2*alpha)*(B1 - B2) + 2*(-4*B1 + 3*B2)*alp3*sig2 + ((-7 + 2*alpha)*B1 + 3*B2)*alp2*alp2*sig2*sig2)*
           sin(2*(-1 + alpha)*sigma))));
  }
  else {
    return -prefactor*0.03125*((exp(- 2.*II*sigma)*((-3*(B2 + B2*alp2*sig2) + B1*(3 + 3*alp2*sig2 - 2.*II*alp3*sig3))*
          (-3.*II*B3*alp3*(1 + bet2*sig2) + 2*B1*bet3*bet3*sig3 - 2*B2*bet3*bet3*sig3 + 
            B2*alp3*(3.*II + 3.*II*bet2*sig2 + 2*bet3*sig3)) + 
         3*(B1 - B2)*exp(-2.*II*(-1 + alpha)*sigma)*(3.*II*B3*alp3*(1 + bet2*sig2) + 2*B1*bet3*bet3*sig3 - 
            2*B2*bet3*bet3*sig3 + B2*alp3*(-3.*II - 3.*II*bet2*sig2 + 2*bet3*sig3))*
          (-II + alpha*sigma)*(-II + alpha*sigma) + 3*(B2 - B3)*alp3*exp(-2.*II*(-1 + beta)*sigma)*
          (-3.*II*(B2 + B2*alp2*sig2) + B1*(3.*II + 3.*II*alp2*sig2 + 2*alp3*sig3))*
          (-II + beta*sigma)*(-II + beta*sigma) - 9.*II*(B1 - B2)*(B2 - B3)*alp3*exp(-2.*II*(alpha - beta)*sigma)*(-II + alpha*sigma)*(-II + alpha*sigma)*
          (II + beta*sigma)*(II + beta*sigma))*(-((1. - II*sigma)*exp(2.*II*sigma)*
            ((3*(B2 + B2*alp2*sig2) + B1*(-3 - 3*alp2*sig2 - 2.*II*alp3*sig3))*
               (-3*B3*alp3*(1 + bet2*sig2) + 2.*II*B1*bet3*bet3*sig3 - 2.*II*B2*bet3*bet3*sig3 + 
                 B2*alp3*(3 + 3*bet2*sig2 + 2.*II*bet3*sig3)) + 
              9*(B1 - B2)*(B2 - B3)*alp3*exp(2.*II*(alpha - beta)*sigma)*(II + alpha*sigma)*(II + alpha*sigma)*(-II + beta*sigma)*(-II + beta*sigma))) + 
         3.*II*(1. + II*sigma)*exp(2.*II*beta*sigma)*
          ((B1 - B2)*exp(2.*II*(alpha - beta)*sigma)*(-3.*II*B3*alp3*(1 + bet2*sig2) + 2*B1*bet3*bet3*sig3 - 
               2*B2*bet3*bet3*sig3 + B2*alp3*(3.*II + 3.*II*bet2*sig2 + 2*bet3*sig3))*
             (II + alpha*sigma)*(II + alpha*sigma) + (B2 - B3)*alp3*
             (3.*II*(-B1 + B2) + 3.*II*(-B1 + B2)*alp2*sig2 + 2*B1*alp3*sig3)*(II + beta*sigma)*(II + beta*sigma)))/
       pw2(B2*(alp3 - bet3) + B1*bet3)).imag()/pw2(alp3*B1*bet3*sig3*sig3)/sigma);
  }
}

// The condition at the end of inflation
inline double EoI(const state_type &phi) {
  return phi[0] - phif;
}

inline double EoN(const state_type &phi) {
  return phi[0] - phiN;
}

#elif MODEL==3

// Model parameters
const std::string model = "hybrid"; // Name of the model
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
const double PSI_INIT = 0.;
const double DPSI_INIT = 0.;

// Conditions of EoI and noise adding time
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

// power spectrum of fields
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

#endif


// EoM
void dphidN(const state_type &phi, state_type &dxdt, const double t) {
  double p1 = phi[1]; // pi
  double HH = hubble(phi);

  dxdt[0] = p1/HH;
  dxdt[1] = -3.*p1 - Vphi(phi)/HH;

  #if MODEL==3
  NFLOOP{
    double p2 = phi[2*nf+1]; // pi
    dxdt[2*nf] = p2/HH;
    dxdt[2*nf+1] = -3.*p2 - Vpsi(phi,nf)/HH;
  }
  #endif
  
}

#endif