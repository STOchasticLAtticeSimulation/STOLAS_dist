#include "STOLAS.hpp"
#include "vec_op.hpp"
#include "fft.hpp"

// useful macro
#define LOOP for(int i = 0; i < NL; i++) for(int j = 0; j < NL; j++) for(int k = 0; k < NL; k++)


STOLAS::STOLAS(std::string Model, double DN, std::string sourcedir, int NoisefileNo, std::vector<double> Phii, double Bias, double NBias, double DNbias) {

#ifdef _OPENMP
  std::cout << "OpenMP : Enabled (Max # of threads = " << omp_get_max_threads() << ")" << std::endl;
#endif
  
  model = Model;
  dN = DN;
  noisefileNo = NoisefileNo;
  phii = Phii;
  bias = Bias;
  Nbias = NBias;
  dNbias = DNbias;

  std::cout << "Noise file No. : " << noisefileNo << std::endl;

  noisefile.open(sourcedir + std::string("/") + noisefilename + std::to_string(noisefileNo) + std::string(".dat"));
  noisefilefail = noisefile.fail();
  biasfile.open(sourcedir + std::string("/") + biasfilename);
  biasfilefail = biasfile.fail();
  
  if (!noisefile.fail() && !biasfile.fail()) {
    std::cout << "model : " << model << std::endl;
    
    std::string str;
    double dd;
    while (std::getline(noisefile, str)) {
      std::vector<double> vv;
      std::stringstream ss(str);
      while (!ss.eof()) {
	ss >> dd;
	vv.push_back(dd);
      }
      vv.pop_back();
      noisedata.push_back(vv);
    }

    while (std::getline(biasfile, str)) {
      std::vector<double> vv;
      std::stringstream ss(str);
      while (!ss.eof()) {
	ss >> dd;
	vv.push_back(dd);
      }
      vv.pop_back();
      biasdata.push_back(vv);
    }

    NL = cbrt(noisedata.size());
    std::cout << "Noise/Bias data imported. Box size is " << NL << "." << std::endl;
    Nfile.open(Nfileprefix + std::to_string(NL) + std::string("_") + std::to_string(noisefileNo) + std::string(".dat"));

    Hdata = std::vector<std::vector<double>>(noisedata[0].size(), std::vector<double>(NL*NL*NL,0));
    pidata = std::vector<std::vector<double>>(noisedata[0].size(), std::vector<double>(NL*NL*NL,0));
    Ndata = std::vector<double>(NL*NL*NL,0);
    Nmap3D = std::vector<std::vector<std::vector<std::complex<double>>>>(NL, std::vector<std::vector<std::complex<double>>>(NL, std::vector<std::complex<double>>(NL, 0)));
  }
}


bool STOLAS::checknoisefile() {
  return !noisefilefail;
}

bool STOLAS::checkbiasfile() {
  return !biasfilefail;
}

bool STOLAS::noisebiassize() {
  return (noisedata.size() == biasdata.size());
}

bool STOLAS::Nfilefail() {
  return Nfile.fail();
}


void STOLAS::dNmap() {
  Nfile << std::setprecision(10);
  int complete = 0;
  
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=0; i<NL*NL*NL; i++) {
    double N=0;
    #if BrokenPoint==1
      double N0;
      bool broken = false;
    #endif
    std::vector<double> phi = phii;
    for (size_t n=0; n<noisedata[i].size(); n++) {
      if (sanimation) {
        Hdata[n][i] = pow(hubble(phi[0],phi[1]),2);
        pidata[n][i] = phi[1]*phi[1];
      }
      
      #if BrokenPoint==0
        RK4Mbias(N,phi,dN,noisedata[i][n],biasdata[i][n]);
      #elif BrokenPoint==1
        RK4Mbias(N,phi,dN,noisedata[i][n],biasdata[i][n],N0,broken);
        if (!broken && phi[0] < 0) {
          N0 = N;
          broken = true;
        }
      #endif

    }

    double dN1 = dN;
    std::vector<double> prephi(2);
    while (dN1 >= Nprec) {
      while (EoI(phi)) {
	prephi[0] = phi[0];
	prephi[1] = phi[1];
	RK4(N,phi,dN1);
      }
      N -= dN1;

      phi[0] = prephi[0];
      phi[1] = prephi[1];
      dN1 *= 0.1;
    }

#ifdef _OPENMP
#pragma omp critical
#endif
    {
      Ndata[i] = N;
      Nfile << i << ' ' << N << std::endl;
      complete++;
      std::cout << "\rLatticeSimulation : " << std::setw(3) << 100*complete/NL/NL/NL << "%" << std::flush;
    }
  }
  std::cout << std::endl;

}

// Calculate power spectrum
void STOLAS::spectrum() {
  powfile.open(powfileprefix + std::to_string(NL) + std::string("_") + std::to_string(noisefileNo) + std::string(".dat"));
  powfile << std::setprecision(10);
  for (int i=0; i<NL*NL*NL; i++) {
    int x=i/NL/NL ,y=(i%(NL*NL))/NL, z=i%NL;
    Nmap3D[x][y][z] = Ndata[i];
  }
  std::vector<std::vector<std::vector<std::complex<double>>>> Nk=fft(Nmap3D);
  
  powsfile.open(powsfileprefix + std::string(".dat"), std::ios::app);
  powsfile << std::setprecision(10);
  int imax = ceil(log(NL/2)/dlnk);
  std::vector<double> disc_power(imax, 0);

  LOOP{
    int nxt, nyt, nzt; // shifted index
    if (i<=NL/2) {
      nxt = i;
    } else {
      nxt = i-NL;
    }

    if (j<=NL/2) {
      nyt = j;
    } else {
      nyt = j-NL;
    }

    if (k<=NL/2) {
      nzt = k;
    } else {
      nzt = k-NL;
    }
    
    double rk=nxt*nxt+nyt*nyt+nzt*nzt;
    powfile << sqrt(rk) << " " << norm(Nk[i][j][k])/NL/NL/NL/NL/NL/NL << std::endl;

    double LogNk = log(sqrt(rk));
    double calPk = norm(Nk[i][j][k])/NL/NL/NL/NL/NL/NL;
    for (size_t ii = 0; ii < imax; ii++) {
      if (std::abs(dlnk*ii-LogNk)<=dlnk/2.) {
        disc_power[ii] += calPk/dlnk;
        break;
      }
    }
  }
  powsfile << noisefileNo << " ";
  for (size_t ii = 0; ii < imax; ii++) {
    powsfile << disc_power[ii] << " " ;
  }
  powsfile << std::endl;
  std::cout << "ExportPowerSpectrum" << std::endl;
}

// Calculate compaction function
void STOLAS::compaction() {
  prbfile.open(prbfileprefix + std::string(".dat"), std::ios::app);
  cmpfile.open(cmpfileprefix + std::to_string(NL) + std::string("_") + std::to_string(noisefileNo) + std::string(".dat"));
  prbfile << std::setprecision(10);
  cmpfile << std::setprecision(10);

  //calculation of weight
  double logw = 0.;
  for (size_t n=0; n<noisedata[0].size(); n++) {
    double N = n*dN;
    double Bias = bias/dNbias/sqrt(2*M_PI)*exp(-(N-Nbias)*(N-Nbias)/2./dNbias/dNbias);
    logw -= Bias*noisedata[0][n]*sqrt(dN) + (Bias*Bias*dN)/2;
  }

  double Naverage = 0;
  int dr = 1;
  for (size_t n = 0; n < Ndata.size(); n++) {
    Naverage += Ndata[n];
  }
  Naverage /= NL*NL*NL;

  // zeta map
  for (size_t n = 0; n < Ndata.size(); n++) {
    Ndata[n] -= Naverage;
  }

  // radial profile
  std::vector<std::vector<double>> zetar(2, std::vector<double>(NL/2,0));
  for (size_t i=0; i<NL*NL*NL; i++) {
    int nx=i/NL/NL ,ny=(i%(NL*NL))/NL, nz=i%NL;

    // centering
    if (nx<=NL/2) {
      nx = nx;
    }
    else {
      nx = nx-NL;
    }
    if (ny<=NL/2) {
      ny = ny;
    }
    else {
      ny = ny-NL;
    }
    if (nz<=NL/2) {
      nz = nz;
    }
    else {
      nz = nz-NL;
    }

    for (size_t ri=0; ri<NL/2; ri++) {
      double norm = std::abs(sqrt(nx*nx+ny*ny+nz*nz)-ri);
      if (norm<=dr/2.) {
        zetar[0][ri]++;
        zetar[1][ri]+=Ndata[i];
        break;
      }
    }
  }
  for (size_t ri=0; ri<NL/2; ri++) {
    zetar[1][ri] /= zetar[0][ri]; // average
  }

  // derivative zeta
  std::vector<double> dzetar(NL/2,0);
  for (size_t ri=1; ri<NL/2-1; ri++) {
    dzetar[ri] = (zetar[1][ri+1] - zetar[1][ri-1])/(2.*dr);
  }

  // compaction function
  double CompactionMax=0, CompactionInt=0, rmax=0, Rmax=0, IntTemp=0;
  bool Cnegative = false;
  for (size_t ri=0; ri<NL/2; ri++) {
    double CompactionTemp = 2./3.*(1. - pow(1 + ri*dzetar[ri], 2));
    IntTemp += ri*ri*CompactionTemp*exp(3.*zetar[1][ri])*(1 + ri*dzetar[ri]);

    if (!Cnegative) {
      if (CompactionTemp < -0.2) {
	Cnegative = true;
      }
      
      if (CompactionMax<CompactionTemp) {
	CompactionMax = CompactionTemp;
	rmax = ri;
	Rmax = exp(zetar[1][ri])*ri;
	CompactionInt += IntTemp;
	IntTemp = 0;
      }
    }
    
    cmpfile << ri << ' ' << CompactionTemp << std::endl;
  }
  CompactionInt /= pow(Rmax, 3)/3.;

  prbfile << noisefileNo << ' ' << logw << ' ' << CompactionInt << ' ' << CompactionMax << ' ' << Rmax << ' ' << rmax << std::endl;
  std::cout << "ExportCompactionFunction" << std::endl;
}

// Export animation
void STOLAS::animation() {
  Hfile.open(Hfileprefix + std::to_string(NL) + std::string("_") + std::to_string(noisefileNo) + std::string(".dat"));
  pifile.open(pifileprefix + std::to_string(NL) + std::string("_") + std::to_string(noisefileNo) + std::string(".dat"));
  Hfile << std::setprecision(14);
  pifile << std::setprecision(14);

  for (size_t n=0; n<Hdata.size(); n++) {
    for (size_t i=0; i<Hdata[n].size(); i++) {
      Hfile << Hdata[n][i] << ' ';
      pifile << pidata[n][i] << ' ';
    }
    Hfile << std::endl;
    pifile << std::endl;
    std::cout << "\rAnimeDataExporting : " << std::setw(3) << 100*(n+1)/Hdata.size() << "%" << std::flush;
  }
  std::cout << std::endl;
}


double STOLAS::ep(double phi, double pi) {
  double HH = hubble(phi,pi);
  return pi*pi/2./HH/HH;
}

double STOLAS::hubble(double phi, double pi) {
  return sqrt((pi*pi/2. + VV(phi))/3.);
}


std::vector<double> STOLAS::dphidN(double N, std::vector<double> phi) {
  std::vector<double> dphidN(2);

  double xx = phi[0]; // phi
  double pp = phi[1]; // pi
  double HH = hubble(xx,pp);

  dphidN[0] = pp/HH;
  dphidN[1] = -3*pp - Vp(xx)/HH;
  
  return dphidN;
}

void STOLAS::RK4(double &t, std::vector<double> &x, double dt) {
  std::vector<double> kx[4]; // 4-stage slopes kx
  double a[4][4],b[4],c[4]; // Butcher

  // -------------- initialise kx, a, b, c --------------- //
  for (int i=0;i<=3;i++) {
    kx[i] = x;
    vec_op::init(kx[i]);
    
    for (int j=0;j<=3;j++) {
      a[i][j]=0.;
    }
  }
  
  a[1][0]=1./2;  a[2][1]=1./2;  a[3][2]=1.;
  b[0]=1./6;     b[1]=1./3;     b[2]=1./3;    b[3]=1./6; 
  c[0]=0;        c[1]=1./2;     c[2]=1./2;    c[3]=1;
  // ----------------------------------------------------- //
  

  std::vector<double> X = x; // position at i-stage
    
  for (int i=0;i<=3;i++) {
    X = x; // initialise X
    
    for (int j=0;j<=3;j++) {
      X += dt * a[i][j] * kx[j];
      kx[i] = dphidN(t,X);
    }
  }

  t += dt;
  x += dt*(b[0]*kx[0] + b[1]*kx[1] + b[2]*kx[2] + b[3]*kx[3]);
}

#if BrokenPoint==0
void STOLAS::RK4Mbias(double &N, std::vector<double> &phi, double dN, double dw, double Bias) {
  double phiamp = sqrt(calPphi(phi));

  RK4(N,phi,dN);
  phi[0] += phiamp * dw * sqrt(dN);

  double GaussianFactor = 1./dNbias/sqrt(2*M_PI) * exp(-(N-Nbias)*(N-Nbias)/2./dNbias/dNbias);
  phi[0] += phiamp * bias * Bias * GaussianFactor * dN;
}
#elif BrokenPoint==1
void STOLAS::RK4Mbias(double &N, std::vector<double> &phi, double dN, double dw, double Bias,
                      double N0, bool broken) {
  
  double phiamp = sqrt(calPphi(N,phi,N0,broken));
  double piamp = sqrt(calPpi(N,phi,N0,broken));
  double crosscor = RecalPphipi(N,phi,N0,broken);

  RK4(N,phi,dN);
  phi[0] += phiamp * dw * sqrt(dN);

  if (crosscor > 0) {
    phi[1] += piamp * dw * sqrt(dN);
  } else {
    phi[1] -= piamp * dw * sqrt(dN);
  }

  double GaussianFactor = 1./dNbias/sqrt(2*M_PI) * exp(-(N-Nbias)*(N-Nbias)/2./dNbias/dNbias);
  phi[0] += phiamp * bias * Bias * GaussianFactor * dN;

  if (crosscor > 0) {
    phi[1] += piamp * bias * Bias * GaussianFactor * dN;
  } else {
    phi[1] -= piamp * bias * Bias * GaussianFactor * dN;
  }
}
#endif