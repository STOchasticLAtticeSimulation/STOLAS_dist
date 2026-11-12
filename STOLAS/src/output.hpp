#ifndef INCLUDED_output_hpp_
#define INCLUDED_output_hpp_


// Common file name suffix of the per-simulation outputs:
// <NLnoise>[_<NFIELDS>]_<NoisefiledirNo>[_<Interpolatingnumber>][_<calPzeta0>]
//   NFIELDS             : only for MODEL==3 (hybrid)
//   Interpolatingnumber : only when internumber != 0 (zoom-in is used)
//   calPzeta0 (SIGMA)   : only for MODEL==2 (USR)
std::string SimFileName(int NoisefiledirNo, int Interpolatingnumber){
  std::string name = std::to_string(NLnoise);
#if MODEL==3
  name += std::string("_") + std::to_string(NFIELDS);
#endif
  name += std::string("_") + std::to_string(NoisefiledirNo);
  if (internumber != 0) name += std::string("_") + std::to_string(Interpolatingnumber);
#if MODEL==2
  name += std::string("_") + std::to_string(int(calPzeta0));
#endif
  return name;
}

// Write a '#' header line at the top of a text output file.
// Files opened in append mode (probabilities.dat) get the header only
// when they are empty, so a header is written exactly once.
inline void write_header_if_empty(std::ofstream& ofs, const std::string& header) {
  ofs.seekp(0, std::ios_base::end);
  if (ofs.tellp() == std::streampos(0)) ofs << header << std::endl;
}

// column headers shared by compaction() and compaction_recenter() (functions.cpp)
const std::string prbfileheader = "# seed  mu2hat  k3hat  k_star*r_m  zeta_m  Cbar_m  lnW";
const std::string cmpfileheader = "# r  C(r)   (r: radius in lattice units, C: compaction function)";

void OpenFiles(int NoisefiledirNo, int Interpolatingnumber){
  std::string FileName = SimFileName(NoisefiledirNo, Interpolatingnumber);

  Nfile.open(Nfileprefix + FileName + std::string(".bin"), std::ios::binary);
  Nfile << std::setprecision(10);
  Nfilefail = Nfile.fail();

  if (sfield) {
    fieldfile.open(fieldfileprefix + FileName + std::string(".dat"));
    fieldfile << std::setprecision(10);
  }
  
  if (strajectory) {
    trajectoryfile.open(trajectoryfileprefix + FileName + std::string(".dat"));
    trajectoryfile << std::setprecision(10);
    #if MODEL==3
      write_header_if_empty(trajectoryfile, "# N  phi  sum(psi^2)  EoI-2  H");
    #else
      write_header_if_empty(trajectoryfile, "# N  phi");
    #endif
  }

  if (doWeight) {
    logwfile.open(logwfileprefix + FileName + std::string(".dat"));
  }

  // if (snoisemap) {
  //   Noisefile.open(sdatadir + "/" + model + "/noisedata/map_" + FileName + std::string(".bin"), std::ios::binary);
  //   Noisefile << std::setprecision(10);
  // }
  
}

#if MODEL==2
void USRLength(int NoisefiledirNo){
  std::ofstream N1N2fileA;
  std::string FileName = SimFileName(NoisefiledirNo, 0);
  N1N2fileA.open(sdatadir + "/" + model + "/N1N2Length_" + FileName + std::string(".bin"), std::ios::binary);
  N1N2fileA << std::setprecision(15);

  static std::array<double, 2 * NLnoiseAll> buffer;
  for (int i=0; i<NLnoiseAll; i++) {
      buffer[2*i] = N1list[i];
      buffer[2*i + 1] =  N2list[i];
  }
  N1N2fileA.write(reinterpret_cast<const char*>(buffer.data()), sizeof(double) * buffer.size());
  std::cout << "Export N1 and N2 Length" << std::endl;
  N1N2fileA.close();
}
#endif

void save_zeta(){
  Nfile.write(reinterpret_cast<const char*>(&Ndata), sizeof(double) * NLnoiseAll);
  std::cout << "Export delta N map" << std::endl;
}

void save_field(){
  for (int i=0; i<NLnoiseAll; i++) {
    #if MODEL==3
      NFLOOP{
        fieldfile << phievol[i][2*nf] << ' ';
      }
      fieldfile << std::endl;
    #else
      fieldfile << phievol[i][0] << std::endl;
    #endif
  }
}

void save_trajectory(state_type PHI, double Ntime){
  trajectoryfile << Ntime << ' ' << PHI[0];
  #if MODEL==3
    double PSISUM = 0;
    NFLOOP PSISUM += pw2(PHI[2*nf]);
    trajectoryfile << ' ' << PSISUM << ' ' << EoI(PHI)-2 << ' ' << hubble(PHI);
  #endif
  trajectoryfile << std::endl;
}

void animation(std::array<state_type,NLnoiseAll>& phievol, int NoisefiledirNo, int Ntime, int Interpolatingnumber){
  std::string FileName = SimFileName(NoisefiledirNo, Interpolatingnumber) + std::string("_") + std::to_string(Ntime);
  fieldfileA.open(animationfileprefix + FileName + std::string(".bin"), std::ios::binary);
  fieldfileA << std::setprecision(10);

  static std::array<double, 2 * NLnoiseAll> buffer;

  for (int i=0; i<NLnoiseAll; i++) {
      buffer[2*i] = phievol[i][1];
      buffer[2*i + 1] = hubble(phievol[i]);
  }
  fieldfileA.write(reinterpret_cast<const char*>(buffer.data()), sizeof(double) * buffer.size());

  // for (int i=0; i<NLnoiseAll; i++){
  //   fieldfileA << phievol[i][1] << ' ' << hubble(phievol[i]) << std::endl;

  //   #if MODEL==3
  //     double PSISUM = 0;
  //     NFLOOP PSISUM += pw2(phievol[i][2*nf]);
  //     fieldfileA << sqrt(PSISUM) << std::endl;

  //   #endif
  // }
  fieldfileA.close();
}

// Calculate power spectrum
// Output: power_<SimFileName>.dat with 2 columns, n and P(n)
void spectrum(std::array<double,NLnoiseAll>& Ndata, int noisefiledirNo, int Interpolatingnumber) {
  disc_power.fill(0.); // each call (each zoom-in step) gives its own spectrum
  powfile.open(powfileprefix + SimFileName(noisefiledirNo, Interpolatingnumber) + std::string(".dat"));
  powfile << std::setprecision(10);
  {
    std::ostringstream hdr;
    hdr << "# n  P   (n = |k| in units of the fundamental mode 2pi/L; P: power spectrum of zeta per ln n;"
        << " bin centre ln n = dlogn*ii, dlogn = " << dlogn << ")";
    write_header_if_empty(powfile, hdr.str());
  }

  fft_1D_real(Ndata);

  for (int i = 0; i < NLnoise; i++) {
    int nxt = (i<=NLnoise/2 ? i : i-NLnoise);
    for (int j = 0; j < NLnoise; j++) {
      int nyt = (j<=NLnoise/2 ? j : j-NLnoise);
      for (int k = 0; k < NLnoiseHalf; k++) {
        int nzt = k;
        int idx = i*NLnoise*NLnoiseHalf + j*NLnoiseHalf + k;
        double weight = (k==0 || k==NLnoise/2) ? 1. : 2.;

        double rk=nxt*nxt+nyt*nyt+nzt*nzt;

        double LogNk = log(sqrt(rk));
        double calPk = weight*norm(bkspectrum[idx])/NLnoise/NLnoise/NLnoise/NLnoise/NLnoise/NLnoise;
        for (size_t ii = 0; ii < imax; ii++) {
          if ((dlogn*ii-LogNk)<dlogn/2. && -dlogn/2.<=(dlogn*ii-LogNk) && rk!=0) {
            disc_power[ii] += calPk/dlogn;
            break;
          }
        }
      }
    }
  }
  for (size_t ii = 0; ii < imax; ii++) {
    powfile << std::exp(dlogn*ii) << ' ' << disc_power[ii] << std::endl;
  }
  powfile.close();
  std::cout << "ExportPowerSpectrum" << std::endl;
}

// calculation of weight
void weight(int seed) {
  double logw = 0.;
  for (size_t n=0; n<weightlist.size(); n++) {
    double N = n*dN;
    double Bias = bias/dNbias/sqrt(2*M_PI)*exp(-(N-Nbias)*(N-Nbias)/2./dNbias/dNbias);
    if(weightbool[n]) logw -= Bias*weightlist[n]*sqrt(dN) + (Bias*Bias*dN)/2;
  }
  logwfile << seed << ' ' << logw << std::endl;
}

// Calculate mu2
void ComputeLaplacian(std::array<double,NLnoiseAll>& Ndata, int noisefiledirNo) {
  std::string FileName = SimFileName(noisefiledirNo, 0);

  mu2file.open(mu2fileprefix + FileName + std::string(".bin"), std::ios::binary);
  k3file.open(k3fileprefix + FileName + std::string(".bin"), std::ios::binary);
  mu2file << std::setprecision(10);
  k3file << std::setprecision(10);

  computeLaplacian(Ndata);
  for (size_t n = 0; n < NLnoiseAll; n++) {
    mutwo[n] = -1.*laplacian[n];
  }

  mu2file.write(reinterpret_cast<const char*>(&mutwo), sizeof(double) * NLnoiseAll);
  std::cout << "Export mu2 map" << std::endl;

  computeLaplacian(laplacian);
  for (size_t n = 0; n < NLnoiseAll; n++) {
    if(mutwo[n]<1.e-12) {
      kthree[n] = 0;
    }
    else {
      kthree[n] = laplacian[n]/mutwo[n];
    }
  }

  k3file.write(reinterpret_cast<const char*>(&kthree), sizeof(double) * NLnoiseAll);
  std::cout << "Export k3 map" << std::endl;
}

// Lattice Laplacian in lattice units (dx=1) at a single point
double LatticeLaplacianAt(const std::array<double,NLnoiseAll>& x, int i, int j, int k) {
  return x[index(INCREMENT(i),j,k)] + x[index(DECREMENT(i),j,k)]
       + x[index(i,INCREMENT(j),k)] + x[index(i,DECREMENT(j),k)]
       + x[index(i,j,INCREMENT(k))] + x[index(i,j,DECREMENT(k))]
       - 6.*x[index(i,j,k)];
}

// Lattice bi-Laplacian in lattice units (dx=1) at a single point
double LatticeBiLaplacianAt(const std::array<double,NLnoiseAll>& x, int i, int j, int k) {
  return LatticeLaplacianAt(x,INCREMENT(i),j,k) + LatticeLaplacianAt(x,DECREMENT(i),j,k)
       + LatticeLaplacianAt(x,i,INCREMENT(j),k) + LatticeLaplacianAt(x,i,DECREMENT(j),k)
       + LatticeLaplacianAt(x,i,j,INCREMENT(k)) + LatticeLaplacianAt(x,i,j,DECREMENT(k))
       - 6.*LatticeLaplacianAt(x,i,j,k);
}

// Calculate compaction function
void compaction(std::array<double,NLnoiseAll>& Ndata, int noisefiledirNo) {
  prbfile.open(prbfileprefix + std::string(".dat"), std::ios::app);
  cmpfile.open(cmpfileprefix + SimFileName(noisefiledirNo, 0) + std::string(".dat"));
  prbfile << std::setprecision(10);
  cmpfile << std::setprecision(10);
  write_header_if_empty(prbfile, prbfileheader);
  write_header_if_empty(cmpfile, cmpfileheader);
  
  double Naverage = 0;
  double dr = 1;
  for (size_t n = 0; n < NLnoiseAll; n++) {
    Naverage += Ndata[n];
  }
  Naverage /= NLnoise*NLnoise*NLnoise;

  // zeta map
  for (size_t n = 0; n < NLnoiseAll; n++) {
    Ndata[n] -= Naverage;
  }

  int xmax = 0, ymax = 0, zmax = 0;
  // {
  //   constexpr int centerSearchRadius = 5; // lattice units; tune to taste
  //   double bestVal = Ndata[0];
  //   for (int dx = -centerSearchRadius; dx <= centerSearchRadius; dx++) {
  //     int nx = ((dx % NLnoise) + NLnoise) % NLnoise;
  //     for (int dy = -centerSearchRadius; dy <= centerSearchRadius; dy++) {
  //       int ny = ((dy % NLnoise) + NLnoise) % NLnoise;
  //       for (int dz = -centerSearchRadius; dz <= centerSearchRadius; dz++) {
  //         int nz = ((dz % NLnoise) + NLnoise) % NLnoise;
  //         double val = Ndata[nx*NLnoise*NLnoise + ny*NLnoise + nz];
  //         if (val > bestVal) {
  //           bestVal = val;
  //           xmax = dx; ymax = dy; zmax = dz;
  //         }
  //       }
  //     }
  //   }
  // }

  // radial profile
  for (size_t i=0; i<NLnoiseAll; i++) {
    int nx=i/NLnoise/NLnoise ,ny=(i%(NLnoise*NLnoise))/NLnoise, nz=i%NLnoise;
    
    // centering
    int nxt, nyt, nzt; // shifted index
    nxt = (nx<=NLnoise/2+xmax ? nx-xmax : nx-xmax-NLnoise);
    nyt = (ny<=NLnoise/2+ymax ? ny-ymax : ny-ymax-NLnoise);
    nzt = (nz<=NLnoise/2+zmax ? nz-zmax : nz-zmax-NLnoise);

    for (size_t ri=0; ri<NLnoise/2; ri++) {
      double norm = std::abs(sqrt(nxt*nxt+nyt*nyt+nzt*nzt)-ri);
      if (norm<=dr/2.) {
        zetar[0][ri]++;
        zetar[1][ri]+=Ndata[i];
        break;
      }
    }
  }
  for (size_t ri=0; ri<NLnoise/2; ri++) {
    zetar[1][ri] /= zetar[0][ri]; // average
  }

  // derivative zeta
  for (size_t ri=1; ri<NLnoise/2-1; ri++) {
    dzetar[ri] = (zetar[1][ri+1] - zetar[1][ri-1])/(2.*dr);
  }

  // compaction function
  double CompactionMax=0, CompactionInt=0, rmax=0, Rmax=0, IntTemp=0;
  bool Cnegative = false;
  for (size_t ri=0; ri<NLnoise/2; ri++) {
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

  // compute weight
  double logw = 0.;
  for (size_t n=0; n<weightlist.size(); n++) {
    double Nt = n*dN;
    double Bias = bias/dNbias/sqrt(2*M_PI)*exp(-(Nt-Nbias)*(Nt-Nbias)/2./dNbias/dNbias);
    if(weightbool[n]) logw -= Bias*weightlist[n]*sqrt_dN + (Bias*Bias*dN)/2;
  }

  // peak curvature in lattice units: mu2hat = -lap(zeta_pk), k3hat^2 = lap^2(zeta_pk)/(-lap(zeta_pk))
  int ipk = (xmax+NLnoise)%NLnoise, jpk = (ymax+NLnoise)%NLnoise, kpk = (zmax+NLnoise)%NLnoise;
  double mu2hat = -LatticeLaplacianAt(Ndata, ipk, jpk, kpk);
  double k3hat2 = LatticeBiLaplacianAt(Ndata, ipk, jpk, kpk)/mu2hat;
  double k3hat = (k3hat2 >= 0 ? sqrt(k3hat2) : -sqrt(-k3hat2)); // signed if k3hat^2<0

  double kstar = 2.*M_PI*nbias/NLnoise; // bias wavenumber in lattice units
  double zetam = zetar[1][(size_t)rmax];

  // seed, mu2hat, k3hat, k_* r_m, zeta_m, Cbar_m, ln W
  prbfile << noisefiledirNo << ' ' << mu2hat << ' ' << k3hat << ' ' << kstar*rmax
  << ' ' << zetam << ' ' << CompactionInt << ' ' << logw << std::endl;
  std::cout << "ExportCompactionFunction" << std::endl;
}

#endif