#include "STOLAS.hpp"

std::ifstream Ninfile;

// input files are stored in subdirectories: <sdatadir>/<model>/Nfile/, <sdatadir>/<model>/logwfile/
const std::string Ninfileprefix = sdatadir + "/" + model + "/Nfile/Nmap_";
const std::string logwinfileprefix = sdatadir + "/" + model + "/logwfile/logw_";
// recentered compaction profiles go to <sdatadir>/<model>/cmpfile_recenter/
const std::string cmprecenterprefix = sdatadir + "/" + model + "/cmpfile_recenter/compaction_";

// compaction() recentered on the mu2 (= -lap zeta) peak near the origin.
// ln W is read from the logw file saved by STOLAS, since weightlist is empty here.
// Outputs: cmpfile_recenter/compaction_<...>_recenter.dat, probabilities_recenter.dat
void compaction_recenter(std::array<double,NLnoiseAll>& Ndata, int noisefiledirNo, double logw) {
  {
    std::error_code ec;
    std::filesystem::create_directories(sdatadir + "/" + model + "/cmpfile_recenter", ec);
  }
  prbfile.open(prbfileprefix + std::string("_recenter.dat"), std::ios::app);
  cmpfile.open(cmprecenterprefix + SimFileName(noisefiledirNo, 0) + std::string("_recenter.dat"));
  if (cmpfile.fail()) std::cout << "cmpfile fail: mkdir " << sdatadir + "/" + model + "/cmpfile_recenter" << std::endl;
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

  // center: mu2 peak
  int xmax = 0, ymax = 0, zmax = 0;
  {
    constexpr int centerSearchRadius = 5; // lattice units; tune to taste
    double bestVal = -LatticeLaplacianAt(Ndata, 0, 0, 0);
    for (int dx = -centerSearchRadius; dx <= centerSearchRadius; dx++) {
      int nx = ((dx % NLnoise) + NLnoise) % NLnoise;
      for (int dy = -centerSearchRadius; dy <= centerSearchRadius; dy++) {
        int ny = ((dy % NLnoise) + NLnoise) % NLnoise;
        for (int dz = -centerSearchRadius; dz <= centerSearchRadius; dz++) {
          int nz = ((dz % NLnoise) + NLnoise) % NLnoise;
          double val = -LatticeLaplacianAt(Ndata, nx, ny, nz);
          if (val > bestVal) {
            bestVal = val;
            xmax = dx; ymax = dy; zmax = dz;
          }
        }
      }
    }
  }

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
  prbfile.close();
  cmpfile.close();
  std::cout << "ExportCompactionFunction (recentered)" << std::endl;
}

int main(int argc, char* argv[])
{
  if (argc!=2) {
    std::cout << "Specify the noise file number correctly." << std::endl;
    return -1;
  }
  if (!sbias) {
    std::cout << "functions computes the compaction function, which needs the importance sampling. Set sbias = true in parameters.hpp." << std::endl;
    return -1;
  }
  
  // ---------- start timer ----------
  struct timeval Nv;
  struct timezone Nz;
  double before, after;
  
  gettimeofday(&Nv, &Nz);
  before = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  // --------------------------------------

  int noisefiledirNo = atoi(argv[1]);
  std::string FileName = SimFileName(noisefiledirNo, 0);
  
  std::cout << Ninfileprefix + FileName << std::endl;
  Ninfile.open(Ninfileprefix + FileName + std::string(".bin"), std::ios::binary);
  if (Ninfile.fail()) {
    std::cout << "Nfile fail" << std::endl;
    return -1;
  }

  Ninfile.read(reinterpret_cast<char*>(Ndata.data()), sizeof(double) * NLnoiseAll);

  // ln W saved by STOLAS
  double logw = 0., seedtmp;
  std::ifstream logwinfile(logwinfileprefix + FileName + std::string(".dat"));
  if (!(logwinfile >> seedtmp >> logw)) std::cout << "logw file fail (ln W = 0)" << std::endl;

  // if(spower) spectrum(Ndata,noisefiledirNo);
  // mu2(Ndata, noisefiledirNo);
  // if(scompaction) compaction(Ndata,noisefiledirNo);
  compaction_recenter(Ndata,noisefiledirNo,logw);


  // ---------- stop timer ----------
  gettimeofday(&Nv, &Nz);
  after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  std::cout << after - before << " sec." << std::endl;
  std::cout << std::endl;
  // -------------------------------------
}
