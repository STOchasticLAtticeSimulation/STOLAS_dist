#include "bias_noise.hpp"

std::vector<double> dwlist(double N);

int main(int argc, char* argv[]) 
{
  if (argc!=2) {
    std::cout << "Specify the noise file number correctly." << std::endl;
    return -1;
  }

  std::ofstream ofs(noisefilename + std::string(argv[1]) + std::string(".dat"));
  if (ofs.fail()) {
    std::cout << "The noise file couldn't be opened. 'mkdir noisedata'" << std::endl;
    return -1;
  }

  
  // ---------- start timer ----------
  struct timeval Nv;
  struct timezone Nz;
  double before, after;
  
  gettimeofday(&Nv, &Nz);
  before = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  // --------------------------------------

#ifdef _OPENMP
  std::cout << "OpenMP : Enabled (Max # of threads = " << omp_get_max_threads() << ")" << std::endl;
#endif

  std::cout << "Box size : " << NL << std::endl;
  
  int totalstep = ceil(log((NL/2-1)/sigma)/dN), count = 0;
  std::vector<std::vector<double>> noisedata(totalstep, std::vector<double>(NL*NL*NL,0));
  
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=0; i<totalstep; i++) {
    noisedata[i] = dwlist(i*dN);
#ifdef _OPENMP
#pragma omp critical
#endif
    {
      count++;
      std::cout << "\rNoiseGenerating : " << std::setw(3) << 100*count/totalstep << "%" << std::flush;
    }
  }
  std::cout << std::endl;
  
  for (size_t i=0; i<noisedata[0].size(); i++) {
    for (size_t n=0; n<noisedata.size(); n++) {
      ofs << noisedata[n][i] << ' ';
    }
    ofs << std::endl;
    std::cout << "\rExporting : " << std::setw(3) << 100*i/noisedata[0].size() << "%" << std::flush;
  }
  std::cout << "\rExporting : 100%" << std::endl;

  // ---------- stop timer ----------
  gettimeofday(&Nv, &Nz);
  after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  std::cout << after - before << " sec." << std::endl;
  // -------------------------------------
}


std::vector<double> dwlist(double N) {
  std::vector<std::vector<std::vector<std::complex<double>>>> dwk(NL, std::vector<std::vector<std::complex<double>>>(NL, std::vector<std::complex<double>>(NL, 0)));
  int count = 0;
  double nsigma = sigma*exp(N);
  std::vector<double> dwlist(NL*NL*NL,0);
  
  LOOP{
    if (innsigma(i,j,k,NL,nsigma,dn)) {
      if (realpoint(i,j,k,NL)) {
	dwk[i][j][k] = dist(engine);
	count++;
      } else if (complexpoint(i,j,k,NL)) {
	dwk[i][j][k] = (dist(engine) + II*dist(engine))/sqrt(2);
	count++;
      }
    }
  }

  // reflection
  int ip, jp, kp; // reflected index
  LOOP{
    if (innsigma(i,j,k,NL,nsigma,dn)) {
      if (!(realpoint(i,j,k,NL)||complexpoint(i,j,k,NL))) {
	if (i==0) {
	  ip = 0;
	} else {
	  ip = NL-i;
	}
	
	if (j==0) {
	  jp = 0;
	} else {
	  jp = NL-j;
	}
	
	if (k==0) {
	  kp = 0;
	} else {
	  kp = NL-k;
	}
	
	dwk[i][j][k] = conj(dwk[ip][jp][kp]);
	count++;
      }
    }
  }

  if (count==0) {
    return dwlist;
  }
  dwk /= sqrt(count);

  std::vector<std::vector<std::vector<std::complex<double>>>> dwlattice = fft(dwk);
  LOOP{
    dwlist[i*NL*NL + j*NL + k] = dwlattice[i][j][k].real();
  }

  return dwlist;
}
