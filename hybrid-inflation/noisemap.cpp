#include "src/bias_noise.hpp"

std::vector<double> dwlist_fftw(double N);

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

  for (int n = 0; n < internumber+2; n++) {
    int noisestepnum = firststep;
    if (n>0) noisestepnum = itpstep;
    if (n==internumber+1) noisestepnum = totalstep-firststep;

    std::vector<std::unique_ptr<std::ofstream>> ofs_vector;
    for (int i = 0; i < totalnoiseNo; i++) {
      std::string filename = noisefiledir + std::string(argv[1]) + noisefilenamediv + std::to_string(i) + "_" + std::to_string(n) + ".bin";
      if (n==internumber+1) filename = noisefiledir + std::string(argv[1]) + noisefilenamediv + std::to_string(i) + "_" + std::to_string(999) + ".bin";

      ofs_vector.emplace_back(std::make_unique<std::ofstream>(filename, std::ios::out | std::ios::binary));
          
      if (!ofs_vector.back()->is_open()) {
          std::cerr << "Failed to open file: " << filename << std::endl;
          return 1;
      }
    }

fftw_init_threads();
#ifdef _OPENMP
  std::cout << "OpenMP : Enabled (Max # of threads = " << omp_get_max_threads() << ")" << std::endl;
  fftw_plan_with_nthreads(omp_get_max_threads());
#endif

    std::cout << "Box size : " << NLnoise << "  Step number : " << noisestepnum << std::endl;
    
    int divstep = int((noisestepnum)/totalnoiseNo);
    int modstep = int((noisestepnum)%totalnoiseNo);
    
    for (int l=0; l<totalnoiseNo+1; l++) {
      std::vector<std::vector<double>> noisedata;
      if (l<totalnoiseNo) {
        noisedata = std::vector<std::vector<double>>(divstep, std::vector<double>(NLnoise*NLnoise*NLnoise,0));
      } else {
        noisedata = std::vector<std::vector<double>>(modstep, std::vector<double>(NLnoise*NLnoise*NLnoise,0));
      }

#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (int i=0; i<divstep; i++) {
        int Nstep = i+l*divstep;
        if (n>0) Nstep += firststep-itpstep;
        if (n==internumber+1) Nstep += itpstep;
        if (l<totalnoiseNo || i<modstep) noisedata[i] = dwlist_fftw(Nstep*dN);
      }

      for (size_t n=0; n<noisedata.size(); n++) {
        int divcount=0;
        for (size_t m=0; m<noisedata[0].size(); m++) {
          ofs_vector[divcount]->write(reinterpret_cast<const char*>(&noisedata[n][m]), sizeof(double));
          if ((m+1)%NL==0) {
            divcount++;
          }
        }
      }

      // std::cout << "\rExporting :       100%" << std::endl;
      std::cout << "\rNoiseGenerating : " << std::setw(3) << 100*l/(totalnoiseNo) << "%" << std::flush;
    }
    std::cout << "\rNoiseGenerating : 100%" << std::endl;
  }

  // ---------- stop timer ----------
  gettimeofday(&Nv, &Nz);
  after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  std::cout << after - before << " sec." << std::endl;
  // -------------------------------------
}



std::vector<double> dwlist_fftw(double N) {
  std::vector<std::vector<std::vector<std::complex<double>>>> dwk(NLnoise, std::vector<std::vector<std::complex<double>>>(NLnoise, std::vector<std::complex<double>>(NLnoise, 0)));
  int count = 0;
  double nsigma = sigma*exp(N);
  std::vector<double> dwlist(NLnoise*NLnoise*NLnoise,0);
  
  LOOP{
    if (innsigma(i,j,k,NLnoise,nsigma,dn)) {
      if (realpoint(i,j,k,NLnoise)) {
	dwk[i][j][k] = dist(engine);
	count++;
      } else if (complexpoint(i,j,k,NLnoise)) {
	dwk[i][j][k] = (dist(engine) + II*dist(engine))/sqrt(2);
	count++;
      }
    }
  }

  // reflection
  int ip, jp, kp; // reflected index
  LOOP{
    if (innsigma(i,j,k,NLnoise,nsigma,dn)) {
      if (!(realpoint(i,j,k,NLnoise)||complexpoint(i,j,k,NLnoise))) {
	if (i==0) {
	  ip = 0;
	} else {
	  ip = NLnoise-i;
	}
	
	if (j==0) {
	  jp = 0;
	} else {
	  jp = NLnoise-j;
	}
	
	if (k==0) {
	  kp = 0;
	} else {
	  kp = NLnoise-k;
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

  std::vector<std::vector<std::vector<std::complex<double>>>> dwlattice = fft_fftw(dwk);
  LOOP{
    dwlist[i*NLnoise*NLnoise + j*NLnoise + k] = dwlattice[i][j][k].real();
  }

  return dwlist;
}
