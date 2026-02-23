#include "STOLAS.hpp"

// main function
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
  // ---------------------------------

#ifdef _OPENMP
  std::cout << "OpenMP : Enabled (Max # of threads = " << omp_get_max_threads() << ")" << std::endl;
#endif

  int noisefiledirNo = atoi(argv[1]);

  std::cout << "Noise file No. : " << noisefiledirNo << std::endl;
  std::cout << "model : " << model << std::endl;

  // initialize fields
  initialize();

  // Solve fields
  for (size_t interpolatingnumber = 0; interpolatingnumber < internumber+1; interpolatingnumber++) {
    OpenFiles(noisefiledirNo, interpolatingnumber);
    
    if (checkNfilefail()) {
      std::cout << "The export file couldn't be opened. 'mkdir data'" << std::endl;
      return -1;
    }

    // fields evolution
    for (int noiseNo = 0; noiseNo < totalnoiseNo; noiseNo++) {
      STOLAS(noisefiledirNo,noiseNo,interpolatingnumber);
      if (checknoisefile()) {
        std::cout << "The noise file couldn't be opened." << std::endl;
        return -1;
      }

      evolution(noiseNo);

      std::cout << "\rLatticeSimulation   : " << std::setw(3) << 100*noiseNo/(totalnoiseNo-1) << "%" << std::flush;
    }
    std::cout << std::endl;


    // keep fields map
    Phidata = phievol;

    // Solve untill superH
    sanisuperH = true;
    if (sallslice||interpolatingnumber==internumber) {
      for (int noiseNo = 0; noiseNo < totalnoiseNo; noiseNo++) {
        STOLAS(noisefiledirNo,noiseNo,999);
        if (checknoisefile()) {
          std::cout << "The noise file couldn't be opened." << std::endl;
          return -1;
        }
    
        evolution(noiseNo);

        // Take an average for each patches
        if(EoI_noise && MeanNumber!=0) {
          std::vector<std::vector<double>> PhidataAv = phievol;
          for (int av = 0; av < MeanNumber; av++) {
            evolutionNoise(noiseNo); // Adding the noise until EoN()
            dNmap(noiseNo,interpolatingnumber);
            Ntotal += Ndata;
            Ndata.assign(Ndata.size(), 0.0); // reset vector
            phievol = PhidataAv;
          }

          for (int i=0; i<NL; i++) Naverage[i + NL*noiseNo] = Ntotal[i + NL*noiseNo]/(double)MeanNumber;
        }
        else {
          // evolutionNoise(noiseNo); // Adding the noise until EoN()
          dNmap(noiseNo,interpolatingnumber);
          Naverage = Ndata;
        }

        save_zeta(noiseNo,interpolatingnumber); // svave delta N map
        if(sfield) save_field(noiseNo);

        std::cout << "\rCompute delta N map : " << std::setw(3) << 100*noiseNo/(totalnoiseNo-1) << "%" << std::flush;
      }
      std::cout << std::endl;
      if(spower) spectrum(Naverage,interpolatingnumber);
      if(sanimation && interpolatingnumber==0) animation(noisefiledirNo,0);
    }
    Nfile.close();
    fieldfile.close();

    

    // Interpolation
    if (interpolatingnumber!=internumber) {
      std::vector<int> shift{NLnoise/4,NLnoise/4,NLnoise/4}; // std::vector<int> shift = findNMaxBox(Ndata);

      phievol = Phidata; // reset field values
      Ntotal.assign(Ntotal.size(), 0.0); // reset vector
      phievol = InterpolatingPhi(shift);
      std::cout << "Number of interpolation: " << interpolatingnumber+1 << std::endl;
    }
  }

  // ---------- stop timer ----------
  gettimeofday(&Nv, &Nz);
  after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  std::cout << after - before << " sec." << std::endl;
  before = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  // --------------------------------
}
