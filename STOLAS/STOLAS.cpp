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

  init_fftw_global();
  if(FFTwisdomFirst) return 0;

  int seed_val = atoi(argv[1]);
  std::mt19937 engine_int(seed_val);

  std::cout << "================ STOLAS ================" << std::endl;
  std::cout << "model          : " << model << " (MODEL=" << MODEL << ")" << std::endl;
  std::cout << "noise seed No. : " << seed_val << std::endl;
  std::cout << "SIGMA (calPzeta0) : " << calPzeta0 << std::endl;
  std::cout << "lattice        : " << NLnoise << "^3, sigma = " << sigma << ", dN = " << dN << std::endl;
  std::cout << "time steps     : firststep = " << firststep << ", totalstep = " << totalstep
            << ", zoom-in number = " << internumber << std::endl;
  std::cout << "output dir     : " << sdatadir << "/" << model << std::endl;
  if (sbias) {
    std::cout << "importance sampling : ON (nbias = " << nbias << ", bias = " << bias
              << ", weight output: " << (doWeight ? "yes" : "no")
              << ", compaction output: " << (doCompaction ? "yes" : "no")
              << ", power spectrum output: no)" << std::endl;
  } else {
    std::cout << "importance sampling : OFF (weight and compaction are not computed, power spectrum output: "
              << (doPower ? "yes" : "no") << ")" << std::endl;
  }
  std::cout << "========================================" << std::endl;

  // create the output directories if they do not exist yet
  {
    std::error_code ec;
    std::filesystem::create_directories(sdatadir + "/" + model, ec);
    if (sanimation) std::filesystem::create_directories(sdatadir + "/" + model + "/animation", ec);
    if (snoisemap) std::filesystem::create_directories(sdatadir + "/" + model + "/noisedata", ec);
  }

  // initialize fields
  initialize();

  // Solve fields
  for (size_t interpolatingnumber = 0; interpolatingnumber < internumber+1; interpolatingnumber++) {
    if (internumber > 0) std::cout << "--- Zoom-in step " << interpolatingnumber << " / " << internumber << " ---" << std::endl;
    OpenFiles(seed_val, interpolatingnumber);
    
    if (checkNfilefail()) {
      std::cout << "The export file couldn't be opened. 'mkdir data'" << std::endl;
      return -1;
    }

    if(interpolatingnumber==0) {
      evolution(seed_val,0,firststep,interpolatingnumber);
    }
    else {
      evolution(seed_val,firststep-itpstep,firststep,interpolatingnumber);
    }

    Phidata = phievol; // save field values for zoom
    evolution(seed_val,firststep,totalstep,interpolatingnumber);

    if(EoI_noise && MeanNumber!=0 && interpolatingnumber==internumber) {
      PhidataAv = phievol;
      for (int av = 0; av < MeanNumber; av++) {
        evolutionNoise(seed_val,av); // Adding the noise until EoN()
        dNmap(interpolatingnumber);
        Ntotal += Ndata;
        Ndata.fill(0.0); // reset vector
        phievol = PhidataAv;
        std::cout << "\rAveraging           : " << std::setw(1) << int(100.*(av+1)/(double)MeanNumber) << "%" << std::flush;
      }
      std::cout << std::endl;
      for (int i=0; i<NLnoiseAll; i++) Ndata[i] = Ntotal[i]/(double)MeanNumber;
    }
    else if(EoI_noise && MeanNumber==0 && interpolatingnumber==internumber) {
      evolutionNoise(seed_val,0); // Adding the noise until EoN()
      dNmap(interpolatingnumber);
    }
    else {
      dNmap(interpolatingnumber);
    }

    save_zeta(); // save delta N map
    if(doPower) spectrum(Ndata,seed_val,interpolatingnumber);
    if(sfield) save_field();
    if(doWeight) weight(seed_val);

    // ComputeLaplacian(Ndata, seed_val);

    if(doCompaction) compaction(Ndata,seed_val);
    #if MODEL==2
      if(susrlength) USRLength(seed_val);
    #endif
    Nfile.close();
    fieldfile.close();


    // Interpolation
    if (interpolatingnumber!=internumber) {
      int shiftx = (int)(NLnoise*0.25*(1.0+std::erf(dist(engine_int)*inv_sqrt2)));// NLnoise/4; //
      int shifty = (int)(NLnoise*0.25*(1.0+std::erf(dist(engine_int)*inv_sqrt2)));
      int shiftz = (int)(NLnoise*0.25*(1.0+std::erf(dist(engine_int)*inv_sqrt2)));
      shift = {shiftx,shifty,shiftz};
      std::cout << shiftx << ' ' << shifty << ' ' << shiftz << std::endl;
      // findNMaxBox(Ndata);
      Ntotal.fill(0.0); // reset vector
      InterpolatingPhi(shift);

      for (int i=0; i<NLnoiseAll; i++){
        #if MODEL==1
          if(phievol[i][0]>0) brokenlist[i] = false;
        #elif MODEL==2
          if (phievol[i][0] > phi1) broken1list[i] = false;
          if (phievol[i][0] > phi2) broken2list[i] = false;
        #endif
      }

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
