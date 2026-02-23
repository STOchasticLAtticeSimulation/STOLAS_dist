#ifndef INCLUDED_output_hpp_
#define INCLUDED_output_hpp_


void OpenFiles(int NoisefiledirNo, int Interpolatingnumber){
  std::string NLfilename = std::to_string(NLnoise) + std::string("_") + std::to_string(NFIELDS) + std::string("_") + std::to_string(NoisefiledirNo);
  std::string InterFileName = NLfilename + std::string("_") + std::to_string(Interpolatingnumber) + std::string(".dat");

  Nfile.open(Nfileprefix + InterFileName);
  Nfile << std::setprecision(10);
  Nfilefail = Nfile.fail();

  if (sfield) {
    fieldfile.open(fieldfileprefix + InterFileName);
    fieldfile << std::setprecision(10);
  }
  
  if (strajectory) {
    trajectoryfile.open(trajectoryfileprefix + InterFileName);
    trajectoryfile << std::setprecision(10);
  }

  if (sanimation) {
    phianimation=std::vector<std::vector<std::vector<double>>> (int(totalstep/aninum),std::vector<std::vector<double>>(NLnoise*NLnoise*NLnoise, std::vector<double>(2*NFIELDS+2,0.)));
  }
  
}

void save_zeta(int NoiseNo, int Interpolatingnumber){
  for (int i=0; i<NL; i++){
    double Nsave = 0.;
    Nsave = Naverage[i + NL*NoiseNo];
    Nfile << i + NL*NoiseNo << ' ' << Nsave << std::endl;
    }
}

void save_field(int NoiseNo){
  for (int i=0; i<NL; i++) {
    fieldfile << i + NL*NoiseNo << ' ';
    NFLOOP{
      fieldfile << phievol[i + NL*NoiseNo][2*nf] << ' ';
    }
    fieldfile << std::endl;
  }
}

void save_trajectory(state_type PHI, double Ntime){
  trajectoryfile << Ntime << ' ' << PHI[0] << ' ';
  double PSISUM = 0;
  NFLOOP PSISUM += pw2(PHI[2*nf]);
  trajectoryfile << PSISUM << ' ' << EoI(PHI)-2 << ' ' << hubble(PHI) << std::endl;
}

void animation(int NoisefiledirNo, int EachField){
  std::cout << "ExportAnimation" << std::endl;

  if (EachField==0){
    for (int n = 0; n < phianimation.size(); n++) {
      std::string fFileName = std::to_string(NLnoise) + std::string("_") + std::to_string(NFIELDS) + std::string("_") + std::to_string(NoisefiledirNo) + std::string("_") + std::to_string(n);
      fieldfileA.open(animationfileprefix + fFileName + std::string(".dat"));
      fieldfileA << std::setprecision(10);

      for (int i=0; i<phianimation[n].size(); i++) {
        double PSISUM = 0;
        NFLOOP{
          PSISUM += pw2(phianimation[n][i][2*nf]);
        }
        fieldfileA << sqrt(PSISUM) << ' ';
      }
      fieldfileA << std::endl;
      fieldfileA.close();
    }
  }
  else{ // export each field
    NFLOOP{
      std::string fFileName = std::to_string(NLnoise) + std::string("_") + std::to_string(NFIELDS) + std::string("_") + std::to_string(noisefiledirNo);
      fieldfileA.open(animationfileprefix + fFileName + std::string("_") + std::to_string(nf) + std::string(".dat"));
      fieldfileA << std::setprecision(10);

      for (int n = 0; n < phianimation.size(); n++) {
        for (int i=0; i<phianimation[n].size(); i++) {
          fieldfileA << phianimation[n][i][2*nf] << ' ';
        }
        fieldfileA << std::endl;
      }
      fieldfileA.close();
    }
  }
}


// Calculate power spectrum
void spectrum(std::vector<double> Ndata, int noisefiledirNo) {
  std::vector<std::vector<std::vector<std::complex<double>>>> Nmap3D = std::vector<std::vector<std::vector<std::complex<double>>>>(NLnoise, std::vector<std::vector<std::complex<double>>>(NLnoise, std::vector<std::complex<double>>(NLnoise, 0)));

  for (int i=0; i<NLnoise*NLnoise*NLnoise; i++) {
    int x=i/NLnoise/NLnoise ,y=(i%(NLnoise*NLnoise))/NLnoise, z=i%NLnoise;
    if(EoI_noise) {
      Nmap3D[x][y][z] = Ndata[i];
    }
    else {
      Nmap3D[x][y][z] = Ndata[i];
    }
  }
  std::vector<std::vector<std::vector<std::complex<double>>>> Nk=fft_fftw(Nmap3D);
  
  powsfile.open(powsfileprefix + std::string(".dat"), std::ios::app);
  powsfile << std::setprecision(10);
  int imax = ceil(log(NLnoise/2)/dlogn);
  std::vector<double> disc_power(imax, 0);

  LOOPLONG{
    int nxt, nyt, nzt; // shifted index
    if (i<=NLnoise/2) {
      nxt = i;
    } else {
      nxt = i-NLnoise;
    }

    if (j<=NLnoise/2) {
      nyt = j;
    } else {
      nyt = j-NLnoise;
    }

    if (k<=NLnoise/2) {
      nzt = k;
    } else {
      nzt = k-NLnoise;
    }
    
    double rk=nxt*nxt+nyt*nyt+nzt*nzt;

    double LogNk = log(sqrt(rk));
    double calPk = norm(Nk[i][j][k])/NLnoise/NLnoise/NLnoise/NLnoise/NLnoise/NLnoise;
    for (size_t ii = 0; ii < imax; ii++) {
      if ((dlogn*ii-LogNk)<dlogn/2. && -dlogn/2.<=(dlogn*ii-LogNk) && rk!=0) {
        disc_power[ii] += calPk/dlogn;
        break;
      }
    }
  }
  powsfile << noisefiledirNo << " ";
  for (size_t ii = 0; ii < imax; ii++) {
    powsfile << disc_power[ii] << " " ;
  }
  powsfile << std::endl;
  powsfile.close();
  std::cout << "ExportPowerSpectrum" << std::endl;
}




#endif