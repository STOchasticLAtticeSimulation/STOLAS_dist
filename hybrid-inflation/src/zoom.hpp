#ifndef INCLUDED_zoom_hpp_
#define INCLUDED_zoom_hpp_


inline int index(int i, int j, int k) {
  return i * NLnoise * NLnoise + j * NLnoise + k;
}

inline int indexsmall(int i, int j, int k) {
  return i * NLnoise * NLnoise/4 + j * NLnoise/2 + k;
}

inline int INCREMENT(int i) {
  // return (i == NLnoise - 1) ? 0 : i + 1;
  return (i == NLnoise - 1) ? i : i + 1;
}

inline int DECREMENT(int i) {
  // return (i == 0) ? NLnoise - 1 : i - 1;
  return i - 1;
}

bool compareByFirstColumn(const std::vector<double>& a, const std::vector<double>& b) {
  return a[0] < b[0];
}


double LinearInterplation(std::vector<std::vector<std::vector<double>>> LatticeField, int x, int y, int z){
  int xp = INCREMENT(x);
  int xm = DECREMENT(x);
  int yp = INCREMENT(y);
  int ym = DECREMENT(y);
  int zp = INCREMENT(z);
  int zm = DECREMENT(z);
  int oddMask = (x % 2) * 4 + (y % 2) * 2 + (z % 2);

  switch (oddMask) {
    case 0:
      return LatticeField[x/2][y/2][z/2];
    
    case 4: // on the x-axis
      return 0.5*(LatticeField[xp/2][y/2][z/2] + LatticeField[xm/2][y/2][z/2]);
      
    case 2: // on the y-axis
      return 0.5*(LatticeField[x/2][yp/2][z/2] + LatticeField[x/2][ym/2][z/2]);
    
    case 1: // on the z-axis
      return 0.5*(LatticeField[x/2][y/2][zp/2] + LatticeField[x/2][y/2][zm/2]);

    case 6: // on the xy-plane
      return 0.25*( LatticeField[xp/2][yp/2][z/2] + LatticeField[xp/2][ym/2][z/2]
                  + LatticeField[xm/2][yp/2][z/2] + LatticeField[xm/2][ym/2][z/2]);

    case 5: // on the xz-plane
      return 0.25*( LatticeField[xp/2][y/2][zp/2] + LatticeField[xp/2][y/2][zm/2]
                  + LatticeField[xm/2][y/2][zp/2] + LatticeField[xm/2][y/2][zm/2]);

    case 3: // on the yz-plane
      return 0.25*( LatticeField[x/2][yp/2][zp/2] + LatticeField[x/2][yp/2][zm/2]
                   + LatticeField[x/2][ym/2][zp/2] + LatticeField[x/2][ym/2][zm/2]);
    
    case 7: // on the mid-point
      return 0.125*(LatticeField[xp/2][yp/2][zp/2] + LatticeField[xp/2][yp/2][zm/2]
                  + LatticeField[xp/2][ym/2][zp/2] + LatticeField[xp/2][ym/2][zm/2]
                  + LatticeField[xm/2][yp/2][zp/2] + LatticeField[xm/2][yp/2][zm/2]
                  + LatticeField[xm/2][ym/2][zp/2] + LatticeField[xm/2][ym/2][zm/2]);

    default:
      std::cout << "fail interpolation" << std::endl;
      return -1;
  }
}

int findMaxZeta(){
  std::vector<double> Nblock(8,0);

  LOOPLONG{
    int idxsum = (i < NLnoise/2 ? 0 : 1) + (j < NLnoise/2 ? 0 : 2) + (k < NLnoise/2 ? 0 : 4);
    int idxf = index(i,j,k);
    NFLOOP{
      Nblock[idxsum] += Ndata[idxf];
    }
  }
  
  return std::distance(Nblock.begin(), std::max_element(Nblock.begin(), Nblock.end()));
}


std::vector<int> findNMaxBox(std::vector<double> NData){
  int maxNpoint = std::distance(NData.begin(), std::max_element(NData.begin(), NData.end()));

  int x = maxNpoint/NLnoise/NLnoise, y = (maxNpoint%(NLnoise*NLnoise))/NLnoise, z = maxNpoint%NLnoise;
  int xb = x - NLnoise/4, yb = y - NLnoise/4, zb = z - NLnoise/4;

  if(x > NLnoise/2) xb = NLnoise/2;
  if(xb<0) xb = 0;
  if(y > NLnoise/2) yb = NLnoise/2;
  if(yb<0) yb = 0;
  if(z > NLnoise/2) zb = NLnoise/2;
  if(zb<0) zb = 0;

  std::vector<int> ShiftVector{xb,yb,zb};

  return ShiftVector;
}


std::vector<std::vector<double>> InterpolatingPhi(std::vector<int> Shift){
  int NLsmall = NLnoise/2;

  std::vector<std::vector<std::vector<std::vector<double>>>> LatticeField(2*NFIELDS+2,std::vector<std::vector<std::vector<double>>>(NLsmall,std::vector<std::vector<double>>(NLsmall,std::vector<double>(NLsmall,0.))));
  std::vector<std::vector<double>> phinew(NLnoise*NLnoise*NLnoise, std::vector<double>(2*NFIELDS+2,0));

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int x = 0; x < NLsmall; x++) {
    for (int y = 0; y < NLsmall; y++) {
      for (int z = 0; z < NLsmall; z++) {
        int idx = index(x+Shift[0],y+Shift[1],z+Shift[2]);
        for (int nf = 0; nf < 2*NFIELDS+2; nf++) {
          LatticeField[nf][x][y][z] = phievol[idx][nf];
        }
      }
    }
  }

  // interpolation
  int complete=0;
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int x = 0; x < NLnoise; x++) {
    for (int y = 0; y < NLnoise; y++) {
      for (int z = 0; z < NLnoise; z++) {
        int idx = index(x,y,z);
        for (int nf = 0; nf < 2*NFIELDS+2; nf++) {
          phinew[idx][nf] = LinearInterplation(LatticeField[nf],x,y,z);
        }
      }
    }
#ifdef _OPENMP
#pragma omp critical
#endif
    {
      complete++;
      std::cout << "\rInterpolation       : " << std::setw(3) << 100*complete/NLnoise << "%" << std::flush;
    }
  }
  std::cout << std::endl;

  return phinew;
}


std::ifstream fieldinitfile;
std::vector<std::vector<double>> fieldindata;
std::vector<std::vector<double>> getphifile(const int order, int NoisefiledirNo){
  std::vector<std::vector<double>> FieldinData;
  int NLsmall = NLnoise/2;
  std::string FieldFileName = fieldfileprefix + std::to_string(NLnoise) + "_" + std::to_string(NoisefiledirNo) + "_" + std::to_string(order-1) + std::string(".dat");
  std::cout << FieldFileName << std::endl;

  fieldinitfile.open(FieldFileName);
  if (fieldinitfile.fail()) std::cout << "field init file fail" << std::endl;

  double f1, f2, f3, f4, f5;
  while (fieldinitfile >> f1 >> f2 >> f3 >> f4 >> f5) {
    FieldinData.push_back({f1, f2, f3, f4, f5});
  }
  fieldinitfile.close();

  std::sort(FieldinData.begin(), FieldinData.end(), compareByFirstColumn);

  return FieldinData;
}


#endif
