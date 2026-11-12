#ifndef INCLUDED_zoom_hpp_
#define INCLUDED_zoom_hpp_

std::array<int,3> shift{};

inline int indexsmall(int i, int j, int k) {
  return i * NLnoise * NLnoise/4 + j * NLnoise/2 + k;
}

// zoom
double LinearInterplation(std::vector<std::vector<std::vector<double>>> LatticeField, int x, int y, int z){
  int xp = (x == NLnoise - 1) ? x : x + 1;
  int xm = (x == 0) ? x : x - 1;
  int yp = (y == NLnoise - 1) ? y : y + 1;
  int ym = (y == 0) ? y : y - 1;
  int zp = (z == NLnoise - 1) ? z : z + 1;
  int zm = (z == 0) ? z : z - 1;
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

  LOOP{
    int idxsum = (i < NLnoise/2 ? 0 : 1) + (j < NLnoise/2 ? 0 : 2) + (k < NLnoise/2 ? 0 : 4);
    int idxf = index(i,j,k);
    NFLOOP{
      Nblock[idxsum] += Ndata[idxf];
    }
  }
  
  return std::distance(Nblock.begin(), std::max_element(Nblock.begin(), Nblock.end()));
}


// std::vector<int> 
void findNMaxBox(std::array<double,NLnoiseAll>& Ndata){
  int maxNpoint = std::distance(Ndata.begin(), std::max_element(Ndata.begin(), Ndata.end()));

  int x = maxNpoint/NLnoise/NLnoise, y = (maxNpoint%(NLnoise*NLnoise))/NLnoise, z = maxNpoint%NLnoise;
  int xb = x - NLnoise/4, yb = y - NLnoise/4, zb = z - NLnoise/4;

  if(x > NLnoise/2) xb = NLnoise/2;
  if(xb<0) xb = 0;
  if(y > NLnoise/2) yb = NLnoise/2;
  if(yb<0) yb = 0;
  if(z > NLnoise/2) zb = NLnoise/2;
  if(zb<0) zb = 0;

  shift={xb,yb,zb};
}


// std::vector<std::vector<std::vector<std::vector<double>>>> LatticeField(2*NFIELDS+2,std::vector<std::vector<std::vector<double>>>(NLnoise/2,std::vector<std::vector<double>>(NLnoise/2,std::vector<double>(NLnoise/2,0.))));

void InterpolatingPhi(std::array<int,3> shift){
  int NLhalf = NLnoise/2;

  std::vector<std::vector<std::vector<std::vector<double>>>> LatticeField(2*NFIELDS+2,std::vector<std::vector<std::vector<double>>>(NLhalf,std::vector<std::vector<double>>(NLhalf,std::vector<double>(NLhalf,0.))));

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int x = 0; x < NLhalf; x++) {
    for (int y = 0; y < NLhalf; y++) {
      for (int z = 0; z < NLhalf; z++) {
        int idx = index(x+shift[0],y+shift[1],z+shift[2]);
        for (int nf = 0; nf < 2*NFIELDS+2; nf++) {
          LatticeField[nf][x][y][z] = Phidata[idx][nf];
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
          phievol[idx][nf] = LinearInterplation(LatticeField[nf],x,y,z);
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
}


#endif
