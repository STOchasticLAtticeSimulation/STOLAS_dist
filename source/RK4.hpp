#ifndef INCLUDED_RK4_
#define INCLUDED_RK4_

#define _USR_MATH_DEFINES
#include <cmath>
#include <functional>
#include "vec_op.hpp"

template <class T>
void RK4(std::function<T(double, T)> &dxdt, double &t, T &x, double dt) {
  T kx[4]; // 4-stage slopes kx
  double a[4][4],b[4],c[4]; // Butcher

  // -------------- initialise kx, a, b, c --------------- //
  for(int i=0;i<=3;i++){
    kx[i] = x;
    vec_op::init(kx[i]);
    
    for(int j=0;j<=3;j++){
      a[i][j]=0.;
    }
  }
  
  a[1][0]=1./2;  a[2][1]=1./2;  a[3][2]=1.;
  b[0]=1./6;     b[1]=1./3;     b[2]=1./3;    b[3]=1./6; 
  c[0]=0;        c[1]=1./2;     c[2]=1./2;    c[3]=1;
  // ----------------------------------------------------- //
  

  T X = x; // position at i-stage
    
  for(int i=0;i<=3;i++){
    X = x; // initialise X
    
    for(int j=0;j<=3;j++){
      X += dt * a[i][j] * kx[j];
      kx[i] = dxdt(t,X);
    }
  }

  t += dt;
  x += dt*(b[0]*kx[0] + b[1]*kx[1] + b[2]*kx[2] + b[3]*kx[3]);
}

template <class T>
void RK4(std::function<T(double, T, std::vector<double>)> &dxdt, double &t, T &x, double dt, std::vector<double> params) {
  T kx[4]; // 4-stage slopes kx
  double a[4][4],b[4],c[4]; // Butcher

  // -------------- initialise kx, a, b, c --------------- //
  for(int i=0;i<=3;i++){
    kx[i] = x;
    vec_op::init(kx[i]);
    
    for(int j=0;j<=3;j++){
      a[i][j]=0.;
    }
  }
  
  a[1][0]=1./2;  a[2][1]=1./2;  a[3][2]=1.;
  b[0]=1./6;     b[1]=1./3;     b[2]=1./3;    b[3]=1./6; 
  c[0]=0;        c[1]=1./2;     c[2]=1./2;    c[3]=1;
  // ----------------------------------------------------- //
  

  T X = x; // position at i-stage
    
  for(int i=0;i<=3;i++){
    X = x; // initialise X
    
    for(int j=0;j<=3;j++){
      X += dt * a[i][j] * kx[j];
      kx[i] = dxdt(t,X,params);
    }
  }

  t += dt;
  x += dt*(b[0]*kx[0] + b[1]*kx[1] + b[2]*kx[2] + b[3]*kx[3]);
}


#endif
