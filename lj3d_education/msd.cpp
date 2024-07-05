#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <sstream>
#include "msd.h"

#define PI 3.1415926
#define DIM 3

MSD::MSD(std::string xyz, int na, int nframe, double dt){
  this->na = na;
  this->nframe = nframe;
  this->dt = dt;
  fin.open(xyz, std::ios::in);
  r = new double**[na];
  for (int i = 0; i < na; ++i) {
    r[i] = new double*[DIM];
    for (int j = 0; j < DIM; ++j)
      r[i][j] = new double[nframe];
  }
  msd = new double[nframe];
  for (int i = 0; i < nframe; ++i)
    msd[i] = 0.0;
}
MSD::~MSD(){
  delete [] msd;
  for (int i = 0; i < na; ++i) {
    for (int j = 0; j < DIM; ++j)
      delete [] r[i][j];
    delete [] r[i];
  }
  delete [] r;
  fin.close();
}
// do MSD calculation
void MSD::run(){
  std::string str, dummy;
  int istep; 
  double boxl;
  double rii[3];
  for (int k = 0; k < nframe; ++k){
    // read 1st line in a xyz file
    std::getline(fin, str);
    std::stringstream s1(str);
    s1 >> dummy;
    // read 2nd line
    std::getline(fin, str);
    std::stringstream s2(str);
    s2 >> istep >> boxl;
    // read xyz coordinate 
    for (int i = 0; i < na; ++i){
      std::getline(fin, str);
      std::stringstream ss(str);
      ss << str;
      ss >> dummy >> r[i][0][k] >> r[i][1][k] >> r[i][2][k];
    }
  }
  for (int k = 0; k < nframe; ++k) {
    for (int l = k; l < nframe; ++l) {
      for (int i = 0; i < na; ++i) {
        // |ri(l) - ri(k)|^2 calculation
        double rsq = 0.0;
        for (int d = 0; d < DIM; d++) {
          rii[d] = r[i][d][l] - r[i][d][k];
          // minimum image convention
          if (rii[d] > 0.5*boxl)
            rii[d] -= boxl;
          if (rii[d] < -0.5*boxl)
            rii[d] += boxl;
          rsq += rii[d]*rii[d];
        }
        msd[l-k] += rsq;
      }
    }
  }
  // average msd over # of particles and frames
  for (int k = 0; k < nframe; ++k) {
    msd[k] /= double(na)*double(nframe-k);
  }
}
// save the result in a text file
void MSD::save(std::string fname){
  std::ofstream fout(fname);
  for (int i = 0; i < nframe; i++)
    fout << i*dt << "\t" << msd[i] <<std::endl;
  fout.close();
}
