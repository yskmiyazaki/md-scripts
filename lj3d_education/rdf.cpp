#include "rdf.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <sstream>

#define PI 3.1415926
#define DIM 3

RDF::RDF(std::string xyz, int na, int nframe, int nbin = 40, double rmax = 2.0){
  this->na = na;
  this->nframe = nframe;
  this->nbin = nbin;
  this->rmax = rmax;
  fin.open(xyz, std::ios::in);
  r = new double*[na];
  for (int i = 0; i < na; ++i)
    r[i] = new double[DIM];
  rdf = new double[nbin];
  hist = new int[nbin];
  bins = new double[nbin];
  for (int i = 0; i < nbin; ++i)
    hist[i] = 0;
}
RDF::~RDF(){
  delete [] rdf;
  delete [] hist;
  delete [] bins;
  for (int i = 0; i < na; ++i)
    delete [] r[i];
  delete [] r;
  fin.close();
}
void RDF::run(){
  std::string str, dummy;
  int istep; 
  double boxl;
  double rij[3];
  double dr = rmax/nbin;
  double V = 0;
  for (int k = 0; k < nframe; ++k){
    std::getline(fin, str);
    std::stringstream s1(str);
    s1 >> dummy;
    std::getline(fin, str);
    std::stringstream s2(str);
    s2 >> istep >> boxl;
    V += boxl*boxl*boxl;
    for (int i = 0; i < na; ++i){
      std::getline(fin, str);
      std::stringstream ss(str);
      ss << str;
      ss >> dummy >> r[i][0] >> r[i][1] >> r[i][2];
    }
    for (int i = 0; i < na - 1; ++i){
      for (int j = i + 1; j < na; ++j){
        double rsq = 0.0;
        for (int k = 0; k < DIM; k++){
          rij[k] = r[i][k] - r[j][k];
          if (rij[k] > 0.5*boxl)
            rij[k] -= boxl;
          if (rij[k] < -0.5*boxl)
            rij[k] += boxl;
          rsq += rij[k]*rij[k];
        }
        double dist = std::sqrt(rsq);
        if (dist > rmax)
          continue;
        int ibin = int(dist/dr) + 1;
        if (ibin < nbin)
          hist[ibin] += 2;
      }
    }
  }
  V /= nframe; 
  double rho = na/V;
  for (int k = 0; k < nbin; ++k) {
    double ri = (k + 0.5)*dr;
    double coef = 1.0/(4.0*PI*ri*ri*dr*(na - 1)*rho);
    rdf[k] = coef*double(hist[k])/nframe;
    bins[k] = ri;
  }
}
void RDF::save(std::string fname){
  std::ofstream fout(fname);
  for (int i = 0; i < nbin; i++)
    fout << bins[i] << "\t" << rdf[i] <<std::endl;
  fout.close();
}
