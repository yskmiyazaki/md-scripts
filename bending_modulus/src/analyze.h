#ifndef ANALYZE_H
#define ANALYZE_H

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <string.h>
#include <string>
#include <vector>
#include <map>
#include "complex.h"

// constant
const int NDIM = 3;
const int X = 0;
const int Y = 1;
const int Z = 2;
const double PI = 3.1415926535;
const double DOUBLE_ERR = 1e-10;
// paramters
const int Nf = 32;  // grid size 2^n

class Average {
 private:
 public:
  int n;
  double ave;
  double ave2;
  double sum;
  double sum2;
  double std;

  Average() { init(); }
  ~Average() {}
  void init() 
  { 
    n = 0; 
    sum = sum2 = ave = ave2 = 0.0;
  }

  void add(const double x)
  {
    sum += x;
    sum2 += x*x;
    n++;
    ave = sum/n;
    ave2 = sum2/n;
    std = sqrt((ave2 - ave*ave)*n/(n - 1.0));
  }
};

class File {
 public:
  FILE *p;

  File() { p = nullptr; }
  File(const char *file_name, const char *mode) { p = nullptr; init(file_name, mode); }
  ~File() { if (p) fclose(p); }

  void init(const char *file_name, const char *mode)
  {
    if ((p = fopen(file_name, mode)) == NULL) {
      printf("Error: can't open file \"%s\" with \"%s\"\n", file_name, mode);
      exit(0);
    }
  }
};

class Atom {
 public:
  double r[NDIM];
  char atmname[5];
  char resname[5];

  Atom() {}
  ~Atom() {}
};

class Coordinate {
 public:
  Atom *atom;
  int natom;
  int istep;
  double box[6];

  Coordinate() {}
  ~Coordinate();
  int get_coordinate(const char *dcd_file_name);
  int get_atom_name(const char *psf_file_name);
};

class LipidBilayer : public Coordinate {
 private:
  Average av_lx;
  Average **av_qnq2;
  Complex **vnx;
  Complex **vny;
  Complex **vnHx;
  Complex **vnHy;
  Complex **vnLx;
  Complex **vnLy;
  int **count;

  void init();
  void get_orientation(Complex **vnLx, Complex **vnLy, const int layer);
  void calc_spectrum();
  void output(const char *out_file_name);
  void debug_output();
  void fft(Complex **a);

  double abs(double r[NDIM]) {
    return sqrt(r[X]*r[X] + r[Y]*r[Y] + r[Z]*r[Z]);
  }

  double shift_inbox(double x, const double box) {
    for (int i = 0; i < 1000; i++) {
      if (x < 0.0) x += box - DOUBLE_ERR;
      else if (box <= x) x -= box;
      else return x;
    }
    return 1e300; // for error
  }

  double* shift_inbox(double r[NDIM]) {
    for (int d = 0; d < NDIM; d++){
      for (int i = 0; i < 1000; i++) {
        if (r[d] < 0.0) r[d] += box[d] - DOUBLE_ERR;
        else if (box[d] <= r[d]) r[d] -= box[d];
        else break;
      }
    }
    return r;
  }

  double* shift_nearest(double r0[NDIM], const double r1[NDIM]) {
    // move r0 to where r0 is closer to r1
    for (int d = 0; d < NDIM; d++) {
      for (int i = 0; i < 1000; i++) {
        const double r10 = r1[d] - r0[d];
        if (-box[d]/2.0 <= r10 && r10 < box[d]/2.0) {
          break;
        } else if (r10 < -box[d]/2.0) {
          r0[d] -= box[d];
        } else {
          r0[d] += box[d] - DOUBLE_ERR;
        }
      }
    }
    return r0;
  }

 public:
  int nmol;
  std::string lipid_name;
  std::vector<std::string> heads, tails;
  std::map<std::string, std::vector<int>> idx_map;
  LipidBilayer(std::string lipid_name, std::vector<std::string> heads, std::vector<std::string> tails)
  { 
    vnx = nullptr; 
    this->lipid_name = lipid_name;
    this->heads = heads;
    this->tails = tails;
  }
  ~LipidBilayer() {}
  void do_analysis(const char *psf_file_name, const char *dcd_file_name, const char *out_file_name, int beg, int end);
};

#endif // ANALYZE_H
