#include <iostream>
#include <string>
#include "rdf.h"
#include "msd.h"

int main(){
  std::string xyzf, outf;
  int na = 256;
  int nframe = 100;
  int nbin = 75;
  double rmax;
  double dt = 1e-3;

  rmax = 1.5;
  xyzf = "liq.xyz";
  outf = "rdf_liq.dat";
  RDF rdf_l(xyzf, na, nframe, nbin, rmax);  
  rdf_l.run();
  rdf_l.save(outf);

  outf = "msd_liq.dat";
  MSD msd_l(xyzf, na, nframe, dt);
  msd_l.run();
  msd_l.save(outf);

  rmax = 15;
  xyzf = "gas.xyz";
  outf = "rdf_gas.dat";
  RDF rdf_g(xyzf, na, nframe, nbin, rmax);  
  rdf_g.run();
  rdf_g.save(outf);

  outf = "msd_gas.dat";
  MSD msd_g(xyzf, na, nframe, dt);
  msd_g.run();
  msd_g.save(outf);

  return 0;
}
