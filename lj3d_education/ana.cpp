#include "rdf.h"
#include <string>

int main(){
  std::string xyzf, outf;
  int nbin = 75;
  int na = 256;
  int nframe = 100;
  double rmax;

  rmax = 1.5;
  xyzf = "liq.xyz";
  outf = "rdf_liq.dat";
  RDF rdf_l(xyzf, na, nframe, nbin, rmax);  
  rdf_l.run();
  rdf_l.save(outf);

  rmax = 15;
  xyzf = "gas.xyz";
  outf = "rdf_gas.dat";
  RDF rdf_g(xyzf, na, nframe, nbin, rmax);  
  rdf_g.run();
  rdf_g.save(outf);

  return 0;
}
