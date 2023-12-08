#include <fstream>

class RDF {
  private:
    int *hist;
    double **r;
    int na, nframe, nbin;
    double rmax;
    std::ifstream fin;
  public:
    double *rdf;
    double *bins;
    RDF(std::string, int, int, int, double);
    ~RDF();
    void run();
    void save(std::string);
};
