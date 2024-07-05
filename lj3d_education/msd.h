#include <fstream>

class MSD {
  private:
    double ***r;
    int na, nframe, nbin;
    double dt;
    std::ifstream fin;
  public:
    double *msd;
    MSD(std::string, int, int, double);
    ~MSD();
    void run();
    void save(std::string);
};
