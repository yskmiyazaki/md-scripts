#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <sstream>

#define DIM 3

// Define a strucuture for simulation
struct MDsys_t {
  int num, nstep;
  double kb, nfree;
  double t, dt, sig, eps, kin, pot;
  double m;
  double **r, **v, **f;
  double boxl;
  double temp, temp0, rho, rho0;
};

// Zeros all components of a 2D array (n x m)
void zero2d(double **a, const int n, const int m)
{
   for (int i = 0; i < n; ++i){
     for (int j = 0; j < m; ++j)
       a[i][j] = 0.0;
   }
}

// Periodic Boundary Condition
void pbc(MDsys_t *mdsys)
{
   for (int i = 0; i < mdsys->num; i++){
     for (int j = 0; j < DIM; j++)
       mdsys->r[i][j] -= mdsys->boxl*std::round(mdsys->r[i][j]/mdsys->boxl);
   }
}

// Minimumc Image Convention
double mic(double x, const double boxl)
{
  if (x > 0.5*boxl)
    x -= boxl;
  if (x < -0.5*boxl)
    x += boxl;
  return x;
}

// Calculate forces
void force(MDsys_t *mdsys)
{
   double pot;
   double dx, dy, dz, xi, xj, yi, yj, zi, zj;
   double c12, c6, rsq;
   double **r, **f;
   double r2inv, r6, fmag;
   c12 = 4.0*mdsys->eps*pow(mdsys->sig, 12.0);
   c6  = 4.0*mdsys->eps*pow(mdsys->sig, 6.0);
   pot = 0.0;
   zero2d(mdsys->f, mdsys->num, DIM);
   for (int i = 0; i < mdsys->num - 1; i++){
     xi = mdsys->r[i][0];
     yi = mdsys->r[i][1];
     zi = mdsys->r[i][2];
     for (int j = i + 1; j < mdsys->num; j++){
       xj = mdsys->r[j][0];
       yj = mdsys->r[j][1];
       zj = mdsys->r[j][2];
       dx = mic(xi - xj, mdsys->boxl);
       dy = mic(yi - yj, mdsys->boxl);
       dz = mic(zi - zj, mdsys->boxl);
       rsq = dx*dx + dy*dy + dz*dz;
       if (rsq > 0.5*mdsys->boxl * 0.5*mdsys->boxl)
         continue;
       r2inv = 1.0/rsq;
       r6 = r2inv*r2inv*r2inv;
       fmag = (12.0*c12*r6 - 6.0*c6)*r6*r2inv;
       pot += (c12*r6 - c6)*r6;
       mdsys->f[i][0] += fmag*dx;
       mdsys->f[i][1] += fmag*dy;
       mdsys->f[i][2] += fmag*dz;
       mdsys->f[j][0] -= fmag*dx;
       mdsys->f[j][1] -= fmag*dy;
       mdsys->f[j][2] -= fmag*dz;
     }
   }
   mdsys->pot = pot;
}

void compute_temp(MDsys_t *mdsys)
{
    mdsys->kin = 0.0;
    for (int i = 0; i < mdsys->num; i++){
      for (int j = 0; j < DIM; j++) {
        mdsys->kin += mdsys->v[i][j]*mdsys->v[i][j];
      }
    }
    mdsys->kin *= 0.5*mdsys->m;
    mdsys->temp = 2.0*mdsys->kin/(mdsys->nfree*mdsys->kb);
}

// Do velocity Verlet
void vverlet(MDsys_t *mdsys)
{
    for (int i = 0; i < mdsys->num; i++){
      for (int j = 0; j < DIM; j++){
        mdsys->v[i][j] += 0.5*mdsys->dt*mdsys->f[i][j]/mdsys->m;
        mdsys->r[i][j] += mdsys->dt*mdsys->v[i][j];
      }
    }
    pbc(mdsys);
    force(mdsys);
    for (int i = 0; i < mdsys->num; i++){
      for (int j = 0; j < DIM; j++) {
        mdsys->v[i][j] += 0.5*mdsys->dt*mdsys->f[i][j]/mdsys->m;
      }
    }
    compute_temp(mdsys);
}

// Initialize file for output
std::ofstream init_file(const std::string fname, const std::string header = "#")
{
  std::ofstream os(fname);
  if (!os.is_open()){
    std::cout << "ERROR: Cannot open" << fname << std::endl;
    exit(1);
  }
  os << std::fixed << std::setprecision(6) << std::scientific;
  os << header << std::endl;
  return os;
}

// Write MD data and current timestep
void write_out(const int istep, const MDsys_t *mdsys, std::ofstream &os)
{
  std::cout << "Step.\t" << istep << "/" << mdsys->nstep << "\r" << std::flush;
  os << istep << "\t" 
     << mdsys->pot << "\t" 
     << mdsys->kin << "\t" 
     << mdsys->pot + mdsys->kin  << "\t" 
     << mdsys->temp << "\t" 
     << mdsys->rho << std::endl;
}

// Write MD xyz file
void write_xyz(const int istep, const MDsys_t *mdsys, std::ofstream &os)
{
  os << mdsys->num << std::endl;
  os << istep << "\t" << mdsys->boxl << std::endl; 
  for (int i = 0; i < mdsys->num; ++i)
    os << "NAME" << "\t"
       << mdsys->r[i][0] << "\t" 
       << mdsys->r[i][1] << "\t" 
       << mdsys->r[i][2] << std::endl;
}

// Write MD xyz file
void write_rst(const int istep, const MDsys_t *mdsys, std::ofstream &os)
{
  os << mdsys->boxl << std::endl;
  for (int i = 0; i < mdsys->num; ++i)
    os << mdsys->r[i][0] << "\t" 
       << mdsys->r[i][1] << "\t" 
       << mdsys->r[i][2] << "\t" 
       << mdsys->v[i][0] << "\t" 
       << mdsys->v[i][1] << "\t" 
       << mdsys->v[i][2] << std::endl;
}

void read_ini(const std::string fname, MDsys_t *mdsys)
{
  std::ifstream is(fname, std::ios::in);
  std::string str;
  for (int i = 0; i < mdsys->num; ++i){
    std::getline(is, str);
    std::stringstream ss(str);
    ss >> mdsys->r[i][0] 
       >> mdsys->r[i][1] 
       >> mdsys->r[i][2];
  }
}

void read_rst(const std::string fname, MDsys_t *mdsys)
{
  std::ifstream is(fname, std::ios::in);
  std::string str;

  std::getline(is, str);
  std::stringstream ss(str);
  ss >> mdsys->boxl;
  for (int i = 0; i < mdsys->num; ++i){
    std::getline(is, str);
    std::stringstream ss(str);
    ss >> mdsys->r[i][0] 
       >> mdsys->r[i][1] 
       >> mdsys->r[i][2] 
       >> mdsys->v[i][0] 
       >> mdsys->v[i][1] 
       >> mdsys->v[i][2];
  }
}

// Initialize a MDsys_t structure 
void initialize(MDsys_t *mdsys)
{
   mdsys->r = new double*[mdsys->num];
   mdsys->v = new double*[mdsys->num];
   mdsys->f = new double*[mdsys->num];
   for (int i = 0; i < mdsys->num; ++i){
     mdsys->r[i] = new double[DIM];
     mdsys->v[i] = new double[DIM];
     mdsys->f[i] = new double[DIM];
   }
   zero2d(mdsys->r, mdsys->num, DIM);
   zero2d(mdsys->v, mdsys->num, DIM);
}

// Finalize a MDsys_t structure 
void finalize(MDsys_t *mdsys)
{
  for (int i = 0; i < mdsys->num; ++i){
    delete [] mdsys->r[i];
    delete [] mdsys->v[i];
    delete [] mdsys->f[i];
  }
  delete [] mdsys->r;
  delete [] mdsys->v;
  delete [] mdsys->f;
}

void length_scaling(MDsys_t *mdsys)
{
   double alpha = std::pow(mdsys->rho / mdsys->rho0, 1.0/3.0);
   mdsys->boxl *= alpha;
   for (int i = 0; i < mdsys->num; ++i){
     for (int j = 0; j < DIM; ++j)
       mdsys->r[i][j] *= alpha;
   }
   mdsys->rho = mdsys->num / (mdsys->boxl*mdsys->boxl*mdsys->boxl);
}

void velocity_scaling(MDsys_t *mdsys)
{
   double alpha = std::sqrt(mdsys->temp0 / mdsys->temp);
   for (int i = 0; i < mdsys->num; ++i){
     for (int j = 0; j < DIM; ++j)
       mdsys->v[i][j] *= alpha;
   }
}

int main()
{
  // Declare a structure variable
  MDsys_t mdsys;

  // Set simulation conditions
  bool scale_b = false;
  std::string str;
  const int nfreq = 100;
  const int vfreq = 10;
  mdsys.num   = 256; 
  mdsys.dt    = 0.001; 
  mdsys.nstep = 10000;
  mdsys.kb    = 1.0; 
  mdsys.sig   = 0.5; 
  mdsys.eps   = 0.5;
  mdsys.m     = 1.0;
  mdsys.nfree = 3.0*mdsys.num;
  mdsys.boxl  = 4.0;

  // Open an output file
  std::string outf = "monitor.dat";
  std::string header = "#t\tU\tK\tH\tT\tr";
  std::ofstream output = init_file(outf, header);
  std::ofstream outxyz("coord.xyz", std::ios::out);

  // Initialize and set positions
  initialize(&mdsys);

  std::cin >> str;
  if (str == "restart")
    read_rst("start.dat", &mdsys);
  else
    read_ini("input256.txt", &mdsys);
  std::cin >> str;
  if (str == "scale"){
    scale_b = true;
    std::cin >> mdsys.temp0;
    std::cin >> mdsys.rho0;
  }
  // Calculate rho
  mdsys.rho = mdsys.num / (mdsys.boxl*mdsys.boxl*mdsys.boxl);
  // Calculate forces 
  pbc(&mdsys);
  if (scale_b)
      length_scaling(&mdsys);
  force(&mdsys);
  compute_temp(&mdsys);

  // Simulation loop
  int istep = 0;
  write_out(istep, &mdsys, output);
  while (istep < mdsys.nstep){
    vverlet(&mdsys);
    istep++;
    if (istep%vfreq == 0 && scale_b)
        velocity_scaling(&mdsys);
    if (istep%nfreq == 0){
      // Output current timestep and write simulation data
      write_out(istep, &mdsys, output);
      write_xyz(istep, &mdsys, outxyz);
    }
  }

  // Finalize and close the output file
  std::ofstream outrst("restart.dat", std::ios::out);
  write_rst(istep, &mdsys, outrst);
  finalize(&mdsys);
  output.close();
  outxyz.close();
  outrst.close();

  return 0;
}
