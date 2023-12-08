#include "analyze.h"

#define min(a,b)  (((a)<(b))?(a):(b))
#define max(a,b)  (((a)>(b))?(a):(b))
void cdft(int n, int isgn, double *a, int *ip, double *w);

void LipidBilayer::init()
{
  // lipid bilayer z-coordinate(rx,ry)
  if (vnx == nullptr) {
    count = new int* [Nf];
    av_qnq2 = new Average* [Nf];
    vnx = new Complex* [Nf];
    vny = new Complex* [Nf];
    vnHx = new Complex* [Nf];
    vnHy = new Complex* [Nf];
    vnLx = new Complex* [Nf];
    vnLy = new Complex* [Nf];
    for (int i = 0; i < Nf; i++) {
      count[i] = new int [Nf];
      av_qnq2[i] = new Average [Nf];
      vnx[i] = new Complex [Nf];
      vny[i] = new Complex [Nf];
      vnHx[i] = new Complex [Nf];
      vnHy[i] = new Complex [Nf];
      vnLx[i] = new Complex [Nf];
      vnLy[i] = new Complex [Nf];
    }
  }
  // zero variables
  for (int i = 0; i < Nf; i++) {
    for (int j = 0; j < Nf; j++) {
      av_qnq2[i][j].init();
    }
  }
  for (int i = 0; i < natom; i++) {
    if (strcmp(atom[i].resname, lipid_name.c_str()) != 0) continue;
    for (auto head : heads)
      if (strcmp(atom[i].atmname, head.c_str()) == 0) idx_map[head].push_back(i);
    for (auto tail : tails)
      if (strcmp(atom[i].atmname, tail.c_str()) == 0) idx_map[tail].push_back(i);
  }
  nmol = idx_map[heads[0]].size();
  for (auto head : heads) {
    if (idx_map[head].size() != nmol) {
      printf("ERROR: Selected head names are not unique.\n");
      return;
    }
  }
  for (auto tail : tails) {
    if (idx_map[tail].size() != nmol) {
      printf("ERROR: Selected tail names are not unique.\n");
      return;
    }
  }
}

void LipidBilayer::calc_spectrum()
{
  // get vnx, vny
  get_orientation(vnHx, vnHy, 0);
  get_orientation(vnLx, vnLy, 1);
  for (int iy = 0; iy < Nf; iy++) {
    for (int ix = 0; ix < Nf; ix++) {
      vnx[iy][ix] = 0.5*(vnHx[iy][ix] - vnLx[iy][ix]);
      vny[iy][ix] = 0.5*(vnHy[iy][ix] - vnLy[iy][ix]);
    }
  }
  // FFT
  fft(vnx);
  fft(vny);
  // <|nq//|^2> = <|1/|q|(q.n)|^2> = <|(q.n)/|^2 1/q^2>
  const double l = box[X];
  for (int iy = 0; iy < Nf; iy++) {
    for (int ix = 0; ix < Nf; ix++) {
      if (ix == 0 && iy == 0) continue;
      double qx, qy;
      if (ix < Nf/2) {
        qx = 2.0*PI/l*ix;
      } else {
        qx = 2.0*PI/l*(ix - Nf);
      }
      if (iy < Nf/2) {
        qy = 2.0*PI/l*iy;
      } else {
        qy = 2.0*PI/l*(iy - Nf);
      }
      Complex qn = qx*vnx[iy][ix] + qy*vny[iy][ix];
      av_qnq2[iy][ix].add( (qn.real*qn.real + qn.imag*qn.imag) / (qx*qx + qy*qy) );
    }
  }
}

void LipidBilayer::get_orientation(Complex **vnLx, Complex **vnLy, const int layer)
{
  // set 0
  for (int i = 0; i < Nf; i++) {
    for (int j = 0; j < Nf; j++) {
      vnLx[i][j].real = vnLx[i][j].imag = 0.0;
      vnLy[i][j].real = vnLy[i][j].imag = 0.0;
      count[i][j] = 0;
    }
  }
  // vnx[iy][ix].real averaged z-position
  for (int i = 0; i < nmol; i++) {
    double vn[NDIM];
    double vn_h[NDIM] = {0., 0., 0.};
    double vn_t[NDIM] = {0., 0., 0.};
    double *r_ref = shift_inbox(atom[idx_map[heads[0]][i]].r);

    for (auto head : heads) {
      int idx = idx_map[head][i];
      for (int d = 0; d < NDIM; d++) {
        vn_h[d] += shift_nearest(atom[idx].r, r_ref)[d];
      }
    }
    for (auto tail : tails) {
      int idx = idx_map[tail][i];
      for (int d = 0; d < NDIM; d++) {
        vn_t[d] += shift_nearest(atom[idx].r, r_ref)[d];
      }
    }
    for (int d = 0; d < NDIM; d++) { 
      vn_h[d] /= heads.size();
      vn_t[d] /= tails.size();
    }
    for (int d = 0; d < NDIM; d++) { 
      vn[d] = vn_h[d] - vn_t[d];
    }
    double absvn = abs(vn);
    for (int d = 0; d < NDIM; d++) {
      vn[d] /= absvn;
    }
    int ix = max( 0, min( Nf - 1, (int)(Nf*vn_h[X]/box[X]) ) );
    int iy = max( 0, min( Nf - 1, (int)(Nf*vn_h[Y]/box[Y]) ) );
    int layer_wk = vn[Z] > 0 ? 0 : 1;
    if (layer_wk == layer) {
      vnLx[iy][ix].real += vn[X];
      vnLy[iy][ix].real += vn[Y];
      count[iy][ix]++;
    }
  }
  for (int iy = 0; iy < Nf; iy++) {
  	for (int ix = 0; ix < Nf; ix++) {
      if (0 < count[iy][ix]) {
        vnLx[iy][ix].real /= count[iy][ix];
        vnLy[iy][ix].real /= count[iy][ix];
      }
  	}
  }
  // interpolate zero values with neighbors
  for (int i = 0; i < 2; i++) {
  	for (int iy = 0; iy < Nf; iy++) { 	
  	  for (int ix = 0; ix < Nf; ix++) {
  	    if (count[iy][ix] == 0) {
  	      const int N = Nf;
  	      int n = 0;
  	      int nxp = count[iy][(ix + 1) % N];
  	      int nxm = count[iy][(ix + N - 1) % N];
  	      int nyp = count[(iy + 1) % N][ix];
  	      int nym = count[(iy + N - 1) % N][ix];
  	      double xxp = vnLx[iy][(ix + 1) % N].real;
  	      double xxm = vnLx[iy][(ix + N - 1) % N].real;
  	      double xyp = vnLx[(iy + 1) % N][ix].real;
  	      double xym = vnLx[(iy + N - 1) % N][ix].real;
  	      double yxp = vnLy[iy][(ix + 1) % N].real;
  	      double yxm = vnLy[iy][(ix + N - 1) % N].real;
  	      double yyp = vnLy[(iy + 1) % N][ix].real;
  	      double yym = vnLy[(iy + N - 1) % N][ix].real;
  	      if (nxp > 0) n++;
  	      if (nxm > 0) n++;
  	      if (nyp > 0) n++;
  	      if (nym > 0) n++;
  	      if (n > 0) {
  	        vnLx[iy][ix].real = (xxp + xxm + xyp + xym)/n;
  	        vnLy[iy][ix].real = (yxp + yxm + yyp + yym)/n;
  	        count[iy][ix] = 1;
  	      }
  	    }
  	  }
  	}
  }
}

void LipidBilayer::fft(Complex **a)
{
  // a(qx,qy) = 1/L ÅÁdxdy a(x,y)
  static int ip[Nf] = {0};
  static double w[Nf/2];
  // fft rx to kx
  for (int iy = 0; iy < Nf; iy++) {
    cdft(2*Nf, 1, (double*)a[iy], ip, w);
  }
  // fft ry to ky
  for (int ix = 0; ix < Nf; ix++) {
    Complex temp[Nf];
    for (int iy = 0; iy < Nf; iy++) {
      temp[iy] = a[iy][ix];
    }
    cdft(2*Nf, 1, (double*)temp, ip, w);
    for (int iy = 0; iy < Nf; iy++) {
      a[iy][ix] = temp[iy];
    }
  }
  // factor
  for (int iy = 0; iy < Nf; iy++) {
    for (int ix = 0; ix < Nf; ix++) {
      a[iy][ix] *= box[X]/(Nf*Nf);
    }
  }
}

void LipidBilayer::debug_output()
{
  File fo("out.debug", "w+");
  for (int iy = 0; iy < Nf; iy++) {
    for (int ix = 0; ix < Nf; ix++) {
      fprintf(fo.p, "%f	", vnHx[iy][ix].real);
    }
    fprintf(fo.p, "\n");
  }
}

void LipidBilayer::output(const char *out_file_name)
{
  File fo(out_file_name, "w+");
  fprintf(fo.p, "#|q|=2pi|n|/L\t<|nq//|^2>\tstd<|nq//|^2>\n");
  for (int i = 0; i < Nf/2; i++) {
    for (int j = i; j < Nf/2; j++) {
      if (i == 0 && j == 0) continue;
      fprintf(fo.p, "%f\t%e\t%e\n", sqrt( (double)(i*i + j*j) )*2.0*PI/av_lx.ave,
        ( av_qnq2[i][j].ave+av_qnq2[i][(Nf - j) % Nf].ave + av_qnq2[(Nf - i) % Nf][j].ave + av_qnq2[(Nf - i) % Nf][(Nf - j) % Nf].ave
        + av_qnq2[j][i].ave+av_qnq2[j][(Nf - i) % Nf].ave + av_qnq2[(Nf - j) % Nf][i].ave + av_qnq2[(Nf - j) % Nf][(Nf - i) % Nf].ave )/8.0,
        ( av_qnq2[i][j].std+av_qnq2[i][(Nf - j) % Nf].std + av_qnq2[(Nf - i) % Nf][j].std + av_qnq2[(Nf - i) % Nf][(Nf - j) % Nf].std
        + av_qnq2[j][i].std+av_qnq2[j][(Nf - i) % Nf].std + av_qnq2[(Nf - j) % Nf][i].std + av_qnq2[(Nf - j) % Nf][(Nf - i) % Nf].std )/8.0 
      );
    }
  }
}

void LipidBilayer::do_analysis(const char *psf_file_name, const char *dcd_file_name, const char *out_file_name, int beg, int end)
{
  get_atom_name(psf_file_name);
  init();
  av_lx.init();
  for (int i = 0; i <=end; i++) {
    if (i % 10 == 0) {
      printf("\r%d frame", i);
      fflush(stdout);
    }
    if (!get_coordinate(dcd_file_name)) break;
    if (i >= beg){
      calc_spectrum();
      av_lx.add(box[X]);
    }
  }
  output(out_file_name);
  printf("\nDone.\n");
}
