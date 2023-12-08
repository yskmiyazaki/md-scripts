#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "analyze.h"

//#define USE_MOLFILE_PLUGIN
#define USE_MY_READ_DCD

#ifdef USE_MY_READ_DCD

int Coordinate::get_atom_name(const char *psf_file_name)
{
  File psffile;
  psffile.init(psf_file_name, "rb");
  const int Nbuf = 1024;
  char buf[Nbuf];

  // read !NATOM
  while (fgets(buf, Nbuf, psffile.p)) {
    if (strstr(buf, "!NATOM") != NULL) break;
  }
  sscanf(buf, "%d", &natom);
  atom = new Atom[natom];

  // get atom information
  char dummy1[20], dummy3[20], dummy4[20];
  for (int i = 0; i < natom; i++) {
    if (fgets(buf, Nbuf, psffile.p) == NULL) {
      printf("ERROR: get_atom_name() fgets()\n"); 
      return false;
    }
    sscanf(buf, "%s %s %s %s %s", dummy1, atom[i].resname, dummy3, dummy4, atom[i].atmname);
  }
  return true;
}

int Coordinate::get_coordinate(const char *dcd_file_name)
{
  size_t count;
  static File dcdfile;
  if (dcdfile.p == nullptr) {
    dcdfile.init(dcd_file_name, "rb");
    const int pos_natom = 68;
    int buf[pos_natom];
    count = fread(buf, sizeof(int), pos_natom-1, dcdfile.p);
    int natom_;
    count = fread(&natom_, sizeof(int), 1, dcdfile.p);
    if (natom_ != natom) { 
      printf("NAtom=%d in out.dcd is different NAtom=%d in out.psf\n", natom_, natom); 
      exit(0);
    }
    count = fread(buf, sizeof(int), 1, dcdfile.p);
    istep = 0;
  }

  int buf[2];
  double d[7];
  if (fread(buf, sizeof(int), 1, dcdfile.p) == 0) return false;
  if (fread(d, sizeof(double), 7, dcdfile.p) == 0) return false;
  box[0] = d[0];
  box[1] = d[2];
  box[2] = d[5];

  float f;
  for (int i = 0; i < natom; i++) {
    if (fread(&f, sizeof(float), 1, dcdfile.p) == 0) return false;
    atom[i].r[X] = f;
  }
  if (fread(buf, sizeof(int), 2, dcdfile.p) == 0) return false;

  for (int i = 0; i < natom; i++) {
    if (fread(&f, sizeof(float), 1, dcdfile.p) == 0) return false;
    atom[i].r[Y] = f;
  }
  if (fread(buf, sizeof(int), 2, dcdfile.p) == 0) return false;

  for (int i = 0; i < natom; i++) {
    if (fread(&f, sizeof(float), 1, dcdfile.p) == 0) return false;
    atom[i].r[Z] = f;
  }
  if (fread(buf, sizeof(int), 1, dcdfile.p) == 0) return false;

  istep++;
  return true;
}

Coordinate::~Coordinate()
{
  if (atom) delete [] atom;
}

#else 

typedef int int4;
void f77_molfile_open_read(int4 *handle, int4 *NAtoms, const char *infile, const char *intype, const int len_if, const int len_it);
void f77_molfile_init();
void f77_molfile_read_next(int4 *handle, int4 *NAtoms, float *xyz, float *box, int4 *status);
void f77_molfile_close_read(int4 *handle);
void f77_molfile_finish();
int status;
int i, j, handle;

Coordinate::Coordinate()
{
  char *file_name = "out.dcd";
  char *intype = "auto";
  f77_molfile_init();
  f77_molfile_open_read(&handle, &natom, file_name, intype, strlen(file_name), strlen(intype));
  if (handle < 0) {
    printf("ERROR: cannot open file \"%s\"\n", file_name);
  } else {
    printf("File opened successfully\n");
    printf("NATOM = %d\n", natom);
  }
  xyz = new float[3*natom];
  istep = 0;
  status = 1;
}

Coordinate::~Coordinate()
{
  f77_molfile_close_read(&Handle);
  f77_molfile_finish();
  if (xyz) delete [] xyz;
}

int Coordinate::get_coordinate()
{
  f77_molfile_read_next(&Handle, &NAtom, xyz, box, &status);
  if (status == false) return false;
  istep++;
  return true;
}

#endif 
