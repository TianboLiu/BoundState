#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>

#include "LHAPDF/LHAPDF.h"

using namespace std;

extern "C"{
  extern struct { int FINI;} fraginid_;
}
extern "C"{
  void fdss_(int * IH, int * IC, int * IO, double * Z, double * Q2, double * zU, double * zUBAR, double * zD, double * zDBAR, double * zS, double * zSBAR, double * zC, double * zB, double * zG);
}

int main(int argc, char * argv[]){

  if (argc < 3) return 0;

  const LHAPDF::PDF * ffs = LHAPDF::mkPDF("DSSFFlo", 12);

  double z = atof(argv[1]);
  double Q2 = atof(argv[2]);

  int ih = 2;
  int ic = 0;
  int io = 0;
  double u, ubar, d, dbar, s, sbar, c, b, g;

  fdss_ (&ih, &ic, &io, &z, &Q2, &u, &ubar, &d, &dbar, &s, &sbar, &c, &b, &g);
  
  printf("%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
	 b, c, sbar, ubar, dbar, d, u, s, c, b, g);
  
  printf("%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
	 ffs->xfxQ2(-5, z, Q2),
	 ffs->xfxQ2(-4, z, Q2),
	 ffs->xfxQ2(-3, z, Q2),
	 ffs->xfxQ2(-2, z, Q2),
	 ffs->xfxQ2(-1, z, Q2),
	 ffs->xfxQ2(1, z, Q2),
	 ffs->xfxQ2(2, z, Q2),
	 ffs->xfxQ2(3, z, Q2),
	 ffs->xfxQ2(4, z, Q2),
	 ffs->xfxQ2(5, z, Q2),
	 ffs->xfxQ2(0, z, Q2));

  return 0;
}
