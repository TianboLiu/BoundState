#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>

using namespace std;

extern "C"{
  extern struct { int FINI;} fraginid_;
}
extern "C"{
  void fdss_(int * IH, int * IC, int * IO, double * Z, double * Q2, double * zU, double * zUBAR, double * zD, double * zDBAR, double * zS, double * zSBAR, double * zC, double * zB, double * zG);
}

int main(){

  int ih = 4;
  int ic = -1;
  int io = 0;
  
  double u, ubar, d, dbar, s, sbar, c, b, g;
  double z, Q2;

  FILE * fp;
  fp = fopen("DSSFFlo_1001.dat", "w");
  fprintf(fp, "PdfType: central\n");
  fprintf(fp, "Format: lhagrid1\n");
  fprintf(fp, "---\n");

  for (int i = 0; i < 100; i++){
    fprintf(fp, "%.6e ", 0.01 + 0.01 * i);
  }
  fprintf(fp, "\n");
  for (int i = 0; i < 90; i++){
    fprintf(fp, "%.6e ", 1.0 + 0.1 * i);
  }
  fprintf(fp, "\n");
  fprintf(fp, "-5 -4 -3 -2 -1 1 2 3 4 5 21\n");

  for (int i = 0; i < 100; i++){
    for (int j = 0; j < 90; j++){
      z = 0.01 + 0.01 * i;
      Q2 = pow(1.0 + 0.1 * j, 2);
      fdss_ (&ih, &ic, &io, &z, &Q2, &u, &ubar, &d, &dbar, &s, &sbar, &c, &b, &g);
      if (b < 0) b = 0;
      if (c < 0) c = 0;
      fprintf(fp, "%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
	      b, c, sbar, ubar, dbar, d, u, s, c, b, g);
    }
  }
  fprintf(fp, "---\n");

  fclose(fp);

  return 0;
}
