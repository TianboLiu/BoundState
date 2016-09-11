#include "model.h"

using namespace std;

int main(){
  readgrid();

  Mnuclear = MN;
  Mmeson = Mphi;
  Mreduced = Mnuclear * Mmeson / (Mnuclear + Mmeson);
  NA = 12.0;

  // double a = one();
  // cout << a << endl;

  // FILE * fw;
  // fw = fopen("Veff.dat","w");
  // for (double r = 0.1; r < 15.0; r+=0.1){
  //   fprintf(fw, "%E  %E\n", r, ur_r_2(r) / (2.0 * Mreduced * ur(r)));
  // }
  // fclose(fw);



  // for (double E = -1.000; E <= -0.0; E+= 0.1){
  //   shooting(E, 30);
  // }

  // shooting(0.0415,30,true);

  return 0;
}



