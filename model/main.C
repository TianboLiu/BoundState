#include "model.h"

using namespace std;

int main(){
  readgrid();

  Mnuclear = MN;
  Mmeson = Mphi;
  Mreduced = Mnuclear * Mmeson / (Mnuclear + Mmeson);
  
  if (false){//write potential file
    FILE * fw;
    fw = fopen("a.dat","w");
    for (double x = 0; x < 18.0; x +=0.01){
      fprintf(fw, "%E  %E\n", x, Veff(x));
    }
    fclose(fw);
  }

  if (true){//shooting Schroedinger equation
    for (double E0 = 0.0; E0 <= 0.02; E0+=1.0e-3) 
      shooting(E0, 10.0);
  }

  if (false){
    shooting(-0.016202, 800.0);
  }

  allfree();
  return 0;
}
