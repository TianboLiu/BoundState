#include "bound.h"

using namespace std;

int main(){
  Ebeam = 1.45;

  bool Prepare = false;
  if (Prepare){
    //Get wave function and potential grid, and setup interpolation
    readgrid();
    //Wave function renormalization
    urone();
    //Generate FQ grid
    makeFQ();
  }

  //makeDS("ds2.dat");
  //double Pd[2] = {0.3, M_PI/6.0};
  //cout << dsigma(Pd) << endl;

  return 0;
}


