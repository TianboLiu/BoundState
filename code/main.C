#include "bound.h"

using namespace std;

int main(){
  Lambda = 0;
  bool Prepare = true;
  if (Prepare){
    //Get wave function and potential grid, and setup interpolation
    readgrid();
    //Wave function renormalization
    urone();
    //Generate FQ grid
    makeFQ();
  }

  double Pd[2] = {0.3, M_PI/6.0};

  for (double th = 0.0; th < M_PI; th += M_PI/100){
    Pd[1] = th;
    cout << th << " " << dsigma(Pd) << endl;
  }

  return 0;
}


