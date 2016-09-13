#include "bound.h"

using namespace std;

int main(){
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

  double one = sigma();
  cout << one << endl;

  return 0;
}


