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

  makeDS();
  

  return 0;
}


