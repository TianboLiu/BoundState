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

  double one = FQ(0.05);

  cout << one << endl;

  return 0;
}
