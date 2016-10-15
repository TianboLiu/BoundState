#include "phifit.h"

using namespace std;

int main(int argc, char * argv[]){
  double s = atof(argv[1]);
  creatset("dataset.dat", "clasdb_E63M1.txt", s);
  addset("dataset.dat", "clasdb_E63M2.txt", s);

  //plotdata("dataset.dat");
  f0Fit();

  return 0;
}
