#include "phifit.h"

using namespace std;

int main(int argc, char * argv[]){
  //double s = atof(argv[1]);
  double s = 1.985;
  for (int i = 0; i < 100; i++){
    s = 1.985 + 0.01 * i;
    creatset("dataset.dat", "clasdb_E63M1.txt", s);
    addset("dataset.dat", "clasdb_E63M2.txt", s);
    f0Fit();
    if (s > 2.835) break;
  }

  //plotdata("dataset.dat");
  //f0Fit();

  return 0;
}
