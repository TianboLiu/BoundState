#include "Phifit.h"

using namespace std;

int main(){
  double Eslimit[2] = {1.98, 1.99};
  makeset("dataset.dat", Eslimit);
  plotdata("dataset.dat", Eslimit);

  return 0;
}
