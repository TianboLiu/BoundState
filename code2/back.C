#include "Lsimulation.h"

using namespace std;

int main(){
  double E = 1.019;
  double bw[2] = {1.019455, 0.00426};
  //cout << BreitWigner(&E, bw) << endl;
  double Nw = BreitWignerNormalization(bw);
  cout << Nw << endl;


  return 0;
}
