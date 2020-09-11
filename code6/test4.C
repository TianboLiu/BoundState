#include "Lcore.h"

int main(){

  JPSIPomLQCD::SetModel();

  cout << JPSIPomLQCD::dSigmaJpsi(4.636, 0) * pow(0.197327, 2) / 1e-7 << endl;

  return 0;
}
