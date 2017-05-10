#include "Lcore.h"

using namespace MODEL;

int main(){

  SetMODEL();

  cout << TF_fMass.Integral(TF_fMass.GetXmin(), TF_fMass.GetXmax()) << endl;

  return 0;
}
