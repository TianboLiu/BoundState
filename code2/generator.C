#include "bound.h"

using namespace std;

int main(int argc, char * argv[]){
  bool Prepare = true;
  if (Prepare){
    //Get wave function and potential grid, and setup interpolation
    readgrid("wf.dat", "Veff.dat");
    //Wave function renormalization
    urone();
    //cout << urone() << endl;
    //Generate or load FQ grid
    //makeFQ(fFQ);
    LoadFQ("FQ.dat");
  }
  Long64_t Nsim = (Long64_t) atof(argv[2]);
  genData(argv[1], Nsim);//create root file, Nsim

  return 0;
}


