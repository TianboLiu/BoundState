#include "bound.h"

using namespace std;

int main(int argc, char * argv[]){
  cout << argv[1] << " " << argv[2] << " " << argv[3] << " " << argv[4] << " " << argv[5] << endl;

  Ebeam = atof(argv[1]);//Set beam energy
  char * fwf = argv[2];//Set wf.dat file
  char * fVeff = argv[3];//Set Veff.dat file
  char * fFQ = argv[4];//Set FQ.dat file
  char * fds = argv[5];//Set ds.dat file

  bool Prepare = true;
  if (Prepare){
    //Get wave function and potential grid, and setup interpolation
    readgrid(fwf, fVeff);
    //Wave function renormalization
    cout << urone() << endl;
    //Generate or load FQ grid
    //makeFQ(fFQ);
    LoadFQ(fFQ);
  }

  makeDS(fds);

  return 0;
}


