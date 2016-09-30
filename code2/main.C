#include "bound.h"

using namespace std;

int main(int argc, char * argv[]){
  //cout << argv[1] << " " << argv[2] << " " << argv[3] << endl;
  //Lambda = 2.0 * GeVfm;
  //b = 1.64;

  Ebeam = atof(argv[1]);//Set beam energy
  char * fFQ = argv[2];//Set FQ.dat file
  char * fds = argv[3];//Set ds.dat file

  char fwf[] = "wf.dat";
  char fVeff[] = "Veff.dat";

  bool Prepare = true;
  if (Prepare){
    //Get wave function and potential grid, and setup interpolation
    readgrid(fwf, fVeff);
    //Wave function renormalization
    urone();
    //cout << urone() << endl;
    //Generate or load FQ grid
    //makeFQ(fFQ);
    LoadFQ(fFQ);
  }

  bool Parallel = false;
  if (Parallel){
    makeDS_Parallel(fds);
  }
  else if (false){
    makeDS(fds);
  }

  bool Total = true;
  if (Total){
    //double Pd[2] = {0.3, M_PI/6};
    printf("%.3f     %.6f  nb\n", Ebeam, sigmaT()*3.89379e5);
  }
  return 0;
}


