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

  bool Parallel = true;
  if (Parallel){
    FILE * fp;
    fp = fopen(fds, "w");
    double Pd[2];
    double ds;
    cout << "Generating dsigma grid Parallel..." << endl;
#pragma omp parallel for
    for (int x = 0; x < 70; x++){
      for (int y = 0; y <= 50; y++){
	Pd[0] = 0.02 * x;
	Pd[1] = 0.01 * M_PI * y;
	ds = dsigma(Pd);
	cout << Pd[0] << "   " << Pd[1] << "  " << cos(Pd[1]) << "   " << ds <<endl;
	fprintf(fp, "%.6E  %.6E  %.6E\n", Pd[0], cos(Pd[1]), ds);
      }
    }
    fclose(fp);
  }
  else{
    makeDS(fds);
  }
  return 0;
}


