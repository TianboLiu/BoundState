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
    double Pd[70][51][2];
    double ds[70][51];
    cout << "Generating dsigma grid Parallel..." << endl;
    #pragma omp parallel
    {
      #pragma omp for
      for (int i = 0; i < 70; i++){
	for (int j = 0; j <= 50; j++){
	  Pd[i][j][0] = 0.02 * i;
	  Pd[i][j][1] = 0.01 * M_PI * j;
	  ds[i][j] = dsigma(Pd[i][j]);
	  cout << Pd[i][j][0] << "   " << Pd[i][j][1] << "  " << cos(Pd[i][j][1]) << "   " << ds[i][j] <<endl;
	}
      }
    }
    cout << "Writing File ..." << endl;
    for (int i = 0; i < 70; i++){
      for (int j = 0; j <= 50; j++){
	fprintf(fp, "%.6E  %.6E  %.6E\n", Pd[i][j][0], cos(Pd[i][j][1]), ds[i][j]);
      }
    }   
    fclose(fp);
  }
  else if (false){
    makeDS(fds);
  }
  return 0;
}


