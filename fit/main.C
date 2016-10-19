#include "phifit.h"

using namespace std;

int main(int argc, char * argv[]){
  //double fixE = atof(argv[1]);
  double s = 1.985;
  //double s = fixE;
  double par[5];
  if (false){
    double Chi2 = 0.0;
    FILE * fp = fopen("fitphi.txt", "w");
    fprintf(fp, "sqrt(s)\t a0\t a2\t b1\t\n");
    for (int i = 0; i < 100; i++){
      s = 1.985 + 0.01 * i;
      creatset("dataset.dat", "clasdb_E63M1.txt", s);
      addset("dataset.dat", "clasdb_E63M2.txt", s);
      Chi2 += f0Fit(par);
      //plotdata("dataset.dat", par);
      fprintf(fp, "%.3f\t %.4f\t %.4f\t %.4f\n", s, par[0], par[1], par[2]);
      if (s >= 2.835) break;
    }
    fclose(fp);
    cout << "Chi2 = " << Chi2 << endl;
  }
 
  if (false){
    ifstream infile("fitphi.txt");
    char tmp[256];
    infile.getline(tmp, 256);
    TString filename;
    while (infile >> par[0] >> par[1] >> par[2] >> par[3]){
      filename = Form("gallery/ds_phi_%.4d.pdf", (int) (par[0] * 1000));
      printplot(filename.Data(), par);
      //break;
    }
    infile.close();
  }

  if (false){
    double globalpar[5] = {0.224538, 2.00410, 3.75802, 1.38555, 0.908769};
    sigma2.SetParameters(globalpar);
    TString filename;
    double ss = 1.985;
    while (ss <= 2.835){
      filename = Form("gallery/global_ds_phi_%.4d.pdf", (int) (ss * 1000 + 1.0e-3));
      printplot2(filename.Data(), &ss);
      ss += 0.01;
    }
  }

  if (true){
    double globalpar[5] = {0.232612, 1.95038, 4.02454, 1.52884, 0.525636};
    sigma2.SetParameters(globalpar);
    TString filename;
    double ss = 1.985;
    while (ss <= 2.45){
      filename = Form("gallery/low_ds_phi_%.4d.pdf", (int) (ss * 1000 + 1.0e-3));
      printplot2(filename.Data(), &ss);
      ss += 0.01;
    }
  }

  if (false){
    double globalpar[5] = {0.224538, 2.00410, 3.75802, 1.38555, 0.908769};
    sigma2.SetParameters(globalpar);
    TString filename;
    double ss = 1.985;
    while (ss <= 2.835){
      filename = Form("gallery/global_logds_phi_%.4d.pdf", (int) (ss * 1000 + 1.0e-3));
      printplot2log(filename.Data(), &ss);
      ss += 0.01;
    }
  }

  if (false){
    double Chi2 = 0.0;
    creatset("dataset.dat", "clasdb_E63M1.txt", 1.9, 2.45);
    addset("dataset.dat", "clasdb_E63M2.txt", 1.9, 2.45);
    Chi2 = f1Fit(par);
    cout << Chi2 << endl;
  }
    
  if (false){
    parameterplot("fitphi.txt");
  }

  return 0;
}
