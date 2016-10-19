#include "lambdafit.h"

using namespace std;

int main(int argc, char * argv[]){
  double par[5];

  if (false){
    creatset("dataset.dat", "clasdb_E57M1.txt", 2.0, 3.0);
    t0Fit(par);
    sigma0.SetParameters(par);
    t0plot("sig_lambda.pdf", "dataset.dat");
  }
  
  if (false){
    double Chi2 = 0.0;
    double Wa = 1.95;
    double Wb = 2.05;
    double chi2;
    FILE * fp = fopen("fitlambda.txt", "w");
    fprintf(fp, "sqrt(s)\t a0\t a2\t b1\t\n");
    TString filename;
    for (int i = 0; i < 10; i++){
      //cout <<  i << " " << Wa << " " << Wb << endl;
      creatset("dataset.dat", "clasdb_E57M2.txt", Wa, Wb);
      chi2 = f0Fit(par);
      cout << chi2 << endl;
      Chi2 += chi2;
      filename = Form("plotdata%.2d.pdf", i); 
      plotdata("dataset.dat", par, filename.Data());
      fprintf(fp, "%.3f\t %.4f\t %.4f\t %.4f\n", (Wa+Wb)/2.0, par[0], par[1], par[2]);
      Wa += 0.1;
      Wb += 0.1;
      if (Wa >= 2.8) break;
    }
    fclose(fp);
    cout << "Chi2 = " << Chi2 << endl;
  }

  if (false){
    parameterplot("fitlambda.txt");
  }

  if (false){
    double Chi2 = 0.0;
    creatset("dataset.dat", "clasdb_E57M2.txt", 1.9, 2.45);
    Chi2 = f1Fit(par);
    cout << Chi2 << endl;
  }

  if (true){
    double globalpar[5] = {11.299, 4.60959, 0.835621, 0.54681, 0.185236};
    sigma2.SetParameters(globalpar);
    TString filename;
    double ss = 2.0;
    while (ss <= 2.45){
      filename = Form("gallery/global_ds_lambda1520_%.4d.pdf", (int) (ss * 1000 + 1.0e-3));
      printplot2(filename.Data(), &ss);
      ss += 0.1;
    }
  }


  return 0;
}
