#ifndef _SPECTRA_H_
#define _SPECTRA_H_

#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>

#include "TString.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TF1.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TCanvas.h"

int makeset(){
  double pm, sm, errp, errs;
  std::ifstream cfs1("spm_all/caspmderads.out");
  std::ifstream cfs2("spm_all/cbspmderads.out");
  std::ifstream cfs3("spm_all/ccspmderads.out");
  std::ifstream cfs4("spm_all/cfspmderads.out");
  FILE * cfs = fopen("cfs.dat", "w");
  while (!cfs1.eof()){
    cfs1 >> pm >> sm >> errp >> errs;
    if (sm < 1.0e-12) continue;
    fprintf(cfs, "%.4E  %.4E  %.4E  %.4E\n", pm, sm, errp, errs);
  }
  cfs1.close();
  while (!cfs2.eof()){
    cfs2 >> pm >> sm >> errp >> errs;
    if (sm < 1.0e-12) continue;
    fprintf(cfs, "%.4E  %.4E  %.4E  %.4E\n", pm, sm, errp, errs);
  }
  cfs2.close();
  while (!cfs3.eof()){
    cfs3 >> pm >> sm >> errp >> errs;
    if (sm < 1.0e-12) continue;
    fprintf(cfs, "%.4E  %.4E  %.4E  %.4E\n", pm, sm, errp, errs);
  }
  cfs3.close();
  while (!cfs4.eof()){
    cfs4 >> pm >> sm >> errp >> errs;
    if (sm < 1.0e-12) continue;
    fprintf(cfs, "%.4E  %.4E  %.4E  %.4E\n", pm, sm, errp, errs);
  }
  cfs4.close();
  fclose(cfs);
  std::ifstream cfp1("spm_all/caspmderadp.out");
  std::ifstream cfp2("spm_all/cbspmderadp.out");
  std::ifstream cfp3("spm_all/ccspmderadp.out");
  std::ifstream cfp4("spm_all/cfspmderadp.out");
  FILE * cfp = fopen("cfp.dat", "w");
  while (!cfp1.eof()){
    cfp1 >> pm >> sm >> errp >> errs;
    if (sm < 1.0e-12) continue;
    fprintf(cfp, "%.4E  %.4E  %.4E  %.4E\n", pm*1.0e-3, sm*1.0e8, errp*1.0e-3, errs*1.0e8);
  }
  cfp1.close();
  while (!cfp2.eof()){
    cfp2 >> pm >> sm >> errp >> errs;
    if (sm < 1.0e-12) continue;
    fprintf(cfp, "%.4E  %.4E  %.4E  %.4E\n", pm*1.0e-3, sm*1.0e8, errp*1.0e-3, errs*1.0e8);
  }
  cfp2.close();
  while (!cfp3.eof()){
    cfp3 >> pm >> sm >> errp >> errs;
    if (sm < 1.0e-12) continue;
    fprintf(cfp, "%.4E  %.4E  %.4E  %.4E\n", pm*1.0e-3, sm*1.0e8, errp*1.0e-3, errs*1.0e8);
  }
  cfp3.close();
  while (!cfp4.eof()){
    cfp4 >> pm >> sm >> errp >> errs;
    if (sm < 1.0e-12) continue;
    fprintf(cfp, "%.4E  %.4E  %.4E  %.4E\n", pm*1.0e-3, sm*1.0e8, errp*1.0e-3, errs*1.0e8);
  }
  cfp4.close();
  fclose(cfp);
  return 0;
}

double f0(const double * par){
  double pm, sm, errp, errs;
  std::ifstream fs("cfs.dat");
  std::ifstream fp("cfp.dat");
  double sum = 0.0;
  int Nt = 0;
  while (!fs.eof()){
    fs >> pm >> sm >> errp >> errs;continue;
    sum += pow(par[0] * exp(-par[2]*par[2]*pm*pm) - sm, 2) / (errs * errs);
    Nt++;
  }
  fs.close();
  while (!fp.eof()){
    fp >> pm >> sm >> errp >> errs;//continue;
    sum += pow(par[1] * par[2] * par[2] * pm * pm * exp(-par[2]*par[2]*pm*pm) - sm, 2) / (errs * errs);
    Nt++;
  }
  //std::cout << "Points: " << Nt << std::endl;
  return sum;
}

int fits(const char * minName = "Minuit", const char * algoName = "Migrad"){
  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  min->SetMaxFunctionCalls(1000000);
  min->SetTolerance(1.0e-4);
  min->SetPrintLevel(1);
  ROOT::Math::Functor f(&f0, 3);
  min->SetFunction(f);
  min->SetVariable(0, "A", 1.0, 1.0e-3);
  min->SetVariable(1, "B", 1.0, 1.0e-3);
  min->SetVariable(2, "C", 0.2, 1.0e-3);
  min->Minimize();
  const double * xs = min->X();
  const double chi2 = min->MinValue();
  std::cout << xs[0] << "   " << xs[1] << "   " << xs[2] << "   " << chi2 << std::endl;
  return 0;
}

#endif
