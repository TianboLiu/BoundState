#ifndef _PHIFIT_H_
#define _PHIFIT_H_

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

int makeset(const char * savefile, const double * Eslimit){
  std::ifstream infile1("clasdb_E63M1.txt");
  std::ifstream infile2("clasdb_E63M2.txt");
  char tmp[256];
  for (int i = 0; i < 8; i++){
    infile1.getline(tmp, 256);
    infile2.getline(tmp, 256);
  }
  double Es, cth, ds, err;
  FILE * fp;
  fp = fopen(savefile, "w");
  while (!infile1.eof()){
    infile1 >> Es >> cth >> ds >> err;
    if (Es < Eslimit[0] || Es > Eslimit[1]) continue;
    fprintf(fp, "%.3f  %.2f  %.4f  %.4f\n", Es, cth, ds, err);
  }
  infile1.close();
  while (!infile2.eof()){
    infile2 >> Es >> cth >> ds >> err;
    if (Es < Eslimit[0] || Es > Eslimit[1]) continue;
    fprintf(fp, "%.3f  %.2f  %.4f  %.4f\n", Es, cth, ds, err);
  }
  infile2.close();
  fclose(fp);
  return 0;
}

int plotdata(const char * datafile, const double * Eslimit){
  int Nt = 0;
  double Es, cth, ds, err;
  std::ifstream infile(datafile);
  while (!infile.eof()){
    infile >> Es >> cth >> ds >> err;
    if (Es > Eslimit[0] && Es < Eslimit[1]) Nt++;
  }
  infile.clear();
  infile.seekg(0, std::ios::beg);
  TH1D * h0 = new TH1D("", "", 1, -1.0, 1.0);
  h0->GetYaxis()->SetRangeUser(0.0, 0.15);
  TGraphErrors * g0 = new TGraphErrors(Nt);
  for (int i = 0; i < Nt; ){
    infile >> Es >> cth >> ds >> err;
    if (Es > Eslimit[0] && Es < Eslimit[1]);
    g0->SetPoint(i, cth, ds);
    g0->SetPointError(i, 0, err);
    i++;
  }
  TF1 f0("f0", "[0]*exp(-[1]*(1-x)-[2]*(1-x)*(1-x))", -1, 1);
  g0->Fit(&f0, "Q");
  g0->SetMarkerStyle(8);
  g0->SetMarkerSize(0.5);
  g0->SetMarkerColor(2);
  g0->SetLineWidth(0.02);
  TCanvas * c0 = new TCanvas("", "", 800, 600);
  //h0->Draw();
  g0->Draw("0apesame");
  c0->Print("points.pdf");
  return 0;
}





#endif
