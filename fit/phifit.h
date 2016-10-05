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
  int Nt = 0;
  fp = fopen(savefile, "w");
  while (!infile1.eof()){
    infile1 >> Es >> cth >> ds >> err;
    if (Es < Eslimit[0] || Es > Eslimit[1]) continue;
    fprintf(fp, "%.3f  %.2f  %.4f  %.4f\n", Es, cth, ds, err);
    Nt++;
  }
  infile1.close();
  while (!infile2.eof()){
    infile2 >> Es >> cth >> ds >> err;
    if (Es < Eslimit[0] || Es > Eslimit[1]) continue;
    fprintf(fp, "%.3f  %.2f  %.4f  %.4f\n", Es, cth, ds, err);
    Nt++;
  }
  infile2.close();
  fclose(fp);
  std::cout << Es << "  " << Nt << std::endl;
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
  g0->Fit(&f0);
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

double f0(const double * par){
  std::ifstream infile("dataset.dat");
  double sum = 0.0;
  double Es, cth, ds, err;
  while (!infile.eof()){
    infile >> Es >> cth >> ds >> err;
    sum += pow((par[0] * exp(-par[1]*(1.0-cth)-par[2]*pow(1.0-cth, 2))) - ds, 2) / (err * err);
  }
  infile.close();
  return sum;
}

int f0Fit(const char * minName = "Minuit2", const char * algoName = "Migrad"){
  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  min->SetMaxFunctionCalls(100000);
  min->SetTolerance(1.0e-6);
  min->SetPrintLevel(0);
  ROOT::Math::Functor f(&f0, 3);
  min->SetFunction(f);
  min->SetVariable(0, "A", 0.12, 1.0e-3);
  min->SetVariable(1, "B", 1.03, 1.0e-3);
  min->SetVariable(2, "C", -0.455, 1.0e-3);
  min->Minimize();
  const double * xs = min->X();
  const double chi2 = min->MinValue();
  std::cout << xs[0] << "   " << xs[1] << "   " << xs[2] << "   " << chi2 << std::endl;
  return 0;
}



#endif
