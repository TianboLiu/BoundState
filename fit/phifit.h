#ifndef _PHIFIT_H_
#define _PHIFIT_H_

#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>

#include "TStyle.h"
#include "TString.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TF1.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TCanvas.h"

int creatset(const char * savefile, const char * datafile, const double Emin, const double Emax = 0){
  double Erange[2] = {Emin, Emax};
  if (Emax < Emin){
    Erange[0] = Emin - 1.0e-9;
    Erange[1] = Emin + 1.0e-9;
  }
  std::ifstream infile(datafile);//Load datafile
  char tmp[256];
  for (int i = 0; i < 8; i++){//Pass the header
    infile.getline(tmp, 256);
  }
  double Es, cth, ds, err;
  FILE * fp;
  int Nt = 0;
  fp = fopen(savefile, "w");
  while (infile >> Es >> cth >> ds >> err){
    if (Es < Erange[0] || Es > Erange[1]) continue;
    if (err < 1.0e-9) continue;
    fprintf(fp, "%.3f  %.2f  %.4f  %.4f\n", Es, cth, ds, err);
    Nt++;
  }
  infile.close();
  fclose(fp);
  //std::cout << Es << "  " << Nt << std::endl;
  return 0;
}

int addset(const char * savefile, const char * datafile, const double Emin, const double Emax = 0){
  double Erange[2] = {Emin, Emax};
  if (Emax == 0){
    Erange[0] = Emin - 1e-9;
    Erange[1] = Emin + 1e-9;
  }
  std::ifstream infile(datafile);//Load datafile
  char tmp[256];
  for (int i = 0; i < 8; i++){//Pass the header
    infile.getline(tmp, 256);
  }
  double Es, cth, ds, err;
  FILE * fp;
  int Nt = 0;
  fp = fopen(savefile, "a");
  while (infile >> Es >> cth >> ds >> err){
    if (Es < Erange[0] || Es > Erange[1]) continue;
    if (err == 0) continue;
    fprintf(fp, "%.3f  %.2f  %.4f  %.4f\n", Es, cth, ds, err);
    Nt++;
  }
  infile.close();
  fclose(fp);
  //std::cout << Es << "  " << Nt << std::endl;
  return 0;
}

int plotdata(const char * datafile){
  int Nt = 0;
  double Es, cth, ds, err;
  std::ifstream infile(datafile);
  while (infile >> Es >> cth >> ds >> err) Nt++;//count data points
  std::cout << Nt << std::endl;
  infile.clear();
  infile.seekg(0, std::ios::beg);//go back to the beginning of the datafile
  TGraphErrors * g0 = new TGraphErrors(Nt);
  double upper, lower;
  for (int i = 0; i < Nt; i++){
    infile >> Es >> cth >> ds >> err;
    g0->SetPoint(i, cth, ds);
    g0->SetPointError(i, 0, err);
    if (i == 0) {upper = ds + err; lower = ds - err;}
    if (upper < ds + err) upper = ds + err;
    if (lower > ds - err) lower = ds - err;
  }
  infile.close();
  TH1D * h0 = new TH1D("h0", "", 1, -1.0, 1.0);
  h0->GetYaxis()->SetRangeUser(0.0, lower+upper);
  h0->GetXaxis()->SetTitle("cos#theta_{c.m.}");
  h0->GetXaxis()->CenterTitle();
  h0->GetYaxis()->SetTitle("d#sigma / dcos#theta_{c.m.} (#mub)");
  h0->GetYaxis()->CenterTitle();
  h0->GetYaxis()->SetTitleOffset(1.5);
  TF1 f0("f0", "[0]", -1, 1); 
  f0.SetParameter(0, (lower+upper)/2.0);
  g0->Fit(&f0);
  g0->SetMarkerStyle(8);
  g0->SetMarkerSize(0.6);
  g0->SetMarkerColor(4);
  g0->SetLineWidth(0.5);
  g0->GetFunction("f0")->SetLineWidth(1);
  g0->GetFunction("f0")->SetLineColor(2);
  TCanvas * c0 = new TCanvas("", "", 800, 600);
  gStyle->SetOptStat(0);
  h0->Draw();
  g0->Draw("pesame");
  c0->Print("plotdata.pdf");
  return 0;
}

TF1 sigma0("sigma0", "[0]", -1, 1);
TF1 sigma("sigma", "[0]*exp(-[1]*(1-x)-[2]*(1-x)*(1-x))", -1, 1);
TF1 sigma2("sigma2", "[0]*(1+[1]*x+[2]*x*(x + abs(x)))", -1, 1);
TF1 sigma3("sigma3", "[0]+x*[1]*exp(-[2]*(1.0-x))", -1, 1);
TF1 sigma1("sigma1", "([0] + [1] * x) * exp([2] * x)", -1, 1);




double f0(const double * par){
  std::ifstream infile("dataset.dat");
  double sum = 0.0;
  sigma1.SetParameters(par);
  double Es, cth, ds, err;
  while (!infile.eof()){
    infile >> Es >> cth >> ds >> err;
    sum += pow(sigma1.Eval(cth) - ds, 2) / (err * err);
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
  min->SetVariable(0, "A", 0.1, 1.0e-3);
  min->SetVariable(1, "B", 0.0, 1.0e-3);
  min->SetVariable(2, "C", 0.0, 1.0e-3);
  min->Minimize();
  const double * xs = min->X();
  const double chi2 = min->MinValue();
  sigma1.SetParameters(xs);
  double total = sigma1.Integral(-1,1);
  std::cout << total << " " << chi2 << " " << xs[0] << "   " << xs[1] << "   " << xs[2] << std::endl;
  return 0;
}



#endif
