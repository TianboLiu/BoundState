#ifndef _SPECTRA_H_
#define _SPECTRA_H_

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
#include "TPad.h"
#include "TRandom.h"

int _Np = 0;
double _x[200], _y[200], _ex[200], _ey[200];

int makeset(){
  double pm, sm, errp, errs;
  double scale = 1.0e8;
  TString tfiles[4] = {"spm_all/gaspmderad80.out","spm_all/gbspmderad80.out","spm_all/gcspmderad80.out","spm_all/gdspmderad80.out"};
  FILE * gft = fopen("gft.dat", "w");
  for (int i = 0; i < 4; i++){
    std::ifstream infile(tfiles[i]);
    while (infile >> pm >> sm >> errp >> errs){
      if (sm < 1.0e-12) continue;
      fprintf(gft, "%.4E  %.4E  %.4E  %.4E\n", pm*1.0e-3, sm*scale, errp*1.0e-3, errs*scale);
    }
    infile.close();
  }
  fclose(gft);
  TString gefile[1] = {"sem_all/gbsemderad.out"};
  FILE * get = fopen("get.dat", "w");
  for (int i = 0; i < 1; i++){
    std::ifstream infile(gefile[i]);
    while (infile >> pm >> sm >> errp >> errs){
      if (errs < 1.0e-12) continue;//errs = 1.0e-6;
      fprintf(get, "%.4E  %.4E  %.4E  %.4E\n", pm*1.0e-3, sm, errp*1.0e-3, errs);
    }
    infile.close();
  }
  fclose(get);
  return 0;
}

int readset(char * filename){
  std::ifstream infile(filename);
  int Np = 0;
  while(infile >> _x[Np] >> _y[Np] >> _ex[Np] >> _ey[Np]){
    Np++;
  }
  std::cout << "Points: " << Np << std::endl;
  _Np = Np;
  return 0;
}

double MissingMomentum(double * x, const double * par){
  double pm = x[0];
  double A0 = par[0];
  double A2 = par[1];
  double A4 = par[2];
  double B2 = par[3];
  double B4 = par[4];
  double result = (A0 + pow(A2 * pm, 2) + pow(A4 * pm, 2)) * exp( - pow(B2 * pm, 2) - pow(B4 * pm, 4));
  return result;
}

double Chi2MissingMomentum(const double * par){
  double sum = 0.0;
  for (int i = 0; i < _Np; i++){
    sum += pow( MissingMomentum(&_x[i], par) - _y[i], 2) / pow(_ey[i], 2);
  }
  return sum;
}

int FitMissingMomentum(const char * minName = "Minuit", const char * algoName = "Migrad"){
  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  min->SetMaxFunctionCalls(1000000);
  min->SetTolerance(1.0e-7);
  min->SetPrintLevel(1);
  ROOT::Math::Functor f(&Chi2MissingMomentum, 5);
  min->SetFunction(f);
  min->SetVariable(0, "A0", 58.3382, 1.0e-6);
  min->SetVariable(1, "A2", 69.2938, 1.0e-6);
  min->SetFixedVariable(2, "A4", 0.0);
  min->SetVariable(3, "B2", 7.82756, 1.0e-6);
  min->SetFixedVariable(4, "B4", 0.0);
  min->Minimize();
  std::cout << "chi2: " << min->MinValue() << std::endl;
  //plot
  gStyle->SetOptStat(0);
  TH1D * h0 = new TH1D("", "", 1, -0.31, 0.31);
  h0->GetYaxis()->SetRangeUser(0.0, 80.0);
  h0->GetXaxis()->CenterTitle();
  h0->GetXaxis()->SetTitle("p / GeV");
  h0->GetXaxis()->SetTitleSize(0.04);
  h0->GetXaxis()->SetTitleOffset(0.8);
  TGraphErrors * g0 = new TGraphErrors(_Np, _x, _y, _ex, _ey);
  g0->SetMarkerStyle(8);
  g0->SetMarkerColor(4);
  g0->SetMarkerSize(0.5);
  g0->SetLineWidth(0.04);
  TF1 * f0 = new TF1("f0", MissingMomentum, -0.3, 0.3, 5);
  f0->SetParameters(min->X());
  TCanvas * c0 = new TCanvas("c0", "c0", 800, 600);
  h0->Draw();
  g0->Draw("pesame");
  f0->Draw("same");
  c0->Print("gold_momentum.pdf");
  return 0;
}

double MissingEnergy(double * x, const double * par){
  double E = x[0];
  //double E0 = 0.0125;
  //if (E < E0) return 0;
  double A1 = par[0];
  double a1 = par[1];
  double b1 = par[2];
  double A2 = par[3];
  double a2 = par[4];
  double b2 = par[5];
  double result = A1 * atan(A2 * pow(E/0.01, a1)) * exp(-b1 * pow(E/0.01, a2));
  return result;
}

double Chi2MissingEnergy(const double * par){
  double sum = 0.0;
  for (int i = 0; i < _Np; i++){
    sum += pow( MissingEnergy(&_x[i], par) - _y[i], 2) / pow(_ey[i], 2);
  }
  return sum;
}

int FitMissingEnergy(const char * minName = "Minuit", const char * algoName = "Migrad"){
  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  min->SetMaxFunctionCalls(1000000);
  min->SetTolerance(1.0e-6);
  min->SetPrintLevel(1);
  ROOT::Math::Functor f(&Chi2MissingEnergy, 5);
  min->SetFunction(f);
  min->SetVariable(0, "A1", 1.73622, 1.0e-5);
  min->SetVariable(1, "a1", 3.07375, 1.0e-5);
  min->SetVariable(2, "b1", 0.645561, 1.0e-5);
  min->SetVariable(3, "A2", 14.1433, 1.0e-5);
  min->SetVariable(4, "a2", 0.795058, 1.0e-5);
  min->Minimize();
  std::cout << "chi2: " << min->MinValue() << std::endl;
  //plot
  gStyle->SetOptStat(0);
  TH1D * h0 = new TH1D("", "", 1, 0.0, 0.1);
  h0->GetYaxis()->SetRangeUser(0.0, 2.0);
  h0->GetXaxis()->CenterTitle();
  h0->GetXaxis()->SetTitle("E / GeV");
  h0->GetXaxis()->SetTitleSize(0.04);
  h0->GetXaxis()->SetTitleOffset(0.8);
  TGraphErrors * g0 = new TGraphErrors(_Np, _x, _y, _ex, _ey);
  g0->SetMarkerStyle(8);
  g0->SetMarkerColor(4);
  g0->SetMarkerSize(0.5);
  g0->SetLineWidth(0.04);
  TF1 * f0 = new TF1("f0", MissingEnergy, 0.0, 0.1, 6);
  f0->SetParameters(min->X());
  f0->SetNpx(500);
  TCanvas * c0 = new TCanvas("c0", "c0", 800, 600);
  h0->Draw();
  g0->Draw("pesame");
  f0->Draw("same");
  c0->Print("gold_energy.pdf");
  return 0;
}

#endif
