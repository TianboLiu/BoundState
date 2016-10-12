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
  TString sfiles[4] = {"spm_all/caspmderads.out","spm_all/cbspmderads.out","spm_all/ccspmderads.out","spm_all/cfspmderads.out"};
  FILE * cfs = fopen("cfs.dat", "w");
  for (int i = 0; i < 4; i++){
    std::ifstream infile(sfiles[i]);
    while (infile >> pm >> sm >> errp >> errs){
      if (sm < 1.0e-12) continue;
      fprintf(cfs, "%.4E  %.4E  %.4E  %.4E\n", pm*1.0e-3, sm*scale, errp*1.0e-3, errs*scale);
    }
    infile.close();
  }
  fclose(cfs);
  TString pfiles[4] = {"spm_all/caspmderadp.out","spm_all/cbspmderadp.out","spm_all/ccspmderadp.out","spm_all/cfspmderadp.out"};
  FILE * cfp = fopen("cfp.dat", "w");
  for (int i = 0; i < 4; i++){
    std::ifstream infile(pfiles[i]);
    while (infile >> pm >> sm >> errp >> errs){
      if (sm < 1.0e-12) continue;
      fprintf(cfp, "%.4E  %.4E  %.4E  %.4E\n", pm*1.0e-3, sm*scale, errp*1.0e-3, errs*scale);
    }
    infile.close();
  }
  fclose(cfp);
  TString tfiles[4] = {"spm_all/caspmderad80.out","spm_all/cbspmderad80.out","spm_all/ccspmderad80.out","spm_all/cfspmderad80.out"};
  FILE * cft = fopen("cft.dat", "w");
  for (int i = 0; i < 4; i++){
    std::ifstream infile(tfiles[i]);
    while (infile >> pm >> sm >> errp >> errs){
      if (sm < 1.0e-12) continue;
      fprintf(cft, "%.4E  %.4E  %.4E  %.4E\n", pm*1.0e-3, sm*scale, errp*1.0e-3, errs*scale);
    }
    infile.close();
  }
  fclose(cft);
  TString cefile[1] = {"sem_all/cbsemderad.out"};
  FILE * cet = fopen("cet.dat", "w");
  for (int i = 0; i < 1; i++){
    std::ifstream infile(cefile[i]);
    while (infile >> pm >> sm >> errp >> errs){
      if (errs < 1.0e-12) continue;//errs = 1.0e-6;
      fprintf(cet, "%.4E  %.4E  %.4E  %.4E\n", pm*1.0e-3, sm, errp*1.0e-3, errs);
    }
    infile.close();
  }
  fclose(cet);
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
  min->SetTolerance(1.0e-6);
  min->SetPrintLevel(1);
  ROOT::Math::Functor f(&Chi2MissingMomentum, 5);
  min->SetFunction(f);
  min->SetVariable(0, "A0", 10.6552, 1.0e-5);
  min->SetVariable(1, "A2", 32.3105, 1.0e-5);
  min->SetFixedVariable(2, "A4", 0.0);
  //min->SetVariable(2, "A4", 0.1, 1.0e-5);
  min->SetVariable(3, "B2", 8.43226, 1.0e-5);
  min->SetFixedVariable(4, "B4", 0.0);
  //min->SetVariable(4, "B4", 0.0105249, 1.0e-5);
  min->Minimize();
  std::cout << "chi2: " << min->MinValue() << std::endl;
  //plot
  gStyle->SetOptStat(0);
  TH1D * h0 = new TH1D("", "", 1, -0.31, 0.31);
  h0->GetYaxis()->SetRangeUser(0.0, 14.0);
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
  c0->Print("carbon_momentum.pdf");
  return 0;
}

double MissingEnergy(double * x, const double * par){
  double E = x[0];
  double E0 = 0.0125;
  //if (E < E0) return 0;
  double A1 = par[0];
  double a1 = par[1];
  double b1 = par[2];
  double A2 = par[3];
  double a2 = par[4];
  double b2 = par[5];
  double result = (A1 * exp(-pow((E-a1)/b1, 2)) + A2 * exp(-pow((E-a2)/b2, 2)));
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
  ROOT::Math::Functor f(&Chi2MissingEnergy, 6);
  min->SetFunction(f);
  min->SetVariable(0, "A1", 0.569161, 1.0e-3);
  min->SetFixedVariable(1, "a1", 0.0175);
  //min->SetVariable(1, "a1", 0.0175, 1.0e-5);
  min->SetVariable(2, "b1", 0.00510586, 1.0e-5);
  min->SetVariable(3, "A2", 0.0683607, 1.0e-5);
  //min->SetFixedVariable(4, "a2", 0.0375);
  min->SetVariable(4, "a2", 0.0369238, 1.0e-5);
  min->SetVariable(5, "b2", 0.0173035, 1.0e-5);
  min->Minimize();
  std::cout << "chi2: " << min->MinValue() << std::endl;
  //plot
  gStyle->SetOptStat(0);
  TH1D * h0 = new TH1D("", "", 1, 0.0, 0.1);
  h0->GetYaxis()->SetRangeUser(0.0, 0.8);
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
  c0->Print("carbon_energy.pdf");
  return 0;
}

#endif
