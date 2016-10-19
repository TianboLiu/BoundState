#ifndef _LAMBDAFIT_H_
#define _LAMBDAFIT_H_

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
#include "TF2.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TCanvas.h"

double dsigma(const double * var, const double * par){
  double x = var[0];
  double r0 = exp(par[2] * x) * par[2] / (2.0 * sinh(par[2]));
  double r2 = x * x * exp(par[2] * x) * pow(par[2], 3) / (2.0 * (par[2] * par[2] + 2.0) * sinh(par[2]) - 4.0 * par[2] * cosh(par[2]));
  double result = par[0] * (r0 + par[1] * r2) / (1.0 + par[1]);
  return result;
}

double ds_global(const double * var, const double * par){
  double sr = var[0];//sqrt(s)
  double x = var[1];//cos theta
  double M = 1.5195 + 0.493677;
  double b1 = par[3] * (sr - par[4]);
  double a2 = 0.25;
  double a0 = par[0] * (sr - M) * exp(-par[1] * pow(sr - M, par[2]));
  double r0 = exp(b1 * x) * b1 / (2.0 * sinh(b1));
  double r2 = x * x * exp(b1 * x) * pow(b1, 3) / (2.0 * (b1 * b1 + 2.0) * sinh(b1) - 4.0 * b1 * cosh(b1));
  double result = a0 * (r0 + a2 * r2) / (1.0 + a2);
  return result;
}

double dsigma0(const double * var, const double * par){
  double Eg = var[0];
  double x = Eg - 1.69062;
  if (x < 0) return 0;
  double result = par[0] * pow(x, par[1]) * exp(-par[2] * x);
  return result;
}

double dsigma1(const double * var, const double * par){
  double W = var[0];
  double x = var[1];
  const double Mp = 0.938272;
  double dE = (W*W - Mp*Mp) / (2.0 * Mp) - 1.69062;
  double a0 = par[0] * pow(dE, par[1]) * exp(-par[2] * dE);
  double a2 = par[3] * dE;
  double b1 = par[4] * pow(dE, par[5]);
  double r0 = exp(b1 * x) * b1 / (2.0 * sinh(b1));
  double r2 = x * x * exp(b1 * x) * pow(b1, 3) / (2.0 * (b1 * b1 + 2.0) * sinh(b1) - 4.0 * b1 * cosh(b1));
  double result = a0 * (r0 + a2 * r2) / (1.0 + a2);
  return result;
}


TF1 sigma0("sigma0", dsigma0, 1.5, 4.5, 3);
TF2 sigma4("sigma4", dsigma1, 2.0, 3.0, -1.0, 1.0, 6);
TF1 sigma("sigma", "[0]*exp(-[1]*(1-x)-[2]*(1-x)*(1-x))", -1, 1);
TF1 sigma3("sigma3", "[0]+x*[1]*exp(-[2]*(1.0-x))", -1, 1);
TF1 sigma1("sigma1", dsigma, -1, 1, 3);
TF2 sigma2("sigma1", ds_global, 1.96, 3.0, -1.0, 1.0, 5);


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

int plotdata(const char * datafile, const double * par, const char * savefile){
  int Nt = 0;
  double Es, cth, ds, err;
  std::ifstream infile(datafile);
  while (infile >> Es >> cth >> ds >> err) Nt++;//count data points
  std::cout << Nt << std::endl;
  infile.clear();
  infile.seekg(0, std::ios::beg);//go back to the beginning of the datafile
  TGraphErrors * g0 = new TGraphErrors(Nt);
  double upper = 0;
  double lower = 0;
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
  sigma1.SetParameters(par);
  sigma1.SetLineWidth(1);
  sigma1.SetLineColor(2);
  //g0->Fit(&f0);
  g0->SetMarkerStyle(8);
  g0->SetMarkerSize(0.6);
  g0->SetMarkerColor(4);
  g0->SetLineWidth(0.5);
  //g0->GetFunction("f0")->SetLineWidth(1);
  //g0->GetFunction("f0")->SetLineColor(2);
  TCanvas * c0 = new TCanvas("", "", 800, 600);
  gStyle->SetOptStat(0);
  h0->Draw();
  g0->Draw("pesame");
  sigma1.Draw("same");
  c0->Print(savefile);
  h0->Delete();
  g0->Delete();
  c0->Close();
  return 0;
}

int printplot(const char * savefile, const double * par){
  double s = par[0];
  creatset("plotset.dat", "clasdb_E63M1.txt", s);
  addset("plotset.dat", "clasdb_E63M2.txt", s);
  int Nt = 0;
  double Es, cth, ds, err;
  std::ifstream infile("plotset.dat");
  while (infile >> Es >> cth >> ds >> err) Nt++;//count data points
  std::cout << Nt << std::endl;
  infile.clear();
  infile.seekg(0, std::ios::beg);//go back to the beginning of the datafile
  TGraphErrors * g0 = new TGraphErrors(Nt);
  double upper = 0;
  double lower = 0;
  for (int i = 0; i < Nt; i++){
    infile >> Es >> cth >> ds >> err;
    g0->SetPoint(i, cth, ds);
    g0->SetPointError(i, 0, err);
    if (i == 0) {upper = ds + err; lower = ds - err;}
    if (upper < ds + err) upper = ds + err;
    if (lower > ds - err) lower = ds - err;
  }
  infile.close();
  TString title = Form("#sqrt{s} = %.3f GeV", s);
  TH1D * h0 = new TH1D("h0", title, 1, -1.0, 1.0);
  h0->Fill(0.0, -100.0);
  h0->GetYaxis()->SetRangeUser(-0.02, lower+upper);
  h0->GetXaxis()->SetTitle("cos#theta_{c.m.}");
  h0->GetXaxis()->CenterTitle();
  h0->GetYaxis()->SetTitle("d#sigma / dcos#theta_{c.m.} (#mub)");
  h0->GetYaxis()->CenterTitle();
  h0->GetYaxis()->SetTitleOffset(1.5);
  sigma1.SetParameters(&par[1]);
  sigma1.SetLineWidth(1);
  sigma1.SetLineColor(2);
  g0->SetMarkerStyle(8);
  g0->SetMarkerSize(0.6);
  g0->SetMarkerColor(4);
  g0->SetLineWidth(0.5);
  TCanvas * c0 = new TCanvas("", "", 800, 600);
  gStyle->SetOptStat(0);
  h0->Draw();
  g0->Draw("pesame");
  sigma1.Draw("same");
  c0->Print(savefile);
  h0->Delete();
  g0->Delete();
  c0->Close();
  return 0;
}

int printplot2(const char * savefile, const double * par){
  double s = par[0];
  creatset("plotset.dat", "clasdb_E57M2.txt", s-0.05, s+0.05);
  int Nt = 0;
  double Es, cth, ds, err;
  std::ifstream infile("plotset.dat");
  while (infile >> Es >> cth >> ds >> err) Nt++;//count data points
  std::cout << Nt << std::endl;
  infile.clear();
  infile.seekg(0, std::ios::beg);//go back to the beginning of the datafile
  TGraphErrors * g0 = new TGraphErrors(Nt);
  double upper = 0;
  double lower = 0;
  for (int i = 0; i < Nt; i++){
    infile >> Es >> cth >> ds >> err;
    g0->SetPoint(i, cth, ds);
    g0->SetPointError(i, 0, err);
    if (i == 0) {upper = ds + err; lower = ds - err;}
    if (upper < ds + err) upper = ds + err;
    if (lower > ds - err) lower = ds - err;
  }
  infile.close();
  TString title = Form("#sqrt{s} = %.3f GeV", Es);
  TH1D * h0 = new TH1D("h0", title, 1, -1.0, 1.0);
  h0->Fill(0.0, -100.0);
  h0->GetYaxis()->SetRangeUser(-0.02, lower+upper);
  h0->GetXaxis()->SetTitle("cos#theta_{c.m.}");
  h0->GetXaxis()->CenterTitle();
  h0->GetYaxis()->SetTitle("d#sigma / dcos#theta_{c.m.} (#mub)");
  h0->GetYaxis()->CenterTitle();
  h0->GetYaxis()->SetTitleOffset(1.5);
  TGraph * g1 = new TGraph(101);
  for (int i = 0; i < 101; i++){
    g1->SetPoint(i, -1.0+0.02*i, sigma2.Eval(Es, -1.0+0.02*i));
    //std::cout << Es << " " << -1.0+0.02*i << " " << sigma2.Eval(Es, -1.0+0.02*i) << std::endl;
  }
  g1->SetLineColor(2);
  g1->SetLineWidth(2);
  g0->SetMarkerStyle(8);
  g0->SetMarkerSize(0.6);
  g0->SetMarkerColor(4);
  TCanvas * c0 = new TCanvas("", "", 800, 600);
  gStyle->SetOptStat(0);
  h0->Draw();
  g0->Draw("pesame");
  g1->Draw("csame");
  c0->Print(savefile);
  h0->Delete();
  g0->Delete();
  g1->Delete();
  c0->Close();
  return 0;
}

double t0(const double * par){
  std::ifstream infile("dataset.dat");
  double sum = 0.0;
  double W, Eg, sig, err;
  sigma0.SetParameters(par);
  while(infile >> W >> Eg >> sig >> err){
    sum += pow(sigma0.Eval(Eg) - sig, 2) / (err * err);
  }
  infile.close();
  return sum;
}

double t1(const double * par){
  std::ifstream infile("dataset.dat");
  double sum = 0.0;
  double W, th, sig, err, err2;
  sigma0.SetParameters(par);
  while(infile >> W >> th >> sig >> err){
    err2 = err * err + 0.116 * 0.116 * sig * sig;
    if (W < 2.05) err2 = err2 + 0.36 * 0.36 * sig * sig;
    sum += pow(sigma4.Eval(W, th) - sig, 2) / err2;
  }
  infile.close();
  return sum;
}

double t0Fit(double * par, const char * minName = "Minuit2", const char * algoName = "Migrad"){
  const int NPAR = 3;
  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  min->SetMaxFunctionCalls(100000);
  min->SetTolerance(1.0e-6);
  min->SetPrintLevel(1);
  ROOT::Math::Functor f(&t0, NPAR);
  min->SetFunction(f);
  min->SetVariable(0, "A", 0.1, 1.0e-3);
  min->SetVariable(1, "B", 1.0, 1.0e-3);
  min->SetVariable(2, "C", 0.5, 1.0e-3);
  min->Minimize();
  const double * xs = min->X();
  for (int i = 0; i < NPAR; i++)
    par[i] = xs[i];
  const double chi2 = min->MinValue();
  sigma1.SetParameters(xs);
  //double total = sigma1.Integral(-1,1);
  std::cout << "chi2 = " << chi2 << std::endl;
  return chi2;
}

double t1Fit(double * par, const char * minName = "Minuit2", const char * algoName = "Migrad"){
  const int NPAR = 6;
  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  min->SetMaxFunctionCalls(100000);
  min->SetTolerance(1.0e-6);
  min->SetPrintLevel(1);
  ROOT::Math::Functor f(&t1, NPAR);
  min->SetFunction(f);
  min->SetFixedVariable(0, "A", 1.390);
  min->SetFixedVariable(1, "B", 0.347);
  min->SetFixedVariable(2, "C", 0.986);
  min->SetVariable(3, "D", 3.7, 1.0e-1);
  min->SetVariable(4, "E", 1.4, 1.0e-1);
  min->SetVariable(5, "F", 0.9, 1.0e-1);
  min->Minimize();
  const double * xs = min->X();
  for (int i = 0; i < NPAR; i++)
    par[i] = xs[i];
  const double chi2 = min->MinValue();
  sigma1.SetParameters(xs);
  //double total = sigma1.Integral(-1,1);
  std::cout << "chi2 = " << chi2 << std::endl;
  return chi2;
}

double t0plot(const char * savefile, const char * dataset){
  int Nt = 0;
  double W, Eg, ds, err;
  std::ifstream infile(dataset);
  while (infile >> W >> Eg >> ds >> err) Nt++;//count data points
  std::cout << Nt << std::endl;
  infile.clear();
  infile.seekg(0, std::ios::beg);//go back to the beginning of the datafile
  TGraphErrors * g0 = new TGraphErrors(Nt);
  double upper = 0;
  double lower = 0;
  for (int i = 0; i < Nt; i++){
    infile >> W >> Eg >> ds >> err;
    g0->SetPoint(i, Eg, ds);
    g0->SetPointError(i, 0, err);
    if (i == 0) {upper = ds + err; lower = ds - err;}
    if (upper < ds + err) upper = ds + err;
    if (lower > ds - err) lower = ds - err;
  }
  infile.close();
  TH1D * h0 = new TH1D("h0", "#Lambda(1520) photoproduction", 1, 1.5, 4.5);
  h0->Fill(0.0, -100.0);
  h0->GetYaxis()->SetRangeUser(0.0, lower+upper);
  h0->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  h0->GetXaxis()->CenterTitle();
  h0->GetYaxis()->SetTitle("#sigma (#mub)");
  h0->GetYaxis()->CenterTitle();
  h0->GetYaxis()->SetTitleOffset(1.5);
  TGraph * g1 = new TGraph(101);
  for (int i = 0; i < 101; i++){
    g1->SetPoint(i, 1.5 + 0.03 * i, sigma0.Eval(1.5 + 0.03*i));
    //std::cout << Es << " " << -1.0+0.02*i << " " << sigma2.Eval(Es, -1.0+0.02*i) << std::endl;
  }
  g1->SetLineColor(2);
  g1->SetLineWidth(2);
  g0->SetMarkerStyle(8);
  g0->SetMarkerSize(0.6);
  g0->SetMarkerColor(4);
  TCanvas * c0 = new TCanvas("", "", 800, 600);
  gStyle->SetOptStat(0);
  h0->Draw();
  g0->Draw("pesame");
  g1->Draw("csame");
  c0->Print(savefile);
  h0->Delete();
  g0->Delete();
  g1->Delete();
  c0->Close();
  return 0;
}

double f0(const double * par){
  std::ifstream infile("dataset.dat");
  double sum = 0.0;
  sigma1.SetParameters(par);
  double Es, cth, ds, err, err2;
  while (!infile.eof()){
    infile >> Es >> cth >> ds >> err;
    err2 = err*err + 0.116*0.116*ds*ds;
    if (Es < 2.05) err2 += 0.36*0.36*ds*ds;
    sum += pow(sigma1.Eval(cth) - ds, 2) / err2;
  }
  infile.close();
  return sum;
}

double f1(const double * par){
  std::ifstream infile("dataset.dat");
  double sum = 0.0;
  sigma2.SetParameters(par);
  double Es, cth, ds, err, err2;
  while (!infile.eof()){
    infile >> Es >> cth >> ds >> err;
    err2 = err*err + 0.116*0.116*ds*ds;
    if (Es < 2.05) err2 += 0.36*0.36*ds*ds;
    sum += pow(sigma2.Eval(Es, cth) - ds, 2) / err2;
  }
  infile.close();
  return sum;
}


double f0Fit(double * par, const char * minName = "Minuit2", const char * algoName = "Migrad"){
  const int NPAR = 3;
  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  min->SetMaxFunctionCalls(100000);
  min->SetTolerance(1.0e-6);
  min->SetPrintLevel(0);
  ROOT::Math::Functor f(&f0, NPAR);
  min->SetFunction(f);
  min->SetVariable(0, "A", 0.1, 1.0e-3);
  //min->SetFixedVariable(1, "b", 0.25);
  min->SetVariable(1, "B", 0.0, 1.0e-3);//, 0.0, 3.0);
  min->SetVariable(2, "C", 0.05, 1.0e-3);
  min->Minimize();
  const double * xs = min->X();
  for (int i = 0; i < NPAR; i++)
    par[i] = xs[i];
  const double chi2 = min->MinValue();
  sigma1.SetParameters(xs);
  double total = sigma1.Integral(-1,1);
  //std::cout << total << " " << chi2 << " " << xs[0] << "   " << xs[1] << "   " << xs[2] << std::endl;
  return chi2;
}

double f1Fit(double * par, const char * minName = "Minuit2", const char * algoName = "Migrad"){
  const int NPAR = 5;
  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  min->SetMaxFunctionCalls(100000);
  min->SetTolerance(1.0e-6);
  min->SetPrintLevel(0);
  ROOT::Math::Functor f(&f1, NPAR);
  min->SetFunction(f);
  min->SetLowerLimitedVariable(0, "a0", 0.2267, 1.0e-3, 0.0);
  min->SetLimitedVariable(1, "a1", 2.06, 1.0e-3, 0.0, 10.0);
  min->SetLimitedVariable(2, "a2", 0.2, 1.0e-3, 0.0, 2.0);
  min->SetVariable(3, "b1", 2.0, 1.0e-3);
  min->SetVariable(4, "b2", 1.75, 1.0e-3);
  min->Minimize();
  const double * xs = min->X();
  for (int i = 0; i < NPAR; i++)
    par[i] = xs[i];
  const double chi2 = min->MinValue();
  sigma1.SetParameters(xs);
  double total = sigma1.Integral(-1,1);
  std::cout << total << " " << chi2 << " " << xs[0] << "   " << xs[1] << "   " << xs[2] << "   " << xs[3] << "   " << xs[4] << std::endl;
  return chi2;
}

int parameterplot(const char * parset){
  double par[4];
  int Np = -1;
  std::ifstream infile(parset);
  char tmp[256];
  while(infile.getline(tmp, 256)){
    Np++;
  }
  infile.clear();
  infile.seekg(0, std::ios::beg);
  TH1D * h0 = new TH1D("h0", "", 1, 2.0, 3.0);
  h0->GetYaxis()->SetRangeUser(0.0, 3.0);
  TGraph * g0 = new TGraph(Np);
  infile.getline(tmp, 256);
  for (int i = 0; i < Np - 1; i++){
    infile >> par[0] >> par[1] >> par[2] >> par[3];
    std::cout << par[0] << " " << par[3] << std::endl;
    g0->SetPoint(i, par[0], par[3]);
  }
  g0->SetMarkerStyle(8);
  g0->SetMarkerSize(0.6);
  g0->SetMarkerColor(4);
  gStyle->SetOptStat(0);
  TCanvas * c0 = new TCanvas("c0", "", 800, 600);
  h0->Draw();
  g0->Draw("psame");
  c0->Print("parplot.pdf");
  h0->Delete();
  g0->Delete();
  c0->Close();
  return 0;
}




#endif
