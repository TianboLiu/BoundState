#include <iostream>
#include <fstream>
#include <cmath>

#include "TStyle.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1D.h"
#include "TH2D.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TString.h"

using namespace std;

class Data{
 private:
  double value;
  double cth;
  double E;
  double stat;
  double syst;
 public:
  int SetData(const double Value, const double Cth, const double Ec, const double Stat = 1.0, const double Syst = 0.0){
    value = Value;
    cth = Cth;
    E = Ec;
    stat = Stat;
    syst = Syst;
    return 0;
  }
  
  double Value(){return value;}
  double Cth(){return cth;}
  double Es(){return E;}
  double Stat(){return stat;}
  double Syst(){return syst;}
      
};


double dsigma(const double * var, const double * par){
  double x = var[0];
  double result = (par[0] + par[1] * x + par[2] * x * x) * exp(par[3] * (x - 1.0) + par[4] * pow(x - 1.0, 2));
  return result;
}

Data data[100];
int Nd = 0; 

double FitFunction(const double * par){
  double sum = 0;
  double th = 0;
  double cth = 0;
  for (int i = 0; i < Nd; i++){
    cth = data[i].Cth();
    th = dsigma(&cth, par);
    sum += pow( th - data[i].Value(), 2) / (data[i].Stat() * data[i].Stat() + data[i].Syst() * data[i].Syst());
  }
  return sum;
}

double Fit(double * par, const char * minName = "Minuit2", const char * algoName = "Migrad"){
  if (Nd == 0) return 0;
  const int NPAR = 5;
  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  min->SetMaxFunctionCalls(100000);
  min->SetTolerance(1.0e-3);
  min->SetPrintLevel(0);
  ROOT::Math::Functor f(&FitFunction, NPAR);
  min->SetFunction(f);
  min->SetVariable(0, "N", 0.5, 1e-3);
  min->SetVariable(1, "a1", 1.0, 1e-3);
  //min->SetFixedVariable(1, "a1", 0.0);
  min->SetVariable(2, "a2", 0.0, 1e-3);
  //min->SetFixedVariable(2, "a2", 0.0);
  min->SetVariable(3, "b1", 0.1, 1e-3);
  min->SetVariable(4, "b2", 0.0, 1e-3);
  //min->SetFixedVariable(4, "b2", 0.0);
  min->Minimize();
  const double * xs = min->X();
  for (int i = 0; i < NPAR; i++)
    par[i] = xs[i];
  const double chi2 = min->MinValue();
  return chi2;
}



int CreateSet(const double Emin, const double Emax){
  Nd = 0;
  ifstream infile("clasdb_E63M2.txt");
  char tmp[256];
  for (int i = 0; i < 8; i++) infile.getline(tmp, 256);
  double Es, cth, ds, err;
  while (infile >> Es >> cth >> ds >> err){
    if (Es > Emin && Es < Emax){
      cth = ((int) ((cth + 0.01)/0.025)) * 0.025;
      data[Nd].SetData(ds, cth, Es, err);
      Nd++;
    }
  }
  infile.close();
  cout << "#Data points:  " << Nd << endl;
  return 0;
}
  
int SetStyle(TGraphErrors * g0){
  g0->SetMarkerStyle(20);
  g0->SetMarkerSize(0.5);
  g0->SetLineWidth(1.0);
  g0->SetMarkerColor(4);
  return 0;
}

int main(int argc, char * argv[]){

  double par[5];
  double Es = 1.985;
  double chi2 = 0;
  double cth = 0;
  
  TCanvas * c0 = new TCanvas("c0", "", 800, 600);
  TGraphErrors * g0;
  TF1 f0("ds", dsigma, -1.0, 1.0, 5);
  TH1D * h0 = new TH1D("h0", "", 1, -1.0, 1.0);
  h0->SetMinimum(0.0);
  h0->SetMaximum(1.0);
  h0->SetStats(0);
  c0->Print("dsphi.pdf(", "pdf");

  FILE * f = fopen("thphi.dat", "w");
  FILE * fexp = fopen("exphi.dat", "w");

  while(Es < 2.8351){
    CreateSet(Es - 1e-6, Es + 1e-6);
    chi2 = Fit(par);
    cout << Es << " : " << par[0] << " " <<  par[1] << " " << par[2] << " " << par[3] << " " << par[4] << "  (chi2 = " << chi2 << ")" << endl;

    f0.SetParameters(par);
    f0.SetLineColor(2);
    f0.SetLineWidth(2);
    g0 = new TGraphErrors(Nd);
    for (int i = 0; i < Nd; i++){
      g0->SetPoint(i, data[i].Cth(), data[i].Value());
      g0->SetPointError(i, 0.0, data[i].Stat());
    }
    SetStyle(g0);
    c0->Clear();
    h0->SetTitle(Form("%.3f", Es));
    h0->Draw();
    f0.DrawClone("lsame");
    g0->DrawClone("pesame");
    c0->Print("dsphi.pdf", "pdf");
    
    fprintf(f, "%.6f\t%.6f\t%.6f\n", Es, -1.0, 0.0);
    fprintf(fexp, "%.6f\t%.6f\t%.6f\n", Es, -1.0, 0.0);
    for (int i = 0; ; i++){
      cth = -0.375 + 0.025 * i;
      if (cth > 0.99) break;
      fprintf(f, "%.6f\t%.6f\t%.6f\n", Es, cth, f0.Eval(cth));
    }
    for (int j = 0; ; j++){
      cth = -0.975 + 0.025 * j;
      if (cth > 0.99) break;
      for (int i = 0; i < Nd; i++){
	if (data[i].Cth() > cth - 1.0e-5 && data[i].Cth() < cth + 1.0e-5)
	  fprintf(fexp, "%.6f\t%.6f\t%.6f\n", data[i].Es(), data[i].Cth(), data[i].Value());
      }
    }
    fprintf(f, "%.6f\t%.6f\t%.6f\n", Es, 1.0, f0.Eval(1.0));
    fprintf(fexp, "%.6f\t%.6f\t%.6f\n", Es, 1.0, f0.Eval(1.0));
    Es += 0.01;
  }
  c0->Clear();
  c0->Print("dsphi.pdf)", "pdf");

  fclose(f);
  fclose(fexp);

  return 0;
}
