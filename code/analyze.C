#include <iostream>
#include <fstream>
#include <cmath>

#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TMath.h"

using namespace std;

int analyze(){
  const double degtorad = TMath::Pi() / 180.0;
  const double radtodeg = 180.0 / TMath::Pi();
 
  TString simfile = "/var/phy/project/mepg/tl190/sim2.root";
  TFile * fs = new TFile(simfile, "r");
  TTree * Ts = (TTree *) fs->Get("data");
  double Nt = Ts->GetEntries();
  double ds, Aw, Bw;
  TLorentzVector * Pd = 0;
  TLorentzVector * AP0 = 0;
  TLorentzVector * AP1 = 0;
  TLorentzVector * AP2 = 0;
  TLorentzVector * BP0 = 0;
  TLorentzVector * BP1 = 0;
  TLorentzVector * BP2 = 0;
  Ts->SetBranchAddress("ds", &ds);
  Ts->SetBranchAddress("Aw", &Aw);
  Ts->SetBranchAddress("Bw", &Bw);
  Ts->SetBranchAddress("Pd", &Pd);
  Ts->SetBranchAddress("AP0", &AP0);
  Ts->SetBranchAddress("AP1", &AP1);
  Ts->SetBranchAddress("AP2", &AP2);
  Ts->SetBranchAddress("BP0", &BP0);
  Ts->SetBranchAddress("BP1", &BP1);
  Ts->SetBranchAddress("BP2", &BP2);

  //normalization
  double nAw = 0.390587;
  double nBw = 1.0;

  TH1D * h1 = new TH1D("h1", "h1", 5, 0.0, 5.0);
  double theta[3];
  double thmin = 5.0 * degtorad;
  double thmax1 = 35.0 * degtorad;
  double thmax2 = 125.0 * degtorad;
  double lumi = 5.8e+30 * 1.0e-26 * pow(0.197327, 2);
  double time = 3600.0;
  double volume = 1.0;
  double ctK = 3.712;
  double ctpi = 7.8045;
  double LOF = 2.0;
  for (int i = 0; i < Nt; i++){
    Ts->GetEntry(i);
    theta[0] = AP0->Theta();
    theta[1] = AP1->Theta();
    theta[2] = AP2->Theta();
    h1->Fill(0.5, ds*Aw/nAw*exp(-LOF/(AP1->Beta()*AP1->Gamma()*ctK))*exp(-LOF/(AP2->Beta()*AP2->Gamma()*ctK)));
    if ( (theta[0] > thmin && theta[0] < thmax2) && (theta[1] > thmin && theta[1] < thmax2) && (theta[2] > thmin && theta[2] < thmax2))
      h1->Fill(1.5, ds*Aw/nAw*exp(-LOF/(AP1->Beta()*AP1->Gamma()*ctK))*exp(-LOF/(AP2->Beta()*AP2->Gamma()*ctK)));
    theta[0] = BP0->Theta();
    theta[1] = BP1->Theta();
    theta[2] = BP2->Theta();
    h1->Fill(2.5, ds*Bw/nBw*exp(-LOF/(AP1->Beta()*AP1->Gamma()*ctpi))*exp(-LOF/(AP2->Beta()*AP2->Gamma()*ctK)));
    if ( (theta[0] > thmin && theta[0] < thmax2) && (theta[1] > thmin && theta[1] < thmax1) && (theta[2] > thmin && theta[2] < thmax2))
      h1->Fill(3.5, ds*Bw/nBw*exp(-LOF/(AP1->Beta()*AP1->Gamma()*ctpi))*exp(-LOF/(AP2->Beta()*AP2->Gamma()*ctK)));
    h1->Fill(4.5, ds);
  }
  cout << h1->Integral(5, 5) * volume / Nt *3.89379e5 << endl;
  h1->Scale(0.5*lumi*time*volume/Nt);

  double BrA = 0.904;
  double BrB = 0.067*0.639;
  cout << h1->Integral(1,1) * BrA << " " << h1->Integral(2,2) * BrA << endl;
  cout << h1->Integral(3,3) * BrB << " " << h1->Integral(4,4) * BrB << endl;
  
  return 0;
}

#ifndef __CINT__
int main(){
  return analyze();
}
#endif
