#include "Lsimulation.h"
#include "TH1D.h"
#include "TCanvas.h"

using namespace std;

int main(){
  gRandom->SetSeed(1);
  SetFunctions();
  TLorentzVector q(0.0, 0.0, 1.45, 1.45);
  TLorentzVector p1;
  TLorentzVector p, k;
  TLorentzVector ki[2], kf[4];
  TLorentzVector P0;
  double weight[5];
  ki[0] = q;
  TH1D * h0 = new TH1D("h0", "", 1000, 1.8, 2.5);
  TH1D * h1 = new TH1D("h1", "", 1000, 1.8, 2.5);
  TH1D * h2 = new TH1D("h2", "", 1000, 1.8, 2.5);
  Long64_t Nsim = 20000000;
  TLorentzVector pc;
  for (int i = 0; i < Nsim; i++){
    if (i%100000==0) cout << i << endl;
    //GenerateNucleonInCarbon(&p1);
    //ki[1] = p1;
    //GeneratePhiProduction(ki, kf, &sigma);
    //P0 = kf[0] + kf[1];
    //h0->Fill(P0.M(), sigma);
    GenerateEventNNKKwithBoundState(ki, kf, weight);
    pc = kf[1] + kf[2] + kf[3];
    h0->Fill(pc.M(), weight[0]);
    pc = kf[0] + kf[2] + kf[3];
    h1->Fill(pc.M(), weight[0]);
    GenerateEventNKKwithoutBoundState(ki, kf, weight);
    pc = kf[0] + kf[1] + kf[2];
    h2->Fill(pc.M(), weight[0]/10 );
  }
  gStyle->SetOptStat(0);
  h0->GetXaxis()->SetTitle("M(PK^{+}K^{-}) / GeV");
  h0->GetXaxis()->CenterTitle();
  TCanvas * c0 = new TCanvas("c0", "c0", 800, 600);
  h0->SetLineColor(1);
  h0->Draw();
  h1->SetLineColor(4);
  h1->Draw("same");
  h2->SetLineColor(2);
  h2->Draw("same");

  c0->Print("c0.pdf");

  cout << h0->Integral(1, -1) / Nsim << endl;
  cout << h1->Integral(1, -1) / Nsim << endl;
  cout << h2->Integral(1, -1) / Nsim << endl;

  return 0;
}
