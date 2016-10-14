#include "Lsimulation.h"


using namespace std;

int main(){
  gRandom->SetSeed(1);
  SetFunctions();
  const Long64_t Nsim = 10000000;
  const int Npt = 5000;

  double res[3] = {0.01, 0.001, 0.004};//resolution dp/p, dth, dphi

  TLorentzVector q;
  TLorentzVector ki[2], kf[4];
  double weight[5];

  TH1D * b1D = new TH1D("b1D", "", 1, 1.8, 2.5);
  TH2D * b2D = new TH2D("b2D", "", 1, 0.0, 1.2, 1, 0.0, 1.5);
  gStyle->SetOptStat(0);
  b1D->SetTitle("M(pK^{+}K^{-}) spectra");
  b1D->GetXaxis()->SetTitle("M(pK^{+}K^{-}) / GeV");
  b1D->GetXaxis()->CenterTitle();
  b1D->GetXaxis()->SetTitleOffset(1.0);
  b1D->GetXaxis()->SetTitleSize(0.05);
  b1D->GetYaxis()->SetRangeUser(0.0, 1.0);
  b2D->SetTitle("proton kaon momentum correlation");
  b2D->GetXaxis()->SetTitle("P(K) / GeV");
  b2D->GetXaxis()->CenterTitle();
  b2D->GetXaxis()->SetTitleOffset(1.0);
  b2D->GetXaxis()->SetTitleSize(0.05);
  b2D->GetYaxis()->SetTitle("P(p) / GeV");
  b2D->GetYaxis()->CenterTitle();
  b2D->GetYaxis()->SetTitleOffset(1.0);
  b2D->GetYaxis()->SetTitleSize(0.05); 

  TH1D * h0 = new TH1D("h0", "", 1000, 1.8, 2.5);
  TH1D * h1 = new TH1D("h1", "", 1000, 1.8, 2.5);
  TH1D * h2 = new TH1D("h2", "", 1000, 1.8, 2.5);
  h0->SetLineColor(1);
  h1->SetLineColor(4);
  h2->SetLineColor(2);

  TH2D * d0 = new TH2D("d0", "", 500, 0.0, 1.5, 500, 0.0, 1.5);
  TH2D * d1 = new TH2D("d1", "", 500, 0.0, 1.5, 500, 0.0, 1.5);
  TH2D * d2 = new TH2D("d2", "", 500, 0.0, 1.5, 500, 0.0, 1.5);

  TGraph * g0 = new TGraph(Npt);
  TGraph * g1 = new TGraph(Npt);
  TGraph * g2 = new TGraph(Npt);
  g0->SetMarkerColor(1);
  g0->SetMarkerStyle(20);
  g0->SetMarkerSize(0.3);
  g1->SetMarkerColor(4);
  g1->SetMarkerStyle(24);
  g1->SetMarkerSize(0.2);
  g2->SetMarkerColor(2);
  g2->SetMarkerStyle(26);
  g2->SetMarkerSize(0.2);

  double x, y;

  TLorentzVector l0, l1, l2;
  for (int i = 0; i < Nsim; i++){
    if (i%100000==0) cout << i << endl;
    GenerateBremsstrahlungPhoton(&q);//Generate a photon
    ki[0] = q;//Set initial photon
    if (q.P() < 1.6){
      GenerateEventNNKKwithBoundState(ki, kf, weight);//Generate a bound state event kf: N (N, K+, K-)
      DetectorResolutionSmear(&kf[0], res);//smear final particles
      DetectorResolutionSmear(&kf[1], res);
      DetectorResolutionSmear(&kf[2], res);
      DetectorResolutionSmear(&kf[3], res);
      l0 = kf[1] + kf[2] + kf[3];//pKK from bound
      l1 = kf[0] + kf[2] + kf[3];//KK + recoil p
      h0->Fill(l0.M(), weight[0]);
      h1->Fill(l1.M(), weight[0]);
      d0->Fill(kf[2].P(), kf[1].P(), weight[0]);
      d1->Fill(kf[2].P(), kf[0].P(), weight[0]);
      GenerateEventNKKwithoutBoundState(ki, kf, weight);//Generate a unbound event with phi production
      DetectorResolutionSmear(&kf[0], res);//smear final particles
      DetectorResolutionSmear(&kf[1], res);
      DetectorResolutionSmear(&kf[2], res);
      l2 = kf[0] + kf[1] + kf[2];//pKK unbound
      h2->Fill(l2.M(), weight[0]);
      d2->Fill(kf[1].P(), kf[0].P(), weight[0]);
    }
  }
  double hmax = h2->GetMaximum();
  h0->Scale(0.8/hmax);
  h1->Scale(0.8/hmax);
  h2->Scale(0.8/hmax);

  for (int i = 0; i < Npt; i++){
    d0->GetRandom2(x, y);
    g0->SetPoint(i, x, y);
    d1->GetRandom2(x, y);
    g1->SetPoint(i, x, y);
    d2->GetRandom2(x, y);
    g2->SetPoint(i, x, y);
  }

  TCanvas * c0 = new TCanvas("c0", "c0", 800, 600);
  b1D->Draw();
  h0->Draw("same");
  h1->Draw("same");
  h2->Draw("same");
  c0->Print("c0.pdf");


  hmax = h0->GetMaximum();
  h0->Scale(0.8/hmax);
  h1->Scale(0.8/hmax);
  h2->Scale(0.8/hmax/10);
  TCanvas * c1 = new TCanvas("c1", "c1", 800, 600);
  b1D->Draw();
  h0->Draw("same");
  h1->Draw("same");
  h2->Draw("same");
  c1->Print("c1.pdf");

  TCanvas * c2 = new TCanvas("c2", "c2", 800, 600);
  b2D->Draw();
  g0->Draw("psame");
  g1->Draw("psame");
  g2->Draw("psame");
  c2->Print("c2.pdf");

  return 0;
}
