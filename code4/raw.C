#include "Lcore.h"

int main(const int argc, const char * argv[]){

  Initialize();

  TLorentzVector ki[2], kf[5];
  double weight = 0;
  //ki[0].SetXYZT(0, 0, sqrt(4.4 * 4.4 - PARTICLE::e.M() * PARTICLE::e.M()), 4.4);
  ki[0].SetXYZT(0, 0, 4.4, 4.4);

  TFile * fs = new TFile("result/raw.root", "RECREATE");

  TH1D * h0 = new TH1D("MpKK_BoundStateAll", "", 2200, 1.88, 2.32);
  TH1D * h1 = new TH1D("MpKK_BoundStateKK", "", 2200, 1.88, 2.32);
  TH1D * h2 = new TH1D("MpKK_phi", "", 2200, 1.88, 2.32);
  TH1D * h3 = new TH1D("MpKK_Lambda1520", "", 2200, 1.88, 2.32);
  TH1D * h4 = new TH1D("MpKK_directKK", "", 2200, 1.88, 2.32);

  h0->SetDirectory(fs);
  h1->SetDirectory(fs);
  h2->SetDirectory(fs);
  h3->SetDirectory(fs);
  h4->SetDirectory(fs);

  TH2D * d0 = new TH2D("Momentum_p_Kp_BoundStateAll", "", 200, 0.0, 2.0, 200, 0.0, 2.0);
  TH2D * d1 = new TH2D("Momentum_p_Kp_BoundStateKK", "", 200, 0.0, 2.0, 200, 0.0, 2.0);
  TH2D * d2 = new TH2D("Momentum_p_Kp_phi", "", 200, 0.0, 2.0, 200, 0.0, 2.0);
  TH2D * d3a = new TH2D("Momentum_p_Kp_Lambda1520", "", 200, 0.0, 2.0, 200, 0.0, 2.0);
  TH2D * d3b = new TH2D("Momentum_p_Km_Lambda1520", "", 200, 0.0, 2.0, 200, 0.0, 2.0);
  TH2D * d4 = new TH2D("Momentum_p_Kp_KK", "", 200, 0.0, 2.0, 200, 0.0, 2.0);

  d0->SetDirectory(fs);
  d1->SetDirectory(fs);
  d2->SetDirectory(fs);
  d3a->SetDirectory(fs);
  d3b->SetDirectory(fs);
  d4->SetDirectory(fs);

  

  Long64_t Nsim = 100000000;

  TLorentzVector PP;
  for (Long64_t i = 0; i < Nsim; i++){
    if (i%1000000 == 0) cout << i << endl;
 
    weight = GENERATE::Event_eNKKN_BoundState(ki, kf);
    if (weight > 0){
      PP = kf[2] + kf[3] + kf[4];
      h0->Fill(PP.M(), weight);
      d0->Fill(kf[2].P(), kf[4].P(), weight);
      PP = kf[1] + kf[2] + kf[3];
      h1->Fill(PP.M(), weight);
      d1->Fill(kf[2].P(), kf[1].P(), weight);
    }

    weight = GENERATE::Event_eNKK_Phi(ki, kf);
    if (weight > 0){
      PP = kf[1] + kf[2] + kf[3];
      h2->Fill(PP.M(), weight);
      d2->Fill(kf[2].P(), kf[1].P(), weight);
    }

    weight = GENERATE::Event_eNKK_L1520(ki, kf);
    if (weight > 0){
      PP = kf[1] + kf[2] + kf[3];
      h3->Fill(PP.M(), weight);
      d3a->Fill(kf[2].P(), kf[1].P(), weight);
      d3b->Fill(kf[3].P(), kf[1].P(), weight);
    }
 
    weight = GENERATE::Event_eNKK_KK(ki, kf);
    if (weight > 0){
      PP = kf[1] + kf[2] + kf[3];
      h4->Fill(PP.M(), weight);
      d4->Fill(kf[2].P(), kf[1].P(), weight);
    }

  }

  h0->Scale(1.0/Nsim);
  h1->Scale(1.0/Nsim);
  h2->Scale(1.0/Nsim);
  h3->Scale(1.0/Nsim);
  h4->Scale(1.0/Nsim);

  d0->Scale(1.0/Nsim);
  d1->Scale(1.0/Nsim);
  d2->Scale(1.0/Nsim);
  d3a->Scale(1.0/Nsim);
  d3b->Scale(1.0/Nsim);
  d4->Scale(1.0/Nsim);

  
 
  fs->Write();

  return 0;
}
