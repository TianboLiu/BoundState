#include "Lcore.h"

int main(const int argc, const char * argv[]){

  Initialize();

  TLorentzVector ki[2], kf[5];
  double weight = 0;
  ki[0].SetXYZT(0, 0, sqrt(4.4 * 4.4 - PARTICLE::e.M() * PARTICLE::e.M()), 4.4);
  //ki[0].SetXYZT(0, 0, 4.4, 4.4);

  TFile * fs = new TFile("result/detected.root", "RECREATE");

  TH1D * h0 = new TH1D("MpKK_BoundStateAll", "", 2200, 1.88, 2.32);
  TH1D * h1 = new TH1D("MpKK_BoundStateKK", "", 2200, 1.88, 2.32);
  TH1D * h2 = new TH1D("MpKK_phi", "", 2200, 1.88, 2.32);
  TH1D * h3 = new TH1D("MpKK_Lambda1520", "", 2200, 1.88, 2.32);
  TH1D * h4 = new TH1D("MpKK_directKK", "", 2200, 1.88, 2.32);

  TH1D * h0a = new TH1D("MpKp_BoundStateAll", "", 2200, 1.38, 1.82);
  TH1D * h1a = new TH1D("MpKp_BoundStateKK", "", 2200, 1.38, 1.82);
  TH1D * h2a = new TH1D("MpKp_phi", "", 2200, 1.38, 1.82);
  TH1D * h3a = new TH1D("MpKp_Lambda1520", "", 2200, 1.38, 1.82);
  TH1D * h4a = new TH1D("MpKp_directKK", "", 2200, 1.38, 1.82);

  TH1D * h0b = new TH1D("MpKm_BoundStateAll", "", 2200, 1.38, 1.82);
  TH1D * h1b = new TH1D("MpKm_BoundStateKK", "", 2200, 1.38, 1.82);
  TH1D * h2b = new TH1D("MpKm_phi", "", 2200, 1.38, 1.82);
  TH1D * h3b = new TH1D("MpKm_Lambda1520", "", 2200, 1.38, 1.82);
  TH1D * h4b = new TH1D("MpKm_directKK", "", 2200, 1.38, 1.82);

  TH1D * h0c = new TH1D("MKK_BoundStateAll", "", 2200, 0.98, 1.42);
  TH1D * h1c = new TH1D("MKK_BoundStateKK", "", 2200, 0.98, 1.42);
  TH1D * h2c = new TH1D("MKK_phi", "", 2200, 0.98, 1.42);
  TH1D * h3c = new TH1D("MKK_Lambda1520", "", 2200, 0.98, 1.42);
  TH1D * h4c = new TH1D("MKK_directKK", "", 2200, 0.98, 1.42);

  h0->SetDirectory(fs);
  h1->SetDirectory(fs);
  h2->SetDirectory(fs);
  h3->SetDirectory(fs);
  h4->SetDirectory(fs);

  h0a->SetDirectory(fs);
  h1a->SetDirectory(fs);
  h2a->SetDirectory(fs);
  h3a->SetDirectory(fs);
  h4a->SetDirectory(fs);

  h0b->SetDirectory(fs);
  h1b->SetDirectory(fs);
  h2b->SetDirectory(fs);
  h3b->SetDirectory(fs);
  h4b->SetDirectory(fs);

  h0c->SetDirectory(fs);
  h1c->SetDirectory(fs);
  h2c->SetDirectory(fs);
  h3c->SetDirectory(fs);
  h4c->SetDirectory(fs);

  TH2D * d0a = new TH2D("Momentum_p_Kp_BoundStateAll", "", 200, 0.0, 2.0, 200, 0.0, 2.0);
  TH2D * d0b = new TH2D("Momentum_p_Km_BoundStateAll", "", 200, 0.0, 2.0, 200, 0.0, 2.0);
  TH2D * d0c = new TH2D("Momentum_Kp_Km_BoundStateAll", "", 200, 0.0, 2.0, 200, 0.0, 2.0);
  TH2D * d1a = new TH2D("Momentum_p_Kp_BoundStateKK", "", 200, 0.0, 2.0, 200, 0.0, 2.0);
  TH2D * d1b = new TH2D("Momentum_p_Km_BoundStateKK", "", 200, 0.0, 2.0, 200, 0.0, 2.0);
  TH2D * d1c = new TH2D("Momentum_Kp_Km_BoundStateKK", "", 200, 0.0, 2.0, 200, 0.0, 2.0);
  TH2D * d2a = new TH2D("Momentum_p_Kp_phi", "", 200, 0.0, 2.0, 200, 0.0, 2.0);
  TH2D * d2b = new TH2D("Momentum_p_Km_phi", "", 200, 0.0, 2.0, 200, 0.0, 2.0);
  TH2D * d2c = new TH2D("Momentum_Kp_Km_phi", "", 200, 0.0, 2.0, 200, 0.0, 2.0);
  TH2D * d3a = new TH2D("Momentum_p_Kp_Lambda1520", "", 200, 0.0, 2.0, 200, 0.0, 2.0);
  TH2D * d3b = new TH2D("Momentum_p_Km_Lambda1520", "", 200, 0.0, 2.0, 200, 0.0, 2.0);
  TH2D * d3c = new TH2D("Momentum_Kp_Km_Lambda1520", "", 200, 0.0, 2.0, 200, 0.0, 2.0);
  TH2D * d4a = new TH2D("Momentum_p_Kp_KK", "", 200, 0.0, 2.0, 200, 0.0, 2.0);
  TH2D * d4b = new TH2D("Momentum_p_Km_KK", "", 200, 0.0, 2.0, 200, 0.0, 2.0);
  TH2D * d4c = new TH2D("Momentum_Kp_Km_KK", "", 200, 0.0, 2.0, 200, 0.0, 2.0);

  d0a->SetDirectory(fs);
  d0b->SetDirectory(fs);
  d0c->SetDirectory(fs);
  d1a->SetDirectory(fs);
  d1b->SetDirectory(fs);
  d1c->SetDirectory(fs);
  d2a->SetDirectory(fs);
  d2b->SetDirectory(fs);
  d2c->SetDirectory(fs);
  d3a->SetDirectory(fs);
  d3b->SetDirectory(fs);
  d3c->SetDirectory(fs);
  d4a->SetDirectory(fs);
  d4b->SetDirectory(fs);
  d4c->SetDirectory(fs);
  

  Long64_t Nsim = 100000000;

  TLorentzVector PP, Pa, Pb, Pc;
  for (Long64_t i = 0; i < Nsim; i++){
    if (i%1000000 == 0) cout << i << endl;
 
    weight = GENERATE::Event_eNKKN_BoundState(ki, kf);
    if (weight > 0){
      double factorp = DETECTOR::Acceptance(kf[4], "p");
      weight *= DETECTOR::Acceptance(kf[2], "K+") * DETECTOR::Acceptance(kf[3], "K-");
      PP = kf[2] + kf[3] + kf[4];
      Pa = kf[2] + kf[4];
      Pb = kf[3] + kf[4];
      Pc = kf[2] + kf[3];
      h0->Fill(PP.M(), weight * factorp);
      h0a->Fill(Pa.M(), weight * factorp);
      h0b->Fill(Pb.M(), weight * factorp);
      h0c->Fill(Pc.M(), weight * factorp);
      d0a->Fill(kf[2].P(), kf[4].P(), weight * factorp);
      d0b->Fill(kf[3].P(), kf[4].P(), weight * factorp);
      d0c->Fill(kf[2].P(), kf[3].P(), weight * factorp);
      factorp = DETECTOR::Acceptance(kf[1], "p");
      PP = kf[1] + kf[2] + kf[3];
      Pa = kf[1] + kf[2];
      Pb = kf[1] + kf[3];
      Pc = kf[2] + kf[3];
      h1->Fill(PP.M(), weight * factorp);
      h1a->Fill(Pa.M(), weight * factorp);
      h1b->Fill(Pb.M(), weight * factorp);
      h1c->Fill(Pc.M(), weight * factorp);
      d1a->Fill(kf[2].P(), kf[1].P(), weight * factorp);
      d1b->Fill(kf[3].P(), kf[1].P(), weight * factorp);
      d1c->Fill(kf[2].P(), kf[3].P(), weight * factorp);
    }

    weight = GENERATE::Event_eNKK_Phi(ki, kf);
    if (weight > 0){
      weight *= DETECTOR::Acceptance(kf[1], "p") * DETECTOR::Acceptance(kf[2], "K+") * DETECTOR::Acceptance(kf[3], "K-");
      PP = kf[1] + kf[2] + kf[3];
      Pa = kf[1] + kf[2];
      Pb = kf[1] + kf[3];
      Pc = kf[2] + kf[3];
      h2->Fill(PP.M(), weight);
      h2a->Fill(Pa.M(), weight);
      h2b->Fill(Pb.M(), weight);
      h2c->Fill(Pc.M(), weight);
      d2a->Fill(kf[2].P(), kf[1].P(), weight);
      d2b->Fill(kf[3].P(), kf[1].P(), weight);
      d2c->Fill(kf[2].P(), kf[3].P(), weight);
    }

    weight = GENERATE::Event_eNKK_L1520(ki, kf);
    if (weight > 0){
      weight *= DETECTOR::Acceptance(kf[1], "p") * DETECTOR::Acceptance(kf[2], "K+") * DETECTOR::Acceptance(kf[3], "K-");
      PP = kf[1] + kf[2] + kf[3];
      Pa = kf[1] + kf[2];
      Pb = kf[1] + kf[3];
      Pc = kf[2] + kf[3];
      h3->Fill(PP.M(), weight);
      h3a->Fill(Pa.M(), weight);
      h3b->Fill(Pb.M(), weight);
      h3c->Fill(Pc.M(), weight);
      d3a->Fill(kf[2].P(), kf[1].P(), weight);
      d3b->Fill(kf[3].P(), kf[1].P(), weight);
      d3c->Fill(kf[2].P(), kf[3].P(), weight);
    }
 
    weight = GENERATE::Event_eNKK_KK(ki, kf);
    if (weight > 0){
      weight *= DETECTOR::Acceptance(kf[1], "p") * DETECTOR::Acceptance(kf[2], "K+") * DETECTOR::Acceptance(kf[3], "K-");
      PP = kf[1] + kf[2] + kf[3];
      Pa = kf[1] + kf[2];
      Pb = kf[1] + kf[3];
      Pc = kf[2] + kf[3];
      h4->Fill(PP.M(), weight);
      h4a->Fill(Pa.M(), weight);
      h4b->Fill(Pb.M(), weight);
      h4c->Fill(Pc.M(), weight);
      d4a->Fill(kf[2].P(), kf[1].P(), weight);
      d4b->Fill(kf[3].P(), kf[1].P(), weight);
      d4c->Fill(kf[2].P(), kf[3].P(), weight);
    }

  }

  h0->Scale(1.0/Nsim);
  h1->Scale(1.0/Nsim);
  h2->Scale(1.0/Nsim);
  h3->Scale(1.0/Nsim);
  h4->Scale(1.0/Nsim);

  h0a->Scale(1.0/Nsim);
  h1a->Scale(1.0/Nsim);
  h2a->Scale(1.0/Nsim);
  h3a->Scale(1.0/Nsim);
  h4a->Scale(1.0/Nsim);

  h0b->Scale(1.0/Nsim);
  h1b->Scale(1.0/Nsim);
  h2b->Scale(1.0/Nsim);
  h3b->Scale(1.0/Nsim);
  h4b->Scale(1.0/Nsim);

  h0c->Scale(1.0/Nsim);
  h1c->Scale(1.0/Nsim);
  h2c->Scale(1.0/Nsim);
  h3c->Scale(1.0/Nsim);
  h4c->Scale(1.0/Nsim);

  d0a->Scale(1.0/Nsim);
  d0b->Scale(1.0/Nsim);
  d0c->Scale(1.0/Nsim);
  d1a->Scale(1.0/Nsim);
  d1b->Scale(1.0/Nsim);
  d1c->Scale(1.0/Nsim);
  d2a->Scale(1.0/Nsim);
  d2b->Scale(1.0/Nsim);
  d2c->Scale(1.0/Nsim);
  d3a->Scale(1.0/Nsim);
  d3b->Scale(1.0/Nsim);
  d3c->Scale(1.0/Nsim);
  d4a->Scale(1.0/Nsim);
  d4b->Scale(1.0/Nsim);
  d4c->Scale(1.0/Nsim);

 
  fs->Write();

  return 0;
}
