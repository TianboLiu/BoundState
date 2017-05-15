#include "Lcore.h"

int MisPID(TLorentzVector * k);

int main(const int argc, const char * argv[]){

  Long64_t Nsim;

  if (argc < 2) Nsim = 100000000;
  else Nsim = atoi(argv[1]);

  Initialize();

  TLorentzVector ki[2], kf[5];
  double weight = 0;
  ki[0].SetXYZT(0, 0, sqrt(4.4 * 4.4 - PARTICLE::e.M() * PARTICLE::e.M()), 4.4);

  TFile * fs = new TFile("result/twopi.root", "RECREATE");

  TH1D * h0 = new TH1D("MpPiPi", "", 2200, 1.18, 1.62);
  TH1D * h1 = new TH1D("MpKK_TwoPi", "", 2200, 1.88, 2.32);

  TH1D * h0a = new TH1D("MpPip", "", 2200, 1.03, 1.47);
  TH1D * h1a = new TH1D("MpKp_TwoPi", "", 2200, 1.38, 1.82);

  TH1D * h0b = new TH1D("MpPim", "", 2200, 1.03, 1.47);
  TH1D * h1b = new TH1D("MpKm_TwoPi", "", 2200, 1.38, 1.82);

  TH1D * h0c = new TH1D("MPiPi", "", 2200, 0.26, 1.70);
  TH1D * h1c = new TH1D("MKK_TwoPi", "", 2200, 0.98, 1.42);

  h0->SetDirectory(fs);
  h1->SetDirectory(fs);

  h0a->SetDirectory(fs);
  h1a->SetDirectory(fs);

  h0b->SetDirectory(fs);
  h1b->SetDirectory(fs);

  h0c->SetDirectory(fs);
  h1c->SetDirectory(fs);

  TH2D * d0a = new TH2D("Momentum_p_Pip", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * d0b = new TH2D("Momentum_p_Pim", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * d0c = new TH2D("Momentum_Pip_Pim", "", 100, 0.0, 2.0, 100, 0.0, 2.0);

  d0a->SetDirectory(fs);
  d0b->SetDirectory(fs);
  d0c->SetDirectory(fs);

  TLorentzVector PP, Pa, Pb, Pc;
  for (Long64_t i = 0; i < Nsim; i++){
    if (i%1000000 == 0) cout << i << endl;
 
    weight = GENERATE::Event_eNPiPi_PiPi(ki, kf);
    if (weight > 0){
      weight *= DETECTOR::Smear(&kf[1], "p") * DETECTOR::Smear(&kf[2], "pi+") * DETECTOR::Smear(&kf[3], "pi-");
      PP = kf[1] + kf[2] + kf[3];
      Pa = kf[1] + kf[2];
      Pb = kf[1] + kf[3];
      Pc = kf[2] + kf[3];
      h0->Fill(PP.M(), weight);
      h0a->Fill(Pa.M(), weight);
      h0b->Fill(Pb.M(), weight);
      h0c->Fill(Pc.M(), weight);
      d0a->Fill(kf[2].P(), kf[1].P(), weight);
      d0b->Fill(kf[3].P(), kf[1].P(), weight);
      d0c->Fill(kf[2].P(), kf[3].P(), weight);
      MisPID(&kf[2]);
      MisPID(&kf[3]);
      PP = kf[1] + kf[2] + kf[3];
      Pa = kf[1] + kf[2];
      Pb = kf[1] + kf[3];
      Pc = kf[2] + kf[3];
      h1->Fill(PP.M(), weight);
      h1a->Fill(Pa.M(), weight);
      h1b->Fill(Pb.M(), weight);
      h1c->Fill(Pc.M(), weight);
    }
  }

  h0->Scale(1.0/Nsim);
  h1->Scale(1.0/Nsim);
 
  h0a->Scale(1.0/Nsim);
  h1a->Scale(1.0/Nsim);
 
  h0b->Scale(1.0/Nsim);
  h1b->Scale(1.0/Nsim);
 
  h0c->Scale(1.0/Nsim);
  h1c->Scale(1.0/Nsim);
 
  d0a->Scale(1.0/Nsim);
  d0b->Scale(1.0/Nsim);
  d0c->Scale(1.0/Nsim);
 
  fs->Write();

  return 0;
}

int MisPID(TLorentzVector * k){
  k->SetE(sqrt(PARTICLE::K.M() * PARTICLE::K.M() + k->P() * k->P()));
  return 0;
}
