#include "Lcore.h"

int main(const int argc, const char * argv[]){

  Long64_t Nsim;
  if (argc < 2) Nsim = 100000000;
  else Nsim = atoi(argv[1]);

  Initialize();

  TLorentzVector ki[2], kf[5];
  double weight = 0;
  //ki[0].SetXYZT(0, 0, sqrt(4.4 * 4.4 - PARTICLE::e.M() * PARTICLE::e.M()), 4.4);
  ki[0].SetXYZT(0, 0, 4.4, 4.4);
  

  TFile * fs = new TFile("result/bumpraw.root", "RECREATE");

  TH1D * h0 = new TH1D("Mass_KK", "", 1000, 0.9, 1.9);
  TH1D * h1 = new TH1D("Mass_gN", "", 1000, 1.9, 2.9);
  TH1D * h2 = new TH1D("CosThetaPhi", "", 1000, -1.0, 1.0);
  TH1D * h3 = new TH1D("E_g", "", 1000, 0.0, 4.4);
  
  TH1D * H0 = new TH1D("Mass_KK_Gold", "", 1000, 0.9, 1.9);
  TH1D * H1 = new TH1D("Mass_gN_Gold", "", 1000, 1.9, 2.9);
  TH1D * H2 = new TH1D("CosThetaPhi_Gold", "", 1000, -1.0, 1.0);
  TH1D * H3 = new TH1D("E_g_Gold", "", 1000, 0.0, 4.4);

  TH2D * s0 = new TH2D("ECosTheta", "", 200, 0.0, 4.4, 200, -1.0, 1.0);
  
  TH2D * S0 = new TH2D("ECosTheta_Gold", "", 200, 0.0, 4.4, 200, -1.0, 1.0);

  TLorentzVector P2(0, 0, 0, Mp);//rest nucleon
  TLorentzVector P1, P3, P4;
  TLorentzVector Ptotal;
  //double p1cm, p3cm, t, t0;
  //double Ecm;
  for (Long64_t i = 0; i < Nsim; i++){
    if (i%1000000 == 0) cout << i << endl;
 
    weight = GENERATE::Event_eNKK_Phi_Nucleon(ki, kf);
    if (weight > 0){
      P1 = ki[0] - kf[0];//virtual photon
      P3 = kf[2] + kf[3];//phi
      P4 = kf[1];//N'
      P2.SetXYZT(0, 0, 0, Mp);
      Ptotal = kf[1] + kf[2] + kf[3];

      P1.Boost(-Ptotal.BoostVector());
      P2.Boost(-Ptotal.BoostVector());
      P3.Boost(-Ptotal.BoostVector());
      P4.Boost(-Ptotal.BoostVector());


      h0->Fill(P3.M(), weight);
      h1->Fill(Ptotal.M(), weight);
      h2->Fill(cos(P3.Angle(P1.Vect())), weight);
      h3->Fill(P1.E(), weight);
      
      s0->Fill(P1.E(), cos(P3.Angle(P1.Vect())), weight);
    }

    weight = GENERATE::Event_eNKK_Phi(ki, kf);
    if (weight > 0){
      P1 = ki[0] - kf[0];//virtual photon
      P3 = kf[2] + kf[3];//phi
      P4 = kf[1];//N'
      P2 = P3 + P4 - P1;
      Ptotal = kf[1] + kf[2] + kf[3];

      P1.Boost(-Ptotal.BoostVector());
      P2.Boost(-Ptotal.BoostVector());
      P3.Boost(-Ptotal.BoostVector());
      P4.Boost(-Ptotal.BoostVector());

      H0->Fill(P3.M(), weight);
      H1->Fill(Ptotal.M(), weight);
      H2->Fill(cos(P3.Angle(P1.Vect())), weight);//|t-t0|
      H3->Fill(P1.E(), weight);
      
      S0->Fill(P1.E(), cos(P3.Angle(P1.Vect())), weight);
    }
  }

  h0->Scale(1.0/Nsim);
  h1->Scale(1.0/Nsim);
  h2->Scale(1.0/Nsim);
  h3->Scale(1.0/Nsim);

  s0->Scale(1.0/Nsim);
 
  H0->Scale(1.0/Nsim);
  H1->Scale(1.0/Nsim);
  H2->Scale(1.0/Nsim);
  H3->Scale(1.0/Nsim);

  S0->Scale(1.0/Nsim);
 
  fs->Write();

  return 0;
}
