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
  TH1D * h2 = new TH1D("T_gPhi", "", 1000, 0.0, 3.0);
  TH1D * h3 = new TH1D("E_g", "", 1000, 0.0, 4.4);
  
  TH1D * H0 = new TH1D("Mass_KK_Gold", "", 1000, 0.9, 1.9);
  TH1D * H1 = new TH1D("Mass_gN_Gold", "", 1000, 1.9, 2.9);
  TH1D * H2 = new TH1D("T_gPhi_Gold", "", 1000, 0.0, 3.0);
  TH1D * H3 = new TH1D("E_g_Gold", "", 1000, 0.0, 4.4);

  TH2D * s0 = new TH2D("ET", "", 200, 0.0, 4.4, 200, 0.0, 3.0);
  
  TH2D * S0 = new TH2D("ET_Gold", "", 200, 0.0, 4.4, 200, 0.0, 3.0);

  TLorentzVector PP, Pa, Pb, Pc;
  TLorentzVector P0(0, 0, 0, Mp);//rest nucleon
  for (Long64_t i = 0; i < Nsim; i++){
    if (i%1000000 == 0) cout << i << endl;
 
    weight = GENERATE::Event_eNKK_Phi_Nucleon(ki, kf);
    if (weight > 0){
      PP = ki[0] - kf[0];//virtual photon
      Pa = kf[2] + kf[3];//K+K-
      Pb = kf[1] + kf[2] + kf[3];//NKK
      Pc = PP - Pa;//photon - phi

      h0->Fill(Pa.M(), weight);
      h1->Fill(Pb.M(), weight);
      h2->Fill(-Pc.M2(), weight);//-t
      h3->Fill(PP.E(), weight);
      
      s0->Fill(PP.E(), -Pc.M2(), weight);
    }

    weight = GENERATE::Event_eNKK_Phi(ki, kf);
    if (weight > 0){
      PP = ki[0] - kf[0];//virtual photon
      Pa = kf[2] + kf[3];//K+K-
      Pb = kf[1] + kf[2] + kf[3];//NKK
      Pc = PP - Pa;//photon - phi

      H0->Fill(Pa.M(), weight);
      H1->Fill(Pb.M(), weight);
      H2->Fill(-Pc.M2(), weight);//-t
      H3->Fill(PP.E(), weight);

      S0->Fill(PP.E(), -Pc.M2(), weight);
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
