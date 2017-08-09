#include "Lcore.h"

int main(const int argc, const char * argv[]){

  Long64_t Nsim;
  if (argc < 2) Nsim = 100000000;
  else Nsim = atoi(argv[1]);

  Initialize();

  TLorentzVector ki[2], kf[5];
  double weight = 0;
  //ki[0].SetXYZT(0, 0, sqrt(4.4 * 4.4 - PARTICLE::e.M() * PARTICLE::e.M()), 4.4);
  ki[0].SetXYZM(0, 0, 4.4, PARTICLE::e.M());
  

  TFile * fs = new TFile("result/bumpdetected.root", "RECREATE");

  TH1D * h0 = new TH1D("Mass_KK", "", 1000, 0.9, 1.9);
  TH1D * h1 = new TH1D("Mass_gN", "", 1000, 1.9, 2.9);
  TH1D * h2 = new TH1D("CosThetaPhi", "", 1000, -1.0, 1.0);
  TH1D * h3 = new TH1D("E_g", "", 1000, 0.0, 4.4);
  
  TH1D * H0 = new TH1D("Mass_KK_Gold", "", 1000, 0.9, 1.9);
  TH1D * H1 = new TH1D("Mass_gN_Gold", "", 1000, 1.9, 2.9);
  TH1D * H2 = new TH1D("CosThetaPhi_Gold", "", 1000, -1.0, 1.0);
  TH1D * H3 = new TH1D("E_g_Gold", "", 1000, 0.0, 4.4);

  TH2D * s0 = new TH2D("ECosTheta", "", 100, 1.9, 2.9, 200, -1.0, 1.0);
  TH2D * s1 = new TH2D("ECosTheta_clas12", "", 100, 1.9, 2.9, 200, -1.0, 1.0);
  TH2D * s2 = new TH2D("ECosTheta_bonus12", "", 100, 1.9, 2.9, 200, -1.0, 1.0);
  
  TH2D * S0 = new TH2D("ECosTheta_Gold", "", 100, 1.9, 2.9, 200, -1.0, 1.0);
  TH2D * S1 = new TH2D("ECosTheta_Gold_clas12", "", 100, 1.9, 2.9, 200, -1.0, 1.0);
  TH2D * S2 = new TH2D("ECosTheta_Gold_bonus12", "", 100, 1.9, 2.9, 200, -1.0, 1.0);

  TLorentzVector P2(0, 0, 0, Mp);//rest nucleon
  TLorentzVector P1, P3, P4;
  TLorentzVector Ptotal;
  double acc = 0.0;
  for (Long64_t i = 0; i < Nsim; i++){
    if (i%1000000 == 0) cout << i << endl;
 
    //weight = GENERATE::Event_eNKK_Phi_Nucleon(ki, kf);
    weight = GENERATE::Event_eNKK_Phi_Nucleon(ki, kf) / 0.489;
    
    if (weight > 0){
      //acc = DETECTOR::Acceptance(kf[2], "K+") * DETECTOR::Acceptance(kf[3], "K-");
      acc = DETECTOR::AcceptancePhi(kf[2] + kf[3], "all");
      if (acc > 0){
	P1 = ki[0] - kf[0];//virtual photon
	P2.SetXYZT(0, 0, 0, Mp);//rest proton
	P3 = kf[2] + kf[3];//phi
	P4 = kf[1];//N'

	Ptotal = P3 + P4;
	
	P1.Boost(-Ptotal.BoostVector());
	P2.Boost(-Ptotal.BoostVector());
	P3.Boost(-Ptotal.BoostVector());
	P4.Boost(-Ptotal.BoostVector());
	
	
	h0->Fill(P3.M(), weight * acc);
	h1->Fill(Ptotal.M(), weight * acc);
	h2->Fill(cos(P3.Angle(P1.Vect())), weight * acc);
	h3->Fill(P1.E(), weight * acc);
	
	s0->Fill(Ptotal.M(), cos(P3.Angle(P1.Vect())), weight * acc);

	//acc = DETECTOR::AcceptanceCLAS12(kf[2], "K+") * DETECTOR::AcceptanceCLAS12(kf[3], "K-");
	acc = DETECTOR::AcceptancePhi(kf[2] + kf[3], "clas12");
	s1->Fill(Ptotal.M(), cos(P3.Angle(P1.Vect())), weight * acc);

	//acc = DETECTOR::AcceptanceBONUS12(kf[2], "K+") * DETECTOR::AcceptanceBONUS12(kf[3], "K-");
	acc = DETECTOR::AcceptancePhi(kf[2] + kf[3], "bonus12");
	s2->Fill(Ptotal.M(), cos(P3.Angle(P1.Vect())), weight * acc);	

      }
    }

    weight = GENERATE::Event_eNKK_Phi(ki, kf);
    if (weight > 0){
      weight = weight / 0.489;
      acc = DETECTOR::AcceptancePhi(kf[2] + kf[3], "all") * DETECTOR::Acceptance(kf[1], "p");
      if (acc > 0){
	P1 = ki[0] - kf[0];//virtual photon
	P3 = kf[2] + kf[3];//phi
	P4 = kf[1];//N'
	P2 = P3 + P4 - P1;
	Ptotal = P3 + P4;
	
	P1.Boost(-Ptotal.BoostVector());
	P2.Boost(-Ptotal.BoostVector());
	P3.Boost(-Ptotal.BoostVector());
	P4.Boost(-Ptotal.BoostVector());
	
	H0->Fill(P3.M(), weight * acc);
	H1->Fill(Ptotal.M(), weight * acc);
	H2->Fill(cos(P3.Angle(P1.Vect())), weight * acc);
	H3->Fill(P1.E(), weight * acc);
	
	S0->Fill(Ptotal.M(), cos(P3.Angle(P1.Vect())), weight * acc);

	acc = DETECTOR::AcceptancePhi(kf[1], "clas12") * DETECTOR::AcceptanceCLAS12(kf[2], "p");
	S1->Fill(Ptotal.M(), cos(P3.Angle(P1.Vect())), weight * acc);

	acc = DETECTOR::AcceptancePhi(kf[1], "bonus12") * DETECTOR::AcceptanceBONUS12(kf[2], "p");
	S2->Fill(Ptotal.M(), cos(P3.Angle(P1.Vect())), weight * acc);
      }
    }
  }

  h0->Scale(1.0/Nsim);
  h1->Scale(1.0/Nsim);
  h2->Scale(1.0/Nsim);
  h3->Scale(1.0/Nsim);

  s0->Scale(1.0/Nsim);
  s1->Scale(1.0/Nsim);
  s2->Scale(1.0/Nsim);
 
  H0->Scale(1.0/Nsim);
  H1->Scale(1.0/Nsim);
  H2->Scale(1.0/Nsim);
  H3->Scale(1.0/Nsim);

  S0->Scale(1.0/Nsim);
  S1->Scale(1.0/Nsim);
  S2->Scale(1.0/Nsim);
  
  fs->Write();

  return 0;
}
