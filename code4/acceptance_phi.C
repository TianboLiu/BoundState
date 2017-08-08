#include "Lcore.h"

int main(const int argc, const char * argv[]){
  
  Long64_t Nsim;
  if (argc < 2) Nsim = 100000000;
  else Nsim = atoi(argv[1]);

  Initialize();

  TLorentzVector P(0, 0, 0, 0);
  TLorentzVector P1, P2;
  double p, theta, phi, Mphi;
  double acc = 0.0;
  double branch = 0.489;

  double mass[2] = {PARTICLE::K.M(), PARTICLE::K.M()};
  
  TFile * fs = new TFile("acceptance/acceptance_phi.root", "RECREATE");
  TH3F * hB = new TH3F("generate_PThetaPhi", "", 400, 0.0, 4.0, 90, 0.0, M_PI, 180, -M_PI, M_PI);
  TH3F * h0 = new TH3F("acceptance_PThetaPhi_clas12", "", 400, 0.0, 4.0, 90, 0.0, M_PI, 180, -M_PI, M_PI);
  TH3F * h1 = new TH3F("acceptance_PThetaPhi_bonus12", "", 400, 0.0, 4.0, 90, 0.0, M_PI, 180, -M_PI, M_PI);
  TH3F * h2 = new TH3F("acceptance_PThetaPhi_all", "", 400, 0.0, 4.0, 90, 0.0, M_PI, 180, -M_PI, M_PI);

  TH2D * gB = new TH2D("generate_PTheta", "", 400, 0.0, 4.0, 90, 0.0, M_PI);
  TH2D * g0 = new TH2D("acceptance_PTheta_clas12", "", 400, 0.0, 4.0, 90, 0.0, M_PI);
  TH2D * g1 = new TH2D("acceptance_PTheta_bonus12", "", 400, 0.0, 4.0, 90, 0.0, M_PI);
  TH2D * g2 = new TH2D("acceptance_PTheta_all", "", 400, 0.0, 4.0, 90, 0.0, M_PI);


  TRandom3 ran(0);
  for (Long64_t i = 0; i < Nsim; i++){
    if (i%10000000 == 0) cout << i << endl;
    p = ran.Uniform(0.0, 4.0);
    theta = ran.Uniform(0.0, M_PI);
    phi = ran.Uniform(-M_PI, M_PI);
    Mphi = PARTICLE::phi.RandomM(mass[0] + mass[1], 1.1);

    P.SetXYZM(p * sin(theta) * cos(phi), p * sin(theta) * sin(phi), p * cos(theta), Mphi);
    
    GENERATE::GenPhase.SetDecay(P, 2, mass);
    GENERATE::GenPhase.Generate();
    P1 = * GENERATE::GenPhase.GetDecay(0);
    P2 = * GENERATE::GenPhase.GetDecay(1);

    hB->Fill(p, theta, phi, 1.0);
    gB->Fill(p, theta, 1.0);

    acc = DETECTOR::AcceptanceCLAS12(P1, "K+") * DETECTOR::AcceptanceCLAS12(P2, "K-");
    h0->Fill(p, theta, phi, branch * acc);
    g0->Fill(p, theta, branch * acc);
    
    acc = DETECTOR::AcceptanceBONUS12(P1, "K+") * DETECTOR::AcceptanceBONUS12(P2, "K-");
    h1->Fill(p, theta, phi, branch * acc);
    g1->Fill(p, theta, branch * acc);
   
    
    acc = DETECTOR::Acceptance(P1, "K+") * DETECTOR::Acceptance(P2, "K-");
    h2->Fill(p, theta, phi, branch * acc);
    g2->Fill(p, theta, branch * acc);
  }

  h0->Divide(hB);
  h1->Divide(hB);
  h2->Divide(hB);

  g0->Divide(gB);
  g1->Divide(gB);
  g2->Divide(gB);  

  fs->Write();

  return 0;
}
