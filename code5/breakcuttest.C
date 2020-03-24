#include "Lcore.h"

int Cut(const TLorentzVector p, const TLorentzVector Kp, const TLorentzVector Km, const TLorentzVector kfe);

int main(const int argc, const char * argv[]){

  if (argc < 3) {
    cout << "./breakcuttest <pemin> <savepath> <Nsim>" << endl;
    return 0;
  }

  const double pemin = atof(argv[1]);
  const double pemax = 4.0;

  TString path = argv[2];
  TString filename = path + "breakcuttest.root";
  
  Long64_t Nsim;
  if (argc < 4) Nsim = 100000000;
  else Nsim = atoi(argv[3]);

  Initialize();

  TLorentzVector ki[2], kf[5];
  TLorentzVector kgold;
  double weight = 0;
  ki[0].SetXYZT(0, 0, sqrt(4.4 * 4.4 - PARTICLE::e.M() * PARTICLE::e.M()), 4.4);
  kgold.SetXYZT(0, 0, 0, 183.43);

  TFile * fs = new TFile(filename.Data(), "RECREATE");

  TH1D * h0 = new TH1D("MpKK_BoundStateAll", "", 2200, 1.88, 2.32);
  TH1D * h1 = new TH1D("MpKK_BoundStateKK", "", 2200, 1.88, 2.32);
  TH1D * h2 = new TH1D("MpKK_phi", "", 2200, 1.88, 2.32);
  TH1D * h3 = new TH1D("MpKK_Lambda1520", "", 2200, 1.88, 2.32);
  TH1D * h4 = new TH1D("MpKK_directKK", "", 2200, 1.88, 2.32);
  
  TH1D * h0z = new TH1D("MpKK_z_BoundStateAll", "", 2200, 1.88, 2.32);  //z: center of mass frame
  TH1D * h1z = new TH1D("MpKK_z_BoundStateKK", "", 2200, 1.88, 2.32);
  TH1D * h2z = new TH1D("MpKK_z_phi", "", 2200, 1.88, 2.32);
  TH1D * h3z = new TH1D("MpKK_z_Lambda1520", "", 2200, 1.88, 2.32);
  TH1D * h4z = new TH1D("MpKK_z_directKK", "", 2200, 1.88, 2.32);
  
  
  TH1D * h0a = new TH1D("MpKp_BoundStateAll", "", 2200, 1.38, 1.82);
  TH1D * h1a = new TH1D("MpKp_BoundStateKK", "", 2200, 1.38, 1.82);
  TH1D * h2a = new TH1D("MpKp_phi", "", 2200, 1.38, 1.82);
  TH1D * h3a = new TH1D("MpKp_Lambda1520", "", 2200, 1.38, 1.82);
  TH1D * h4a = new TH1D("MpKp_directKK", "", 2200, 1.38, 1.82);
  
  TH1D * h0az= new TH1D("MpKp_z_BoundStateAll", "", 2200, 1.38, 1.82);
  TH1D * h1az = new TH1D("MpKp_z_BoundStateKK", "", 2200, 1.38, 1.82);
  TH1D * h2az = new TH1D("MpKp_z_phi", "", 2200, 1.38, 1.82);
  TH1D * h3az = new TH1D("MpKp_z_Lambda1520", "", 2200, 1.38, 1.82);
  TH1D * h4az = new TH1D("MpKp_z_directKK", "", 2200, 1.38, 1.82);


  TH1D * h0b = new TH1D("MpKm_BoundStateAll", "", 2200, 1.38, 1.82);
  TH1D * h1b = new TH1D("MpKm_BoundStateKK", "", 2200, 1.38, 1.82);
  TH1D * h2b = new TH1D("MpKm_phi", "", 2200, 1.38, 1.82);
  TH1D * h3b = new TH1D("MpKm_Lambda1520", "", 2200, 1.38, 1.82);
  TH1D * h4b = new TH1D("MpKm_directKK", "", 2200, 1.38, 1.82);

  TH1D * h0bz = new TH1D("MpKm_z_BoundStateAll", "", 2200, 1.38, 1.82);
  TH1D * h1bz = new TH1D("MpKm_z_BoundStateKK", "", 2200, 1.38, 1.82);
  TH1D * h2bz = new TH1D("MpKm_z_phi", "", 2200, 1.38, 1.82);
  TH1D * h3bz = new TH1D("MpKm_z_Lambda1520", "", 2200, 1.38, 1.82);
  TH1D * h4bz = new TH1D("MpKm_z_directKK", "", 2200, 1.38, 1.82);

  TH1D * h0c = new TH1D("MKK_BoundStateAll", "", 2200, 0.98, 1.42);
  TH1D * h1c = new TH1D("MKK_BoundStateKK", "", 2200, 0.98, 1.42);
  TH1D * h2c = new TH1D("MKK_phi", "", 2200, 0.98, 1.42);
  TH1D * h3c = new TH1D("MKK_Lambda1520", "", 2200, 0.98, 1.42);
  TH1D * h4c = new TH1D("MKK_directKK", "", 2200, 0.98, 1.42);
  
  TH1D * h0cz = new TH1D("MKK_z_BoundStateAll", "", 2200, 0.98, 1.42);
  TH1D * h1cz = new TH1D("MKK_z_BoundStateKK", "", 2200, 0.98, 1.42);
  TH1D * h2cz = new TH1D("MKK_z_phi", "", 2200, 0.98, 1.42);
  TH1D * h3cz = new TH1D("MKK_z_Lambda1520", "", 2200, 0.98, 1.42);
  TH1D * h4cz = new TH1D("MKK_z_directKK", "", 2200, 0.98, 1.42);


  TH1D * h0d = new TH1D("Missing_BoundStateAll", "", 2200, 182.41, 184.61);
  TH1D * h1d = new TH1D("Missing_BoundStateKK", "", 2200, 182.41, 184.61);
  TH1D * h2d = new TH1D("Missing_phi", "", 2200, 182.41, 184.61);
  TH1D * h3d = new TH1D("Missing_Lambda1520", "", 2200, 182.41, 184.61);
  TH1D * h4d = new TH1D("Missing_directKK", "", 2200, 182.41, 184.61);
  
  
  
  
  
  TH1D * h0e = new TH1D("Theta_pKp_BoundStateAll", "", 1000, 0.0,  M_PI);
  TH1D * h1e = new TH1D("Theta_pKp_BoundStateKK", "", 1000, 0.0,  M_PI);
  TH1D * h2e = new TH1D("Theta_pKp_phi", "", 1000, 0.0,  M_PI);
  TH1D * h3e = new TH1D("Theta_pKp_Lambda1520", "", 1000, 0.0,  M_PI);
  TH1D * h4e = new TH1D("Theta_pKp_directKK", "", 1000, 0.0,  M_PI);
 
  TH1D * h0ez = new TH1D("Theta_z_pKp_BoundStateAll", "", 1000, 0.0,  M_PI);
  
  
  TH1D * h0f = new TH1D("Theta_pKm_BoundStateAll", "", 1000, 0.0,  M_PI);
  TH1D * h1f = new TH1D("Theta_pKm_BoundStateKK", "", 1000, 0.0,  M_PI);
  TH1D * h2f = new TH1D("Theta_pKm_phi", "", 1000, 0.0,  M_PI);
  TH1D * h3f = new TH1D("Theta_pKm_Lambda1520", "", 1000, 0.0,  M_PI);
  TH1D * h4f = new TH1D("Theta_pKm_directKK", "", 1000, 0.0,  M_PI);
  
  TH1D * h0fz = new TH1D("Theta_z_pKm_BoundStateAll", "", 1000, 0.0,  M_PI);
  
  
  TH1D * h0g = new TH1D("Theta_KpKm_BoundStateAll", "", 1000, 0.0,  M_PI);
  TH1D * h1g = new TH1D("Theta_KpKm_BoundStateKK", "", 1000, 0.0,  M_PI);
  TH1D * h2g = new TH1D("Theta_KpKm_phi", "", 1000, 0.0,  M_PI);
  TH1D * h3g = new TH1D("Theta_KpKm_Lambda1520", "", 1000, 0.0,  M_PI);
  TH1D * h4g = new TH1D("Theta_KpKm_directKK", "", 1000, 0.0,  M_PI);
  
  TH1D * h0gz = new TH1D("Theta_z_KpKm_BoundStateAll", "", 1000, 0.0,  M_PI);
  
  TH1D * h1vp = new TH1D("Theta_vp_BoundStateAll", "", 1000, 0.0, M_PI); //v:virtual photon  p:proton kp: Kaon +   km: Kaon -
  TH1D * h1vkp = new TH1D("Theta_vKp_BoundStateAll", "", 1000, 0.0, M_PI);
  TH1D * h1vkm = new TH1D("Theta_vKm_BoundStateAll", "", 1000, 0.0, M_PI);

  TH1D * h1vpz = new TH1D("Theta_z_vp_BoundStateAll", "", 1000, 0.0, M_PI); //v:virtual photon  p:proton kp: Kaon +   km: Kaon -
  TH1D * h1vkpz = new TH1D("Theta_z_vKp_BoundStateAll", "", 1000, 0.0, M_PI);
  TH1D * h1vkmz = new TH1D("Theta_z_vKm_BoundStateAll", "", 1000, 0.0, M_PI);
  
  h0->SetDirectory(fs);
  h1->SetDirectory(fs);
  h2->SetDirectory(fs);
  h3->SetDirectory(fs);
  h4->SetDirectory(fs);
  h0z->SetDirectory(fs);
  h1z->SetDirectory(fs);
  h2z->SetDirectory(fs);
  h3z->SetDirectory(fs);
  h4z->SetDirectory(fs);
 
  h0a->SetDirectory(fs);
  h1a->SetDirectory(fs);
  h2a->SetDirectory(fs);
  h3a->SetDirectory(fs);
  h4a->SetDirectory(fs);
  h0az->SetDirectory(fs);
  h1az->SetDirectory(fs);
  h2az->SetDirectory(fs);
  h3az->SetDirectory(fs);
  h4az->SetDirectory(fs);
 
  h0b->SetDirectory(fs);
  h1b->SetDirectory(fs);
  h2b->SetDirectory(fs);
  h3b->SetDirectory(fs);
  h4b->SetDirectory(fs);
  h0bz->SetDirectory(fs);
  h1bz->SetDirectory(fs);
  h2bz->SetDirectory(fs);
  h3bz->SetDirectory(fs);
  h4bz->SetDirectory(fs);
 
  h0c->SetDirectory(fs);
  h1c->SetDirectory(fs);
  h2c->SetDirectory(fs);
  h3c->SetDirectory(fs);
  h4c->SetDirectory(fs);
  h0cz->SetDirectory(fs);
  h1cz->SetDirectory(fs);
  h2cz->SetDirectory(fs);
  h3cz->SetDirectory(fs);
  h4cz->SetDirectory(fs);

  h0d->SetDirectory(fs);
  h1d->SetDirectory(fs);
  h2d->SetDirectory(fs);
  h3d->SetDirectory(fs);
  h4d->SetDirectory(fs);
  
  h0e->SetDirectory(fs);
  h1e->SetDirectory(fs);
  h2e->SetDirectory(fs);
  h3e->SetDirectory(fs);
  h4e->SetDirectory(fs);
  h0ez->SetDirectory(fs);
  h0f->SetDirectory(fs);
  h1f->SetDirectory(fs);
  h2f->SetDirectory(fs);
  h3f->SetDirectory(fs);
  h4f->SetDirectory(fs);
  h0fz->SetDirectory(fs);
  h0g->SetDirectory(fs);
  h1g->SetDirectory(fs);
  h2g->SetDirectory(fs);
  h3g->SetDirectory(fs);
  h4g->SetDirectory(fs);
  h0gz->SetDirectory(fs);
 
  h1vp->SetDirectory(fs);
  h1vkp->SetDirectory(fs);
  h1vkm->SetDirectory(fs);
  
  h1vpz->SetDirectory(fs);
  h1vkpz->SetDirectory(fs);
  h1vkmz->SetDirectory(fs);
  
  TH2D * d0a = new TH2D("Momentum_p_Kp_BoundStateAll", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * d0b = new TH2D("Momentum_p_Km_BoundStateAll", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * d0c = new TH2D("Momentum_Kp_Km_BoundStateAll", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * d1a = new TH2D("Momentum_p_Kp_BoundStateKK", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * d1b = new TH2D("Momentum_p_Km_BoundStateKK", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * d1c = new TH2D("Momentum_Kp_Km_BoundStateKK", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * d2a = new TH2D("Momentum_p_Kp_phi", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * d2b = new TH2D("Momentum_p_Km_phi", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * d2c = new TH2D("Momentum_Kp_Km_phi", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * d3a = new TH2D("Momentum_p_Kp_Lambda1520", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * d3b = new TH2D("Momentum_p_Km_Lambda1520", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * d3c = new TH2D("Momentum_Kp_Km_Lambda1520", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * d4a = new TH2D("Momentum_p_Kp_KK", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * d4b = new TH2D("Momentum_p_Km_KK", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * d4c = new TH2D("Momentum_Kp_Km_KK", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  
  TH2D * d0az = new TH2D("Momentum_z_p_Kp_BoundStateAll", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * d0bz = new TH2D("Momentum_z_p_Km_BoundStateAll", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * d0cz = new TH2D("Momentum_z_Kp_Km_BoundStateAll", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * d1az = new TH2D("Momentum_z_p_Kp_BoundStateKK", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * d1bz = new TH2D("Momentum_z_p_Km_BoundStateKK", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * d1cz = new TH2D("Momentum_z_Kp_Km_BoundStateKK", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * d2az = new TH2D("Momentum_z_p_Kp_phi", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * d2bz = new TH2D("Momentum_z_p_Km_phi", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * d2cz = new TH2D("Momentum_z_Kp_Km_phi", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * d3az = new TH2D("Momentum_z_p_Kp_Lambda1520", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * d3bz = new TH2D("Momentum_z_p_Km_Lambda1520", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * d3cz = new TH2D("Momentum_z_Kp_Km_Lambda1520", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * d4az = new TH2D("Momentum_z_p_Kp_KK", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * d4bz = new TH2D("Momentum_z_p_Km_KK", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * d4cz = new TH2D("Momentum_z_Kp_Km_KK", "", 100, 0.0, 2.0, 100, 0.0, 2.0);

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
  
  d0az->SetDirectory(fs);
  d0bz->SetDirectory(fs);
  d0cz->SetDirectory(fs);
  d1az->SetDirectory(fs);
  d1bz->SetDirectory(fs);
  d1cz->SetDirectory(fs);
  d2az->SetDirectory(fs);
  d2bz->SetDirectory(fs);
  d2cz->SetDirectory(fs);
  d3az->SetDirectory(fs);
  d3bz->SetDirectory(fs);
  d3cz->SetDirectory(fs);
  d4az->SetDirectory(fs);
  d4bz->SetDirectory(fs);
  d4cz->SetDirectory(fs);
  
  TH2D * r0a = new TH2D("PTheta_p_BoundStateAll", "", 90, 0.0, 90.0, 100, 0.0, 2.0);
  TH2D * r0b = new TH2D("PTheta_Kp_BoundStateAll", "", 90, 0.0, 90.0, 100, 0.0, 2.0);
  TH2D * r0c = new TH2D("PTheta_Km_BoundStateAll", "", 90, 0.0, 90.0, 100, 0.0, 2.0);
  TH2D * r1a = new TH2D("PTheta_p_BoundStateKK", "", 90, 0.0, 90.0, 100, 0.0, 2.0);
  TH2D * r1b = new TH2D("PTheta_Kp_BoundStateKK", "", 90, 0.0, 90.0, 100, 0.0, 2.0);
  TH2D * r1c = new TH2D("PTheta_Km_BoundStateKK", "", 90, 0.0, 90.0, 100, 0.0, 2.0);
  TH2D * r2a = new TH2D("PTheta_p_phi", "", 90, 0.0, 90.0, 100, 0.0, 2.0);
  TH2D * r2b = new TH2D("PTheta_Kp_phi", "", 90, 0.0, 90.0, 100, 0.0, 2.0);
  TH2D * r2c = new TH2D("PTheta_Km_phi", "", 90, 0.0, 90.0, 100, 0.0, 2.0);
  TH2D * r3a = new TH2D("PTheta_p_Lambda1520", "", 90, 0.0, 90.0, 100, 0.0, 2.0);
  TH2D * r3b = new TH2D("PTheta_Kp_Lambda1520", "", 90, 0.0, 90.0, 100, 0.0, 2.0);
  TH2D * r3c = new TH2D("PTheta_Km_Lambda1520", "", 90, 0.0, 90.0, 100, 0.0, 2.0);
  TH2D * r4a = new TH2D("PTheta_p_KK", "", 90, 0.0, 90.0, 100, 0.0, 2.0);
  TH2D * r4b = new TH2D("PTheta_Kp_KK", "", 90, 0.0, 90.0, 100, 0.0, 2.0);
  TH2D * r4c = new TH2D("PTheta_Km_KK", "", 90, 0.0, 90.0, 100, 0.0, 2.0);
  
  TH2D * r0az = new TH2D("PTheta_z_p_BoundStateAll", "", 90, 0.0, 90.0, 100, 0.0, 2.0);
  TH2D * r0bz = new TH2D("PTheta_z_Kp_BoundStateAll", "", 90, 0.0, 90.0, 100, 0.0, 2.0);
  TH2D * r0cz = new TH2D("PTheta_z_Km_BoundStateAll", "", 90, 0.0, 90.0, 100, 0.0, 2.0);
  TH2D * r1az = new TH2D("PTheta_z_p_BoundStateKK", "", 90, 0.0, 90.0, 100, 0.0, 2.0);
  TH2D * r1bz = new TH2D("PTheta_z_Kp_BoundStateKK", "", 90, 0.0, 90.0, 100, 0.0, 2.0);
  TH2D * r1cz = new TH2D("PTheta_z_Km_BoundStateKK", "", 90, 0.0, 90.0, 100, 0.0, 2.0);
  TH2D * r2az = new TH2D("PTheta_z_p_phi", "", 90, 0.0, 90.0, 100, 0.0, 2.0);
  TH2D * r2bz = new TH2D("PTheta_z_Kp_phi", "", 90, 0.0, 90.0, 100, 0.0, 2.0);
  TH2D * r2cz = new TH2D("PTheta_z_Km_phi", "", 90, 0.0, 90.0, 100, 0.0, 2.0);
  TH2D * r3az = new TH2D("PTheta_z_p_Lambda1520", "", 90, 0.0, 90.0, 100, 0.0, 2.0);
  TH2D * r3bz = new TH2D("PTheta_z_Kp_Lambda1520", "", 90, 0.0, 90.0, 100, 0.0, 2.0);
  TH2D * r3cz = new TH2D("PTheta_z_Km_Lambda1520", "", 90, 0.0, 90.0, 100, 0.0, 2.0);
  TH2D * r4az = new TH2D("PTheta_z_p_KK", "", 90, 0.0, 90.0, 100, 0.0, 2.0);
  TH2D * r4bz = new TH2D("PTheta_z_Kp_KK", "", 90, 0.0, 90.0, 100, 0.0, 2.0);
  TH2D * r4cz = new TH2D("PTheta_z_Km_KK", "", 90, 0.0, 90.0, 100, 0.0, 2.0);
  
  
  
  
  
  TH2D * s0a = new TH2D("Pphi_p_BoundStateAll", "", 180, 0.0, 180.0, 100, 0.0, 2.0);
  TH2D * s0b = new TH2D("Pphi_Kp_BoundStateAll", "", 180, 0.0, 180.0, 100, 0.0, 2.0);
  TH2D * s0c = new TH2D("Pphi_Km_BoundStateAll", "", 180, 0.0, 180.0, 100, 0.0, 2.0);
  
  
  
  const double deg = 180.0 / M_PI;

  r0a->SetDirectory(fs);
  r0b->SetDirectory(fs);
  r0c->SetDirectory(fs);
  r1a->SetDirectory(fs);
  r1b->SetDirectory(fs);
  r1c->SetDirectory(fs);
  r2a->SetDirectory(fs);
  r2b->SetDirectory(fs);
  r2c->SetDirectory(fs);
  r3a->SetDirectory(fs);
  r3b->SetDirectory(fs);
  r3c->SetDirectory(fs);
  r4a->SetDirectory(fs);
  r4b->SetDirectory(fs);
  r4c->SetDirectory(fs);
  
  r0az->SetDirectory(fs);
  r0bz->SetDirectory(fs);
  r0cz->SetDirectory(fs);
  r1az->SetDirectory(fs);
  r1bz->SetDirectory(fs);
  r1cz->SetDirectory(fs);
  r2az->SetDirectory(fs);
  r2bz->SetDirectory(fs);
  r2cz->SetDirectory(fs);
  r3az->SetDirectory(fs);
  r3bz->SetDirectory(fs);
  r3cz->SetDirectory(fs);
  r4az->SetDirectory(fs);
  r4bz->SetDirectory(fs);
  r4cz->SetDirectory(fs);
  
  s0a->SetDirectory(fs);
  s0b->SetDirectory(fs);
  s0c->SetDirectory(fs);
  
  
  TH2D * p0a = new TH2D("phitheta_p_BoundStateAll", "", 150, 0.0, 150.0, 200, 0.0, 200.0);
  TH2D * p0b = new TH2D("phitheta_Kp_BoundStateAll", "", 150, 0.0, 150.0, 200, 0.0, 200.0);
  TH2D * p0c = new TH2D("phitheta_Km_BoundStateAll", "", 150, 0.0, 150.0, 200, 0.0, 200.0);
  
  p0a->SetDirectory(fs);
  p0b->SetDirectory(fs);
  p0c->SetDirectory(fs);
  
  
  TLorentzVector PP, Pa, Pb, Pc, Pd, pbreak;
  TLorentzVector Pv; // virtual photon
  
  double angle1, angle2, angle3;
  double angle4, angle5, angle6;
  
  for (Long64_t i = 0; i < Nsim; i++){
    if (i%1000000 == 0) cout << i << endl;

    GENERATE::NucleonGold(&pbreak);
    pbreak.SetXYZM(pbreak.Px(), pbreak.Py(), pbreak.Pz(), Mp);
    
    weight = GENERATE::Event_eNKKN_BoundState(ki, kf);
    if (weight > 0 && kf[0].P() > pemin && kf[0].P() < pemax){
      weight *= DETECTOR::Smear(&kf[2], "K+") * DETECTOR::Smear(&kf[3], "K-") * DETECTOR::Smear(&kf[4], "p");
      if ( Cut(kf[4], kf[2], kf[3],kf[0]) ){
	PP = kf[2] + kf[3] + kf[4];
	Pa = kf[2] + kf[4];
	Pb = kf[3] + kf[4];
	Pc = kf[2] + kf[3];
	Pd = ki[0] + kgold - kf[0] - PP;
        angle1 = kf[4].Angle(kf[2].Vect());
        angle2 = kf[4].Angle(kf[3].Vect());
        angle3 = kf[2].Angle(kf[3].Vect());
	Pv = ki[0] - kf[0];
	angle4 = Pv.Angle(kf[4].Vect());
        angle5 = Pv.Angle(kf[2].Vect());
        angle6 = Pv.Angle(kf[3].Vect());
      
	h0->Fill(PP.M(), weight);
	h0a->Fill(Pa.M(), weight);
	h0b->Fill(Pb.M(), weight);
	h0c->Fill(Pc.M(), weight);
	h0d->Fill(Pd.M(), weight);
	d0a->Fill(kf[2].P(), kf[4].P(), weight);
	d0b->Fill(kf[3].P(), kf[4].P(), weight);
	d0c->Fill(kf[2].P(), kf[3].P(), weight);
	r0a->Fill(kf[4].Theta() * deg, kf[4].P(), weight);
	r0b->Fill(kf[2].Theta() * deg, kf[2].P(), weight);
	r0c->Fill(kf[3].Theta() * deg, kf[3].P(), weight);
	
	h0e->Fill(angle1 , weight);
        h0f->Fill(angle2 , weight);
        h0g->Fill(angle3 , weight);
	
	h1vp->Fill(angle4 , weight);
        h1vkp->Fill(angle5 , weight);
        h1vkm->Fill(angle6 , weight);
	

	
	p0a->Fill(kf[4].Theta() * deg, kf[4].Phi() * deg, weight);
	p0b->Fill(kf[2].Theta() * deg, kf[2].Phi() * deg, weight);
	p0c->Fill(kf[3].Theta() * deg, kf[3].Phi() * deg, weight);
	
	s0a->Fill(kf[4].Phi() * deg, kf[4].P(), weight);
	s0b->Fill(kf[2].Phi() * deg, kf[2].P(), weight);
	s0c->Fill(kf[3].Phi() * deg, kf[3].P(), weight);
	
	
	kf[4].Boost(-PP.BoostVector());    //change the lab-frame to center of mass frame
	kf[2].Boost(-PP.BoostVector());
	kf[3].Boost(-PP.BoostVector());
	
	PP = kf[2] + kf[3] + kf[4];
	Pa = kf[2] + kf[4];
	Pb = kf[3] + kf[4];
	Pc = kf[2] + kf[3];
	angle1 = kf[4].Angle(kf[2].Vect());
        angle2 = kf[4].Angle(kf[3].Vect());
        angle3 = kf[2].Angle(kf[3].Vect());
	angle4 = Pv.Angle(kf[4].Vect());
        angle5 = Pv.Angle(kf[2].Vect());
        angle6 = Pv.Angle(kf[3].Vect());
	
	h0z->Fill(PP.M(), weight);
	h0az->Fill(Pa.M(), weight);
	h0bz->Fill(Pb.M(), weight);
	h0cz->Fill(Pc.M(), weight);
	
	d0az->Fill(kf[2].P(), kf[4].P(), weight);
	d0bz->Fill(kf[3].P(), kf[4].P(), weight);
	d0cz->Fill(kf[2].P(), kf[3].P(), weight);
	r0az->Fill(kf[4].Theta() * deg, kf[4].P(), weight);
	r0bz->Fill(kf[2].Theta() * deg, kf[2].P(), weight);
	r0cz->Fill(kf[3].Theta() * deg, kf[3].P(), weight);
	
	h0ez->Fill(angle1 , weight);
        h0fz->Fill(angle2 , weight);
        h0gz->Fill(angle3 , weight);
	h1vpz->Fill(angle4 , weight);
        h1vkpz->Fill(angle5 , weight);
        h1vkmz->Fill(angle6 , weight);
      }
    }
    
    weight = GENERATE::Event_eNKKN_BoundState(ki, kf);
    if (weight > 0 && kf[0].P() > pemin && kf[0].P() < pemax){
      kf[4] = pbreak;
      weight *= DETECTOR::Smear(&kf[2], "K+") * DETECTOR::Smear(&kf[3], "K-") * DETECTOR::Smear(&kf[4], "p");
      if ( Cut(kf[4], kf[2], kf[3], kf[0]) ){
	PP = kf[2] + kf[3] + kf[4];
	Pa = kf[2] + kf[4];
	Pb = kf[3] + kf[4];
	Pc = kf[2] + kf[3];
	Pd = ki[0] + kgold - kf[0] - PP;
      
	angle1 = kf[4].Angle(kf[2].Vect());
        angle2 = kf[4].Angle(kf[3].Vect());
        angle3 = kf[2].Angle(kf[3].Vect());
      
      
	h1->Fill(PP.M(), weight);
	h1a->Fill(Pa.M(), weight);
	h1b->Fill(Pb.M(), weight);
	h1c->Fill(Pc.M(), weight);
	h1d->Fill(Pd.M(), weight);
	d1a->Fill(kf[2].P(), kf[4].P(), weight);
	d1b->Fill(kf[3].P(), kf[4].P(), weight);
	d1c->Fill(kf[2].P(), kf[3].P(), weight);
	r1a->Fill(kf[4].Theta() * deg, kf[4].P(), weight);
	r1b->Fill(kf[2].Theta() * deg, kf[2].P(), weight);
	r1c->Fill(kf[3].Theta() * deg, kf[3].P(), weight);
	
	h1e->Fill(angle1 , weight);
        h1f->Fill(angle2 , weight);
        h1g->Fill(angle3 , weight);
	
	
	
	
	
	kf[4].Boost(-PP.BoostVector());
	kf[2].Boost(-PP.BoostVector());
	kf[3].Boost(-PP.BoostVector());
	
	PP = kf[2] + kf[3] + kf[4];
	Pa = kf[2] + kf[4];
	Pb = kf[3] + kf[4];
	Pc = kf[2] + kf[3];
	
	h1z->Fill(PP.M(), weight);
	h1az->Fill(Pa.M(), weight);
	h1bz->Fill(Pb.M(), weight);
	h1cz->Fill(Pc.M(), weight);
	d1az->Fill(kf[2].P(), kf[4].P(), weight);
	d1bz->Fill(kf[3].P(), kf[4].P(), weight);
	d1cz->Fill(kf[2].P(), kf[3].P(), weight);
	r1az->Fill(kf[4].Theta() * deg, kf[4].P(), weight);
	r1bz->Fill(kf[2].Theta() * deg, kf[2].P(), weight);
	r1cz->Fill(kf[3].Theta() * deg, kf[3].P(), weight);
      }
    }


    weight = GENERATE::Event_eNKK_Phi(ki, kf);
    if (weight > 0 && kf[0].P() > pemin && kf[0].P() < pemax){
      kf[1] = pbreak;
      weight *= DETECTOR::Smear(&kf[1], "p") * DETECTOR::Smear(&kf[2], "K+") * DETECTOR::Smear(&kf[3], "K-");
      if ( Cut(kf[1], kf[2], kf[3],kf[0]) ){
	PP = kf[1] + kf[2] + kf[3];
	Pa = kf[1] + kf[2];
	Pb = kf[1] + kf[3];
	Pc = kf[2] + kf[3];
	Pd = ki[0] + kgold - kf[0] - PP;
      
	angle1 = kf[1].Angle(kf[2].Vect());
        angle2 = kf[1].Angle(kf[3].Vect());
        angle3 = kf[2].Angle(kf[3].Vect());
      
      
	h2->Fill(PP.M(), weight);
	h2a->Fill(Pa.M(), weight);
	h2b->Fill(Pb.M(), weight);
	h2c->Fill(Pc.M(), weight);
	h2d->Fill(Pd.M(), weight);
	d2a->Fill(kf[2].P(), kf[1].P(), weight);
	d2b->Fill(kf[3].P(), kf[1].P(), weight);
	d2c->Fill(kf[2].P(), kf[3].P(), weight);
	r2a->Fill(kf[1].Theta() * deg, kf[1].P(), weight);
	r2b->Fill(kf[2].Theta() * deg, kf[2].P(), weight);
	r2c->Fill(kf[3].Theta() * deg, kf[3].P(), weight);
	
	h2e->Fill(angle1 , weight);
        h2f->Fill(angle2 , weight);
        h2g->Fill(angle3 , weight);
	
	kf[1].Boost(-PP.BoostVector());
	kf[2].Boost(-PP.BoostVector());
	kf[3].Boost(-PP.BoostVector());
	
	PP = kf[1] + kf[2] + kf[3];
	Pa = kf[1] + kf[2];
	Pb = kf[1] + kf[3];
	Pc = kf[2] + kf[3];
	
	h2z->Fill(PP.M(), weight);
	h2az->Fill(Pa.M(), weight);
	h2bz->Fill(Pb.M(), weight);
	h2cz->Fill(Pc.M(), weight);
	d2az->Fill(kf[2].P(), kf[1].P(), weight);
	d2bz->Fill(kf[3].P(), kf[1].P(), weight);
	d2cz->Fill(kf[2].P(), kf[3].P(), weight);
	r2az->Fill(kf[1].Theta() * deg, kf[1].P(), weight);
	r2bz->Fill(kf[2].Theta() * deg, kf[2].P(), weight);
	r2cz->Fill(kf[3].Theta() * deg, kf[3].P(), weight);
      }
    }
      
    weight = GENERATE::Event_eNKK_L1520(ki, kf);
    if (weight > 0 && kf[0].P() > pemin && kf[0].P() < pemax){
      kf[1] = pbreak;
      weight *= DETECTOR::Smear(&kf[1], "p") * DETECTOR::Smear(&kf[2], "K+") * DETECTOR::Smear(&kf[3], "K-");
      if ( Cut(kf[1], kf[2], kf[3],kf[0]) ){
	PP = kf[1] + kf[2] + kf[3];
	Pa = kf[1] + kf[2];
	Pb = kf[1] + kf[3];
	Pc = kf[2] + kf[3];
	Pd = ki[0] + kgold - kf[0] - PP;
      
	angle1 = kf[1].Angle(kf[2].Vect());
        angle2 = kf[1].Angle(kf[3].Vect());
        angle3 = kf[2].Angle(kf[3].Vect());
	
	h3->Fill(PP.M(), weight);
	h3a->Fill(Pa.M(), weight);
	h3b->Fill(Pb.M(), weight);
	h3c->Fill(Pc.M(), weight);
	h3d->Fill(Pd.M(), weight);
	d3a->Fill(kf[2].P(), kf[1].P(), weight);
	d3b->Fill(kf[3].P(), kf[1].P(), weight);
	d3c->Fill(kf[2].P(), kf[3].P(), weight);
	r3a->Fill(kf[1].Theta() * deg, kf[1].P(), weight);
	r3b->Fill(kf[2].Theta() * deg, kf[2].P(), weight);
	r3c->Fill(kf[3].Theta() * deg, kf[3].P(), weight);
	
	h3e->Fill(angle1 , weight);
        h3f->Fill(angle2 , weight);
        h3g->Fill(angle3 , weight);
	
	kf[1].Boost(-PP.BoostVector());
	kf[2].Boost(-PP.BoostVector());
	kf[3].Boost(-PP.BoostVector());
	
	PP = kf[1] + kf[2] + kf[3];
	Pa = kf[1] + kf[2];
	Pb = kf[1] + kf[3];
	Pc = kf[2] + kf[3];
	
	h3z->Fill(PP.M(), weight);
	h3az->Fill(Pa.M(), weight);
	h3bz->Fill(Pb.M(), weight);
	h3cz->Fill(Pc.M(), weight);
	d3az->Fill(kf[2].P(), kf[1].P(), weight);
	d3bz->Fill(kf[3].P(), kf[1].P(), weight);
	d3cz->Fill(kf[2].P(), kf[3].P(), weight);
	r3az->Fill(kf[1].Theta() * deg, kf[1].P(), weight);
	r3bz->Fill(kf[2].Theta() * deg, kf[2].P(), weight);
	r3cz->Fill(kf[3].Theta() * deg, kf[3].P(), weight);
      }
    }
      
    weight = GENERATE::Event_eNKK_KK(ki, kf);
    if (weight > 0 && kf[0].P() > pemin && kf[0].P() < pemax){
      kf[1] = pbreak;
      weight *= DETECTOR::Smear(&kf[1], "p") * DETECTOR::Smear(&kf[2], "K+") * DETECTOR::Smear(&kf[3], "K-");
      if ( Cut(kf[1], kf[2], kf[3], kf[0]) ){
	PP = kf[1] + kf[2] + kf[3];
	Pa = kf[1] + kf[2];
	Pb = kf[1] + kf[3];
	Pc = kf[2] + kf[3];
	Pd = ki[0] + kgold - kf[0] - PP;
      
	angle1 = kf[1].Angle(kf[2].Vect());
        angle2 = kf[1].Angle(kf[3].Vect());
        angle3 = kf[2].Angle(kf[3].Vect());
      
	h4->Fill(PP.M(), weight);
	h4a->Fill(Pa.M(), weight);
	h4b->Fill(Pb.M(), weight);
	h4c->Fill(Pc.M(), weight);
	h4d->Fill(Pd.M(), weight);
	d4a->Fill(kf[2].P(), kf[1].P(), weight);
	d4b->Fill(kf[3].P(), kf[1].P(), weight);
	d4c->Fill(kf[2].P(), kf[3].P(), weight);
	r4a->Fill(kf[1].Theta() * deg, kf[1].P(), weight);
	r4b->Fill(kf[2].Theta() * deg, kf[2].P(), weight);
	r4c->Fill(kf[3].Theta() * deg, kf[3].P(), weight);
	
	h4e->Fill(angle1 , weight);
        h4f->Fill(angle2 , weight);
        h4g->Fill(angle3 , weight);
	
	kf[1].Boost(-PP.BoostVector());
	kf[2].Boost(-PP.BoostVector());
	kf[3].Boost(-PP.BoostVector());
	
	PP = kf[1] + kf[2] + kf[3];
	Pa = kf[1] + kf[2];
	Pb = kf[1] + kf[3];
	Pc = kf[2] + kf[3];
	
	h4z->Fill(PP.M(), weight);
	h4az->Fill(Pa.M(), weight);
	h4bz->Fill(Pb.M(), weight);
	h4cz->Fill(Pc.M(), weight);
	d4az->Fill(kf[2].P(), kf[1].P(), weight);
	d4bz->Fill(kf[3].P(), kf[1].P(), weight);
	d4cz->Fill(kf[2].P(), kf[3].P(), weight);
	r4az->Fill(kf[1].Theta() * deg, kf[1].P(), weight);
	r4bz->Fill(kf[2].Theta() * deg, kf[2].P(), weight);
	r4cz->Fill(kf[3].Theta() * deg, kf[3].P(), weight);
      }
    }

  }

  h0->Scale(1.0/Nsim);
  h1->Scale(1.0/Nsim);
  h2->Scale(1.0/Nsim);
  h3->Scale(1.0/Nsim);
  h4->Scale(1.0/Nsim);
  h0z->Scale(1.0/Nsim);
  h1z->Scale(1.0/Nsim);
  h2z->Scale(1.0/Nsim);
  h3z->Scale(1.0/Nsim);
  h4z->Scale(1.0/Nsim);
 
  h0a->Scale(1.0/Nsim);
  h1a->Scale(1.0/Nsim);
  h2a->Scale(1.0/Nsim);
  h3a->Scale(1.0/Nsim);
  h4a->Scale(1.0/Nsim);
  h0az->Scale(1.0/Nsim);
  h1az->Scale(1.0/Nsim);
  h2az->Scale(1.0/Nsim);
  h3az->Scale(1.0/Nsim);
  h4az->Scale(1.0/Nsim);
 
  h0b->Scale(1.0/Nsim);
  h1b->Scale(1.0/Nsim);
  h2b->Scale(1.0/Nsim);
  h3b->Scale(1.0/Nsim);
  h4b->Scale(1.0/Nsim);
  h0bz->Scale(1.0/Nsim);
  h1bz->Scale(1.0/Nsim);
  h2bz->Scale(1.0/Nsim);
  h3bz->Scale(1.0/Nsim);
  h4bz->Scale(1.0/Nsim);
 
  h0c->Scale(1.0/Nsim);
  h1c->Scale(1.0/Nsim);
  h2c->Scale(1.0/Nsim);
  h3c->Scale(1.0/Nsim);
  h4c->Scale(1.0/Nsim);
  h0cz->Scale(1.0/Nsim);
  h1cz->Scale(1.0/Nsim);
  h2cz->Scale(1.0/Nsim);
  h3cz->Scale(1.0/Nsim);
  h4cz->Scale(1.0/Nsim);

  h0d->Scale(1.0/Nsim);
  h1d->Scale(1.0/Nsim);
  h2d->Scale(1.0/Nsim);
  h3d->Scale(1.0/Nsim);
  h4d->Scale(1.0/Nsim);
 
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
  
  d0az->Scale(1.0/Nsim);
  d0bz->Scale(1.0/Nsim);
  d0cz->Scale(1.0/Nsim);
  d1az->Scale(1.0/Nsim);
  d1bz->Scale(1.0/Nsim);
  d1cz->Scale(1.0/Nsim);
  d2az->Scale(1.0/Nsim);
  d2bz->Scale(1.0/Nsim);   
  d2cz->Scale(1.0/Nsim);
  d3az->Scale(1.0/Nsim);
  d3bz->Scale(1.0/Nsim);
  d3cz->Scale(1.0/Nsim);
  d4az->Scale(1.0/Nsim);
  d4bz->Scale(1.0/Nsim);
  d4cz->Scale(1.0/Nsim);
 
  r0a->Scale(1.0/Nsim);
  r0b->Scale(1.0/Nsim);
  r0c->Scale(1.0/Nsim);
  r1a->Scale(1.0/Nsim);
  r1b->Scale(1.0/Nsim);
  r1c->Scale(1.0/Nsim);
  r2a->Scale(1.0/Nsim);
  r2b->Scale(1.0/Nsim);
  r2c->Scale(1.0/Nsim);
  r3a->Scale(1.0/Nsim);
  r3b->Scale(1.0/Nsim);
  r3c->Scale(1.0/Nsim);
  r4a->Scale(1.0/Nsim);
  r4b->Scale(1.0/Nsim);
  r4c->Scale(1.0/Nsim);
  r0az->Scale(1.0/Nsim);
  r0bz->Scale(1.0/Nsim);
  r0cz->Scale(1.0/Nsim);
  r1az->Scale(1.0/Nsim);
  r1bz->Scale(1.0/Nsim);
  r1cz->Scale(1.0/Nsim);
  r2az->Scale(1.0/Nsim);
  r2bz->Scale(1.0/Nsim);
  r2cz->Scale(1.0/Nsim);
  r3az->Scale(1.0/Nsim);
  r3bz->Scale(1.0/Nsim);
  r3cz->Scale(1.0/Nsim);
  r4az->Scale(1.0/Nsim);
  r4bz->Scale(1.0/Nsim);
  r4cz->Scale(1.0/Nsim);
  
  h0e->Scale(1.0/Nsim);
  h1e->Scale(1.0/Nsim);
  h2e->Scale(1.0/Nsim);
  h3e->Scale(1.0/Nsim);
  h4e->Scale(1.0/Nsim);
  
  h0f->Scale(1.0/Nsim);
  h1f->Scale(1.0/Nsim);
  h2f->Scale(1.0/Nsim);
  h3f->Scale(1.0/Nsim);
  h4f->Scale(1.0/Nsim);
  
  h0g->Scale(1.0/Nsim);
  h1g->Scale(1.0/Nsim);
  h2g->Scale(1.0/Nsim);
  h3g->Scale(1.0/Nsim);
  h4g->Scale(1.0/Nsim);
  
  h0ez->Scale(1.0/Nsim);
  h0fz->Scale(1.0/Nsim);
  h0gz->Scale(1.0/Nsim);
  
  h1vp->Scale(1.0/Nsim);
  h1vkp->Scale(1.0/Nsim);
  h1vkm->Scale(1.0/Nsim);
  h1vpz->Scale(1.0/Nsim);
  h1vkpz->Scale(1.0/Nsim);
  h1vkmz->Scale(1.0/Nsim);
  
  p0a->Scale(1.0/Nsim);
  p0b->Scale(1.0/Nsim);
  p0c->Scale(1.0/Nsim);
  
  s0a->Scale(1.0/Nsim);
  s0b->Scale(1.0/Nsim);
  s0c->Scale(1.0/Nsim);
  
  
  fs->Write();

  return 0;
}

int Cut(const TLorentzVector p, const TLorentzVector Kp, const TLorentzVector Km, const TLorentzVector kfe){
  TLorentzVector ki;
  TLorentzVector kgold;
  ki.SetXYZT(0, 0, sqrt(4.4 * 4.4 - PARTICLE::e.M() * PARTICLE::e.M()), 4.4);
  kgold.SetXYZT(0, 0, 0, 183.43);
  
  
  TLorentzVector kk = p + Kp;
  if (kk.M() > 1.48) return 0;
  kk = p + Km;
  if (kk.M() > 1.48) return 0;
  kk = Kp + Km;
  if (kk.M() > 1.04) return 0;
  
  
  TLorentzVector miss= ki + kgold - kfe - p - Kp - Km;
  if (miss.M() < 182.49) return 0;
  
  
  if (p.P() > 0.8) return 0;
  if (Kp.P() > 0.5) return 0;
  if (Km.P() >0.5) return 0;
  if (p.Theta() > M_PI / 3.0) return 0; 
  
  TLorentzVector PP = p + Kp + Km;    // change the frame to cm frame
  TLorentzVector p1 = p; 
  TLorentzVector Kp1 = Kp;
  TLorentzVector Km1 = Km;
  
  p1.Boost(-PP.BoostVector());
  Kp1.Boost(-PP.BoostVector());
  Km1.Boost(-PP.BoostVector());
  
  if(Km1.P()+p1.P()-0.26 > 0) return 0;
  if(Kp1.P()+p1.P()-0.26 > 0) return 0;
  if(Km1.P()+Kp1.P()-0.26 > 0) return 0;
  return 1;
}
