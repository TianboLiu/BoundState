#include "Lcore.h"

int MisPID(TLorentzVector * k);
int Cut(const TLorentzVector p, const TLorentzVector Kp, const TLorentzVector Km, const TLorentzVector kfe);

int main(const int argc, const char * argv[]){

  if (argc < 3){
    cout << "./twopinewcuttest <pemin> <savepath> <filename> <Nsim>" << endl;
    return 0;
  }

  const double pemin = atof(argv[1]);
  const double pemax = 4.0;

  TString path = argv[2];
  TString filename = path + "modifiedtwopinew.root";
  
  double Nsim = atof(argv[4]);

  Initialize();

  double convert = 1e-4 / pow(Phys::hbar, 2) * 79.0;//convert unit to GeV^-2 per gold

  TLorentzVector ki[2], kf[5];
  TLorentzVector kgold;
  ki[0].SetXYZT(0, 0, sqrt(4.4 * 4.4 - PARTICLE::e.M() * PARTICLE::e.M()), 4.4);
  kgold.SetXYZT(0, 0, 0, 183.43);

  TFile * fs = new TFile(filename.Data(), "RECREATE");

  TH1D * h0 = new TH1D("MpPiPi", "", 3000, 1.24, 2.6);    //
  TH1D * h1 = new TH1D("MpKK_PiPi", "", 3000, 1.93, 2.85);
  TH1D * h1z = new TH1D("MpKK_z_PiPi", "", 3000, 1.93, 2.85); // z: cm frame 
  
  
  TH1D * h0a = new TH1D("MpPip", "", 3000, 1.08, 2.0);   //
  TH1D * h1a = new TH1D("MpKp_PiPi", "", 3000, 1.43, 2.1);
  TH1D * h1az = new TH1D("MpKp_z_PiPi", "", 3000, 1.43, 2.1);
 
  
  TH1D * h0b = new TH1D("MpPim", "", 3000, 1.08, 1.80);   //
  TH1D * h1b = new TH1D("MpKm_PiPi", "", 3000, 1.43, 1.87);
  TH1D * h1bz = new TH1D("MpKm_z_PiPi", "", 3000, 1.43, 1.87);
  
 
  TH1D * h0c = new TH1D("MPiPi", "", 3000, 0.28, 1.10);    //
  TH1D * h1c = new TH1D("MKK_PiPi", "", 3000, 0.98, 1.50); //
  TH1D * h1cz = new TH1D("MKK_z_PiPi", "", 3000, 0.98, 1.50); //
  
  
  TH1D * h0d = new TH1D("Missing", "", 2200, 182.45, 182.52);
  TH1D * h1d = new TH1D("Missing_PiPi", "", 2200, 181.80, 182.23);//183.6, 186.0);
  TH1D * h0dp = new TH1D("Missingp", "", 2200, 0.0, 0.6);
  TH1D * h1dp = new TH1D("Missingp_PiPi", "", 2200, 0.0, 0.6);
  
  TH1D * h1e = new TH1D("Theta_pKp_PiPi", "", 1000, 0.0,  M_PI);
  TH1D * h1ez = new TH1D("Theta_z_pKp_PiPi", "", 1000, 0.0,  M_PI);
  
  TH1D * h1f = new TH1D("Theta_pKm_PiPi", "", 1000, 0.0,  M_PI);
  TH1D * h1fz = new TH1D("Theta_z_pKm_PiPi", "", 1000, 0.0,  M_PI);
  
  TH1D * h1g = new TH1D("Theta_KpKm_PiPi", "", 1000, 0.0,  M_PI);
  TH1D * h1gz = new TH1D("Theta_z_KpKm_PiPi", "", 1000, 0.0,  M_PI);
  
  TH1D * h1vp = new TH1D("Theta_vp_PiPi", "", 1000, 0.0, M_PI); //v:virtual photon  p:proton kp: Kaon +   km: Kaon -
  TH1D * h1vkp = new TH1D("Theta_vKp_PiPi", "", 1000, 0.0, M_PI);
  TH1D * h1vkm = new TH1D("Theta_vKm_PiPi", "", 1000, 0.0, M_PI);

  TH1D * h1vpz = new TH1D("Theta_z_vp_PiPi", "", 1000, 0.0, M_PI); //z:center of mass frame
  TH1D * h1vkpz = new TH1D("Theta_z_vKp_PiPi", "", 1000, 0.0, M_PI);//v:virtual photon  p:proton kp: Kaon +   km: Kaon -
  TH1D * h1vkmz = new TH1D("Theta__zvKm_PiPi", "", 1000, 0.0, M_PI);
  
  h0->SetDirectory(fs);
  h1->SetDirectory(fs);
  h1z->SetDirectory(fs);
  
  h0a->SetDirectory(fs);
  h1a->SetDirectory(fs);
  h1az->SetDirectory(fs);

  h0b->SetDirectory(fs);
  h1b->SetDirectory(fs);
  h1bz->SetDirectory(fs);

  h0c->SetDirectory(fs);
  h1c->SetDirectory(fs);
  h1cz->SetDirectory(fs);
  
  h0d->SetDirectory(fs);
  h1d->SetDirectory(fs);
  
  h0dp->SetDirectory(fs);
  h1dp->SetDirectory(fs);
  
  h1e->SetDirectory(fs);
  h1f->SetDirectory(fs);
  h1g->SetDirectory(fs);
  
  h1ez->SetDirectory(fs);
  h1fz->SetDirectory(fs);
  h1gz->SetDirectory(fs);
  
  h1vp->SetDirectory(fs);
  h1vkp->SetDirectory(fs);
  h1vkm->SetDirectory(fs);
  h1vpz->SetDirectory(fs);
  h1vkpz->SetDirectory(fs);
  h1vkmz->SetDirectory(fs);
  
  TH2D * d0a = new TH2D("Momentum_p_Pip", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * d0b = new TH2D("Momentum_p_Pim", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * d0c = new TH2D("Momentum_Pip_Pim", "", 100, 0.0, 2.0, 100, 0.0, 2.0);

  d0a->SetDirectory(fs);
  d0b->SetDirectory(fs);
  d0c->SetDirectory(fs);

  TH2D * d1a = new TH2D("Momentum_p_Pip_PiPi", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * d1b = new TH2D("Momentum_p_Pim_PiPi", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  TH2D * d1c = new TH2D("Momentum_Pip_Pim_PiPi", "", 100, 0.0, 2.0, 100, 0.0, 2.0);
  
  TH2D * d1az = new TH2D("Momentum_z_p_Pip_PiPi", "", 100, 0.0, 1.35, 100, 0.0, 1.35);  //
  TH2D * d1bz = new TH2D("Momentum_z_p_Pim_PiPi", "", 100, 0.0, 1.2, 100, 0.0, 1.35);  //
  TH2D * d1cz = new TH2D("Momentum_z_Pip_Pim_PiPi", "", 100, 0.0, 1.35, 100, 0.0, 1.2);  //

  d1a->SetDirectory(fs);
  d1b->SetDirectory(fs);
  d1c->SetDirectory(fs);
  d1az->SetDirectory(fs);
  d1bz->SetDirectory(fs);
  d1cz->SetDirectory(fs);
  
  TH2D * r0a = new TH2D("PTheta_p", "", 90,0.0, 90.0, 100, 0.0, 2.0);
  TH2D * r0b = new TH2D("PTheta_Pip", "", 90,0.0, 90.0,100, 0.0, 2.0);
  TH2D * r0c = new TH2D("PTheta_Pim", "", 90,0.0, 90.0,100, 0.0, 2.0);
  
  TH2D * r1a = new TH2D("PTheta_p_PiPi", "", 90,0.0, 90.0,100, 0.0, 2.0);
  TH2D * r1b = new TH2D("PTheta_Kp_PiPi", "", 90,0.0, 90.0,100, 0.0, 2.0);
  TH2D * r1c = new TH2D("PTheta_Km_PiPi", "", 90,0.0, 90.0,100, 0.0, 2.0);
  
  TH2D * r1az = new TH2D("PTheta_z_p_PiPi", "", 90,0.0, 90.0,100, 0.0, 1.35);  //
  TH2D * r1bz = new TH2D("PTheta_z_Kp_PiPi", "", 90,0.0, 90.0,100, 0.0, 1.35);  //
  TH2D * r1cz= new TH2D("PTheta_z_Km_PiPi", "", 90,0.0, 90.0,100, 0.0, 1.2);   //
  
  
  
  TH2D * s1a = new TH2D("Pphi_p_PiPi", "", 180,0.0, 180.0,100, 0.0, 2.0);
  TH2D * s1b = new TH2D("Pphi_Kp_PiPi", "", 180,0.0, 180.0,100, 0.0, 2.0);
  TH2D * s1c = new TH2D("Pphi_Km_PiPi", "", 180,0.0, 180.0,100, 0.0, 2.0);
  
  const double deg = 180.0 / M_PI;
  
  r0a->SetDirectory(fs);
  r0b->SetDirectory(fs);
  r0c->SetDirectory(fs);
  r1a->SetDirectory(fs);
  r1b->SetDirectory(fs);
  r1c->SetDirectory(fs);
  r1az->SetDirectory(fs);
  r1bz->SetDirectory(fs);
  r1cz->SetDirectory(fs);
  
  s1a->SetDirectory(fs);
  s1b->SetDirectory(fs);
  s1c->SetDirectory(fs);
 
  
  
  TLorentzVector PP, Pa, Pb, Pc, Pd;
  TLorentzVector Pv; //virtual photon

  ifstream infile(argv[3]);
  double temp, px, py, pz, E, weight;
  double angle1, angle2, angle3; //theta between pkp pkm kpkm
  double angle4, angle5, angle6; //theta between virtual photon and p, kp, km
  
  int countN = 0;

  
  while (infile >> temp >> temp >> temp >> temp >> temp >> temp >> temp >> temp >> temp >> weight){
    infile >> temp >> temp >> temp >> temp >> temp >> temp >> px >> py >> pz >> E >> temp >> temp >> temp >> temp;
    kf[0].SetXYZT(px, py, pz, E);
    infile >> temp >> temp >> temp >> temp >> temp >> temp >> px >> py >> pz >> E >> temp >> temp >> temp >> temp;
    kf[1].SetXYZT(px, py, pz, E);
    infile >> temp >> temp >> temp >> temp >> temp >> temp >> px >> py >> pz >> E >> temp >> temp >> temp >> temp;
    kf[2].SetXYZT(px, py, pz, E);
    infile >> temp >> temp >> temp >> temp >> temp >> temp >> px >> py >> pz >> E >> temp >> temp >> temp >> temp;
    kf[3].SetXYZT(px, py, pz, E);
    if (weight > 0 && kf[0].P() > pemin && kf[0].P() < pemax){
      if (kf[0].Theta() / M_PI * 180.0 > 4.5) continue;
      if (kf[0].Theta() / M_PI * 180.0 < 2.5) continue;
        
      weight *= DETECTOR::Smear(&kf[1], "p") * DETECTOR::Smear(&kf[2], "pi+") * DETECTOR::Smear(&kf[3], "pi-");
      
      
      PP = kf[1] + kf[2] + kf[3];
      Pa = kf[1] + kf[2];
      Pb = kf[1] + kf[3];
      Pc = kf[2] + kf[3];
      Pd = ki[0] + kgold - kf[0] - PP;
      Pv = ki[0] - kf[0];
      
      h0->Fill(PP.M(), weight);
      h0a->Fill(Pa.M(), weight);
      h0b->Fill(Pb.M(), weight);
      h0c->Fill(Pc.M(), weight);
      h0d->Fill(Pd.M(), weight);
      h0dp->Fill(Pd.P(), weight);
      
      d0a->Fill(kf[2].P(), kf[1].P(), weight);
      d0b->Fill(kf[3].P(), kf[1].P(), weight);
      d0c->Fill(kf[2].P(), kf[3].P(), weight);
      r0a->Fill(kf[1].Theta() * deg, kf[1].P(), weight);
      r0b->Fill(kf[2].Theta() * deg, kf[2].P(), weight);
      r0c->Fill(kf[3].Theta() * deg, kf[3].P(), weight);
      
      
      MisPID(&kf[2]);
      MisPID(&kf[3]);
      
      weight *= Cut(kf[1], kf[2], kf[3], kf[0]);
      
      PP = kf[1] + kf[2] + kf[3];
      Pa = kf[1] + kf[2];
      Pb = kf[1] + kf[3];
      Pc = kf[2] + kf[3];
      Pd = ki[0] + kgold - kf[0] - PP;
      angle1 = kf[1].Angle(kf[2].Vect()); //angle between p kp
      angle2 = kf[1].Angle(kf[3].Vect()); //angle between p km
      angle3 = kf[2].Angle(kf[3].Vect()); //angle between kp km
      angle4 = Pv.Angle(kf[1].Vect());  //angle between virtual photon p
      angle5 = Pv.Angle(kf[2].Vect());  //angle between virtual photon kp
      angle6 = Pv.Angle(kf[3].Vect());  //angle between virtual photon km
      
      h1->Fill(PP.M(), weight);
      h1a->Fill(Pa.M(), weight);
      h1b->Fill(Pb.M(), weight);
      h1c->Fill(Pc.M(), weight);
      h1d->Fill(Pd.M(), weight);
      h1dp->Fill(Pd.P(), weight);
      
      d1a->Fill(kf[2].P(), kf[1].P(), weight);
      d1b->Fill(kf[3].P(), kf[1].P(), weight);
      d1c->Fill(kf[2].P(), kf[3].P(), weight);
      r1a->Fill(kf[1].Theta() * deg, kf[1].P(), weight);
      r1b->Fill(kf[2].Theta() * deg, kf[2].P(), weight);
      r1c->Fill(kf[3].Theta() * deg, kf[3].P(), weight);
      h1e->Fill(angle1 , weight);
      h1f->Fill(angle2 , weight);
      h1g->Fill(angle3 , weight);
      
      s1a->Fill(kf[1].Phi() * deg, kf[1].P(), weight);
      s1b->Fill(kf[2].Phi() * deg, kf[2].P(), weight);
      s1c->Fill(kf[3].Phi() * deg, kf[3].P(), weight);
      h1vp->Fill(angle4 , weight);
      h1vkp->Fill(angle5 , weight);
      h1vkm->Fill(angle6 , weight);
      
      kf[1].Boost(-PP.BoostVector());
      kf[2].Boost(-PP.BoostVector());
      kf[3].Boost(-PP.BoostVector());
      
      PP = kf[1] + kf[2] + kf[3];
      Pa = kf[1] + kf[2];
      Pb = kf[1] + kf[3];
      Pc = kf[2] + kf[3];
      angle1 = kf[1].Angle(kf[2].Vect());
      angle2 = kf[1].Angle(kf[3].Vect());
      angle3 = kf[2].Angle(kf[3].Vect());
      angle4 = Pv.Angle(kf[1].Vect());  //angle between virtual photon p
      angle5 = Pv.Angle(kf[2].Vect());  //angle between virtual photon kp
      angle6 = Pv.Angle(kf[3].Vect());  //angle between virtual photon km
      
      h1z->Fill(PP.M(), weight);
      h1az->Fill(Pa.M(), weight);
      h1bz->Fill(Pb.M(), weight);
      h1cz->Fill(Pc.M(), weight);
      d1az->Fill(kf[2].P(), kf[1].P(), weight);
      d1bz->Fill(kf[3].P(), kf[1].P(), weight);
      d1cz->Fill(kf[2].P(), kf[3].P(), weight);
      r1az->Fill(kf[1].Theta() * deg, kf[1].P(), weight);
      r1bz->Fill(kf[2].Theta() * deg, kf[2].P(), weight);
      r1cz->Fill(kf[3].Theta() * deg, kf[3].P(), weight);
      h1ez->Fill(angle1 , weight);
      h1fz->Fill(angle2 , weight);
      h1gz->Fill(angle3 , weight);
      h1vpz->Fill(angle4 , weight);
      h1vkpz->Fill(angle5 , weight);
      h1vkmz->Fill(angle6 , weight);
      }
    countN++;
    if (countN % 100000 == 0) cout<<countN<<endl;
  }
  

  
  double F1 = Nsim * (2.85-1.93)/ 3000;   // F = Range_x / bin_x        and     Scale(1 / F)  
  double F0 = Nsim * (2.6-1.24)/ 3000;
  
  double F1a = Nsim * (2.2-1.43)/ 3000;
  double F0a = Nsim * (2.0-1.08)/ 3000; 
  
  double F1b = Nsim * (1.87-1.43)/ 3000;
  double F0b = Nsim * (1.80-1.08)/ 3000;
  
  double F1c = Nsim * (1.50-0.98)/ 3000;
  double F0c = Nsim * (1.10-0.28)/ 3000;

  double F1d = Nsim * (182.23-181.80)/ 2200;
  double F0d = Nsim * (182.52-182.45)/ 2200;
  double F1dp = Nsim * 0.6 / 2200;
  double F0dp = Nsim * 0.6 / 2200;
  
  double F1e = Nsim * ( M_PI)/ 1000;
  double F1f = Nsim * ( M_PI)/ 1000;
  double F1g = Nsim * ( M_PI)/ 1000;
  
  h0->Scale(convert);
  h0->Scale(1/F0);
  h1->Scale(convert);
  h1->Scale(1/F1);
  h1z->Scale(convert);
  h1z->Scale(1/F1);
 
  h0a->Scale(convert);
  h0a->Scale(1/F0a);
  h1a->Scale(convert);
  h1a->Scale(1/F1a);
  h1az->Scale(convert);
  h1az->Scale(1/F1a);
  
  h0b->Scale(convert);
  h0b->Scale(1/F0b);
  h1b->Scale(convert);
  h1b->Scale(1/F1b);
  h1bz->Scale(convert);
  h1bz->Scale(1/F1b);
 
  h0c->Scale(convert);
  h0c->Scale(1/F0c);
  h1c->Scale(convert);
  h1c->Scale(1/F1c);
  h1cz->Scale(convert);
  h1cz->Scale(1/F1c);
  
  h0d->Scale(convert);
  h0d->Scale(1/F0d);
  h1d->Scale(convert);
  h1d->Scale(1/F1d);
  
  h0dp->Scale(convert);
  h0dp->Scale(1/F0dp);
  h1dp->Scale(convert);
  h1dp->Scale(1/F1dp);
  
  h1e->Scale(convert);
  h1e->Scale(1/F1e);
  h1f->Scale(convert);
  h1f->Scale(1/F1f);
  h1g->Scale(convert);
  h1g->Scale(1/F1g);
  
  h1vp->Scale(convert);
  h1vp->Scale(1/F1e);
  h1vkp->Scale(convert);
  h1vkp->Scale(1/F1e);
  h1vkm->Scale(convert);
  h1vkm->Scale(1/F1e);
  
  h1vpz->Scale(convert);
  h1vpz->Scale(1/F1e);
  h1vkpz->Scale(convert);
  h1vkpz->Scale(1/F1e);
  h1vkmz->Scale(convert);
  h1vkmz->Scale(1/F1e);
  
  h1ez->Scale(convert);
  h1ez->Scale(1/F1e);
  h1fz->Scale(convert);
  h1fz->Scale(1/F1f);
  h1gz->Scale(convert);
  h1gz->Scale(1/F1g);
  
  
  d0a->Scale(convert/Nsim *100*100/2.0/2.0);  // bin_x * bin_y / Range_x / Range_y 
  d0b->Scale(convert/Nsim *100*100/2.0/2.0);
  d0c->Scale(convert/Nsim *100*100/2.0/2);
   
  d1a->Scale(convert/Nsim *100*100/2.0/2.0); 
  d1b->Scale(convert/Nsim *100*100/2.0/2.0);
  d1c->Scale(convert/Nsim *100*100/2.0/2.0);
  d1az->Scale(convert/Nsim *100*100/1.35/1.35);
  d1bz->Scale(convert/Nsim *100*100/1.35/1.2);
  d1cz->Scale(convert/Nsim *100*100/1.35/1.2);
  
  r0a->Scale(convert/Nsim *90*100/90.0/2.0);
  r0b->Scale(convert/Nsim *90*100/90.0/2.0);
  r0c->Scale(convert/Nsim *90*100/90.0/2.0);
  r1a->Scale(convert/Nsim *90*100/90.0/2.0);
  r1b->Scale(convert/Nsim *90*100/90.0/2.0);
  r1c->Scale(convert/Nsim *90*100/90.0/2.0);
  r1az->Scale(convert/Nsim *90*100/90.0/1.35);
  r1bz->Scale(convert/Nsim *90*100/90.0/1.35);
  r1cz->Scale(convert/Nsim *90*100/90.0/1.2);
  
  
  s1a->Scale(convert/Nsim *100/2.0);
  s1b->Scale(convert/Nsim *100/2.0);
  s1c->Scale(convert/Nsim *100/2.0);
  

  infile.close();
  fs->Write();

  h0->Draw();
  h1->Draw();
  
  
  return 0;
} 

int MisPID(TLorentzVector * k){
  k->SetE(sqrt(PARTICLE::K.M() * PARTICLE::K.M() + k->P() * k->P()));
  return 0;
}

  int Cut(const TLorentzVector p, const TLorentzVector Kp, const TLorentzVector Km, const TLorentzVector kfe){
  TLorentzVector ki;
  TLorentzVector kgold;
  ki.SetXYZT(0, 0, sqrt(4.4 * 4.4 - PARTICLE::e.M() * PARTICLE::e.M()), 4.4);
  kgold.SetXYZT(0, 0, 0, Mp * 183.43);
    
  TLorentzVector kk = p + Kp;
  if (kk.M() > 1.48) return 0;
  kk = p + Km;
  if (kk.M() > 1.48) return 0;
  kk = Kp + Km;
  if (kk.M() > 1.04) return 0;
  
  if (p.P() >0.8) return 0;
  if (Kp.P() > 0.5) return 0;
  if (Km.P() >0.5) return 0;
  if (p.Theta() > M_PI / 3.0) return 0;
  
  TLorentzVector miss= ki + kgold - kfe - p - Kp - Km;
  if (miss.M() < 182.49) return 0; 
  
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
