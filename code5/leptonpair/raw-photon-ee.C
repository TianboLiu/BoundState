#include "../Lcore.h"

int main(const int argc, const char * argv[]){

  if (argc < 2){
    cout << "./raw-photon-ee <Nsim>" << endl;
    return 0;
  }

  TString filename = "result/raw-photon-ee.root";

  Long64_t Nsim;
  Nsim = atoi(argv[1]);

  Initialize();

  TLorentzVector ki[2], kf[4];
  double weight = 0;

  TFile * fs = new TFile(filename.Data(), "RECREATE");

  TH2D * PTheta_phi = new TH2D("PTheta_phi", "", 90, 0.0, 180.0, 100, 0.0, 3.0);
  TH2D * PTheta_electron = new TH2D("PTheta_electron", "", 90, 0.0, 180.0, 100, 0.0, 3.0);
  TH2D * PTheta_positron = new TH2D("PTheta_positron", "", 90, 0.0, 180.0, 100, 0.0, 3.0);
  TH2D * ThetaTheta = new TH2D("ThetaTheta", "", 90, 0.0, 180.0, 90, 0.0, 180.0);
  TH2D * MomMom = new TH2D("PP", "", 100, 0.0, 3.0, 100, 0.0, 3.0);
  TH2D * ThetaAngle = new TH2D("ThetaAngle_pair", "", 90, 0.0, 180.0, 90, 0.0, 180.0);
  TH2D * PAngle = new TH2D("PAngle_pair", "", 100, 0.0, 3.0, 90, 0.0, 180.0);

  PTheta_phi->SetDirectory(fs);
  PTheta_electron->SetDirectory(fs);
  PTheta_positron->SetDirectory(fs);
  ThetaTheta->SetDirectory(fs);
  MomMom->SetDirectory(fs);
  ThetaAngle->SetDirectory(fs);
  PAngle->SetDirectory(fs);
  
  const double deg = 180.0 / M_PI;
  double range[2] = {1.30, 1.60};
  TLorentzVector PP;

  for (Long64_t i = 0; i < Nsim; i++){
    if (i%(Nsim/10)==0) cout << i / (Nsim/100) << "%" << endl;

    GENERATE::BremsstrahlungPhoton(ki, range);

    weight = GENERATE::Event_Nee_Phi(ki, kf);
    if (weight > 0){
      PP = kf[1] + kf[2];
      PTheta_phi->Fill(PP.Theta() * deg, PP.P(), weight);
      PTheta_electron->Fill(kf[1].Theta() * deg, kf[1].P(), weight);
      PTheta_positron->Fill(kf[2].Theta() * deg, kf[2].P(), weight);
      ThetaTheta->Fill(kf[1].Theta() * deg, kf[2].Theta() * deg, weight);
      MomMom->Fill(kf[1].P(), kf[2].P(), weight);
      ThetaAngle->Fill(kf[2].Theta() * deg, kf[2].Angle(kf[1].Vect()) * deg, weight);
      PAngle->Fill(kf[2].P(), kf[2].Angle(kf[1].Vect()) * deg, weight);
    }
  }

  PTheta_phi->Scale(1.0/Nsim);
  PTheta_electron->Scale(1.0/Nsim);
  PTheta_positron->Scale(1.0/Nsim);
  ThetaTheta->Scale(1.0/Nsim);
  MomMom->Scale(1.0/Nsim);
  ThetaAngle->Scale(1.0/Nsim);
  PAngle->Scale(1.0/Nsim);

  fs->Write();

  return 0;
}
  
