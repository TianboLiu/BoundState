#include "Lcore.h"

double CalculateInternalMomentum(const double M){
  double M1 = 3.0969;
  double M2 = 0.938272;
  return sqrt((M*M - pow(M1 + M2, 2)) * (M*M - pow(M1 - M2, 2))) / (2.0 * M);
}

int main(const int argc, const char * argv[]){

  //gErrorIgnoreLevel = 5000;
  // Set simulation
  gRandom->SetSeed(0);

  if (argc < 5){
    cout << "./photo-solid-2h-jpsi-4dim <model> <Ebeam> <filename> <Nsim>" << endl;
    return 0;
  }
  Long64_t Nsim = atoi(argv[4]);

  // Electron beam energy and luminosity
  double Ebeam = atof(argv[2]);//GeV
  double lumi = 1.2e37 * 1.0e-26 * pow(0.197327, 2);//GeV^2 s^-1 eN
  double time = 3600.0;//s
    
  // Set Jpsi production model
  JPSID4d::SetModel(argv[1]);

  // Set bremsstrahlung photon
  GENERATE::SetBremsstrahlung();
  double kmin = 7.2;
  double kmax = Ebeam;

  // Set Acceptance
  DETECTOR::SetDetector("SoLID");

  TString filename = argv[3];

  // detected subthreshold 1
  TFile * fsub = new TFile("result-photon/"+filename, "RECREATE");
  TH1D * hMJpsi2 = new TH1D("Mass_e+e-_Jpsi", ";M[e^{+}e^{-}] (GeV);Events / hour", 100, 2.6, 3.6);
  TH2D * hMomentum2 = new TH2D("Pe+Pe-_Jpsi", ";P[e^{+}] (GeV);P[e^{-}] (GeV)", 100, 0.0, 10.0, 100, 0.0, 10.0);
  TH2D * hThetaPelectron2 = new TH2D("ThetaP_e-_Jpsi", ";#theta[e^{-}] (deg);P[e^{-}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPpositron2 = new TH2D("ThetaP_e+_Jpsi", ";#theta[e^{+}] (deg);P[e^{+}] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH2D * hThetaPJpsi2 = new TH2D("ThetaPJpsi", ";#theta[J/#psi] (deg);P[J/#psi] (GeV)", 90, 0.0, 180.0, 100, 0.0, 10.0);
  TH1D * hkappa2 = new TH1D("kappa_Jpsi", ";#kappa_{J/#psi} (GeV);Events / hour", 100, 0.0, 1.0);
  hMJpsi2->SetDirectory(fsub);
  hMomentum2->SetDirectory(fsub);
  hThetaPelectron2->SetDirectory(fsub);
  hThetaPpositron2->SetDirectory(fsub);
  hThetaPJpsi2->SetDirectory(fsub);
  hkappa2->SetDirectory(fsub);

  TLorentzVector ki[2], kf[4];
  ki[0].SetXYZM(0, 0, Ebeam, PARTICLE::e.M());
  double weight = 0.0;
  //double acceptance = 0.0;

  for (Long64_t i = 0; i < Nsim; i++){
    if (i % (Nsim/100) == 0) cout << i/(Nsim/100) << "%" << endl;
    
    weight = GENERATE::BremsstrahlungPhoton(&ki[0], kmin, kmax, Ebeam) * 1.95 / 2;//15cm LD2 target
    weight *= GENERATE::Event_gD2eep_Jpsi(ki, kf);

    weight *= DETECTOR::SmearSoLID(kf[0], "e+") * DETECTOR::SmearSoLID(kf[1], "e-") * DETECTOR::SmearSoLID(kf[2], "p");

    if (weight > 0.0){

      if (ki[0].E() < 8.2){
	hMJpsi2->Fill( (kf[0]+kf[1]).M(), weight);
	hMomentum2->Fill( kf[0].P(), kf[1].P(), weight);
	hThetaPelectron2->Fill( kf[1].Theta()/M_PI*180, kf[1].P(), weight);
	hThetaPpositron2->Fill( kf[0].Theta()/M_PI*180, kf[0].P(), weight);
	hThetaPJpsi2->Fill( (kf[0]+kf[1]).Theta()/M_PI*180, (kf[0]+kf[1]).P(), weight);
	hkappa2->Fill( CalculateInternalMomentum((kf[0]+kf[1]+kf[2]).M()), weight);
      }
    }
    
  }

  hMJpsi2->Scale(lumi*time/Nsim);
  hMomentum2->Scale(lumi*time/Nsim);
  hThetaPelectron2->Scale(lumi*time/Nsim);
  hThetaPpositron2->Scale(lumi*time/Nsim);
  hThetaPJpsi2->Scale(lumi*time/Nsim);
  hkappa2->Scale(lumi*time/Nsim);
  fsub->Write();
  fsub->Close();

  return 0;
}
