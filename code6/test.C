#include "Lcore.h"

int main(const int argc, const char * argv[]){

  gRandom->SetSeed(0);

  // Electron beam energy
  double Ebeam = 10.6;//GeV
  
  // Set nuclear
  NUCLEAR::SetNuclear("C12");  
  GENERATE::TF_fMomentum = new TF1("fp", NUCLEAR::fMomentum, 0.0, 1.0, 0);
  GENERATE::TF_fMomentum->SetNpx(1000);
  GENERATE::TF_fEnergy = new TF1("fE", NUCLEAR::fEnergy, 0.0, 0.3, 0);
  GENERATE::TF_fEnergy->SetNpx(1000);
  
  // Set Jpsi production model
  JPSIMODEL::SetModel("23g");

  // Set scattered electron range
  GENERATE::cthrange[0] = cos(4.5/180.0 * M_PI);
  GENERATE::cthrange[1] = cos(2.5/180.0 * M_PI);
  GENERATE::perange[0] = 0.5;//GeV
  GENERATE::perange[1] = 5.0;//GeV


  TLorentzVector ki[2], kf[4];
  ki[0].SetXYZM(0, 0, Ebeam, PARTICLE::e.M());
  double weight = 0.0;
  double acceptance = 0.0;
  

  for (int i = 0; i < 100; i++){
    weight = GENERATE::GetNucleon(&ki[1]);

    weight *= GENERATE::Event_eN2eNee_Jpsi(ki, kf); 

    if (weight > 0.0){
      cout << kf[2].P() << endl;
    }
  }

  

  return 0;
}
