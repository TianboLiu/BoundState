#include "../Lelectro.h"

using namespace std;

int main(){

  gRandom->SetSeed(0);
  SetFunctions();
  SetTagger(3.5, 4.5, 0.0, M_PI-0.01);

  Long64_t Nsim = 20000000;

  TLorentzVector ki[2];
  ki[0].SetXYZT(0, 0, 5.5, 5.5);
  ki[1].SetXYZT(0, 0, 0, 0.938272);
  TLorentzVector Pout = ki[0] + ki[1];
  //cout << Pout.M() << endl;
  TLorentzVector kf[2];
  double weight;
  double w1, w2;
  
  TH1D * h0 = new TH1D("h0", "", 1, 0.0, 1.0);

  for (Long64_t i = 0; i < Nsim; i++){
    if (i%1000000 == 0) cout << i << endl;
    
    //w1 = GenerateScatteredElectron(ki, kf);
    //w2 = GeneratePhotoproductionPhiNucleon(ki, kf);
    weight = GenerateElectroproductionPhiNucleon(ki, kf);
    //weight = w1 * w2;
    if (weight > 0){
      h0->Fill(0.5, weight);
    }
  }

  cout << h0->Integral(1, -1) / Nsim * 389.379 << endl;

  return 0;
}
