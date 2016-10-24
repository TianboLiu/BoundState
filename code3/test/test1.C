#include "../Lelectro.h"

using namespace std;

int main(){

  gRandom->SetSeed(1);
  SetFunctions();

  Long64_t Nsim = 2000000;

  TLorentzVector ki[2];
  ki[0].SetXYZT(0, 0, 1.7, 1.7);
  ki[1].SetXYZT(0, 0, 0, 0.938272);
  TLorentzVector kf[2];
  double weight;
  
  TH1D * h0 = new TH1D("h0", "", 1, 0.0, 1.0);

  for (Long64_t i = 0; i < Nsim; i++){
    weight = GeneratePhotoproductionPhi(ki, kf);
    if (weight > 0){
      h0->Fill(0.5, weight);
    }
  }

  cout << h0->Integral(1, -1) / Nsim * 389.379 << endl;

  return 0;
}
