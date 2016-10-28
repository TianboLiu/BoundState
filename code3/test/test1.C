#include "../Lelectro.h"

using namespace std;

int main(){

  gRandom->SetSeed(0);
  SetFunctions();
  SetTagger(2.5, 4.0);

  Long64_t Nsim = 10000000;

  TLorentzVector ki[2];
  ki[0].SetXYZT(0, 0, 5.5, 5.5);
  ki[1].SetXYZT(0, 0, 0, 0.938272);
  TLorentzVector Pout = ki[0] + ki[1];
  //cout << Pout.M() << endl;
  TLorentzVector kf[5];

  double weight;
  double w1, w2;

  cout <<TF_fE.Integral(0.0, 0.1) << endl;
  
  TH1D * h0 = new TH1D("h0", "", 1, 0.0, 1.0);

  Long64_t j = 0;
  for (Long64_t i = 0; i < 0; i++){
    if (i%1000000 == 0) cout << i << endl;
    
    //w1 = GenerateScatteredElectron(ki, kf);
    //w2 = GeneratePhotoproductionPhiNucleon(ki, kf);
    //cout << ki[0].E() << endl;
    //weight = GenerateElectroproductionPhiNucleon(ki, kf);
    weight = GenerateElectroproductionPhiCarbon(ki, kf);
    //weight = GenerateScatteredElectron(ki, kf);
    //weight = w1 * w2;
    if (weight > 0){
      h0->Fill(0.5, weight);
      j++;
    }
  }

  cout << h0->Integral(1, -1) / Nsim * 389.379 << endl;
  cout << j << endl;

  return 0;
}
