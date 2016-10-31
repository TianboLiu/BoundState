#include "../Lelectro.h"

using namespace std;

int main(int argc, char * argv[]){

  gRandom->SetSeed(0);
  SetFunctions();
  SetTagger(2.5, 4.5);

  Long64_t Nsim = atoi(argv[2]);

  double eq = atof(argv[1]);

  TLorentzVector ki[2];
  ki[0].SetXYZT(0, 0, eq, eq);
  ki[1].SetXYZT(0, 0, 0, 0.938272);
  TLorentzVector Pout = ki[0] + ki[1];
  //cout << Pout.M() << endl;
  TLorentzVector kf[5];

  double weight;
  double w1, w2;

  TH1D * h0 = new TH1D("h0", "", 1000, 1.8, 2.5);

  Long64_t j = 0;
  for (Long64_t i = 0; i < Nsim; i++){
    
    if (i%1000000 == 0) cout << i << endl;
    
    //w1 = GenerateScatteredElectron(ki, kf);
    //weight = GeneratePhotoproductionPhiNucleon(ki, kf);
    //cout << ki[0].E() << endl;
    //weight = GenerateElectroproductionPhiNucleon(ki, kf);
    //weight = GenerateElectroproductionPhiCarbon(ki, kf);
    weight = GeneratePhotoproductionBoundState(ki, kf);
    weight = GenerateElectroproductionBoundState(ki, kf);
    //weight = GenerateScatteredElectron(ki, kf);
    //weight = w1 * w2;
    if (weight > 0){
      h0->Fill(kf[2].M(), weight);
      j++;
    }
  }

  cout << h0->Integral(1, -1) / Nsim * 389.379 * 132000 << endl;
  cout << j << endl;

  TCanvas * c0 = new TCanvas("c0", "", 800, 600);
  h0->Draw();

  c0->Print("c0.pdf");

  return 0;
}
