#include "Lgold.h"

using namespace std;

int main(int argc, char * argv[]){

  double Eb = atof(argv[1]);
  double Q2 = atof(argv[2]);

  gRandom->SetSeed(0);
  SetFunctions();
  
  TLorentzVector ki[2], kf[6];
  double weight;
  ki[0].SetXYZT(0, 0, sqrt(Eb*Eb+Q2), Eb);

  Long64_t Nsim = 100000000;
  
  TH1D * h0 = new TH1D("h0", "", 1, 0.0, 1.0);
  TH1D * h1 = new TH1D("h1", "", 1, 0.0, 1.0);

  for (Long64_t i = 0; i < Nsim; i++){
    weight = GeneratePhotoproductionBoundStateGold(ki, kf);
    //weight = GenerateEvent_NKKN_withPhotoproductionBoundStateGold(ki, kf);

    if (weight > 0){
      h0->Fill(0.5, weight * 197.0 * 196);
    }
  }

  bool check = true;
  if (check){
    for (Long64_t i = 0; i < Nsim; i++){
      weight = GeneratePhotoproductionBoundStateGold(ki, kf);
      //weight = GenerateEvent_NKKN_withPhotoproductionBoundStateGold(ki, kf);
      
      if (weight > 0){
	h1->Fill(0.5, weight * 197 * 196.0);
      }
    }
  }

  if (check){
    cout << Eb << "    " << Q2 << "    " << h0->Integral(1, -1) * 3.89379e5 / Nsim << "    " << h1->Integral(1, -1) * 3.89379e5 / Nsim << endl;
  }
  else{
    cout << Eb << "    " << Q2 << "    " << h0->Integral(1, -1) * 3.89379e5 / Nsim << endl;
  }

  return 0;
}
  
