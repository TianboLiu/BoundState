#include "Lgold.h"

using namespace std;

int main(int argc, char * argv[]){

  double Eb = atof(argv[1]);
 
  gRandom->SetSeed(0);
  SetFunctions();
  SetTagger(0.5, Eb);
  //SetTagger(0.5, Eb, 0.0, M_PI / 2.0);
  
  TLorentzVector ki[2], kf[6];
  double weight;
  ki[0].SetXYZT(0, 0, Eb, Eb);

  Long64_t Nsim = 100000000;
  
  TH1D * h0 = new TH1D("h0", "", 1, 0.0, 1.0);
  TH1D * h1 = new TH1D("h1", "", 1, 0.0, 1.0);

  Long64_t j0 = 0;
  Long64_t j1 = 0;
  for (Long64_t i = 0; i < Nsim; i++){
    weight = GenerateElectroproductionBoundStateGold(ki, kf) * 197.0 * 196.0;
 
    if (weight > 0){
      h0->Fill(0.5, weight);
      j0++;
    }
  }

  bool check = true;
  if (check){
    for (Long64_t i = 0; i < Nsim; i++){
      weight = GenerateElectroproductionBoundStateGold(ki, kf) * 197.0 * 196.0;
      
      if (weight > 0){
	h1->Fill(0.5, weight);
	j1++;
      }
    }
  }

  if (check){
    printf("%.1f\t %.2E\t %.2E\t %.6f\t %.6f\n",
	   Eb, 
	   h0->Integral(1, -1) * 3.89379e5 / Nsim, 
	   h1->Integral(1, -1) * 3.89379e5 / Nsim,
	   ((double) j0) / Nsim,
	   ((double) j1) / Nsim);
  }
  else{
    cout << Eb << "    " << h0->Integral(1, -1) * 3.89379e5 / Nsim << endl;
  }

  return 0;
}
  
