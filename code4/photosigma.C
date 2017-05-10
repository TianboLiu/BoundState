#include "Lcore.h"

int main(const int argc, const char * argv[]){
  
  if (argc < 2) return 0;
  
  MODEL::SetMODEL();
  GOLD::SetGOLD();
  
  Long64_t Nsim = 10000000;
  
  TH1D * h0 = new TH1D("h0", "", 1, 0.0, 2.0);

  TH1D * hm = new TH1D("hm", "", 100, 4.0, 4.1);

  TLorentzVector ki[2], kf[2];
  double weight;

  double Ephoton = atof(argv[1]);
  ki[0].SetXYZT(0.0, 0.0, Ephoton, Ephoton);

  for (Long64_t i = 0; i < Nsim; i++){
    //if (i%1000000 == 0) cout << i << endl;
    weight = GENERATE::BoundStatePhotoproductionGold(ki, kf);
    //weight = GENERATE::PhiPhotoproductionGold(ki, kf);
    if (weight > 0){
      h0->Fill(1.0, weight);
      //hm->Fill(kf[1].M(), weight);
    }
  }

  cout << Ephoton << "\t" << h0->Integral(1, -1) * 0.389379e6 / Nsim << "  nb"<< endl;
  
  // TCanvas * c0 = new TCanvas("c0", "", 800, 600);
  // hm->SetStats(0);
  // hm->Draw();

  // c0->Print("c0.pdf");

  return 0;
}
