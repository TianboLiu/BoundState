#include "Lcore.h"

int main(){

  Initialize();

  TLorentzVector ki[2], kf[5], PP;
  ki[0].SetXYZT(0, 0, 4.4, 4.4);
  ki[1].SetXYZT(0, 0, 0, Mp);
  double weight = 0;
  for (int i = 0; i < 0; i++){
    weight = GENERATE::KKPhotoproduction(ki, kf);
    if (weight > 0){
      PP = kf[1] + kf[2];
      if (PP.M() < 1.0)
	cout << PP.M() << "  " << weight << endl;	\
    }
  }
  

  if (true){
  TH1D * hh = new TH1D("hh", "", 2200, 1.88, 3.32);
  for (Long64_t i = 0; i < 10000000; i++){
    if (i%1000000 == 0) cout << i << endl;
    weight = GENERATE::Event_eNKK_KK(ki, kf);
    if (weight > 0){
      PP = kf[1] + kf[2] + kf[3];
      hh->Fill(PP.M(), weight); 
    }
  }
  hh->Scale(1.0/1000000);

  TCanvas * c0 = new TCanvas("c0", "", 800, 600);
  hh->Draw();

  c0->Print("c0.pdf");
  }

  return 0;
}
