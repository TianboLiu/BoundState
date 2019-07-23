#include "Lcore.h"

int main(const int argc, const char * argv[]){

  if (argc < 3){
    cout << "./quasi-elastic <Ebeam> <Nsim>" << endl;
    return 0;
  }

  TFile * fs = new TFile("quasi-elastic.root", "RECREATE");
  TH2D * h2 = new TH2D("quasi-elastic", "", 180, 0, M_PI, 200, 0.0, 2.0);
  h2->GetXaxis()->SetTitle("Theta(rad)");
  h2->GetYaxis()->SetTitle("P(GeV)");

  double binx = M_PI / 180;
  double biny = 2.0 / 200;

  double Ebeam = atof(argv[1]);
  Long64_t Nsim = (Long64_t) atof(argv[2]);

  Initialize();
  TLorentzVector ki(0, 0, Ebeam, Ebeam);
  TLorentzVector kf[2];

  double weight = 0;

  for (Long64_t i = 0; i < Nsim; i++){
    if (i%(Nsim/10)==0) cout << i/(Nsim/10)*10 << "%" << endl;
    weight = GENERATE::Event_ep_QE(&ki, kf);
    h2->Fill(kf[1].Theta(), kf[1].P(), weight);
  }

  h2->Scale(1.0/binx/biny/Nsim);

  fs->Write();

  return 0;
}
