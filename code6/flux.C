/* test virtual photon flux */

#include "Lcore.h"

int main(const int argc, const char * argv[]){

  gRandom->SetSeed(0);

  if (argc < 3){
    cout << "./flux <Ebeam> <savefile>" << endl;
    return 0;
  }

  Long64_t nsim = 10000000;

  // Electron beam energy
  double Ebeam = atof(argv[1]);

  // Set scattered electron range
  GENERATE::perange[0] = 0.0;
  GENERATE::perange[1] = Ebeam;
  double degtorad = M_PI / 180.0;
  GENERATE::cthrange[0] = cos(5.0 * degtorad);
  GENERATE::cthrange[1] = cos(0.0 * degtorad);

  TString filename = argv[2];

  TLorentzVector ki[2], kf[2], q;
  ki[0].SetXYZM(0, 0, Ebeam, PARTICLE::e.M());
  ki[1].SetXYZM(0, 0, 0, Mp);

  double weight = 0.0;

  TFile * fs = new TFile(filename.Data(), "RECREATE");
  TH1D * h0 = new TH1D("h0", ";E_{#gamma};", 100, 0.0, 10.0);
  TH1D * h1 = new TH1D("h1", ";Q^{2} (GeV^{2});", 100, 0.0, 10.0);
  TH1D * h2 = new TH1D("h2", ";E_{eff};", 100, 0.0, 10.0);
  
  for (Long64_t i = 0; i < nsim; i++){
    if (i%(nsim/10)==0) cout << i << endl;

    weight = GENERATE::VirtualPhoton(ki, kf);

    if (weight > 0){
      q = ki[0] - kf[0];
      h0->Fill(q.E(), weight);
      h1->Fill(-q * q, weight);
      h2->Fill(q.E() + q * q / (2.0 * 1.8756), weight);
    }
  }

  h0->Scale(1.0/nsim);
  h1->Scale(1.0/nsim);
  h2->Scale(1.0/nsim);
  fs->Write();
  fs->Close();

  return 0;
}
