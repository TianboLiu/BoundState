#include "Lcore.h"

int main(){

  Initialize();

  TLorentzVector k(0.5, 0.1, 1.6, 2.5);
  //cout << DETECTOR::Acceptance(k, "K+") << endl;

  cout << DETECTOR::acc_pip->GetBinContent(2,3,4) << endl;

  if (false){
    TFile * fs =  new TFile("acceptance/clasev_acceptance.root", "r");
    TH3F * ac = (TH3F *) fs->Get("acceptance_PThetaPhi_pip");
    TH3F * tt = ac;
    cout << tt->GetBinContent(tt->GetXaxis()->FindBin(k.Phi() * 180.0 / M_PI), tt->GetYaxis()->FindBin(k.Theta() * 180.0 / M_PI), tt->GetZaxis()->FindBin(k.P())) << endl;
  }

  return 0;
}
