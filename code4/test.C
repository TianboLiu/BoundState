#include "Lcore.h"

int main(){

  Initialize();

  TLorentzVector k(0.1, 0.0, 0.1, 0.5);
  //cout << DETECTOR::Acceptance(k, "K+") << endl;

  //cout << DETECTOR::AcceptanceBONUS12(k, "K+") << endl;

  GENERATE::NPiPiWeight();

  if (false){
    TFile * fs =  new TFile("acceptance/clasev_acceptance.root", "r");
    TH3F * ac = (TH3F *) fs->Get("acceptance_PThetaPhi_pip");
    TH3F * tt = ac;
    cout << tt->GetBinContent(tt->GetXaxis()->FindBin(k.Phi() * 180.0 / M_PI), tt->GetYaxis()->FindBin(k.Theta() * 180.0 / M_PI), tt->GetZaxis()->FindBin(k.P())) << endl;
  }

  return 0;
}
