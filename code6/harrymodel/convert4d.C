#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TH3D.h"
#include "TString.h"

using namespace std;

int main(const int argc, const char * argv[]){

  if (argc < 2){
    cout << "./convert4d <opt>" << endl;
    return 0;
  }

  const int opt = atoi(argv[1]);

  double eps = 1e-10;
  double zero = 1e-200;

  if (opt == 1){//jpsi from d, 4-fold
    double q, theta_k, k, phi_p, theta_p, p, crst, tmp;
    double deg = M_PI / 180.0;
    char dum[200];

    TFile * fs = new TFile("harrymodel-jpsi-d-4D.root", "RECREATE");
    TH3D * hs, * hp;
    double vol = 0.01 * (2.0 * deg) * (2.0 * M_PI) * (1.0 * deg) * (2.0 * M_PI);
    for (int i = 1; i <= 10; i++){
      ifstream infile(Form("7p2-phip/fort.90%.2d", i));
      hs = new TH3D(Form("ds_E7.2_idx%d", i), "", 700, 3.0, 10.0, 6, 0.0, 12.0, 25, 5.5, 30.5);
      hp = new TH3D(Form("p_E7.2_idx%d", i), "", 700, 3.0, 10.0, 6, 0.0, 12.0, 25, 5.5, 30.5);
      infile.getline(dum, 200);
      while (infile >> q >> theta_k >> k >> phi_p >> theta_p >> p >> crst >> tmp >> tmp){
	hs->Fill(k, theta_k + eps, theta_p, crst * sin((theta_k+1.0)*deg) * sin(theta_p*deg) * vol + zero);
	hp->Fill(k, theta_k + eps, theta_p, p);
      }
      infile.close();
      cout << hs->Integral() << endl;
      fs->Write();
      hs->Delete();
      hp->Delete();
    }
    for (int i = 1; i <= 10; i++){
      ifstream infile(Form("8p2-phip/fort.90%.2d", i));
      hs = new TH3D(Form("ds_E8.2_idx%d", i), "", 700, 3.0, 10.0, 6, 0.0, 12.0, 25, 5.5, 30.5);
      hp = new TH3D(Form("p_E8.2_idx%d", i), "", 700, 3.0, 10.0, 6, 0.0, 12.0, 25, 5.5, 30.5);
      infile.getline(dum, 200);
      while (infile >> q >> theta_k >> k >> phi_p >> theta_p >> p >> crst >> tmp >> tmp){
	hs->Fill(k, theta_k + eps, theta_p, crst * sin((theta_k+1.0)*deg) * sin(theta_p*deg) * vol + zero);
	hp->Fill(k, theta_k + eps, theta_p, p);
      }
      infile.close();
      cout << hs->Integral() << endl;
      fs->Write();
      hs->Delete();
      hp->Delete();
    }
    fs->Close();
  }

  return 0;
}
    
    
