/* Generate and same events */

#include "Lcore.h"

int main(const int argc, const char * argv[]){

  gRandom->SetSeed(0);

  if (argc < 6){
    cout << "./electro-solid-2h-jpsi-cut <model> <Ebeam> <filename> <Nsim> <Nfiles>" << endl;
    return 0;
  }
  Long64_t Nsim = atoi(argv[4]);
  int Nfiles = atoi(argv[5]);
  Long64_t nsim = Nsim / Nfiles;

  // Electron beam energy and luminosity
  double Ebeam = atof(argv[2]);//GeV
    
  // Set Jpsi production model
  JPSID4d::SetModel(argv[1]);

  // Set scattered electron range
  double degtorad = M_PI / 180.0;
  GENERATE::cthrange[0] = cos(40.0 * degtorad);
  GENERATE::cthrange[1] = cos(0.0 * degtorad);
  GENERATE::perange[0] = 0.0;
  GENERATE::perange[1] = Ebeam - 7.2;

  TString filename = argv[3];

  TLorentzVector ki[2], kf[4], q;
  ki[0].SetXYZM(0, 0, Ebeam, PARTICLE::e.M());
  double weight = 0.0;

  FILE * f;
  TString name;
  for (int j = 0; j < Nfiles; j++){
    cout << j << " / " << Nfiles << endl;
    name = filename + Form("%.4d.dat", j);
    f = fopen(name.Data(), "w");
    fprintf(f, "Ebeam=%.2fGeV, Nsim=%lld\n", Ebeam, Nsim);
    for (Long64_t i = 0; i < nsim; i++){
        weight = GENERATE::GetNucleon(&ki[1]);
      weight *= GENERATE::Event_eD2eeep_Jpsi(ki, kf);
      if (weight > 0.0){
	fprintf(f, "%.6E\n", weight);
	fprintf(f, "e':\t%.6E\t%.6E\t%.6E\n", kf[0].X(), kf[0].Y(), kf[0].Z());
	fprintf(f, "e+:\t%.6E\t%.6E\t%.6E\n", kf[1].X(), kf[1].Y(), kf[1].Z());
	fprintf(f, "e-:\t%.6E\t%.6E\t%.6E\n", kf[2].X(), kf[2].Y(), kf[2].Z());
	fprintf(f, "p:\t%.6E\t%.6E\t%.6E\n", kf[3].X(), kf[3].Y(), kf[3].Z());
      }
    }
    fclose(f);    
  }

  return 0;
}
