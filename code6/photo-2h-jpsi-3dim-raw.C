/* To generate and same the events */

#include "Lcore.h"

int main(const int argc, const char * argv[]){

  gRandom->SetSeed(0);

  if (argc < 7){
    cout << "./photo-solid-2h-jpsi-4dim <model> <Ebeam> <filename> <nsim> <Nfiles> <idx0>" << endl;
    return 0;
  }
  Long64_t nsim = atoi(argv[4]);
  int Nfiles = atoi(argv[5]);
  int j0 = atoi(argv[6]);

  // Electron beam energy and luminosity
  double Ebeam = atof(argv[2]);//GeV
    
  // Set Jpsi production model
  JPSID3d::SetModel(argv[1]);

  // Set bremsstrahlung photon
  GENERATE::SetBremsstrahlung();
  double kmin = 7.2;
  double kmax = Ebeam;

  TString filename = argv[3];

  TLorentzVector ki[2], kf[4];
  ki[0].SetXYZM(0, 0, Ebeam, PARTICLE::e.M());
  double weight = 0.0;
  
  FILE * f;
  TString name;
  for (int j = j0; j < j0 + Nfiles; j++){
    cout << j << " / " << "[" << j0 << "," << j0+Nfiles << "]" << endl;
    name = filename + Form("%.4d.dat", j);
    f = fopen(name.Data(), "w");
    fprintf(f, "Ebeam=%.2fGeV,Nsim= %lld\n", Ebeam, nsim);
    for (Long64_t i = 0; i < nsim; i++){
      weight = GENERATE::BremsstrahlungPhoton(&ki[0], kmin, kmax, Ebeam) * 1.95 / 2;//15cm LD2 target
      weight *= GENERATE::Event_gD2eep_Jpsi(ki, kf);
      if (weight > 0.0){
	fprintf(f, "%.6E\n", weight);
	fprintf(f, "q:\t%.6E\t%.6E\t%.6E\t%.6E\n", ki[0].X(), ki[0].Y(), ki[0].Z(), ki[0].E());
	fprintf(f, "e+:\t%.6E\t%.6E\t%.6E\t%.6E\n", kf[0].X(), kf[0].Y(), kf[0].Z(), kf[0].E());
	fprintf(f, "e-:\t%.6E\t%.6E\t%.6E\t%.6E\n", kf[1].X(), kf[1].Y(), kf[1].Z(), kf[1].E());
	fprintf(f, "p:\t%.6E\t%.6E\t%.6E\t%.6E\n", kf[2].X(), kf[2].Y(), kf[2].Z(), kf[2].E());
      }
    }
    fclose(f);
  }

  return 0;
}
