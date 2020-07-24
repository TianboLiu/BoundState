#include "Lcore.h"

double Getkappa(const TLorentzVector P){
  double M = P.M();
  double M1 = 3.0969;
  double M2 = 0.938272;
  return sqrt((M*M - pow(M1 + M2, 2)) * (M*M - pow(M1 - M2, 2))) / (2.0 * M);
}

double CalcEg(const TLorentzVector P){
  double Md = 1.8756;
  double Mn = 0.93957;
  double Eg = (pow(P.E() - Md, 2) - pow(P.P(), 2) - Mn * Mn) / (2.0 * (P.E() - P.Pz() - Md));
  return Eg;
}
  
int main(const int argc, const char * argv[]){

  if (argc < 5){
    cout << "./fill <loadfile> <Nfiles> <savefile> <opt>" << endl;
    cout << "1: raw" << endl;
    cout << "2: 3fold detection" << endl;
    cout << "3: 3fold smeared" << endl;
    return 0;
  }

  const int Nfiles = atoi(argv[2]);
  const int opt = atoi(argv[4]);
  
  const double lumi = 0.6e37 * 1e-26 * pow(0.197326,2);
  const double time = 3600.0;
  const double deg = 180.0 / M_PI;
  const double Md = 1.8756;

  DETECTOR::SetDetector("SoLID");

  TString loadfile = argv[1];
  TString name;
  char tmp[200];
  double weight, acc1, acc2;
  TLorentzVector l, lp;
  l.SetXYZM(0, 0, 8.5, PARTICLE::e.M());
  TLorentzVector q, ep, em, p;
  double px, py, pz, E;
  TFile * fs = new TFile(argv[3], "RECREATE");
  TH1D * hMass = new TH1D("Mass_e+e-", ";M[e^{+}e^{-}] (GeV);Events / hour", 100, 2.6, 3.6);
  TH1D * hMassFalse = new TH1D("Mass_e+e-_False", ";M[e^{+}e^{-}] (GeV);Events / hour", 100, 2.6, 3.6);
  hMass->SetDirectory(fs);
  hMassFalse->SetDirectory(fs);
  
  double Nsim = 0.0;
  double nsim;

  if (opt == 2){//detected
    for (int j = 0; j < Nfiles; j++){
      name = loadfile + Form("%.4d.dat",j);
      ifstream infile(name.Data());
      infile >> tmp >> nsim;
      Nsim += nsim;
      while (infile >> weight){
	infile >> tmp >> px >> py >> pz >> E;
	lp.SetXYZM(px, py, pz, PARTICLE::e.M());
	infile >> tmp >> px >> py >> pz >> E;
	ep.SetXYZM(px, py, pz, PARTICLE::e.M());
	infile >> tmp >> px >> py >> pz >> E;
	em.SetXYZM(px, py, pz, PARTICLE::e.M());
	infile >> tmp >> px >> py >> pz >> E;
	p.SetXYZM(px, py, pz, Mp);
	q = l - lp;
	acc1 = DETECTOR::AcceptanceSoLID(ep, "e+") * DETECTOR::AcceptanceSoLID(em, "e-") * DETECTOR::AcceptanceSoLID(p, "p");
	acc2 = DETECTOR::AcceptanceSoLID(ep, "e+") * DETECTOR::AcceptanceSoLID(lp, "e-") * DETECTOR::AcceptanceSoLID(p, "p");
	hMass->Fill( (ep+em).M(), weight * acc1);
	hMassFalse->Fill( (lp+ep).M(), weight * acc2);
      }
      infile.close();
    }
  }

  if (opt == 3){//smeared
    for (int j = 0; j < Nfiles; j++){
      name = loadfile + Form("%.4d.dat",j);
      ifstream infile(name.Data());
      infile >> tmp >> nsim;
      Nsim += nsim;
      while (infile >> weight){
	infile >> tmp >> px >> py >> pz >> E;
	lp.SetXYZM(px, py, pz, PARTICLE::e.M());
	infile >> tmp >> px >> py >> pz >> E;
	ep.SetXYZM(px, py, pz, PARTICLE::e.M());
	infile >> tmp >> px >> py >> pz >> E;
	em.SetXYZM(px, py, pz, PARTICLE::e.M());
	infile >> tmp >> px >> py >> pz >> E;
	p.SetXYZM(px, py, pz, Mp);
	q = l - lp;
	acc1 = DETECTOR::SmearSoLID(ep, "e+") * DETECTOR::SmearSoLID(em, "e-") * DETECTOR::SmearSoLID(p, "p");
	acc2 = DETECTOR::SmearSoLID(ep, "e+") * DETECTOR::SmearSoLID(lp, "e-") * DETECTOR::SmearSoLID(p, "p");
	hMass->Fill( (ep+em).M(), weight * acc1);
	hMassFalse->Fill( (lp+ep).M(), weight * acc2);
      }
      infile.close();
    }
  }

  hMass->Scale(lumi*time/Nsim);
  hMassFalse->Scale(lumi*time/Nsim);
 
  fs->Write();

  return 0;
}
