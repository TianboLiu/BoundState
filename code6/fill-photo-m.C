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

  DETECTOR::SetDetector("SoLID");

  TString loadfile = argv[1];
  TString name;
  char tmp[200];
  double weight, acc;
  TLorentzVector l;
  l.SetXYZM(0, 0, 8.2, PARTICLE::e.M());
  TLorentzVector q, ep, em, p;
  double px, py, pz, E;
  TFile * fs = new TFile(argv[3], "RECREATE");
  TH2D * hEgPJpsi = new TH2D("EgPJpsi", ";E_{#gamma} (GeV);P_{J/#psi} (GeV)", 140, 7.2, 8.6, 50, 5.0, 10.0);
  TH2D * hEgThetaJpsi = new TH2D("EgThetaJpsi", ";E_{#gamma} (GeV);#theta_{J/#psi} (deg)", 140, 7.2, 8.6, 14, 0.0, 14.0);
  TH2D * hEgkappa = new TH2D("Egkappa", ";E_{#gamma} (GeV);#kappa_{J/#psi} (GeV)", 140, 7.2, 8.6, 100, 0.0, 1.0);
  TH2D * hEgPp = new TH2D("EgPp", ";E_{#gamma} (GeV);P_{p*} (GeV)", 140, 7.2, 8.6, 100, 0.0, 1.0); 
  hEgPJpsi->SetDirectory(fs);
  hEgThetaJpsi->SetDirectory(fs);
  hEgkappa->SetDirectory(fs);
  hEgPp->SetDirectory(fs);
  
  double Nsim = 0.0;
  double nsim;
  
  if (opt == 3){//smeared
    for (int j = 0; j < Nfiles; j++){
      name = loadfile + Form("%.4d.dat",j);
      ifstream infile(name.Data());
      infile >> tmp >> nsim;
      Nsim += nsim;
      while (infile >> weight){
	infile >> tmp >> px >> py >> pz >> E;
	q.SetXYZT(px, py, pz, E);
	infile >> tmp >> px >> py >> pz >> E;
	ep.SetXYZM(px, py, pz, PARTICLE::e.M());
	infile >> tmp >> px >> py >> pz >> E;
	em.SetXYZM(px, py, pz, PARTICLE::e.M());
	infile >> tmp >> px >> py >> pz >> E;
	p.SetXYZM(px, py, pz, Mp);
	acc = DETECTOR::SmearSoLID(ep, "e+") * DETECTOR::SmearSoLID(em, "e-") * DETECTOR::SmearSoLID(p, "p");
	if (acc > 0){
	  hEgPJpsi->Fill( CalcEg(ep+em+p), (ep+em).P(), weight * acc);
	  hEgThetaJpsi->Fill( CalcEg(ep+em+p), (ep+em).Theta()*deg, weight * acc);
	  hEgkappa->Fill( CalcEg(ep+em+p), Getkappa(ep+em+p), weight * acc);
	  hEgPp->Fill( CalcEg(ep+em+p), (ep+em+p-q).P(), weight * acc);
	}
      }
      infile.close();
    }
  }

  hEgPJpsi->Scale(lumi*time/Nsim);
  hEgThetaJpsi->Scale(lumi*time/Nsim);
  hEgkappa->Scale(lumi*time/Nsim);
  hEgPp->Scale(lumi*time/Nsim);
  
  fs->Write();

  return 0;
}
