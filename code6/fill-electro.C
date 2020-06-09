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
  TLorentzVector l, lp;
  l.SetXYZM(0, 0, 8.2, PARTICLE::e.M());
  TLorentzVector q, ep, em, p;
  double px, py, pz, E;
  TFile * fs = new TFile(argv[3], "RECREATE");
  TH1D * hMass = new TH1D("Mass_e+e-", ";M[e^{+}e^{-}] (GeV);Events / hour", 100, 2.6, 3.6);
  TH2D * hep = new TH2D("ThetaP_e+", ";e^{+} polar angle (deg);e^{+} momentum (GeV)", 100, 0.0, 40.0, 100, 0.0, 8.0);
  TH2D * hem = new TH2D("ThetaP_e-", ";e^{-} polar angle (deg);e^{-} momentum (GeV)", 100, 0.0, 40.0, 100, 0.0, 8.0);
  TH2D * hp = new TH2D("ThetaP_p", ";proton polar angle (deg);proton momentum (GeV)", 100, 0.0, 40.0, 100, 0.0, 8.0);
  TH2D * hJpsi  = new TH2D("ThetaP_Jpsi", ";J/#psi polar angle (deg);J/#psi moemntum (GeV)", 100, 0.0, 15.0, 100, 0.0, 10.0);
  TH2D * hAngle = new TH2D("Angle", ";e^{+}e^{-} angle (deg);e^{+} momentum (GeV)", 100, 0.0, 180.0, 100, 0.0, 8.0);
  TH1D * hkappa = new TH1D("kappa", ";#kappa_{J/#psi} (GeV);Events / hour", 100, 0.0, 1.0);
  TH2D * hkappaP = new TH2D("kappaP", ";#kappa_{J/#psi} (GeV);proton momentum (GeV)", 100, 0.0, 1.0, 50, 0.0, 5.0);
  TH2D * hkappaTheta = new TH2D("kappaTheta", ";#kappa_{J/#psi} (GeV);proton polar angle (deg)", 100, 0.0, 1.0, 60, 0.0, 30.0);
  TH1D * hRealEg = new TH1D("RealEg", ";E_{#gamma} (GeV);Events / hour", 140, 7.0, 8.4);
  TH2D * hkappaEg  = new TH2D("kappaEg", ";#kappa_{J/#psi} (GeV);E_{#gamma} (GeV)", 100, 0.0, 1.0, 140, 7.0, 8.4);
  TH2D * hkappaEgLow = new TH2D("kappaEgLow", ";#kappa_{J/#psi} (GeV);E_{#gamma} (GeV)", 100, 0.0, 1.0, 140, 7.0, 8.4);
  TH2D * hkappaEgHigh = new TH2D("kappaEgHigh", ";#kappa_{J/#psi} (GeV);E_{#gamma} (GeV)", 100, 0.0, 1.0, 140, 7.0, 8.4);
  hMass->SetDirectory(fs);
  hep->SetDirectory(fs);
  hem->SetDirectory(fs);
  hp->SetDirectory(fs);
  hJpsi->SetDirectory(fs);
  hAngle->SetDirectory(fs);
  hkappa->SetDirectory(fs);
  hkappaP->SetDirectory(fs);
  hkappaTheta->SetDirectory(fs);
  hRealEg->SetDirectory(fs);
  hkappaEg->SetDirectory(fs);
  hkappaEgLow->SetDirectory(fs);
  hkappaEgHigh->SetDirectory(fs);
  
  double Nsim = 0.0;
  double nsim;
  
  if (opt == 1){//raw
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
	acc = 1.0;
	hMass->Fill( (ep+em).M(), weight * acc);
	hep->Fill( ep.Theta() * deg, ep.P(), weight * acc);
	hem->Fill( em.Theta() * deg, em.P(), weight * acc);
	hp->Fill( p.Theta() * deg, p.P(), weight * acc);
	hJpsi->Fill( (ep+em).Theta() * deg, (ep+em).P(), weight * acc);
	hAngle->Fill( ep.Angle(em.Vect()) * deg, ep.P(), weight * acc);
	hkappa->Fill( Getkappa(ep+em+p), weight * acc);
	hkappaP->Fill( Getkappa(ep+em+p), p.P(), weight * acc);
	hkappaTheta->Fill( Getkappa(ep+em+p), p.Theta() * deg, weight * acc);
	hRealEg->Fill( q.E(), weight * acc);
	hkappaEg->Fill( Getkappa(ep+em+p), CalcEg(ep+em+p), weight * acc);
	if (p.P() < 1.0) hkappaEgLow->Fill( Getkappa(ep+em+p), CalcEg(ep+em+p), weight * acc);
	if (p.P() > 1.0) hkappaEgHigh->Fill( Getkappa(ep+em+p), CalcEg(ep+em+p), weight * acc);
      }
      infile.close();
    }
    hMass->Scale(lumi*time/Nsim);
    hep->Scale(lumi*time/Nsim);
    hem->Scale(lumi*time/Nsim);
    hp->Scale(lumi*time/Nsim);
    hJpsi->Scale(lumi*time/Nsim);
    hAngle->Scale(lumi*time/Nsim);
    hkappa->Scale(lumi*time/Nsim);
    hkappaP->Scale(lumi*time/Nsim);
    hkappaTheta->Scale(lumi*time/Nsim);
    hRealEg->Scale(lumi*time/Nsim);
    hkappaEg->Scale(lumi*time/Nsim);
    hkappaEgLow->Scale(lumi*time/Nsim);
    hkappaEgHigh->Scale(lumi*time/Nsim);
  }

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
	acc = DETECTOR::AcceptanceSoLID(ep, "e+") * DETECTOR::AcceptanceSoLID(em, "e-") * DETECTOR::AcceptanceSoLID(p, "p");
	hMass->Fill( (ep+em).M(), weight * acc);
	hep->Fill( ep.Theta() * deg, ep.P(), weight * acc);
	hem->Fill( em.Theta() * deg, em.P(), weight * acc);
	hp->Fill( p.Theta() * deg, p.P(), weight * acc);
	hJpsi->Fill( (ep+em).Theta() * deg, (ep+em).P(), weight * acc);
	hAngle->Fill( ep.Angle(em.Vect()) * deg, ep.P(), weight * acc);
	hkappa->Fill( Getkappa(ep+em+p), weight * acc);
	hkappaP->Fill( Getkappa(ep+em+p), p.P(), weight * acc);
	hkappaTheta->Fill( Getkappa(ep+em+p), p.Theta() * deg, weight * acc);
	hRealEg->Fill( q.E(), weight * acc);
	hkappaEg->Fill( Getkappa(ep+em+p), CalcEg(ep+em+p), weight * acc);
	if (p.P() < 1.0) hkappaEgLow->Fill( Getkappa(ep+em+p), CalcEg(ep+em+p), weight * acc);
	if (p.P() > 1.0) hkappaEgHigh->Fill( Getkappa(ep+em+p), CalcEg(ep+em+p), weight * acc);
      }
      infile.close();
    }
    hMass->Scale(lumi*time/Nsim);
    hep->Scale(lumi*time/Nsim);
    hem->Scale(lumi*time/Nsim);
    hp->Scale(lumi*time/Nsim);
    hJpsi->Scale(lumi*time/Nsim);
    hAngle->Scale(lumi*time/Nsim);
    hkappa->Scale(lumi*time/Nsim);
    hkappaP->Scale(lumi*time/Nsim);
    hkappaTheta->Scale(lumi*time/Nsim);
    hRealEg->Scale(lumi*time/Nsim);
    hkappaEg->Scale(lumi*time/Nsim);
    hkappaEgLow->Scale(lumi*time/Nsim);
    hkappaEgHigh->Scale(lumi*time/Nsim);
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
	acc = DETECTOR::SmearSoLID(ep, "e+") * DETECTOR::SmearSoLID(em, "e-") * DETECTOR::SmearSoLID(p, "p");
	hMass->Fill( (ep+em).M(), weight * acc);
	hep->Fill( ep.Theta() * deg, ep.P(), weight * acc);
	hem->Fill( em.Theta() * deg, em.P(), weight * acc);
	hp->Fill( p.Theta() * deg, p.P(), weight * acc);
	hJpsi->Fill( (ep+em).Theta() * deg, (ep+em).P(), weight * acc);
	hAngle->Fill( ep.Angle(em.Vect()) * deg, ep.P(), weight * acc);
	hkappa->Fill( Getkappa(ep+em+p), weight * acc);
	hkappaP->Fill( Getkappa(ep+em+p), p.P(), weight * acc);
	hkappaTheta->Fill( Getkappa(ep+em+p), p.Theta() * deg, weight * acc);
	hRealEg->Fill( q.E(), weight * acc);
	hkappaEg->Fill( Getkappa(ep+em+p), CalcEg(ep+em+p), weight * acc);
	if (p.P() < 1.0) hkappaEgLow->Fill( Getkappa(ep+em+p), CalcEg(ep+em+p), weight * acc);
	if (p.P() > 1.0) hkappaEgHigh->Fill( Getkappa(ep+em+p), CalcEg(ep+em+p), weight * acc);
      }
      infile.close();
    }
    hMass->Scale(lumi*time/Nsim);
    hep->Scale(lumi*time/Nsim);
    hem->Scale(lumi*time/Nsim);
    hp->Scale(lumi*time/Nsim);
    hJpsi->Scale(lumi*time/Nsim);
    hAngle->Scale(lumi*time/Nsim);
    hkappa->Scale(lumi*time/Nsim);
    hkappaP->Scale(lumi*time/Nsim);
    hkappaTheta->Scale(lumi*time/Nsim);
    hRealEg->Scale(lumi*time/Nsim);
    hkappaEg->Scale(lumi*time/Nsim);
    hkappaEgLow->Scale(lumi*time/Nsim);
    hkappaEgHigh->Scale(lumi*time/Nsim);
  }
 
  fs->Write();

  return 0;
}
