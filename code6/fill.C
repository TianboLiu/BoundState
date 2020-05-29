#include "Lcore.h"

double Getkappa(const double M){
  double M1 = 3.0969;
  double M2 = 0.938272;
  return sqrt((M*M - pow(M1 + M2, 2)) * (M*M - pow(M1 - M2, 2))) / (2.0 * M);
}
  

int main(const int argc, const char * argv[]){

  if (argc < 5){
    cout << "./fill <loadfile> <Nfiles> <savefile> <opt>" << endl;
    cout << "1: raw" << endl;
    cout << "2: 3fold detection" << endl;
    cout << "3: 3fold smeared" << endl;
    cout << "4: 4fold smeared" << endl;
    return 0;
  }

  const int Nfiles = atoi(argv[2]);
  const double Nsim = 1e5 * Nfiles;
  const int opt = atoi(argv[4]);

  const double lumi = 1.2e37 * 1e-26 * pow(0.197326,2);
  const double time = 3600.0;
  const double deg = 180.0 / M_PI;

  DETECTOR::SetDetector("SoLID");

  TString loadfile = argv[1];
  TString name;
  char tmp[200];
  double weight, acc;
  TLorentzVector l, ep, em, p;
  double px, py, pz;
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

  if (opt == 1){//raw
    for (int j = 0; j < Nfiles; j++){
      name = loadfile + Form("%.4d.dat",j);
      ifstream infile(name.Data());
      infile.getline(tmp, 200);
      while (infile >> weight){
	infile >> tmp >> px >> py >> pz;
	l.SetXYZM(px, py, pz, PARTICLE::e.M());
	infile >> tmp >> px >> py >> pz;
	ep.SetXYZM(px, py, pz, PARTICLE::e.M());
	infile >> tmp >> px >> py >> pz;
	em.SetXYZM(px, py, pz, PARTICLE::e.M());
	infile >> tmp >> px >> py >> pz;
	p.SetXYZM(px, py, pz, Mp);
	acc = 1.0;
	hMass->Fill( (ep+em).M(), weight * acc);
	hep->Fill( ep.Theta() * deg, ep.P(), weight * acc);
	hem->Fill( em.Theta() * deg, em.P(), weight * acc);
	hp->Fill( p.Theta() * deg, p.P(), weight * acc);
	hJpsi->Fill( (ep+em).Theta() * deg, (ep+em).P(), weight * acc);
	hAngle->Fill( ep.Angle(em.Vect()) * deg, ep.P(), weight * acc);
	hkappa->Fill( Getkappa((ep+em+p).M()), weight * acc);
	hkappaP->Fill( Getkappa((ep+em+p).M()), p.P(), weight * acc);
	hkappaTheta->Fill( Getkappa((ep+em+p).M()), p.Theta() * deg, weight * acc);
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
  }

  if (opt == 2){//e+e-p detected
    for (int j = 0; j < Nfiles; j++){
      name = loadfile + Form("%.4d.dat",j);
      ifstream infile(name.Data());
      infile.getline(tmp, 200);
      while (infile >> weight){
	infile >> tmp >> px >> py >> pz;
	l.SetXYZM(px, py, pz, PARTICLE::e.M());
	infile >> tmp >> px >> py >> pz;
	ep.SetXYZM(px, py, pz, PARTICLE::e.M());
	infile >> tmp >> px >> py >> pz;
	em.SetXYZM(px, py, pz, PARTICLE::e.M());
	infile >> tmp >> px >> py >> pz;
	p.SetXYZM(px, py, pz, Mp);
	acc = DETECTOR::AcceptanceSoLID(ep, "e+") * DETECTOR::AcceptanceSoLID(em, "e-") * DETECTOR::AcceptanceSoLID(p, "p");
	hMass->Fill( (ep+em).M(), weight * acc);
	hep->Fill( ep.Theta() * deg, ep.P(), weight * acc);
	hem->Fill( em.Theta() * deg, em.P(), weight * acc);
	hp->Fill( p.Theta() * deg, p.P(), weight * acc);
	hJpsi->Fill( (ep+em).Theta() * deg, (ep+em).P(), weight * acc);
	hAngle->Fill( ep.Angle(em.Vect()) * deg, ep.P(), weight * acc);
	hkappa->Fill( Getkappa((ep+em+p).M()), weight * acc);
	hkappaP->Fill( Getkappa((ep+em+p).M()), p.P(), weight * acc);
	hkappaTheta->Fill( Getkappa((ep+em+p).M()), p.Theta() * deg, weight * acc);
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
  }

  if (opt == 3){//e+e-p detected and smeared
    for (int j = 0; j < Nfiles; j++){
      name = loadfile + Form("%.4d.dat",j);
      ifstream infile(name.Data());
      infile.getline(tmp, 200);
      while (infile >> weight){
	infile >> tmp >> px >> py >> pz;
	l.SetXYZM(px, py, pz, PARTICLE::e.M());
	infile >> tmp >> px >> py >> pz;
	ep.SetXYZM(px, py, pz, PARTICLE::e.M());
	infile >> tmp >> px >> py >> pz;
	em.SetXYZM(px, py, pz, PARTICLE::e.M());
	infile >> tmp >> px >> py >> pz;
	p.SetXYZM(px, py, pz, Mp);
	acc = DETECTOR::SmearSoLID(ep, "e+") * DETECTOR::SmearSoLID(em, "e-") * DETECTOR::SmearSoLID(p, "p");
	hMass->Fill( (ep+em).M(), weight * acc);
	hep->Fill( ep.Theta() * deg, ep.P(), weight * acc);
	hem->Fill( em.Theta() * deg, em.P(), weight * acc);
	hp->Fill( p.Theta() * deg, p.P(), weight * acc);
	hJpsi->Fill( (ep+em).Theta() * deg, (ep+em).P(), weight * acc);
	hAngle->Fill( ep.Angle(em.Vect()) * deg, ep.P(), weight * acc);
	hkappa->Fill( Getkappa((ep+em+p).M()), weight * acc);
	hkappaP->Fill( Getkappa((ep+em+p).M()), p.P(), weight * acc);
	hkappaTheta->Fill( Getkappa((ep+em+p).M()), p.Theta() * deg, weight * acc);
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
  }

  if (opt == 4){//scatterred e and e+e-p
    for (int j = 0; j < Nfiles; j++){
      name = loadfile + Form("%.4d.dat",j);
      ifstream infile(name.Data());
      infile.getline(tmp, 200);
      while (infile >> weight){
	infile >> tmp >> px >> py >> pz;
	l.SetXYZM(px, py, pz, PARTICLE::e.M());
	infile >> tmp >> px >> py >> pz;
	ep.SetXYZM(px, py, pz, PARTICLE::e.M());
	infile >> tmp >> px >> py >> pz;
	em.SetXYZM(px, py, pz, PARTICLE::e.M());
	infile >> tmp >> px >> py >> pz;
	p.SetXYZM(px, py, pz, Mp);
	acc = DETECTOR::SmearSoLID(l, "e-") * DETECTOR::SmearSoLID(ep, "e+") * DETECTOR::SmearSoLID(em, "e-") * DETECTOR::SmearSoLID(p, "p");
	hMass->Fill( (ep+em).M(), weight * acc);
	hep->Fill( ep.Theta() * deg, ep.P(), weight * acc);
	hem->Fill( em.Theta() * deg, em.P(), weight * acc);
	hp->Fill( p.Theta() * deg, p.P(), weight * acc);
	hJpsi->Fill( (ep+em).Theta() * deg, (ep+em).P(), weight * acc);
	hAngle->Fill( ep.Angle(em.Vect()) * deg, ep.P(), weight * acc);
	hkappa->Fill( Getkappa((ep+em+p).M()), weight * acc);
	hkappaP->Fill( Getkappa((ep+em+p).M()), p.P(), weight * acc);
	hkappaTheta->Fill( Getkappa((ep+em+p).M()), p.Theta() * deg, weight * acc);
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
  }

  if (opt == 5){//only e+e
    for (int j = 0; j < Nfiles; j++){
      name = loadfile + Form("%.4d.dat",j);
      ifstream infile(name.Data());
      infile.getline(tmp, 200);
      while (infile >> weight){
        infile >> tmp >> px >> py >> pz;
        l.SetXYZM(px, py, pz, PARTICLE::e.M());
        infile >> tmp >> px >> py >> pz;
        ep.SetXYZM(px, py, pz, PARTICLE::e.M());
        infile >> tmp >> px >> py >> pz;
        em.SetXYZM(px, py, pz, PARTICLE::e.M());
        infile >> tmp >> px >> py >> pz;
        p.SetXYZM(px, py, pz, Mp);
        acc = DETECTOR::SmearSoLID(ep, "e+") * DETECTOR::SmearSoLID(em, "e-");
        hMass->Fill( (ep+em).M(), weight * acc);
        hep->Fill( ep.Theta() * deg, ep.P(), weight * acc);
        hem->Fill( em.Theta() * deg, em.P(), weight * acc);
        hp->Fill( p.Theta() * deg, p.P(), weight * acc);
        hJpsi->Fill( (ep+em).Theta() * deg, (ep+em).P(), weight * acc);
        hAngle->Fill( ep.Angle(em.Vect()) * deg, ep.P(), weight * acc);
        hkappa->Fill( Getkappa((ep+em+p).M()), weight * acc);
	hkappaP->Fill( Getkappa((ep+em+p).M()), p.P(), weight * acc);
	hkappaTheta->Fill( Getkappa((ep+em+p).M()), p.Theta() * deg, weight * acc);
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
  }

  if (opt == 6){//scatterred e and e+e
    for (int j = 0; j < Nfiles; j++){
      name = loadfile + Form("%.4d.dat",j);
      ifstream infile(name.Data());
      infile.getline(tmp, 200);
      while (infile >> weight){
        infile >> tmp >> px >> py >> pz;
        l.SetXYZM(px, py, pz, PARTICLE::e.M());
        infile >> tmp >> px >> py >> pz;
        ep.SetXYZM(px, py, pz, PARTICLE::e.M());
        infile >> tmp >> px >> py >> pz;
        em.SetXYZM(px, py, pz, PARTICLE::e.M());
        infile >> tmp >> px >> py >> pz;
        p.SetXYZM(px, py, pz, Mp);
        acc = DETECTOR::SmearSoLID(l, "e-") * DETECTOR::SmearSoLID(ep, "e+") * DETECTOR::SmearSoLID(em, "e-");
        hMass->Fill( (ep+em).M(), weight * acc);
        hep->Fill( ep.Theta() * deg, ep.P(), weight * acc);
        hem->Fill( em.Theta() * deg, em.P(), weight * acc);
        hp->Fill( p.Theta() * deg, p.P(), weight * acc);
        hJpsi->Fill( (ep+em).Theta() * deg, (ep+em).P(), weight * acc);
        hAngle->Fill( ep.Angle(em.Vect()) * deg, ep.P(), weight * acc);
        hkappa->Fill( Getkappa((ep+em+p).M()), weight * acc);
	hkappaP->Fill( Getkappa((ep+em+p).M()), p.P(), weight * acc);
	hkappaTheta->Fill( Getkappa((ep+em+p).M()), p.Theta() * deg, weight * acc);
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
  }
 
  fs->Write();

  return 0;
}
