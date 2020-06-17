#include "Lcore.h"

double CalcEg(const TLorentzVector P){
  double Md = 1.8756;
  double Mn = 0.93957;
  double Eg = (pow(P.E() - Md, 2) - pow(P.P(), 2) - Mn * Mn) / (2.0 * (P.E() - P.Pz() - Md));
  return Eg;
}

int main(const int argc, const char * argv[]){

  if (argc < 1){
    return 0;
  }

  const double lumi = 0.6e37 * 1e-26 * pow(0.197326,2);
  const double time = 3600.0;
  const double Nsim = 1e6;

  DETECTOR::SetDetector("SoLID");
  
  TLorentzVector q, ep, em, p;
  double px, py, pz, E;
  double weight, acc;
  double Eg0, Eg, dE;

  TFile * fs = new TFile("gammares.root", "RECREATE");
  TH2D * h2 = new TH2D("Eres", ";E_{#gamma};#deltaE_{#gamma}", 100, 7.2, 8.2, 100, -0.5, 0.5);
  TH1D * ha = new TH1D("a", "", 100, -0.2, 0.2);
  TH1D * hb = new TH1D("b", "", 100, -0.2, 0.2);
  TH1D * hc = new TH1D("c", "", 100, -0.2, 0.2);
  TH1D * hd = new TH1D("d", "", 100, -0.2, 0.2);
  TH1D * he = new TH1D("e", "", 100, -0.2, 0.2);
  

  char tmp[200];
  ifstream infile("testres0000.dat");
  infile.getline(tmp, 200);

  while (infile >> weight){
    infile >> tmp >> px >> py >> pz >> E;
    q.SetXYZT(px, py, pz, E);
    infile >> tmp >> px >> py >> pz >> E;
    ep.SetXYZM(px, py, pz, PARTICLE::e.M());
    infile >> tmp >> px >> py >> pz >> E;
    em.SetXYZM(px, py, pz, PARTICLE::e.M());
    infile >> tmp >> px >> py >> pz >> E;
    p.SetXYZM(px, py, pz, Mp);
    Eg0 = CalcEg(ep+em+p);
    acc = DETECTOR::SmearSoLID(ep, "e+") * DETECTOR::SmearSoLID(em, "e-") * DETECTOR::SmearSoLID(p, "p");
    Eg = CalcEg(ep+em+p);
    dE = Eg - Eg0;
    h2->Fill(Eg0, dE, weight * acc);
    if (Eg0 > 7.2 && Eg0 < 7.4)
      ha->Fill(dE, weight * acc);
    else if (Eg0 > 7.4 && Eg0 < 7.6)
      hb->Fill(dE, weight * acc);
    else if (Eg0 > 7.6 && Eg0 < 7.8)
      hc->Fill(dE, weight * acc);
    else if (Eg0 > 7.8 && Eg0 < 8.0)
      hd->Fill(dE, weight * acc);
    else if (Eg0 > 8.0 && Eg0 < 8.2)
      he->Fill(dE, weight * acc);
  }

  h2->Scale(lumi*time/Nsim);
  ha->Scale(lumi*time/Nsim);
  hb->Scale(lumi*time/Nsim);
  hc->Scale(lumi*time/Nsim);
  hd->Scale(lumi*time/Nsim);
  he->Scale(lumi*time/Nsim);
  

  fs->Write();
  fs->Close();

  return 0;
}
    
  
