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
  TH2D * h2 = new TH2D("Eres", ";E_{#gamma};#deltaE_{#gamma}", 130, 7.2, 8.5, 100, -0.5, 0.5);
  TH1D * ha = new TH1D("7.2-7.5", "", 500, -0.5, 0.5);
  TH1D * hb = new TH1D("7.5-7.7", "", 500, -0.5, 0.5);
  TH1D * hc = new TH1D("7.7-7.9", "", 500, -0.5, 0.5);
  TH1D * hd = new TH1D("7.9-8.1", "", 500, -0.5, 0.5);
  TH1D * he = new TH1D("8.1-8.3", "", 500, -0.5, 0.5);
  TH1D * hf = new TH1D("8.3-8.5", "", 500, -0.5, 0.5);
  

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
    if (Eg0 > 7.2 && Eg0 < 7.5)
      ha->Fill(dE, weight * acc);
    else if (Eg0 > 7.5 && Eg0 < 7.7)
      hb->Fill(dE, weight * acc);
    else if (Eg0 > 7.7 && Eg0 < 7.9)
      hc->Fill(dE, weight * acc);
    else if (Eg0 > 7.9 && Eg0 < 8.1)
      hd->Fill(dE, weight * acc);
    else if (Eg0 > 8.1 && Eg0 < 8.3)
      he->Fill(dE, weight * acc);
    else if (Eg0 > 8.3 && Eg0 < 8.5)
      hf->Fill(dE, weight * acc);
  }

  h2->Scale(lumi*time/Nsim);
  ha->Scale(lumi*time/Nsim);
  hb->Scale(lumi*time/Nsim);
  hc->Scale(lumi*time/Nsim);
  hd->Scale(lumi*time/Nsim);
  he->Scale(lumi*time/Nsim);
  hf->Scale(lumi*time/Nsim);
  

  fs->Write();
  fs->Close();

  return 0;
}
    
  
