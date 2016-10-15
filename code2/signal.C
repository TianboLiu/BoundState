#include "Lsimulation.h"


using namespace std;

double Kfactor(const TLorentzVector * k){
  double theta = k->Theta() / M_PI * 180.0;//get polar angle
  if (theta < 5.0 || theta > 125.0) return 0;
  double beta = k->Beta();
  double df = exp(- 2.0 / beta / 3.712);//kaon decay
  return df * 0.8;
}

double pfactor(const TLorentzVector * p){
  double theta = p->Theta() / M_PI * 180.0;//get polar angle
  if (theta < 5.0 || theta > 125.0) return 0;
  return 1.0 * 0.8;
}

int main(int argc, char * argv[]){
  double Egamma = atof(argv[1]);
  gRandom->SetSeed(0);
  SetFunctions();
  const Long64_t Nsim = 2000000;

  double res[3] = {0.01, 0.001, 0.004};//resolution dp/p, dth, dphi
  
  double lumi = 1.0e35 * pow(1.97e-14, 2);//luminosity limit in GeV^2 / s
  double A = 12.0;
  double fL = 0.06;
  double gpar[2] = {fL, A};
  double gk[3] = {11.0, 1.2, 1.6};
  double fgamma = BremsstrahlungPhotonLumi(gk, gpar);
  double time = 3600;
  double tonb = 3.89379e5;

  TLorentzVector q;
  TLorentzVector ki[2], kf[4];
  double weight[5];

  TH1D * h0 = new TH1D("hh0", "", 1000, 1.8, 2.5);
  TH1D * h1 = new TH1D("hh1", "", 1000, 1.8, 2.5);
  TH1D * h2 = new TH1D("hh2", "", 1000, 1.8, 2.5);
  TH1D * h3 = new TH1D("hh3", "", 1000, 1.8, 2.5);
  TH1D * h4 = new TH1D("hh4", "", 1000, 1.8, 2.5);
  TH1D * h5 = new TH1D("hh5", "", 1000, 1.8, 2.5);

  TH1D * a0 = new TH1D("a0", "", 1000, 1.8, 2.5);
  TH1D * a1 = new TH1D("a1", "", 1000, 1.8, 2.5);
  TH1D * a2 = new TH1D("a2", "", 1000, 1.8, 2.5);
  TH1D * a3 = new TH1D("a3", "", 1000, 1.8, 2.5);
  TH1D * a4 = new TH1D("a4", "", 1000, 1.8, 2.5);
  TH1D * a5 = new TH1D("a5", "", 1000, 1.8, 2.5);
  
  TH1D * ht = new TH1D("ht", "", 1, 0.0, 1.0);
  for (int i = 0; i < 1000000; i++){
    q.SetXYZT(0.0, 0.0, Egamma, Egamma);
    ki[0] = q;
    GenerateEventNKKwithPhi(ki, kf, weight);
    ht->Fill(0.5, weight[0]);
  }
  cout << ht->Integral(1, -1) / 1000000 * 0.389379e3 << endl;

    
  return 0;

  double af1, af2;

  TLorentzVector l0, l1, l2;
  for (int i = 0; i < Nsim; i++){
    if (i%100000==0) cout << i << endl;
    GenerateBremsstrahlungPhoton(&q, &gk[1]);//Generate a photon
    q.SetXYZT(0.0, 0.0, Egamma, Egamma);
    ki[0] = q;//Set initial photon
    if (q.E() < 2.0){
      GenerateEventNNKKwithBoundState(ki, kf, weight);//Generate a bound state event kf: N (N, K+, K-)
      if (weight[0] > 0.0){
	af1 = pfactor(&kf[1]) * Kfactor(&kf[2]) * Kfactor(&kf[3]);
	af2 = pfactor(&kf[0]) * Kfactor(&kf[2]) * Kfactor(&kf[3]);
	DetectorResolutionSmear(&kf[0], res);//smear final particles
	DetectorResolutionSmear(&kf[1], res);
	DetectorResolutionSmear(&kf[2], res);
	DetectorResolutionSmear(&kf[3], res);
	l0 = kf[1] + kf[2] + kf[3];//pKK from bound
	l1 = kf[0] + kf[2] + kf[3];//KK + recoil p
	h0->Fill(l0.M(), weight[0] * A * (A - 1) / 2.0);
	h1->Fill(l1.M(), weight[0] * A * (A - 1) / 2.0);
	a0->Fill(l0.M(), weight[0] * A * (A - 1) / 2.0 * af1);
	a1->Fill(l1.M(), weight[0] * A * (A - 1) / 2.0 * af2);
      }

      GenerateEventNKKwithPhi(ki, kf, weight);//Generate a unbound event with phi production
      if (weight[0] > 0){
	af1 = pfactor(&kf[0]) * Kfactor(&kf[1]) * Kfactor(&kf[2]);
	DetectorResolutionSmear(&kf[0], res);//smear final particles
	DetectorResolutionSmear(&kf[1], res);
	DetectorResolutionSmear(&kf[2], res);
	l2 = kf[0] + kf[1] + kf[2];//pKK unbound
	h2->Fill(l2.M(), weight[0] * A / 2.0);
	a2->Fill(l2.M(), weight[0] * A / 2.0 * af1);
      }

      GenerateEventNKKwithLambda1520(ki, kf, weight);//Generate a unbound event with Lambda1520 production
      if (weight[0] > 0){
	af1 = pfactor(&kf[0]) * Kfactor(&kf[1]) * Kfactor(&kf[2]);
	DetectorResolutionSmear(&kf[0], res);//smear final particles
	DetectorResolutionSmear(&kf[1], res);
	DetectorResolutionSmear(&kf[2], res);
	l2 = kf[0] + kf[1] + kf[2];//pKK unbound
	h3->Fill(l2.M(), weight[0] * A / 2.0);
	a3->Fill(l2.M(), weight[0] * A / 2.0 * af1);
      }
      
      if (weight[0] > 0){
	GenerateEventNKKwithout(ki, kf, weight);//Generate a unbound event with direct KK production
	af1 = pfactor(&kf[0]) * Kfactor(&kf[1]) * Kfactor(&kf[2]);
	DetectorResolutionSmear(&kf[0], res);//smear final particles
	DetectorResolutionSmear(&kf[1], res);
	DetectorResolutionSmear(&kf[2], res);
	l2 = kf[0] + kf[1] + kf[2];//pKK unbound
	h4->Fill(l2.M(), weight[0] * A / 2.0);
	a4->Fill(l2.M(), weight[0] * A / 2.0 * af1);
      }
    }
  }

  h5->Add(h0); h5->Add(h1); h5->Add(h2); h5->Add(h3); h5->Add(h4);

  //cout << "cross section" << endl;
  cout << "cross section: " << Egamma << "  " << h0->Integral(1, -1) / Nsim * tonb << " nb" << endl;
  if (false){
    cout << "number / h" << endl;
    cout << "Signal:     " << h0->Integral(1, -1) / Nsim * lumi * fgamma * time << endl;
    cout << "accepted number / h" << endl;
    cout << "Signal:     " << a0->Integral(1, -1) / Nsim * lumi * fgamma * time << endl;
  }

  return 0;
}
