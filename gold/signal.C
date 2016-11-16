#include "Lgold.h"

using namespace std;

int main(int argc, char * argv[]){

  double Ebeam = atof(argv[1]);//Get electron beam energy
  cout << "Ebeam: " << Ebeam << endl;

  gRandom->SetSeed(1);
  SetFunctions();
  SetStoppingPower();

  SetTagger(0.5, Ebeam);
  //SetTagger(0.5, 4.5);
  
  Long64_t Nsim = 100000000;

  double lumi = 1.0e35 * 1.0e-26 * pow(0.197327, 2) / 197.0;//eA GeV^2 s^-1
  double time = 3600.0;

  double res_clas12fa[3] = {0.01, 0.001, 0.004};//resolution dp/p, dth, dphi
  double res_clas12la[3] = {0.05, 0.020, 0.005};//resolution dp/p, dth, dphi
  double res_bonus[3] = {0.1, 0.020, 0.020};//resolution dp/p, dth, dphi

  TLorentzVector ki[2], kf[6];
  double weight;

  ki[0].SetXYZT(0, 0, Ebeam, Ebeam);//Set electron beam 4-momentum

  TLorentzVector P0_p, P0_Kp, P0_Km, P0_all;
  TLorentzVector P_p, P_Kp, P_Km, P_all;

  double fbonus, fclas12fa, fclas12la;
  double factorp, factorKp, factorKm, ftotal;

  double z;
  double hmax;

  TH1D * hB = new TH1D("hB", "", 1000, 1.88, 2.32);
  hB->SetTitle("M(pK^{+}K^{-}) spectra");
  hB->GetXaxis()->SetTitle("M(pK^{+}K^{-}) / GeV");
  hB->GetXaxis()->CenterTitle();
  hB->GetXaxis()->SetTitleOffset(1.0);
  hB->GetXaxis()->SetTitleSize(0.05);
  hB->GetYaxis()->SetRangeUser(0.0, 1.0);
  TH1D * h0 = new TH1D("h0", "", 1000, 1.88, 2.32);
  TH1D * hh0 = new TH1D("hh0", "", 1000, 1.88, 2.32);
  TH1D * h1 = new TH1D("h1", "", 1000, 1.88, 2.32);
  TH1D * h2 = new TH1D("h2", "", 1000, 1.88, 2.32);
  TH1D * h3 = new TH1D("h3", "", 1000, 1.88, 2.32);
  TH1D * h4 = new TH1D("h4", "", 1000, 1.88, 2.32);
  h0->SetLineColor(1);
  hh0->SetLineColor(4);
  h1->SetLineColor(2);
  h2->SetLineColor(6);
  h3->SetLineColor(7);
  h4->SetLineColor(5);

  TH1D * hBs = new TH1D("hBs", "", 1000, 1.88, 2.32);
  hBs->SetTitle("M(pK^{+}K^{-}) spectra");
  hBs->GetXaxis()->SetTitle("M(pK^{+}K^{-}) / GeV");
  hBs->GetXaxis()->CenterTitle();
  hBs->GetXaxis()->SetTitleOffset(1.0);
  hBs->GetXaxis()->SetTitleSize(0.05);
  hBs->GetYaxis()->SetRangeUser(0.0, 1.0);
  TH1D * h0s = new TH1D("h0s", "", 1000, 1.88, 2.32);
  TH1D * hh0s = new TH1D("hh0s", "", 1000, 1.88, 2.32);
  TH1D * h1s = new TH1D("h1s", "", 1000, 1.88, 2.32);
  TH1D * h2s = new TH1D("h2s", "", 1000, 1.88, 2.32);
  TH1D * h3s = new TH1D("h3s", "", 1000, 1.88, 2.32);
  TH1D * h4s = new TH1D("h4s", "", 1000, 1.88, 2.32);
  h0s->SetLineColor(1);
  hh0s->SetLineColor(4);
  h1s->SetLineColor(2);
  h2s->SetLineColor(6);
  h3s->SetLineColor(7);
  h4s->SetLineColor(5);

  TH1D * hBp = new TH1D("hBp", "", 1000, 1.88, 2.32);
  hBp->SetTitle("M(pK^{+}K^{-}) spectra");
  hBp->GetXaxis()->SetTitle("M(pK^{+}K^{-}) / GeV");
  hBp->GetXaxis()->CenterTitle();
  hBp->GetXaxis()->SetTitleOffset(1.0);
  hBp->GetXaxis()->SetTitleSize(0.05);
  hBp->GetYaxis()->SetRangeUser(0.0, 1.0);
  TH1D * h0p = new TH1D("h0p", "", 1000, 1.88, 2.32);
  TH1D * hh0p = new TH1D("hh0p", "", 1000, 1.88, 2.32);
  TH1D * h1p = new TH1D("h1p", "", 1000, 1.88, 2.32);
  TH1D * h2p = new TH1D("h2p", "", 1000, 1.88, 2.32);
  TH1D * h3p = new TH1D("h3p", "", 1000, 1.88, 2.32);
  TH1D * h4p = new TH1D("h4p", "", 1000, 1.88, 2.32);
  h0p->SetLineColor(1);
  hh0p->SetLineColor(4);
  h1p->SetLineColor(2);
  h2p->SetLineColor(6);
  h3p->SetLineColor(7);
  h4p->SetLineColor(5);

  TH1D * hBsp = new TH1D("hBsp", "", 1000, 1.88, 2.32);
  hBsp->SetTitle("M(pK^{+}K^{-}) spectra");
  hBsp->GetXaxis()->SetTitle("M(pK^{+}K^{-}) / GeV");
  hBsp->GetXaxis()->CenterTitle();
  hBsp->GetXaxis()->SetTitleOffset(1.0);
  hBsp->GetXaxis()->SetTitleSize(0.05);
  hBsp->GetYaxis()->SetRangeUser(0.0, 1.0);
  TH1D * h0sp = new TH1D("h0sp", "", 1000, 1.88, 2.32);
  TH1D * hh0sp = new TH1D("hh0sp", "", 1000, 1.88, 2.32);
  TH1D * h1sp = new TH1D("h1sp", "", 1000, 1.88, 2.32);
  TH1D * h2sp = new TH1D("h2sp", "", 1000, 1.88, 2.32);
  TH1D * h3sp = new TH1D("h3sp", "", 1000, 1.88, 2.32);
  TH1D * h4sp = new TH1D("h4sp", "", 1000, 1.88, 2.32);
  h0sp->SetLineColor(1);
  hh0sp->SetLineColor(4);
  h1sp->SetLineColor(2);
  h2sp->SetLineColor(6);
  h3sp->SetLineColor(7);
  h4sp->SetLineColor(5);
 
  for (Long64_t i = 0; i < Nsim; i++){
    //if (i%100000 == 0) cout << i << endl;

    z = gRandom->Uniform(0.0, 1.0);

    //channel 0
    weight = GenerateEvent_eNKKN_withElectroproductionBoundStateGold(ki, kf);//Generate a bound state event kf: N (N, K+, K-)
    if (weight > 0){
      P0_p = kf[4];
      P0_Kp = kf[2];
      P0_Km = kf[3];
      P_p = kf[4];
      P_Kp = kf[2];
      P_Km = kf[3];
      //proton detection
      fbonus = CheckBONUS(&P_p, "p", z);
      fclas12fa = CheckCLAS12FA(&P_p, "p", z);
      fclas12la = CheckCLAS12LA(&P_p, "p", z);
      if (fbonus > 0.0)
	DetectorResolutionSmear(&P_p, res_bonus);
      else if (fclas12fa > 0.0)
	DetectorResolutionSmear(&P_p, res_clas12fa);
      else if (fclas12la > 0.0)
	DetectorResolutionSmear(&P_p, res_clas12la);
      factorp = fbonus + fclas12fa + fclas12la;
      //K+ detection
      fbonus = CheckBONUS(&P_Kp, "K+", z);
      fclas12fa = CheckCLAS12FA(&P_Kp, "K+", z);
      fclas12la = CheckCLAS12LA(&P_Kp, "K+", z);
      if (fbonus > 0.0)
	DetectorResolutionSmear(&P_Kp, res_bonus);
      else if (fclas12fa > 0.0)
	DetectorResolutionSmear(&P_Kp, res_clas12fa);
      else if (fclas12la > 0.0)
	DetectorResolutionSmear(&P_Kp, res_clas12la);
      factorKp = fbonus + fclas12fa + fclas12la;
      //K- detection
      fbonus = CheckBONUS(&P_Km, "K-", z);
      fclas12fa = CheckCLAS12FA(&P_Km, "K-", z);
      fclas12la = CheckCLAS12LA(&P_Km, "K-", z);
      if (fbonus > 0.0)
	DetectorResolutionSmear(&P_Km, res_bonus);
      else if (fclas12fa > 0.0)
	DetectorResolutionSmear(&P_Km, res_clas12fa);
      else if (fclas12la > 0.0)
	DetectorResolutionSmear(&P_Km, res_clas12la);
      factorKm = fbonus + fclas12fa + fclas12la;

      ftotal = factorp * factorKp * factorKm;
      //ftotal = 1.0;
      if (ftotal > 0.0){
	P0_all = P0_p + P0_Kp + P0_Km;

	P_all = P_p + P_Kp + P_Km;

	h0->Fill(P0_all.M(), weight * ftotal / 2.0);

	if (MomentumCut(&P0_p, 0.0, 0.5) && MomentumCut(&P0_Kp, 0.0, 0.35) && MomentumCut(&P0_Km, 0.0, 0.35))
	  h0p->Fill(P0_all.M(), weight * ftotal / 2.0);

	h0s->Fill(P_all.M(), weight * ftotal / 2.0);

	if (MomentumCut(&P_p, 0.0, 0.5) && MomentumCut(&P_Kp, 0.0, 0.35) && MomentumCut(&P_Km, 0.0, 0.35))
	  h0sp->Fill(P_all.M(), weight * ftotal / 2.0);
      }
    }

    //channel 0'
    if (weight > 0){
      P0_p = kf[1];
      P0_Kp = kf[2];
      P0_Km = kf[3];
      P_p = kf[1];
      P_Kp = kf[2];
      P_Km = kf[3];
      //proton detection
      fbonus = CheckBONUS(&P_p, "p", z);
      fclas12fa = CheckCLAS12FA(&P_p, "p", z);
      fclas12la = CheckCLAS12LA(&P_p, "p", z);
      if (fbonus > 0.0)
	DetectorResolutionSmear(&P_p, res_bonus);
      else if (fclas12fa > 0.0)
	DetectorResolutionSmear(&P_p, res_clas12fa);
      else if (fclas12la > 0.0)
	DetectorResolutionSmear(&P_p, res_clas12la);
      factorp = fbonus + fclas12fa + fclas12la;
      //K+ detection
      fbonus = CheckBONUS(&P_Kp, "K+", z);
      fclas12fa = CheckCLAS12FA(&P_Kp, "K+", z);
      fclas12la = CheckCLAS12LA(&P_Kp, "K+", z);
      if (fbonus > 0.0)
	DetectorResolutionSmear(&P_Kp, res_bonus);
      else if (fclas12fa > 0.0)
	DetectorResolutionSmear(&P_Kp, res_clas12fa);
      else if (fclas12la > 0.0)
	DetectorResolutionSmear(&P_Kp, res_clas12la);
      factorKp = fbonus + fclas12fa + fclas12la;
      //K- detection
      fbonus = CheckBONUS(&P_Km, "K-", z);
      fclas12fa = CheckCLAS12FA(&P_Km, "K-", z);
      fclas12la = CheckCLAS12LA(&P_Km, "K-", z);
      if (fbonus > 0.0)
	DetectorResolutionSmear(&P_Km, res_bonus);
      else if (fclas12fa > 0.0)
	DetectorResolutionSmear(&P_Km, res_clas12fa);
      else if (fclas12la > 0.0)
	DetectorResolutionSmear(&P_Km, res_clas12la);
      factorKm = fbonus + fclas12fa + fclas12la;
      
      ftotal = factorp * factorKp * factorKm;
      //ftotal = 1.0;
      if (ftotal > 0.0){
	P0_all = P0_p + P0_Kp + P0_Km;
	P_all = P_p + P_Kp + P_Km;

	hh0->Fill(P0_all.M(), weight * ftotal / 2.0);

	if (MomentumCut(&P0_p, 0.0, 0.5) && MomentumCut(&P0_Kp, 0.0, 0.35) && MomentumCut(&P0_Km, 0.0, 0.35))
	  hh0p->Fill(P0_all.M(), weight * ftotal / 2.0);

	hh0s->Fill(P_all.M(), weight * ftotal / 2.0);

	if (MomentumCut(&P_p, 0.0, 0.5) && MomentumCut(&P_Kp, 0.0, 0.35) && MomentumCut(&P_Km, 0.0, 0.35))
	  hh0sp->Fill(P_all.M(), weight * ftotal / 2.0);
      }
    }

    //channel 1
    weight = GenerateEvent_eNKK_withElectroproductionPhiGold(ki, kf);
    if (weight > 0){
      P0_p = kf[1];
      P0_Kp = kf[2];
      P0_Km = kf[3];
      P_p = kf[1];
      P_Kp = kf[2];
      P_Km = kf[3];
      //proton detection
      fbonus = CheckBONUS(&P_p, "p", z);
      fclas12fa = CheckCLAS12FA(&P_p, "p", z);
      fclas12la = CheckCLAS12LA(&P_p, "p", z);
      if (fbonus > 0.0)
	DetectorResolutionSmear(&P_p, res_bonus);
      else if (fclas12fa > 0.0)
	DetectorResolutionSmear(&P_p, res_clas12fa);
      else if (fclas12la > 0.0)
	DetectorResolutionSmear(&P_p, res_clas12la);
      factorp = fbonus + fclas12fa + fclas12la;
      //K+ detection
      fbonus = CheckBONUS(&P_Kp, "K+", z);
      fclas12fa = CheckCLAS12FA(&P_Kp, "K+", z);
      fclas12la = CheckCLAS12LA(&P_Kp, "K+", z);
      if (fbonus > 0.0)
	DetectorResolutionSmear(&P_Kp, res_bonus);
      else if (fclas12fa > 0.0)
	DetectorResolutionSmear(&P_Kp, res_clas12fa);
      else if (fclas12la > 0.0)
	DetectorResolutionSmear(&P_Kp, res_clas12la);
      factorKp = fbonus + fclas12fa + fclas12la;
      //K- detection
      fbonus = CheckBONUS(&P_Km, "K-", z);
      fclas12fa = CheckCLAS12FA(&P_Km, "K-", z);
      fclas12la = CheckCLAS12LA(&P_Km, "K-", z);
      if (fbonus > 0.0)
	DetectorResolutionSmear(&P_Km, res_bonus);
      else if (fclas12fa > 0.0)
	DetectorResolutionSmear(&P_Km, res_clas12fa);
      else if (fclas12la > 0.0)
	DetectorResolutionSmear(&P_Km, res_clas12la);
      factorKm = fbonus + fclas12fa + fclas12la;

      ftotal = factorp * factorKp * factorKm;
      //ftotal = 1.0;
      if (ftotal > 0.0){
	P0_all = P0_p + P0_Kp + P0_Km;
	P_all = P_p + P_Kp + P_Km;

	h1->Fill(P0_all.M(), weight * ftotal / 2.0);

	if (MomentumCut(&P0_p, 0.0, 0.5) && MomentumCut(&P0_Kp, 0.0, 0.35) && MomentumCut(&P0_Km, 0.0, 0.35))
	  h1p->Fill(P0_all.M(), weight * ftotal / 2.0);

	h1s->Fill(P_all.M(), weight * ftotal / 2.0);

	if (MomentumCut(&P_p, 0.0, 0.5) && MomentumCut(&P_Kp, 0.0, 0.35) && MomentumCut(&P_Km, 0.0, 0.35))
	  h1sp->Fill(P_all.M(), weight * ftotal / 2.0);
      }
    }

    //channel 2
    weight = GenerateEvent_eNKK_withElectroproductionLambda1520Gold(ki, kf);
    if (weight > 0){
      P0_p = kf[1];
      P0_Kp = kf[2];
      P0_Km = kf[3];
      P_p = kf[1];
      P_Kp = kf[2];
      P_Km = kf[3];
      //proton detection
      fbonus = CheckBONUS(&P_p, "p", z);
      fclas12fa = CheckCLAS12FA(&P_p, "p", z);
      fclas12la = CheckCLAS12LA(&P_p, "p", z);
      if (fbonus > 0.0)
	DetectorResolutionSmear(&P_p, res_bonus);
      else if (fclas12fa > 0.0)
	DetectorResolutionSmear(&P_p, res_clas12fa);
      else if (fclas12la > 0.0)
	DetectorResolutionSmear(&P_p, res_clas12la);
      factorp = fbonus + fclas12fa + fclas12la;
      //K+ detection
      fbonus = CheckBONUS(&P_Kp, "K+", z);
      fclas12fa = CheckCLAS12FA(&P_Kp, "K+", z);
      fclas12la = CheckCLAS12LA(&P_Kp, "K+", z);
      if (fbonus > 0.0)
	DetectorResolutionSmear(&P_Kp, res_bonus);
      else if (fclas12fa > 0.0)
	DetectorResolutionSmear(&P_Kp, res_clas12fa);
      else if (fclas12la > 0.0)
	DetectorResolutionSmear(&P_Kp, res_clas12la);
      factorKp = fbonus + fclas12fa + fclas12la;
      //K- detection
      fbonus = CheckBONUS(&P_Km, "K-", z);
      fclas12fa = CheckCLAS12FA(&P_Km, "K-", z);
      fclas12la = CheckCLAS12LA(&P_Km, "K-", z);
      if (fbonus > 0.0)
	DetectorResolutionSmear(&P_Km, res_bonus);
      else if (fclas12fa > 0.0)
	DetectorResolutionSmear(&P_Km, res_clas12fa);
      else if (fclas12la > 0.0)
	DetectorResolutionSmear(&P_Km, res_clas12la);
      factorKm = fbonus + fclas12fa + fclas12la;

      ftotal = factorp * factorKp * factorKm;
      //ftotal = 1.0;
      if (ftotal > 0.0){
	P0_all = P0_p + P0_Kp + P0_Km;
	P_all = P_p + P_Kp + P_Km;

	h2->Fill(P0_all.M(), weight * ftotal / 2.0);

	if (MomentumCut(&P0_p, 0.0, 0.5) && MomentumCut(&P0_Kp, 0.0, 0.35) && MomentumCut(&P0_Km, 0.0, 0.35))
	  h2p->Fill(P0_all.M(), weight * ftotal / 2.0);

	h2s->Fill(P_all.M(), weight * ftotal / 2.0);

	if (MomentumCut(&P_p, 0.0, 0.5) && MomentumCut(&P_Kp, 0.0, 0.35) && MomentumCut(&P_Km, 0.0, 0.35))
	  h2sp->Fill(P_all.M(), weight * ftotal / 2.0);
      }
    }

    //channel 3
    weight = GenerateEvent_eNKK_withElectroproductionKKGold(ki, kf);
    if (weight > 0){
      P0_p = kf[1];
      P0_Kp = kf[2];
      P0_Km = kf[3];
      P_p = kf[1];
      P_Kp = kf[2];
      P_Km = kf[3];
      //proton detection
      fbonus = CheckBONUS(&P_p, "p", z);
      fclas12fa = CheckCLAS12FA(&P_p, "p", z);
      fclas12la = CheckCLAS12LA(&P_p, "p", z);
      if (fbonus > 0.0)
	DetectorResolutionSmear(&P_p, res_bonus);
      else if (fclas12fa > 0.0)
	DetectorResolutionSmear(&P_p, res_clas12fa);
      else if (fclas12la > 0.0)
	DetectorResolutionSmear(&P_p, res_clas12la);
      factorp = fbonus + fclas12fa + fclas12la;
      //K+ detection
      fbonus = CheckBONUS(&P_Kp, "K+", z);
      fclas12fa = CheckCLAS12FA(&P_Kp, "K+", z);
      fclas12la = CheckCLAS12LA(&P_Kp, "K+", z);
      if (fbonus > 0.0)
	DetectorResolutionSmear(&P_Kp, res_bonus);
      else if (fclas12fa > 0.0)
	DetectorResolutionSmear(&P_Kp, res_clas12fa);
      else if (fclas12la > 0.0)
	DetectorResolutionSmear(&P_Kp, res_clas12la);
      factorKp = fbonus + fclas12fa + fclas12la;
      //K- detection
      fbonus = CheckBONUS(&P_Km, "K-", z);
      fclas12fa = CheckCLAS12FA(&P_Km, "K-", z);
      fclas12la = CheckCLAS12LA(&P_Km, "K-", z);
      if (fbonus > 0.0)
	DetectorResolutionSmear(&P_Km, res_bonus);
      else if (fclas12fa > 0.0)
	DetectorResolutionSmear(&P_Km, res_clas12fa);
      else if (fclas12la > 0.0)
	DetectorResolutionSmear(&P_Km, res_clas12la);
      factorKm = fbonus + fclas12fa + fclas12la;
      
      ftotal = factorp * factorKp * factorKm;
      //ftotal = 1.0;
      if (ftotal > 0.0){
	P0_all = P0_p + P0_Kp + P0_Km;
	P_all = P_p + P_Kp + P_Km;

	h3->Fill(P0_all.M(), weight * ftotal / 2.0);

	if (MomentumCut(&P0_p, 0.0, 0.5) && MomentumCut(&P0_Kp, 0.0, 0.35) && MomentumCut(&P0_Km, 0.0, 0.35))
	  h3p->Fill(P0_all.M(), weight * ftotal / 2.0);

	h3s->Fill(P_all.M(), weight * ftotal / 2.0);

	if (MomentumCut(&P_p, 0.0, 0.5) && MomentumCut(&P_Kp, 0.0, 0.35) && MomentumCut(&P_Km, 0.0, 0.35))
	  h3sp->Fill(P_all.M(), weight * ftotal / 2.0);
      }
    }
    
  }
  
  h4->Add(h0); h4->Add(hh0); h4->Add(h1); h4->Add(h2); h4->Add(h3); 
  h4p->Add(h0p); h4p->Add(hh0p); h4p->Add(h1p); h4p->Add(h2p); h4p->Add(h3p); 
  h4s->Add(h0s); h4s->Add(hh0s); h4s->Add(h1s); h4s->Add(h2s); h4s->Add(h3s); 
  h4sp->Add(h0sp); h4sp->Add(hh0sp); h4sp->Add(h1sp); h4sp->Add(h2sp); h4sp->Add(h3sp);

  gStyle->SetOptStat(0);  

  TCanvas * cc = new TCanvas("cc", "", 800, 600);
  hmax = h1->GetMaximum();
  hB->GetYaxis()->SetRangeUser(0.0, 1.05 * hmax);
  TLegend * leg = new TLegend(0.62, 0.7, 0.9, 0.9);
  leg->AddEntry(h0, "with bound state", "l");
  leg->AddEntry(hh0, "with bound state", "l");
  leg->AddEntry(h1, "with #phi production", "l");
  leg->AddEntry(h2, "with #Lambda(1520) prod.", "l");
  leg->AddEntry(h3, "direct KK production", "l");
  //leg->AddEntry(h4, "total", "l");
  hB->Draw();
  h0->Draw("same");
  hh0->Draw("same");
  h1->Draw("same");
  h2->Draw("same");
  h3->Draw("same");
  //h4->Draw("same");
  leg->Draw("same");
  cc->Print("signal_background_detected.pdf");

  TCanvas * ccp = new TCanvas("ccp", "", 800, 600);
  hmax = h4p->GetMaximum();
  hBp->GetYaxis()->SetRangeUser(0.0, 1.1 * hmax);
  TLegend * legp = new TLegend(0.62, 0.7, 0.9, 0.9);
  legp->AddEntry(h0p, "with bound state", "l");
  legp->AddEntry(hh0p, "with bound state", "l");
  legp->AddEntry(h1p, "with #phi production", "l");
  legp->AddEntry(h2p, "with #Lambda(1520) prod.", "l");
  legp->AddEntry(h3p, "direct KK production", "l");
  //legp->AddEntry(h4p, "total", "l");
  hBp->Draw();
  h0p->Draw("same");
  hh0p->Draw("same");
  h1p->Draw("same");
  h2p->Draw("same");
  h3p->Draw("same");
  //h4p->Draw("same");
  legp->Draw("same");
  ccp->Print("signal_background_momentemcut_detected.pdf");

  TCanvas * ccs = new TCanvas("ccs", "", 800, 600);
  hmax = h1s->GetMaximum();
  hBs->GetYaxis()->SetRangeUser(0.0, 1.05 * hmax);
  TLegend * legs = new TLegend(0.62, 0.7, 0.9, 0.9);
  legs->AddEntry(h0s, "with bound state", "l");
  legs->AddEntry(hh0s, "with bound state", "l");
  legs->AddEntry(h1s, "with #phi production", "l");
  legs->AddEntry(h2s, "with #Lambda(1520) prod.", "l");
  legs->AddEntry(h3s, "direct KK production", "l");
  //legs->AddEntry(h4s, "total", "l");
  hBs->Draw();
  h0s->Draw("same");
  hh0s->Draw("same");
  h1s->Draw("same");
  h2s->Draw("same");
  h3s->Draw("same");
  //h4s->Draw("same");
  legs->Draw("same");
  ccs->Print("signal_background_smear_detected.pdf");

  TCanvas * ccsp = new TCanvas("ccsp", "", 800, 600);
  hmax = h4sp->GetMaximum();
  hBsp->GetYaxis()->SetRangeUser(0.0, 1.1 * hmax);
  TLegend * legsp = new TLegend(0.62, 0.7, 0.9, 0.9);
  legsp->AddEntry(h0sp, "with bound state", "l");
  legsp->AddEntry(hh0sp, "with bound state", "l");
  legsp->AddEntry(h1sp, "with #phi production", "l");
  legsp->AddEntry(h2sp, "with #Lambda(1520) prod.", "l");
  legsp->AddEntry(h3sp, "direct KK production", "l");
  //legsp->AddEntry(h4sp, "total", "l");
  hBsp->Draw();
  h0sp->Draw("same");
  hh0sp->Draw("same");
  h1sp->Draw("same");
  h2sp->Draw("same");
  h3sp->Draw("same");
  //h4sp->Draw("same");
  legsp->Draw("same");
  ccsp->Print("signal_background_smear_momentumcut_detected.pdf");


  int bin_low = hB->FindBin(1.94);
  int bin_high = hB->FindBin(1.96);
  cout << "Signal/h:   " << h0->Integral(1, -1) * lumi * time / Nsim << endl;
  cout << "Smeared:    " << h0s->Integral(1, -1) * lumi * time / Nsim << endl;
  cout << "After cut:  " << h0p->Integral(1, -1) * lumi * time / Nsim << endl;
  cout << "Smear & cut:" << h0sp->Integral(1, -1) * lumi * time / Nsim << endl;
  cout << endl;
  cout << "Total/h:    " << h4->Integral(1, -1) * lumi * time / Nsim << endl;
  cout << "Smeared:    " << h4s->Integral(1, -1) * lumi * time / Nsim << endl;
  cout << "After cut:  " << h4p->Integral(1, -1) * lumi * time / Nsim << endl;
  cout << "Smear & cut:" << h4sp->Integral(1, -1) * lumi * time / Nsim << endl;
  cout << endl;
  cout << "Signal/h:   " << h0->Integral(bin_low, bin_high) * lumi * time / Nsim << endl;
  cout << "Smeared:    " << h0s->Integral(bin_low, bin_high) * lumi * time / Nsim << endl;
  cout << "After cut:  " << h0p->Integral(bin_low, bin_high) * lumi * time / Nsim << endl;
  cout << "Smear & cut:" << h0sp->Integral(bin_low, bin_high) * lumi * time / Nsim << endl;
  cout << endl;
  cout << "Total/h:    " << h4->Integral(bin_low, bin_high) * lumi * time / Nsim << endl;
  cout << "Smeared:    " << h4s->Integral(bin_low, bin_high) * lumi * time / Nsim << endl;
  cout << "After cut:  " << h4p->Integral(bin_low, bin_high) * lumi * time / Nsim << endl;
  cout << "Smear & cut:" << h4sp->Integral(bin_low, bin_high) * lumi * time / Nsim << endl;
  cout << "=======================================" << endl;
  cout << endl;


  return 0;
}
