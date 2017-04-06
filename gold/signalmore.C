#include "Lgoldmore.h"

using namespace std;

int main(int argc, char * argv[]){

  double Ebeam = atof(argv[1]);//Get electron beam energy
  cout << "Ebeam: " << Ebeam << endl;

  gRandom->SetSeed(1);
  SetFunctions();
  SetStoppingPower();
  SetSigmaProtonPionpPionm();

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
  hB->GetXaxis()->SetTitle("M(pK^{+}K^{-}) (GeV)");
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
  hBs->SetTitle("before momentum cuts");
  hBs->GetXaxis()->SetTitle("M(pK^{+}K^{-}) / GeV");
  hBs->GetXaxis()->CenterTitle(true);
  hBs->GetXaxis()->SetTitleOffset(1.15);
  hBs->GetXaxis()->SetTitleSize(0.06);
  hBs->GetXaxis()->SetLabelSize(0.055);
  hBs->GetYaxis()->SetTitle("counts / hour");
  hBs->GetYaxis()->CenterTitle(true);
  hBs->GetYaxis()->SetLabelSize(0.055);
  hBs->GetYaxis()->SetTitleSize(0.06);
  hBs->GetYaxis()->SetTitleOffset(1.15);
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
  hBp->GetXaxis()->SetTitle("M(pK^{+}K^{-}) (GeV)");
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
  hBsp->SetTitle("after momentum cuts");
  hBsp->GetXaxis()->SetTitle("M(pK^{+}K^{-}) (GeV)");
  hBsp->GetXaxis()->CenterTitle(true);
  hBsp->GetXaxis()->SetLabelSize(0.055);
  hBsp->GetXaxis()->SetTitleOffset(1.15);
  hBsp->GetXaxis()->SetTitleOffset(1.15);
  hBsp->GetXaxis()->SetTitleSize(0.06);
  hBsp->GetYaxis()->SetTitle("counts / hour");
  hBsp->GetYaxis()->CenterTitle(true);
  hBsp->GetYaxis()->SetLabelFont(42);
  hBsp->GetYaxis()->SetLabelSize(0.045);
  hBsp->GetYaxis()->SetTitleSize(0.06);
  hBsp->GetYaxis()->SetTitleFont(42);
  hBsp->GetYaxis()->SetTitleOffset(1.15);
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
    if (i%1000000 == 0) cout << i << endl;

    z = gRandom->Uniform(0.0, 1.0);

    //channel 0
    weight = GenerateEvent_eNpipi_withElectroproductionpipiGold(ki, kf);
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
      fbonus = CheckBONUS(&P_Kp, "p", z);
      fclas12fa = CheckCLAS12FA(&P_Kp, "p", z);
      fclas12la = CheckCLAS12LA(&P_Kp, "p", z);
      if (fbonus > 0.0)
	DetectorResolutionSmear(&P_Kp, res_bonus);
      else if (fclas12fa > 0.0)
	DetectorResolutionSmear(&P_Kp, res_clas12fa);
      else if (fclas12la > 0.0)
	DetectorResolutionSmear(&P_Kp, res_clas12la);
      factorKp = fbonus + fclas12fa + fclas12la;
      //K- detection
      fbonus = CheckBONUS(&P_Km, "p", z);
      fclas12fa = CheckCLAS12FA(&P_Km, "p", z);
      fclas12la = CheckCLAS12LA(&P_Km, "p", z);
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
	MisPIDtoKaon(&P0_Kp);
	MisPIDtoKaon(&P0_Km);
	MisPIDtoKaon(&P_Kp);
	MisPIDtoKaon(&P_Km);

	P0_all = P0_p + P0_Kp + P0_Km;
	P_all = P_p + P_Kp + P_Km;

	h0->Fill(P0_all.M(), weight * ftotal * 79.0/197.0);

	if (MomentumCut(&P0_p, 0.0, 0.5) && MomentumCut(&P0_Kp, 0.0, 0.35) && MomentumCut(&P0_Km, 0.0, 0.35))
	  h0p->Fill(P0_all.M(), weight * ftotal * 79.0/197.0);

	h0s->Fill(P_all.M(), weight * ftotal * 79.0/197.0);

	if (MomentumCut(&P_p, 0.0, 0.5) && MomentumCut(&P_Kp, 0.0, 0.35) && MomentumCut(&P_Km, 0.0, 0.35))
	  h0sp->Fill(P_all.M(), weight * ftotal * 79.0/197.0);
      }
    }

  }
  
  gStyle->SetOptStat(0);  

  TCanvas * cc = new TCanvas("cc", "", 800, 600);
  hmax = h0->GetMaximum();
  hB->GetYaxis()->SetRangeUser(0.0, 1.05 * hmax);

  hB->Draw();
  h0->Draw("same");
  cc->Print("twopi_detected.pdf");
  cc->Print("Plot/twopi_detected.C");

  TCanvas * ccp = new TCanvas("ccp", "", 800, 600);
  hmax = h0p->GetMaximum();
  hBp->GetYaxis()->SetRangeUser(0.0, 1.1 * hmax);

  hBp->Draw();
  h0p->Draw("same");

  ccp->Print("twopi_momentemcut_detected.pdf");
  ccp->Print("Plot/twopi_momentemcut_detected.C");

  TCanvas * ccs = new TCanvas("ccs", "", 800, 600);
  ccs->SetHighLightColor(2);
  ccs->Range(1.825,-0.02166037,2.375,0.1949433);
  ccs->SetFillColor(0);
  ccs->SetBorderMode(0);
  ccs->SetBorderSize(2);
  ccs->SetFrameBorderMode(0);
  ccs->SetFrameBorderMode(0);
  ccs->SetLeftMargin(0.15);
  ccs->SetBottomMargin(0.15);
  hmax = h0s->GetMaximum();
  hBs->GetYaxis()->SetRangeUser(0.0, 1.05 * hmax);

  hBs->Draw();
  h0s->Draw("same");

  ccs->Print("twopi_smear_detected.pdf");
  ccs->Print("Plot/twopi_smear_detected.C");
 
  TCanvas * ccsp = new TCanvas("ccsp", "", 800, 600);
  ccsp->SetHighLightColor(2);
  ccsp->Range(1.825,-0.002231509,2.375,0.02008358);
  ccsp->SetFillColor(0);
  ccsp->SetBorderMode(0);
  ccsp->SetBorderSize(2);
  ccsp->SetFrameBorderMode(0);
  ccsp->SetFrameBorderMode(0);
  ccsp->SetLeftMargin(0.15);
  ccsp->SetBottomMargin(0.15);
  hmax = h0sp->GetMaximum();
  hBsp->GetYaxis()->SetRangeUser(0.0, 1.1 * hmax);

  hBsp->Draw();
  h0sp->Draw("same");

  ccsp->Print("twopi_smear_momentumcut_detected.pdf");
  ccsp->Print("Plot/twopi_smear_momentumcut_detected.C");

  int bin_low = hB->FindBin(1.94);
  int bin_high = hB->FindBin(1.96);
  cout << "Signal/h:   " << h0->Integral(1, -1) * lumi * time / Nsim << endl;
  cout << "Smeared:    " << h0s->Integral(1, -1) * lumi * time / Nsim << endl;
  cout << "After cut:  " << h0p->Integral(1, -1) * lumi * time / Nsim << endl;
  cout << "Smear & cut:" << h0sp->Integral(1, -1) * lumi * time / Nsim << endl;
  cout << endl;
  cout << "Signal/h:   " << h0->Integral(bin_low, bin_high) * lumi * time / Nsim << endl;
  cout << "Smeared:    " << h0s->Integral(bin_low, bin_high) * lumi * time / Nsim << endl;
  cout << "After cut:  " << h0p->Integral(bin_low, bin_high) * lumi * time / Nsim << endl;
  cout << "Smear & cut:" << h0sp->Integral(bin_low, bin_high) * lumi * time / Nsim << endl;
  cout << "=======================================" << endl;
  cout << endl;


  return 0;
}
