#include "Lsimulation.h"

using namespace std;

int main(int argc, char * argv[]){
  gRandom->SetSeed(1);
  SetFunctions();
  
  Long64_t Nsim = 2000000;
  double Eq[2] = {1.37, 1.57};
  double A = 12.0;
  double lumi_e_limit = 1.0e35 * 1.0e-26 * pow(0.197327, 2);//GeV^2 / s
  double rad_length_factor = 0.06;
  double YFactor = BremsstrahlungYFactor(11.0, Eq);
  double lumi = lumi_e_limit * YFactor * rad_length_factor / (2.0 * A);
  double time = 3600;
  double res[3] = {0.01, 0.001, 0.004};//resolution dp/p, dth, dphi

  TLorentzVector ki[2];
  TLorentzVector kf0[4], kf1[4], kf2[4], kf3[4];
  double weight0, weight1, weight2, weight3;
  double decay[4];

  TLorentzVector tmp1, tmp2;
  double hmax;

  TH1D * hB = new TH1D("hB", "", 1000, 1.88, 2.32);
  hB->SetTitle("M(pK^{+}K^{-}) spectra");
  hB->GetXaxis()->SetTitle("M(pK^{+}K^{-}) / GeV");
  hB->GetXaxis()->CenterTitle();
  hB->GetXaxis()->SetTitleOffset(1.0);
  hB->GetXaxis()->SetTitleSize(0.05);
  hB->GetYaxis()->SetRangeUser(0.0, 1.0);
  TH1D * h0 = new TH1D("h0", "", 1000, 1.88, 2.32);
  TH1D * h1 = new TH1D("h1", "", 1000, 1.88, 2.32);
  TH1D * h2 = new TH1D("h2", "", 1000, 1.88, 2.32);
  TH1D * h3 = new TH1D("h3", "", 1000, 1.88, 2.32);
  TH1D * h4 = new TH1D("h4", "", 1000, 1.88, 2.32);
  h0->SetLineColor(1);
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
  TH1D * h1s = new TH1D("h1s", "", 1000, 1.88, 2.32);
  TH1D * h2s = new TH1D("h2s", "", 1000, 1.88, 2.32);
  TH1D * h3s = new TH1D("h3s", "", 1000, 1.88, 2.32);
  TH1D * h4s = new TH1D("h4s", "", 1000, 1.88, 2.32);
  h0s->SetLineColor(1);
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
  TH1D * h1p = new TH1D("h1p", "", 1000, 1.88, 2.32);
  TH1D * h2p = new TH1D("h2p", "", 1000, 1.88, 2.32);
  TH1D * h3p = new TH1D("h3p", "", 1000, 1.88, 2.32);
  TH1D * h4p = new TH1D("h4p", "", 1000, 1.88, 2.32);
  h0p->SetLineColor(1);
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
  TH1D * h1sp = new TH1D("h1sp", "", 1000, 1.88, 2.32);
  TH1D * h2sp = new TH1D("h2sp", "", 1000, 1.88, 2.32);
  TH1D * h3sp = new TH1D("h3sp", "", 1000, 1.88, 2.32);
  TH1D * h4sp = new TH1D("h4sp", "", 1000, 1.88, 2.32);
  h0sp->SetLineColor(1);
  h1sp->SetLineColor(2);
  h2sp->SetLineColor(6);
  h3sp->SetLineColor(7);
  h4sp->SetLineColor(5);
 
  for (Long64_t i = 0; i < Nsim; i++){
    if (i%100000 == 0) cout << i << endl;

    GenerateBremsstrahlungPhoton(ki, Eq);//Generate a photon
    
    GenerateEventNNKKwithBoundState(ki, kf0, &weight0);
    GenerateEventNKKwithPhi(ki, kf1, &weight1);
    GenerateEventNKKwithLambda1520(ki, kf2, &weight2);
    GenerateEventNKKwithout(ki, kf3, &weight3);
  
    //channel 0
    if (weight0 > 0){
      tmp1 = kf0[1] + kf0[2] + kf0[3];
      decay[1] = 1.0;
      decay[2] = LongLifetimeDecayFactor(&kf0[2], 3.712);
      decay[3] = LongLifetimeDecayFactor(&kf0[3], 3.712);
      h0->Fill(tmp1.M(), weight0 * decay[1] * decay[2] * decay[3] * A * (A - 1.0) / 2.0);
      if (MomentumCut(&kf0[1], 0.0, 0.6) && MomentumCut(&kf0[2], 0.0, 0.4) && MomentumCut(&kf0[3], 0.0, 0.4)){
	h0p->Fill(tmp1.M(), weight0 * decay[1] * decay[2] * decay[3] * A * (A - 1.0) / 2.0);
      }
      DetectorResolutionSmear(&kf0[1], res);
      DetectorResolutionSmear(&kf0[2], res);
      DetectorResolutionSmear(&kf0[3], res);
      tmp2 = kf0[1] + kf0[2] + kf0[3];
      h0s->Fill(tmp2.M(), weight0 * decay[1] * decay[2] * decay[3] * A * (A - 1.0) / 2.0);
      if (MomentumCut(&kf0[1], 0.0, 0.6) && MomentumCut(&kf0[2], 0.0, 0.4) && MomentumCut(&kf0[3], 0.0, 0.4)){
	h0sp->Fill(tmp1.M(), weight0 * decay[1] * decay[2] * decay[3] * A * (A - 1.0) / 2.0);
      }
    }

    //channel 1
    if (weight1 > 0){
      tmp1 = kf1[0] + kf1[1] + kf1[2];
      decay[0] = 1.0;
      decay[1] = LongLifetimeDecayFactor(&kf1[1], 3.712);
      decay[2] = LongLifetimeDecayFactor(&kf1[2], 3.712);
      h1->Fill(tmp1.M(), weight1 * decay[0] * decay[1] * decay[2] * A / 2.0);
      if (MomentumCut(&kf1[0], 0.0, 0.6) && MomentumCut(&kf1[1], 0.0, 0.4) && MomentumCut(&kf1[2], 0.0, 0.4)){
	h1p->Fill(tmp1.M(), weight1 * decay[0] * decay[1] * decay[2] * A / 2.0);
      }
      DetectorResolutionSmear(&kf1[0], res);
      DetectorResolutionSmear(&kf1[1], res);
      DetectorResolutionSmear(&kf1[2], res);
      tmp2 = kf1[0] + kf1[1] + kf1[2];
      h1s->Fill(tmp2.M(), weight1 * decay[0] * decay[1] * decay[2] * A / 2.0);
      if (MomentumCut(&kf1[0], 0.0, 0.6) && MomentumCut(&kf1[1], 0.0, 0.4) && MomentumCut(&kf1[2], 0.0, 0.4)){
	h1sp->Fill(tmp1.M(), weight1 * decay[0] * decay[1] * decay[2] * A / 2.0);
      }
    }

    //channel 2
    if (weight2 > 0){
      tmp1 = kf2[0] + kf2[1] + kf2[2];
      decay[0] = 1.0;
      decay[1] = LongLifetimeDecayFactor(&kf2[1], 3.712);
      decay[2] = LongLifetimeDecayFactor(&kf2[2], 3.712);
      h2->Fill(tmp1.M(), weight2 * decay[0] * decay[1] * decay[2] * A / 2.0);
      if (MomentumCut(&kf2[0], 0.0, 0.6) && MomentumCut(&kf2[1], 0.0, 0.4) && MomentumCut(&kf2[2], 0.0, 0.4)){
	h2p->Fill(tmp1.M(), weight2 * decay[0] * decay[1] * decay[2] * A / 2.0);
      }
      DetectorResolutionSmear(&kf2[0], res);
      DetectorResolutionSmear(&kf2[1], res);
      DetectorResolutionSmear(&kf2[2], res);
      tmp2 = kf2[0] + kf2[1] + kf2[2];
      h2s->Fill(tmp2.M(), weight2 * decay[0] * decay[1] * decay[2] * A / 2.0);
      if (MomentumCut(&kf2[0], 0.0, 0.6) && MomentumCut(&kf2[1], 0.0, 0.4) && MomentumCut(&kf2[2], 0.0, 0.4)){
	h2sp->Fill(tmp1.M(), weight2 * decay[0] * decay[1] * decay[2] * A / 2.0);
      }
    }

    //channel 3
    if (weight3 > 0){
      tmp1 = kf3[0] + kf3[1] + kf3[2];
      decay[0] = 1.0;
      decay[1] = LongLifetimeDecayFactor(&kf3[1], 3.712);
      decay[2] = LongLifetimeDecayFactor(&kf3[2], 3.712);
      h3->Fill(tmp1.M(), weight3 * decay[0] * decay[1] * decay[2] * A / 2.0);
      if (MomentumCut(&kf3[0], 0.0, 0.6) && MomentumCut(&kf3[1], 0.0, 0.4) && MomentumCut(&kf3[2], 0.0, 0.4)){
	h3p->Fill(tmp1.M(), weight3 * decay[0] * decay[1] * decay[2] * A / 2.0);
      }
      DetectorResolutionSmear(&kf3[0], res);
      DetectorResolutionSmear(&kf3[1], res);
      DetectorResolutionSmear(&kf3[2], res);
      tmp2 = kf3[0] + kf3[1] + kf3[2];
      h3s->Fill(tmp2.M(), weight3 * decay[0] * decay[1] * decay[2] * A / 2.0);
      if (MomentumCut(&kf3[0], 0.0, 0.6) && MomentumCut(&kf3[1], 0.0, 0.4) && MomentumCut(&kf3[2], 0.0, 0.4)){
	h3sp->Fill(tmp1.M(), weight3 * decay[0] * decay[1] * decay[2] * A / 2.0);
      }
    }
  }

  h4->Add(h0); h4->Add(h1); h4->Add(h2); h4->Add(h3); 
  h4p->Add(h0p); h4p->Add(h1p); h4p->Add(h2p); h4p->Add(h3p); 
  h4s->Add(h0s); h4s->Add(h1s); h4s->Add(h2s); h4s->Add(h3s); 
  h4sp->Add(h0sp); h4sp->Add(h1sp); h4sp->Add(h2sp); h4sp->Add(h3sp);

  gStyle->SetOptStat(0);
  
  

  TCanvas * cc = new TCanvas("cc", "", 800, 600);
  hmax = h4->GetMaximum();
  hB->GetYaxis()->SetRangeUser(0.0, 1.2 * hmax);
  TLegend * leg = new TLegend(0.62, 0.7, 0.9, 0.9);
  leg->AddEntry(h0, "with bound state", "l");
  leg->AddEntry(h1, "with #phi production", "l");
  leg->AddEntry(h2, "with #Lambda(1520) prod.", "l");
  leg->AddEntry(h3, "direct KK production", "l");
  leg->AddEntry(h4, "total", "l");
  hB->Draw();
  h0->Draw("same");
  h1->Draw("same");
  h2->Draw("same");
  h3->Draw("same");
  h4->Draw("same");
  leg->Draw("same");
  cc->Print("signal_background.pdf");

  TCanvas * ccp = new TCanvas("ccp", "", 800, 600);
  hmax = h4p->GetMaximum();
  hBp->GetYaxis()->SetRangeUser(0.0, 1.2 * hmax);
  TLegend * legp = new TLegend(0.62, 0.7, 0.9, 0.9);
  legp->AddEntry(h0p, "with bound state", "l");
  legp->AddEntry(h1p, "with #phi production", "l");
  legp->AddEntry(h2p, "with #Lambda(1520) prod.", "l");
  legp->AddEntry(h3p, "direct KK production", "l");
  legp->AddEntry(h4p, "total", "l");
  hBp->Draw();
  h0p->Draw("same");
  h1p->Draw("same");
  h2p->Draw("same");
  h3p->Draw("same");
  h4p->Draw("same");
  legp->Draw("same");
  ccp->Print("signal_background_momentemcut.pdf");

  TCanvas * ccs = new TCanvas("ccs", "", 800, 600);
  hmax = h4s->GetMaximum();
  hBs->GetYaxis()->SetRangeUser(0.0, 1.2 * hmax);
  TLegend * legs = new TLegend(0.62, 0.7, 0.9, 0.9);
  legs->AddEntry(h0s, "with bound state", "l");
  legs->AddEntry(h1s, "with #phi production", "l");
  legs->AddEntry(h2s, "with #Lambda(1520) prod.", "l");
  legs->AddEntry(h3s, "direct KK production", "l");
  legs->AddEntry(h4s, "total", "l");
  hBs->Draw();
  h0s->Draw("same");
  h1s->Draw("same");
  h2s->Draw("same");
  h3s->Draw("same");
  h4s->Draw("same");
  legs->Draw("same");
  ccs->Print("signal_background_smear.pdf");

  TCanvas * ccsp = new TCanvas("ccsp", "", 800, 600);
  hmax = h4sp->GetMaximum();
  hBsp->GetYaxis()->SetRangeUser(0.0, 1.2 * hmax);
  TLegend * legsp = new TLegend(0.62, 0.7, 0.9, 0.9);
  legsp->AddEntry(h0sp, "with bound state", "l");
  legsp->AddEntry(h1sp, "with #phi production", "l");
  legsp->AddEntry(h2sp, "with #Lambda(1520) prod.", "l");
  legsp->AddEntry(h3sp, "direct KK production", "l");
  legsp->AddEntry(h4sp, "total", "l");
  hBsp->Draw();
  h0sp->Draw("same");
  h1sp->Draw("same");
  h2sp->Draw("same");
  h3sp->Draw("same");
  h4sp->Draw("same");
  legsp->Draw("same");
  ccsp->Print("signal_background_smear_momentumcut.pdf");

  cout << "Signal/h:   " << h0->Integral(1, -1) * lumi * time / Nsim << endl;
  cout << "Smeared:    " << h0s->Integral(1, -1) * lumi * time / Nsim << endl;
  cout << "After cut:  " << h0p->Integral(1, -1) * lumi * time / Nsim << endl;
  cout << "Smear & cut:" << h0sp->Integral(1, -1) * lumi * time / Nsim << endl;

  return 0;
}
