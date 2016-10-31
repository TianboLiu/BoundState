#include "Lelectro.h"

using namespace std;

int main(int argc, char * argv[]){
  double Ebeam = atof(argv[1]);//Get electron beam energy
  const double Eqmin = 1.0;//Set virtual photon energy minimum
  const double Eqmax = 4.0;//Set virtual photon energy maximum

  gRandom->SetSeed(0);
  SetFunctions();
  //SetTagger(Ebeam - Eqmax, Ebeam - Eqmin);//Set scattered electron energy range
  SetTagger(0.5, 4.5);

  const Long64_t Nsim = 1000000;
  const int Npt = 1000;

  double lumi = 1.0e35 * 1.0e-26 * pow(0.197327, 2) / 12.0;//eA GeV^2 s^-1

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

  TH1D * b1D = new TH1D("b1D", "", 1, 1.88, 2.32);
  TH2D * b2D = new TH2D("b2D", "", 1, 0.0, 1.0, 1, 0.0, 1.4);
  TH2D * b3D = new TH2D("b3D", "", 1, 0.0, 0.6, 1, 0.0, M_PI/2.0);
  gStyle->SetOptStat(0);
  b1D->SetTitle("M(pK^{+}K^{-}) spectra");
  b1D->GetXaxis()->SetTitle("M(pK^{+}K^{-}) / GeV");
  b1D->GetXaxis()->CenterTitle();
  b1D->GetXaxis()->SetTitleOffset(1.0);
  b1D->GetXaxis()->SetTitleSize(0.05);
  b1D->GetYaxis()->SetRangeUser(0.0, 1.0);
  b2D->SetTitle("proton kaon momentum correlation");
  b2D->GetXaxis()->SetTitle("P(K) / GeV");
  b2D->GetXaxis()->CenterTitle();
  b2D->GetXaxis()->SetTitleOffset(1.0);
  b2D->GetXaxis()->SetTitleSize(0.05);
  b2D->GetYaxis()->SetTitle("P(p) / GeV");
  b2D->GetYaxis()->CenterTitle();
  b2D->GetYaxis()->SetTitleOffset(1.0);
  b2D->GetYaxis()->SetTitleSize(0.05);
  b3D->SetTitle("Momentum : Angel distribution");
  b3D->GetXaxis()->SetTitle("Momentum (GeV)");
  b3D->GetXaxis()->CenterTitle();
  b3D->GetXaxis()->SetTitleOffset(1.0);
  b3D->GetXaxis()->SetTitleSize(0.05);
  b3D->GetYaxis()->SetTitle("Polar angle (rad)");
  b3D->GetYaxis()->CenterTitle();
  b3D->GetYaxis()->SetTitleOffset(1.0);
  b3D->GetYaxis()->SetTitleSize(0.05); 

  //invariant mass
  TH1D * h0 = new TH1D("h0", "", 1000, 1.88, 2.32);
  TH1D * h1 = new TH1D("h1", "", 1000, 1.88, 2.32);
  TH1D * h2 = new TH1D("h2", "", 1000, 1.88, 2.32);
  TH1D * h3 = new TH1D("h3", "", 1000, 1.88, 2.32);
  TH1D * h4 = new TH1D("h4", "", 1000, 1.88, 2.32);
  TH1D * h5 = new TH1D("h5", "", 1000, 1.88, 2.32);
  h0->SetLineColor(1);
  h1->SetLineColor(4);
  h2->SetLineColor(2);
  h3->SetLineColor(6);
  h4->SetLineColor(7);
  h5->SetLineColor(5);

  //momentum correlation
  TH2D * d0 = new TH2D("d0", "", 500, 0.0, 1.5, 500, 0.0, 1.5);
  TH2D * d1 = new TH2D("d1", "", 500, 0.0, 1.5, 500, 0.0, 1.5);
  TH2D * d2 = new TH2D("d2", "", 500, 0.0, 1.5, 500, 0.0, 1.5);
  TH2D * d3 = new TH2D("d3", "", 500, 0.0, 1.5, 500, 0.0, 1.5);
  TH2D * d3p = new TH2D("d3p", "", 500, 0.0, 1.5, 500, 0.0, 1.5);
  TH2D * d4 = new TH2D("d4", "", 500, 0.0, 1.5, 500, 0.0, 1.5);
  //momentum correlation graph
  TGraph * g0 = new TGraph(Npt);
  TGraph * g1 = new TGraph(Npt);
  TGraph * g2 = new TGraph(Npt);
  TGraph * g3 = new TGraph(Npt);
  TGraph * g3p = new TGraph(Npt);
  TGraph * g4 = new TGraph(Npt);
  g0->SetMarkerColor(1);
  g0->SetMarkerStyle(20);
  g0->SetMarkerSize(0.3);
  g1->SetMarkerColor(4);
  g1->SetMarkerStyle(24);
  g1->SetMarkerSize(0.2);
  g2->SetMarkerColor(2);
  g2->SetMarkerStyle(26);
  g2->SetMarkerSize(0.2);
  g3->SetMarkerColor(6);
  g3->SetMarkerStyle(25);
  g3->SetMarkerSize(0.2);
  g3p->SetMarkerColor(3);
  g3p->SetMarkerStyle(25);
  g3p->SetMarkerSize(0.2);
  g4->SetMarkerColor(7);
  g4->SetMarkerStyle(27);
  g4->SetMarkerSize(0.2);

  double x, y;

  for (int i = 0; i < Nsim; i++){
    if (i%100000==0) cout << i << endl;
    
    z = gRandom->Uniform(0.0, 1.0);//Get vertex position

    if (true){
      weight = GenerateEvent_eNKKN_withElectroproductionBoundStateCarbon(ki, kf);//Generate a bound state event kf: N (N, K+, K-)
      P0_p = kf[4];//p from bound state
      P0_Kp = kf[2];
      P0_Km = kf[3];
      P0_all = P0_p + P0_Kp + P0_Km;
      h0->Fill(P0_all.M(), weight);
      d0->Fill(P0_Kp.P(), P0_p.P(), weight);

      P0_p = kf[1];//produced p
      P0_all = P0_p + P0_Kp + P0_Km;
      h1->Fill(P0_all.M(), weight);
      d1->Fill(P0_Kp.P(), P0_p.P(), weight);
 
      weight = GenerateEvent_eNKK_withElectroproductionPhiCarbon(ki, kf);//Generate a unbound event with phi production
      P0_p = kf[1];
      P0_Kp = kf[2];
      P0_Km = kf[3];
      P0_all = P0_p + P0_Kp + P0_Km;
      h2->Fill(P0_all.M(), weight);
      d2->Fill(P0_Kp.P(), P0_p.P(), weight);

      weight = GenerateEvent_eNKK_withElectroproductionLambda1520Carbon(ki, kf);//Generate a unbound event with Lambda1520 production
      P0_p = kf[1];
      P0_Kp = kf[2];
      P0_Km = kf[3];
      P0_all = P0_p + P0_Kp + P0_Km;
      h3->Fill(P0_all.M(), weight);
      d3->Fill(P0_Kp.P(), P0_p.P(), weight);
      d3p->Fill(P0_Km.P(), P0_p.P(), weight);

      weight = GenerateEvent_eNKK_withElectroproductionKKCarbon(ki, kf);//Generate a unbound event with direct KK production
      P0_p = kf[1];
      P0_Kp = kf[2];
      P0_Km = kf[3];
      P0_all = P0_p + P0_Kp + P0_Km;
      h4->Fill(P0_all.M(), weight);
      d4->Fill(P0_Kp.P(), P0_p.P(), weight);
    }
  }

  cout << "==============================" << endl;
  cout << h0->Integral(1, -1) * lumi / Nsim * 3600.0 << endl;
  cout << "==============================" << endl;


  for (int i = 0; i < Npt; i++){
    d0->GetRandom2(x, y);
    g0->SetPoint(i, x, y);
    d1->GetRandom2(x, y);
    g1->SetPoint(i, x, y);
    d2->GetRandom2(x, y);
    g2->SetPoint(i, x, y);
    d3->GetRandom2(x, y);
    g3->SetPoint(i, x, y);
    d3p->GetRandom2(x, y);
    g3p->SetPoint(i, x, y);
    d4->GetRandom2(x, y);
    g4->SetPoint(i, x, y);
  }

  h5->Add(h0); h5->Add(h1); h5->Add(h2); h5->Add(h3); h5->Add(h4);
  double hmax = h4->GetMaximum();
  h0->Scale(0.85/hmax);
  h1->Scale(0.85/hmax);
  h2->Scale(0.85/hmax);
  h3->Scale(0.85/hmax);
  h4->Scale(0.85/hmax);
  h5->Scale(0.85/hmax);
  TLegend * leg0 = new TLegend(0.62, 0.7, 0.9, 0.9);
  leg0->AddEntry(h0, "with bound state", "l");
  leg0->AddEntry(h1, "with bound state", "l");
  leg0->AddEntry(h2, "with #phi production", "l");
  leg0->AddEntry(h3, "with #Lambda(1520) prod.", "l");
  leg0->AddEntry(h4, "direct production", "l");
  //leg0->AddEntry(h5, "total", "l");
  TCanvas * c0 = new TCanvas("c0", "c0", 800, 600);
  b1D->Draw();
  h0->Draw("same");
  h1->Draw("same");
  h2->Draw("same");
  h3->Draw("same");
  h4->Draw("same");
  //h5->Draw("same");
  leg0->Draw("same");
  c0->Print("background0.pdf");


  hmax = h5->GetMaximum();
  h0->Scale(0.85/hmax);
  h1->Scale(0.85/hmax);
  h2->Scale(0.85/hmax);
  h3->Scale(0.85/hmax);
  h4->Scale(0.85/hmax);
  h5->Scale(0.85/hmax);
  TLegend * leg1 = new TLegend(0.62, 0.7, 0.9, 0.9);
  leg1->AddEntry(h0, "with bound state", "l");
  leg1->AddEntry(h1, "with bound state", "l");
  leg1->AddEntry(h2, "with #phi production", "l");
  leg1->AddEntry(h3, "with #Lambda(1520) prod.", "l");
  leg1->AddEntry(h4, "direct production", "l");
  leg1->AddEntry(h5, "total", "l");
  TCanvas * c1 = new TCanvas("c1", "c1", 800, 600);
  b1D->Draw();
  h0->Draw("same");
  h1->Draw("same");
  h2->Draw("same");
  h3->Draw("same");
  h4->Draw("same");
  h5->Draw("same");
  leg1->Draw("same");
  c1->Print("background1.pdf");

  TLegend * leg2 = new TLegend(0.62, 0.7, 0.9, 0.9);
  leg2->AddEntry(g0, "with bound state", "p");
  leg2->AddEntry(g1, "with bound state", "p");
  leg2->AddEntry(g2, "with #phi production", "p");
  leg2->AddEntry(g3, "with #Lambda(1520) prod. (K^{+})", "p");
  leg2->AddEntry(g3p, "with #Lambda(1520) prod. (K^{-})", "p");
  leg2->AddEntry(g4, "direct production", "p");
  TCanvas * c2 = new TCanvas("c2", "c2", 800, 600);
  b2D->Draw();
  g0->Draw("psame");
  g1->Draw("psame");
  g2->Draw("psame");
  g3->Draw("psame");
  g3p->Draw("psame");
  g4->Draw("psame");
  leg2->Draw("same");
  c2->Print("momentumcorrelation.pdf");

  return 0;
}
