#include "Lcore.h"

int main(const int argc, const char * argv[]){

  if (argc < 2) return 1;

  const int opt = atoi(argv[1]);

  if (opt == 1){//wave function Ur
    MODEL::SetUr();
    double x[100], y[100];
    for (int i = 0; i < 100; i++){
      x[i] = 0.05 * i;
      y[i] = MODEL::Ur(x[i] / Phys::hbar);
    }
    TGraph * g0 = new TGraph(100, x, y);
    g0->SetLineColor(4);
    g0->SetLineWidth(2);
    g0->SetTitle("");
    g0->GetXaxis()->SetTitle("r (fm)");
    g0->GetXaxis()->CenterTitle(true);
    g0->GetXaxis()->SetTitleSize(0.06);
    g0->GetXaxis()->SetTitleOffset(1.15);
    g0->GetXaxis()->SetLabelSize(0.06);
    g0->GetXaxis()->SetRangeUser(0.0, 5.0);
    g0->GetXaxis()->SetNdivisions(6, 5, 0);
    g0->GetYaxis()->SetTitle("u(r) (GeV^{1/2})");
    g0->GetYaxis()->CenterTitle(true);
    g0->GetYaxis()->SetTitleSize(0.06);
    g0->GetYaxis()->SetTitleOffset(1.15);
    g0->GetYaxis()->SetLabelSize(0.06);
    g0->GetYaxis()->SetRangeUser(0.0, 0.5);
    g0->GetYaxis()->SetNdivisions(5, 5, 0);
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetBottomMargin(0.15);
    c0->SetLeftMargin(0.15);
    g0->DrawClone("al");
    c0->Print("gallary/wf.pdf");
  }

  
  if (opt == 2){//Veff
    MODEL::SetVeff();
    double x[100], y[100];
    for (int i = 0; i < 100; i++){
      x[i] = 0.05 * i;
      y[i] = MODEL::Veff(x[i] / Phys::hbar);
    }
    TGraph * g0 = new TGraph(100, x, y);
    g0->SetLineColor(4);
    g0->SetLineWidth(2);
    g0->SetTitle("");
    g0->GetXaxis()->SetTitle("r (fm)");
    g0->GetXaxis()->CenterTitle(true);
    g0->GetXaxis()->SetTitleSize(0.06);
    g0->GetXaxis()->SetTitleOffset(1.15);
    g0->GetXaxis()->SetLabelSize(0.06);
    g0->GetXaxis()->SetRangeUser(0.0, 5.0);
    g0->GetXaxis()->SetNdivisions(6, 5, 0);
    g0->GetYaxis()->SetTitle("V_{eff} (GeV)");
    g0->GetYaxis()->CenterTitle(true);
    g0->GetYaxis()->SetTitleSize(0.06);
    g0->GetYaxis()->SetTitleOffset(1.15);
    g0->GetYaxis()->SetLabelSize(0.06);
    g0->GetYaxis()->SetRangeUser(-0.6, 0.2);
    g0->GetYaxis()->SetNdivisions(5, 5, 0);
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetBottomMargin(0.15);
    c0->SetLeftMargin(0.15);
    g0->DrawClone("al");
    c0->Print("gallary/Veff.pdf"); 
  }


  if (opt == 3){//F(Q)
    MODEL::SetFQ();
    double x[100], y[100];
    for (int i = 0; i < 100; i++){
      x[i] = 0.02 * i;
      y[i] = MODEL::FQ(x[i]);
    }
    TGraph * g0 = new TGraph(100, x, y);
    g0->SetLineColor(4);
    g0->SetLineWidth(2);
    g0->SetTitle("");
    g0->GetXaxis()->SetTitle("Q (GeV)");
    g0->GetXaxis()->CenterTitle(true);
    g0->GetXaxis()->SetTitleSize(0.06);
    g0->GetXaxis()->SetTitleOffset(1.15);
    g0->GetXaxis()->SetLabelSize(0.06);
    g0->GetXaxis()->SetRangeUser(0.0, 2.0);
    g0->GetXaxis()->SetNdivisions(6, 5, 0);
    g0->GetYaxis()->SetTitle("F (GeV^{-1/2})");
    g0->GetYaxis()->CenterTitle(true);
    g0->GetYaxis()->SetTitleSize(0.06);
    g0->GetYaxis()->SetTitleOffset(1.15);
    g0->GetYaxis()->SetLabelSize(0.06);
    g0->GetYaxis()->SetRangeUser(0.0, 0.06);
    g0->GetYaxis()->SetNdivisions(5, 5, 0);
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetBottomMargin(0.15);
    c0->SetLeftMargin(0.15);
    g0->DrawClone("al");
    c0->Print("gallary/FQ.pdf"); 
  }

  if (opt == 4){//invariant mass vertex
    TFile * fs = new TFile("result/raw.root", "r");
    TH1D * h0 = (TH1D *) fs->Get("MpKK_BoundStateAll");
    TH1D * h1 = (TH1D *) fs->Get("MpKK_BoundStateKK");
    TH1D * h2 = (TH1D *) fs->Get("MpKK_phi");
    TH1D * h3 = (TH1D *) fs->Get("MpKK_Lambda1520");
    TH1D * h4 = (TH1D *) fs->Get("MpKK_directKK");
    double lumi = 1.0e35 * 1.0e-26 * pow(Phys::hbar, 2) / GOLD::NA;//eA GeV^2 s^-1
    double time = 3600.0;
    double binsize = 0.2;//MeV
    h0->GetXaxis()->SetTitle("M(pKK) (GeV)");
    h0->GetXaxis()->CenterTitle(true);
    h0->GetXaxis()->SetTitleSize(0.06);
    h0->GetXaxis()->SetTitleOffset(1.15);
    h0->GetXaxis()->SetLabelSize(0.06);
    h0->GetXaxis()->SetRangeUser(1.88, 2.32);
    h0->GetXaxis()->SetNdivisions(6, 5, 0);
    h0->GetYaxis()->SetTitle("counts / hour / MeV");
    h0->GetYaxis()->CenterTitle(true);
    h0->GetYaxis()->SetTitleSize(0.06);
    h0->GetYaxis()->SetTitleOffset(1.15);
    h0->GetYaxis()->SetLabelSize(0.06);
    h0->GetYaxis()->SetRangeUser(0.0, 15.0);
    h0->GetYaxis()->SetNdivisions(5, 5, 0);
    //
    h0->Scale(lumi*time/binsize);
    h1->Scale(lumi*time/binsize);
    h2->Scale(lumi*time/binsize);
    h3->Scale(lumi*time/binsize);
    h4->Scale(lumi*time/binsize);\
    //
    h0->SetStats(0);
    h0->SetLineColor(1);
    h1->SetLineColor(4);
    h2->SetLineColor(2);
    h3->SetLineColor(6);
    h4->SetLineColor(7);
    TLegend * leg = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg->AddEntry(h0, "pKK from N-#phi", "l");
    leg->AddEntry(h1, "only KK from N-#phi", "l");
    leg->AddEntry(h2, "#phi production", "l");
    leg->AddEntry(h3, "#Lambda(1520) production", "l");
    leg->AddEntry(h4, "direct KK production", "l");
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);
    h0->SetMaximum(15.0);
    h0->DrawClone("axis");
    h0->DrawClone("same");
    h1->DrawClone("same");
    h2->DrawClone("same");
    h3->DrawClone("same");
    h4->DrawClone("same");
    leg->DrawClone("same");
    c0->Print("gallary/MpKKraw.pdf");
  }

  if (opt == 5){//Momentum vertex
    TFile * fs = new TFile("result/raw.root", "r");
    TH2D * d0 = (TH2D *) fs->Get("Momentum_p_Kp_BoundStateAll");
    TH2D * d1 = (TH2D *) fs->Get("Momentum_p_Kp_BoundStateKK");
    TH2D * d2 = (TH2D *) fs->Get("Momentum_p_Kp_phi");
    TH2D * d3a = (TH2D *) fs->Get("Momentum_p_Kp_Lambda1520");
    TH2D * d3b = (TH2D *) fs->Get("Momentum_p_Km_Lambda1520");
    TH2D * d4 = (TH2D *) fs->Get("Momentum_p_Kp_KK");
    TGraph * g0 = new TGraph(1000);
    TGraph * g1 = new TGraph(1000);
    TGraph * g2 = new TGraph(1000);
    TGraph * g3a = new TGraph(1000);
    TGraph * g3b = new TGraph(1000);
    TGraph * g4 = new TGraph(1000);
    double x, y;
    for (int i = 0; i < 1000; i++){
      d0->GetRandom2(x, y);
      g0->SetPoint(i, x, y);
      d1->GetRandom2(x, y);
      g1->SetPoint(i, x, y);
      d2->GetRandom2(x, y);
      g2->SetPoint(i, x, y);
      d3a->GetRandom2(x, y);
      g3a->SetPoint(i, x, y);
      d3b->GetRandom2(x, y);
      g3b->SetPoint(i, x, y);
      d4->GetRandom2(x, y);
      g4->SetPoint(i, x, y);
    }
    d0->GetXaxis()->SetTitle("P(K) (GeV)");
    d0->GetXaxis()->CenterTitle(true);
    d0->GetXaxis()->SetTitleSize(0.06);
    d0->GetXaxis()->SetTitleOffset(1.15);
    d0->GetXaxis()->SetLabelSize(0.06);
    d0->GetXaxis()->SetRangeUser(0.0, 2.0);
    d0->GetXaxis()->SetNdivisions(6, 5, 0);
    d0->GetYaxis()->SetTitle("P(p) (GeV)");
    d0->GetYaxis()->CenterTitle(true);
    d0->GetYaxis()->SetTitleSize(0.06);
    d0->GetYaxis()->SetTitleOffset(1.15);
    d0->GetYaxis()->SetLabelSize(0.06);
    d0->GetYaxis()->SetRangeUser(0.0, 2.0);
    d0->GetYaxis()->SetNdivisions(5, 5, 0);
    d0->SetStats(0);
    //
    g0->SetMarkerStyle(20);
    g0->SetMarkerColor(1);
    g0->SetMarkerSize(0.6);
    g1->SetMarkerStyle(21);
    g1->SetMarkerColor(4);
    g1->SetMarkerSize(0.6);
    g2->SetMarkerStyle(22);
    g2->SetMarkerColor(2);
    g2->SetMarkerSize(0.6);
    g3a->SetMarkerStyle(23);
    g3a->SetMarkerColor(6);
    g3a->SetMarkerSize(0.6);
    g3b->SetMarkerStyle(24);
    g3b->SetMarkerColor(3);
    g3b->SetMarkerSize(0.6);
    g4->SetMarkerStyle(25);
    g4->SetMarkerColor(7);
    g4->SetMarkerSize(0.6);
    TLegend * leg = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg->AddEntry(g0, "pKK from N-#phi", "p");
    leg->AddEntry(g1, "only KK from N-#phi", "p");
    leg->AddEntry(g2, "#phi production", "p");
    leg->AddEntry(g3a, "#Lambda(1520) production (K^{+})", "p");
    leg->AddEntry(g3b, "#Lambda(1520) production (K^{-})", "p");
    leg->AddEntry(g4, "direct KK production", "p");
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);
    d0->DrawClone("axis");
    g0->DrawClone("psame");
    g1->DrawClone("psame");
    g2->DrawClone("psame");
    g3a->DrawClone("psame");
    g3b->DrawClone("psame");
    g4->DrawClone("psame");
    leg->DrawClone("same");
    c0->Print("gallary/Momentumraw.pdf");
  }

  if (opt == 6){//invariant mass vertex
    TFile * fs = new TFile("result/detected.root", "r");
    TH1D * h0 = (TH1D *) fs->Get("MpKK_BoundStateAll");
    TH1D * h1 = (TH1D *) fs->Get("MpKK_BoundStateKK");
    TH1D * h2 = (TH1D *) fs->Get("MpKK_phi");
    TH1D * h3 = (TH1D *) fs->Get("MpKK_Lambda1520");
    TH1D * h4 = (TH1D *) fs->Get("MpKK_directKK");
    double lumi = 1.0e35 * 1.0e-26 * pow(Phys::hbar, 2) / GOLD::NA;//eA GeV^2 s^-1
    double time = 3600.0;
    double binsize = 0.2;//MeV
    h0->GetXaxis()->SetTitle("M(pKK) (GeV)");
    h0->GetXaxis()->CenterTitle(true);
    h0->GetXaxis()->SetTitleSize(0.06);
    h0->GetXaxis()->SetTitleOffset(1.15);
    h0->GetXaxis()->SetLabelSize(0.06);
    h0->GetXaxis()->SetRangeUser(1.88, 2.32);
    h0->GetXaxis()->SetNdivisions(6, 5, 0);
    h0->GetYaxis()->SetTitle("counts / hour / MeV");
    h0->GetYaxis()->CenterTitle(true);
    h0->GetYaxis()->SetTitleSize(0.06);
    h0->GetYaxis()->SetTitleOffset(1.15);
    h0->GetYaxis()->SetLabelSize(0.06);
    h0->GetYaxis()->SetRangeUser(0.0, 15.0);
    h0->GetYaxis()->SetNdivisions(5, 5, 0);
    //
    h0->Scale(lumi*time/binsize);
    h1->Scale(lumi*time/binsize);
    h2->Scale(lumi*time/binsize);
    h3->Scale(lumi*time/binsize);
    h4->Scale(lumi*time/binsize);\
    //
    h0->SetStats(0);
    h0->SetLineColor(1);
    h1->SetLineColor(4);
    h2->SetLineColor(2);
    h3->SetLineColor(6);
    h4->SetLineColor(7);
    TLegend * leg = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg->AddEntry(h0, "pKK from N-#phi", "l");
    leg->AddEntry(h1, "only KK from N-#phi", "l");
    leg->AddEntry(h2, "#phi production", "l");
    leg->AddEntry(h3, "#Lambda(1520) production", "l");
    leg->AddEntry(h4, "direct KK production", "l");
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);
    h0->SetMaximum(1.5);
    h0->DrawClone("axis");
    h0->DrawClone("same");
    h1->DrawClone("same");
    h2->DrawClone("same");
    h3->DrawClone("same");
    h4->DrawClone("same");
    leg->DrawClone("same");
    c0->Print("gallary/MpKKdetected.pdf");
  }

  if (opt == 7){//Momentum vertex
    TFile * fs = new TFile("result/detected.root", "r");
    TH2D * d0 = (TH2D *) fs->Get("Momentum_p_Kp_BoundStateAll");
    TH2D * d1 = (TH2D *) fs->Get("Momentum_p_Kp_BoundStateKK");
    TH2D * d2 = (TH2D *) fs->Get("Momentum_p_Kp_phi");
    TH2D * d3a = (TH2D *) fs->Get("Momentum_p_Kp_Lambda1520");
    TH2D * d3b = (TH2D *) fs->Get("Momentum_p_Km_Lambda1520");
    TH2D * d4 = (TH2D *) fs->Get("Momentum_p_Kp_KK");
    TGraph * g0 = new TGraph(1000);
    TGraph * g1 = new TGraph(1000);
    TGraph * g2 = new TGraph(1000);
    TGraph * g3a = new TGraph(1000);
    TGraph * g3b = new TGraph(1000);
    TGraph * g4 = new TGraph(1000);
    double x, y;
    for (int i = 0; i < 1000; i++){
      d0->GetRandom2(x, y);
      g0->SetPoint(i, x, y);
      d1->GetRandom2(x, y);
      g1->SetPoint(i, x, y);
      d2->GetRandom2(x, y);
      g2->SetPoint(i, x, y);
      d3a->GetRandom2(x, y);
      g3a->SetPoint(i, x, y);
      d3b->GetRandom2(x, y);
      g3b->SetPoint(i, x, y);
      d4->GetRandom2(x, y);
      g4->SetPoint(i, x, y);
    }
    d0->GetXaxis()->SetTitle("P(K) (GeV)");
    d0->GetXaxis()->CenterTitle(true);
    d0->GetXaxis()->SetTitleSize(0.06);
    d0->GetXaxis()->SetTitleOffset(1.15);
    d0->GetXaxis()->SetLabelSize(0.06);
    d0->GetXaxis()->SetRangeUser(0.0, 2.0);
    d0->GetXaxis()->SetNdivisions(6, 5, 0);
    d0->GetYaxis()->SetTitle("P(p) (GeV)");
    d0->GetYaxis()->CenterTitle(true);
    d0->GetYaxis()->SetTitleSize(0.06);
    d0->GetYaxis()->SetTitleOffset(1.15);
    d0->GetYaxis()->SetLabelSize(0.06);
    d0->GetYaxis()->SetRangeUser(0.0, 2.0);
    d0->GetYaxis()->SetNdivisions(5, 5, 0);
    d0->SetStats(0);
    //
    g0->SetMarkerStyle(20);
    g0->SetMarkerColor(1);
    g0->SetMarkerSize(0.6);
    g1->SetMarkerStyle(21);
    g1->SetMarkerColor(4);
    g1->SetMarkerSize(0.6);
    g2->SetMarkerStyle(22);
    g2->SetMarkerColor(2);
    g2->SetMarkerSize(0.6);
    g3a->SetMarkerStyle(23);
    g3a->SetMarkerColor(6);
    g3a->SetMarkerSize(0.6);
    g3b->SetMarkerStyle(24);
    g3b->SetMarkerColor(3);
    g3b->SetMarkerSize(0.6);
    g4->SetMarkerStyle(25);
    g4->SetMarkerColor(7);
    g4->SetMarkerSize(0.6);
    TLegend * leg = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg->AddEntry(g0, "pKK from N-#phi", "p");
    leg->AddEntry(g1, "only KK from N-#phi", "p");
    leg->AddEntry(g2, "#phi production", "p");
    leg->AddEntry(g3a, "#Lambda(1520) production (K^{+})", "p");
    leg->AddEntry(g3b, "#Lambda(1520) production (K^{-})", "p");
    leg->AddEntry(g4, "direct KK production", "p");
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);
    d0->DrawClone("axis");
    g0->DrawClone("psame");
    g1->DrawClone("psame");
    g2->DrawClone("psame");
    g3a->DrawClone("psame");
    g3b->DrawClone("psame");
    g4->DrawClone("psame");
    leg->DrawClone("same");
    c0->Print("gallary/Momentumdetected.pdf");
  }

  if (opt == 8){//Momentum distribution colz raw
    TH2D * hB = new TH2D("hB", "", 1, 0.0, 2.0, 1, 0.0, 2.0);
    hB->GetXaxis()->SetTitle("P(K) (GeV)");
    hB->GetXaxis()->CenterTitle(true);
    hB->GetXaxis()->SetTitleSize(0.06);
    hB->GetXaxis()->SetTitleOffset(1.15);
    hB->GetXaxis()->SetLabelSize(0.06);
    hB->GetXaxis()->SetRangeUser(0.0, 2.0);
    hB->GetXaxis()->SetNdivisions(6, 5, 0);
    hB->GetYaxis()->SetTitle("P(p) (GeV)");
    hB->GetYaxis()->CenterTitle(true);
    hB->GetYaxis()->SetTitleSize(0.06);
    hB->GetYaxis()->SetTitleOffset(1.15);
    hB->GetYaxis()->SetLabelSize(0.06);
    hB->GetYaxis()->SetRangeUser(0.0, 2.0);
    hB->GetYaxis()->SetNdivisions(5, 5, 0);
    hB->SetStats(0);
    TFile * fs = new TFile("result/raw.root", "r");
    TH2D * d0 = (TH2D *) fs->Get("Momentum_p_Kp_BoundStateAll");
    TH2D * d1 = (TH2D *) fs->Get("Momentum_p_Kp_BoundStateKK");
    TH2D * d2 = (TH2D *) fs->Get("Momentum_p_Kp_phi");
    TH2D * d3a = (TH2D *) fs->Get("Momentum_p_Kp_Lambda1520");
    TH2D * d3b = (TH2D *) fs->Get("Momentum_p_Km_Lambda1520");
    TH2D * d4 = (TH2D *) fs->Get("Momentum_p_Kp_KK");
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);
    hB->SetTitle("K and p from bound state");
    hB->DrawClone("axis");
    d0->DrawClone("colsame");
    c0->Print("gallary/distributionraw.pdf(", "pdf");
    hB->SetTitle("K from bound state");
    hB->DrawClone("axis");
    d1->DrawClone("colsame");
    c0->Print("gallary/distributionraw.pdf", "pdf");
    hB->SetTitle("with #phi production");
    hB->DrawClone("axis");
    d2->DrawClone("colsame");
    c0->Print("gallary/distributionraw.pdf", "pdf");
    hB->SetTitle("with #Lambda(1520) production - K^{+}");
    hB->DrawClone("axis");
    d3a->DrawClone("colsame");
    c0->Print("gallary/distributionraw.pdf", "pdf");
    hB->SetTitle("with #Lambda(1520) production - K^{-}");
    hB->DrawClone("axis");
    d3b->DrawClone("colsame");
    c0->Print("gallary/distributionraw.pdf", "pdf");
    hB->SetTitle("with direct KK production");
    hB->DrawClone("axis");
    d4->DrawClone("colsame");
    c0->Print("gallary/distributionraw.pdf)", "pdf");
  }

  if (opt == 9){//Momentum distribution colz raw
    TH2D * hB = new TH2D("hB", "", 1, 0.0, 2.0, 1, 0.0, 2.0);
    hB->GetXaxis()->SetTitle("P(K) (GeV)");
    hB->GetXaxis()->CenterTitle(true);
    hB->GetXaxis()->SetTitleSize(0.06);
    hB->GetXaxis()->SetTitleOffset(1.15);
    hB->GetXaxis()->SetLabelSize(0.06);
    hB->GetXaxis()->SetRangeUser(0.0, 2.0);
    hB->GetXaxis()->SetNdivisions(6, 5, 0);
    hB->GetYaxis()->SetTitle("P(p) (GeV)");
    hB->GetYaxis()->CenterTitle(true);
    hB->GetYaxis()->SetTitleSize(0.06);
    hB->GetYaxis()->SetTitleOffset(1.15);
    hB->GetYaxis()->SetLabelSize(0.06);
    hB->GetYaxis()->SetRangeUser(0.0, 2.0);
    hB->GetYaxis()->SetNdivisions(5, 5, 0);
    hB->SetStats(0);
    TFile * fs = new TFile("result/detected.root", "r");
    TH2D * d0 = (TH2D *) fs->Get("Momentum_p_Kp_BoundStateAll");
    TH2D * d1 = (TH2D *) fs->Get("Momentum_p_Kp_BoundStateKK");
    TH2D * d2 = (TH2D *) fs->Get("Momentum_p_Kp_phi");
    TH2D * d3a = (TH2D *) fs->Get("Momentum_p_Kp_Lambda1520");
    TH2D * d3b = (TH2D *) fs->Get("Momentum_p_Km_Lambda1520");
    TH2D * d4 = (TH2D *) fs->Get("Momentum_p_Kp_KK");
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);
    hB->SetTitle("K and p from bound state");
    hB->DrawClone("axis");
    d0->DrawClone("colsame");
    c0->Print("gallary/distributiondetected.pdf(", "pdf");
    hB->SetTitle("K from bound state");
    hB->DrawClone("axis");
    d1->DrawClone("colsame");
    c0->Print("gallary/distributiondetected.pdf", "pdf");
    hB->SetTitle("with #phi production");
    hB->DrawClone("axis");
    d2->DrawClone("colsame");
    c0->Print("gallary/distributiondetected.pdf", "pdf");
    hB->SetTitle("with #Lambda(1520) production - K^{+}");
    hB->DrawClone("axis");
    d3a->DrawClone("colsame");
    c0->Print("gallary/distributiondetected.pdf", "pdf");
    hB->SetTitle("with #Lambda(1520) production - K^{-}");
    hB->DrawClone("axis");
    d3b->DrawClone("colsame");
    c0->Print("gallary/distributiondetected.pdf", "pdf");
    hB->SetTitle("with direct KK production");
    hB->DrawClone("axis");
    d4->DrawClone("colsame");
    c0->Print("gallary/distributiondetected.pdf)", "pdf");
  }

  return 0;
}
