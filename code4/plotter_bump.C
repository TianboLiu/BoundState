#include "Lcore.h"
#include "TLatex.h"

int main(const int argc, const char * argv[]){
  if (argc < 2) return 1;

  const int opt = atoi(argv[1]);
  gStyle->SetPalette(55);

  cout << opt << endl;

  if (opt == 1){//total Egamma cm distribution
    TFile * fs = new TFile("result/bumpdetected.root", "r");
    TH1D * hEg = (TH1D *) fs->Get("E_g");
    TH1D * HEg = (TH1D *) fs->Get("E_g_Gold");

    double lumi = 1.0e35 * 1.0e-26 * pow(Phys::hbar, 2);
    double time = 3600.0;

    TH1D * h0 = new TH1D("h0", "", 1, 0.7, 1.7);
    h0->SetStats(0);
    h0->GetXaxis()->SetTitle("E_{#gamma (c.m.)} (GeV)");
    h0->GetXaxis()->CenterTitle(true);
    h0->GetXaxis()->SetTitleSize(0.06);
    h0->GetXaxis()->SetTitleOffset(1.15);
    h0->GetXaxis()->SetLabelSize(0.055);
    h0->GetXaxis()->SetRangeUser(0.7, 1.7);
    h0->GetXaxis()->SetNdivisions(6, 5, 0);
    h0->GetYaxis()->SetTitle("counts / hour / MeV");
    h0->GetYaxis()->CenterTitle(true);
    h0->GetYaxis()->SetTitleSize(0.06);
    h0->GetYaxis()->SetTitleOffset(1.15);
    h0->GetYaxis()->SetLabelSize(0.055);

    double binsize = 4.4;//MeV

    hEg->Scale(lumi*time/binsize);
    HEg->Scale(lumi*time/GOLD::NA/binsize);

    h0->SetMinimum(0);
    h0->SetMaximum(hEg->GetMaximum() * 1.2);

    hEg->SetLineColor(4);
    HEg->SetLineColor(2);

    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);
    h0->DrawClone("axis");
    hEg->DrawClone("same");
    HEg->DrawClone("same");

    c0->Print("gallary/bumpdetected.pdf");
  }

  if (opt == 2){//total Egamma cm distribution
    TFile * fs = new TFile("result/bumpdetected.root", "r");
    TH2D * h2 = (TH2D *) fs->Get("ECosTheta");
    TH2D * h2a = (TH2D *) fs->Get("ECosTheta_clas12");
    TH2D * h2b = (TH2D *) fs->Get("ECosTheta_bonus12");

    TH1D * h1[40];
    TH1D * h1a[40];
    TH1D * h1b[40];

    double lumi = 1.0e35 * 1.0e-26 * pow(Phys::hbar, 2);
    double time = 3600.0;

    double binsize = 10.0;

    h2->Scale(lumi*time/binsize);
    h2a->Scale(lumi*time/binsize);
    h2b->Scale(lumi*time/binsize);

    TH1D * h0 = new TH1D("h0", "", 1, 1.9, 2.9);
    h0->SetStats(0);
    h0->GetXaxis()->SetTitle("W (GeV)");
    h0->GetXaxis()->CenterTitle(true);
    h0->GetXaxis()->SetTitleSize(0.06);
    h0->GetXaxis()->SetTitleOffset(1.15);
    h0->GetXaxis()->SetLabelSize(0.055);
    h0->GetXaxis()->SetRangeUser(1.9, 2.9);
    h0->GetXaxis()->SetNdivisions(6, 5, 0);
    h0->GetYaxis()->SetTitle("counts / hour / MeV");
    h0->GetYaxis()->CenterTitle(true);
    h0->GetYaxis()->SetTitleSize(0.06);
    h0->GetYaxis()->SetTitleOffset(1.15);
    h0->GetYaxis()->SetLabelSize(0.055);

    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);

    TLatex latex;
    latex.SetTextSize(0.05);
    latex.SetTextAlign(22);
    latex.SetTextFont(22);
    for (int i = 0; i < 40; i++){
      h1[i] = h2->ProjectionX(Form("hh%d",i), 5 * i + 1, 5 * i + 5);
      h1[i]->SetLineColor(1);
      h1a[i] = h2a->ProjectionX(Form("ta%d",i), 5 * i + 1, 5 * i + 5);
      h1a[i]->SetLineColor(4);
      h1b[i] = h2b->ProjectionX(Form("tb%d",i), 5 * i + 1, 5 * i + 5);
      h1b[i]->SetLineColor(2);
      
      h0->GetXaxis()->SetRangeUser(1.9, 2.9);
      h0->SetMinimum(0);
      h0->SetMaximum(h1[i]->GetMaximum() * 1.2);

      c0->Clear();
      h0->DrawClone("axis");
      h1[i]->DrawClone("same");
      h1a[i]->DrawClone("same");
      h1b[i]->DrawClone("same");
      latex.DrawLatex(2.4, h1[i]->GetMaximum() * 1.1, Form("%.2f < cos#theta^{#phi}_{c.m.} < %.2f", -1.0 + 0.05 * i, -0.95 + 0.05 * i));

      if (i == 0)
	c0->Print("gallary/bumpdetectedcos.pdf(", "pdf");
      else if (i == 39)
	c0->Print("gallary/bumpdetectedcos.pdf)", "pdf");
      else
	c0->Print("gallary/bumpdetectedcos.pdf", "pdf");
    }
  }


  if (opt == 3){//phi acceptance
    TFile * fs = new TFile("acceptance/acceptance_phi.root", "r");
    TH2D * h0 = (TH2D *) fs->Get("acceptance_PTheta_all");
    TH2D * ha = (TH2D *) fs->Get("acceptance_PTheta_clas12");
    TH2D * hb = (TH2D *) fs->Get("acceptance_PTheta_bonus12");
    cout << "aaa" << endl;
    TH2D * hB = new TH2D("hB", "", 1, 0.0, 4.0, 1, 0.0, M_PI);
    hB->SetStats(0);
    hB->GetXaxis()->SetTitle("P (GeV)");
    hB->GetXaxis()->CenterTitle(true);
    hB->GetXaxis()->SetTitleSize(0.06);
    hB->GetXaxis()->SetTitleOffset(1.15);
    hB->GetXaxis()->SetLabelSize(0.055);
    hB->GetXaxis()->SetRangeUser(0.0, 4.0);
    hB->GetXaxis()->SetNdivisions(6, 5, 0);
    hB->GetYaxis()->SetTitle("#theta (rad)");
    hB->GetYaxis()->CenterTitle(true);
    hB->GetYaxis()->SetTitleSize(0.06);
    hB->GetYaxis()->SetTitleOffset(1.15);
    hB->GetYaxis()->SetLabelSize(0.055);
    hB->GetYaxis()->SetRangeUser(0.0, M_PI);
    hB->GetXaxis()->SetNdivisions(6, 5, 0);
    
    TLatex latex;
    latex.SetTextFont(22);
    latex.SetTextAlign(22);
    latex.SetTextSize(0.055);
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);
    hB->DrawClone("axis");
    h0->DrawClone("colzsame");
    latex.DrawLatex(2.2, M_PI * 1.05, "CLAS12 + BONUS12");
    c0->Print("gallary/phiacceptance.pdf(", "pdf");

    c0->Clear();
    hB->DrawClone("axis");
    ha->DrawClone("colzsame");
    latex.DrawLatex(2.2, M_PI * 1.05, "CLAS12 only");
    c0->Print("gallary/phiacceptance.pdf", "pdf");
    
    c0->Clear();
    hB->DrawClone("axis");
    hb->DrawClone("colzsame");
    latex.DrawLatex(2.2, M_PI * 1.05, "BONUS12 only");
    c0->Print("gallary/phiacceptance.pdf)", "pdf");

    

    
  }


  
  return 0;
}


     

