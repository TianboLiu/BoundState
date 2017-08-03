#include "Lcore.h"

int main(const int argc, const char * argv[]){
  if (argc < 2) return 1;

  const int opt = atoi(argv[1]);
  gStyle->SetPalette(55);

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
    TH2D * H2 = (TH2D *) fs->Get("ECosTheta_Gold");
    TH1D * h1[10];
    TH1D * H1[10];

    double lumi = 1.0e35 * 1.0e-26 * pow(Phys::hbar, 2);
    double time = 3600.0;

    double binsize = 22.0;

    h2->Scale(lumi*time/binsize);
    H2->Scale(lumi*time/GOLD::NA/binsize);


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

    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);

    for (int i = 0; i < 10; i++){
      h1[i] = h2->ProjectionX(Form("h%d",i), 20 * i + 1, 20 * i + 20);
      H1[i] = H2->ProjectionX(Form("H%d",i), 20 * i + 1, 20 * i + 20);
      h1[i]->SetLineColor(4);
      H1[i]->SetLineColor(2);
   
      h0->GetXaxis()->SetRangeUser(0.7, 1.7);
      h0->SetMinimum(0);
      h0->SetMaximum(max(h1[i]->GetMaximum(), H1[i]->GetMaximum()) * 1.2);

      c0->Clear();
      h0->DrawClone("axis");
      h1[i]->DrawClone("same");
      H1[i]->DrawClone("same");

      if (i == 0)
	c0->Print("gallary/bumpdetectedcos.pdf(", "pdf");
      else if (i == 9)
	c0->Print("gallary/bumpdetectedcos.pdf)", "pdf");
      else
	c0->Print("gallary/bumpdetectedcos.pdf", "pdf");
    }


  }



  return 0;
}


     

