/* plots for quasi-free production 2020/09 */

#include "Lcore.h"

int SetStyle(TH1D * hB){
  hB->GetXaxis()->CenterTitle(true);
  hB->GetXaxis()->SetTitleSize(0.06);
  hB->GetXaxis()->SetTitleOffset(1.15);
  hB->GetXaxis()->SetLabelSize(0.06);
  hB->GetXaxis()->SetNdivisions(6, 5, 0);
  hB->GetYaxis()->CenterTitle(true);
  hB->GetYaxis()->SetTitleSize(0.06);
  hB->GetYaxis()->SetTitleOffset(1.15);
  hB->GetYaxis()->SetLabelSize(0.06);
  hB->GetYaxis()->SetNdivisions(6, 5, 0);
  hB->SetStats(0);
  return 0;
}

int SetStyle(TH2D * hB){
  hB->GetXaxis()->CenterTitle(true);
  hB->GetXaxis()->SetTitleSize(0.06);
  hB->GetXaxis()->SetTitleOffset(1.15);
  hB->GetXaxis()->SetLabelSize(0.06);
  hB->GetXaxis()->SetNdivisions(6, 5, 0);
  hB->GetYaxis()->CenterTitle(true);
  hB->GetYaxis()->SetTitleSize(0.06);
  hB->GetYaxis()->SetTitleOffset(1.15);
  hB->GetYaxis()->SetLabelSize(0.06);
  hB->GetYaxis()->SetNdivisions(6, 5, 0);
  hB->SetStats(0);
  return 0;
}
  

int main(const int argc, const char * argv[]){
  if (argc < 2){
    cout << "./plot202009-E8.5 <opt>" << endl;
    cout << "1: mass" << endl;
    return 1;
  }

  const int opt = atoi(argv[1]);

  gStyle->SetPalette(55);
  gStyle->SetOptStat(0);

  TFile * fs1 = new TFile("result202009-E8.5/photo-qf-m-smeared.root", "r");
  TFile * fs2 = new TFile("result202009-E8.5/electro-qf-m-smeared.root", "r");
 
  if (opt == 1){
    TH2D * h1 = (TH2D *) fs1->Get("EgPJpsi");
    TH2D * h2 = (TH2D *) fs2->Get("EgPJpsi");
    h1->Add(h2);
    h1->Scale(47*24*0.8);
   
    TH1D * h0 = h1->ProjectionX("Eg", 0, -1);   
    h0->Rebin(10);

    for (int i = 1; i <= 14; i++){
      h0->SetBinError(i, sqrt(h0->GetBinContent(i)));
    }
      
    
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);

    SetStyle(h0);
    h0->GetYaxis()->SetTitle("Counts in 47 days");
    h0->SetLineColor(4);
    h0->SetLineWidth(3);  
    h0->SetFillColorAlpha(4, 0.35); 
    h0->DrawClone("pe2");
    h0->DrawClone("esame");
  
    c0->Print("figures202009-E8.5/gammaproj100MeV.pdf");
  }
  
  if (opt == 2){
    double eps = 1e-10;
    TH2D * h1 = (TH2D *) fs1->Get("EgPJpsi");
    TH2D * h2 = (TH2D *) fs2->Get("EgPJpsi");
    h1->Add(h2);
    h1->Scale(47*24*0.8);

    TH1D * ha = h1->ProjectionY("7.2-7.4", h1->GetXaxis()->FindBin(7.2+eps), h1->GetXaxis()->FindBin(7.4-eps));
    ha->SetTitle("E_{#gamma}: [7.2,7.4]GeV;P_{J/#psi} (GeV);Counts in 47 days");
    TH1D * hb = h1->ProjectionY("7.4-7.6", h1->GetXaxis()->FindBin(7.4+eps), h1->GetXaxis()->FindBin(7.6-eps));
    hb->SetTitle("E_{#gamma}: [7.4,7.6]GeV;P_{J/#psi} (GeV);Counts in 47 days");
    TH1D * hc = h1->ProjectionY("7.6-7.8", h1->GetXaxis()->FindBin(7.6+eps), h1->GetXaxis()->FindBin(7.8-eps));
    hc->SetTitle("E_{#gamma}: [7.6,7.8]GeV;P_{J/#psi} (GeV);Counts in 47 days");
    TH1D * hd = h1->ProjectionY("7.8-8.0", h1->GetXaxis()->FindBin(7.8+eps), h1->GetXaxis()->FindBin(8.0-eps));
    hd->SetTitle("E_{#gamma}: [7.8,8.0]GeV;P_{J/#psi} (GeV);Counts in 47 days");
    TH1D * he = h1->ProjectionY("8.0-8.2", h1->GetXaxis()->FindBin(8.0+eps), h1->GetXaxis()->FindBin(8.2-eps));
    he->SetTitle("E_{#gamma}: [8.0,8.2]GeV;P_{J/#psi} (GeV);Counts in 47 days");
    TH1D * hf = h1->ProjectionY("8.2-8.4", h1->GetXaxis()->FindBin(8.2+eps), h1->GetXaxis()->FindBin(8.4-eps));
    hf->SetTitle("E_{#gamma}: [8.2,8.4]GeV;P_{J/#psi} (GeV);Counts in 47 days");
    TH1D * hg = h1->ProjectionY("8.4-8.6", h1->GetXaxis()->FindBin(8.4+eps), h1->GetXaxis()->FindBin(8.6-eps));
    hg->SetTitle("E_{#gamma} > 8.4GeV;P_{J/#psi} (GeV);Counts in 47 days");
    
    
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);

    SetStyle(ha);
    SetStyle(hb);
    SetStyle(hc);
    SetStyle(hd);
    SetStyle(he);
    SetStyle(hf);
    SetStyle(hg);

 
    ha->DrawClone("");
    c0->Print("figures202009-E8.5/gamma-pjpsi.pdf(", "pdf");

    hb->DrawClone("");
    c0->Print("figures202009-E8.5/gamma-pjpsi.pdf", "pdf");

    hc->DrawClone("");
    c0->Print("figures202009-E8.5/gamma-pjpsi.pdf", "pdf");

    hd->DrawClone("");
    c0->Print("figures202009-E8.5/gamma-pjpsi.pdf", "pdf");

    he->DrawClone("");
    c0->Print("figures202009-E8.5/gamma-pjpsi.pdf", "pdf");
    
    hf->DrawClone("");
    c0->Print("figures202009-E8.5/gamma-pjpsi.pdf", "pdf");

    hg->DrawClone("");
    c0->Print("figures202009-E8.5/gamma-pjpsi.pdf)", "pdf");
  }

  if (opt == 3){
    double eps = 1e-10;
    TH2D * h1 = (TH2D *) fs1->Get("EgThetaJpsi");
    TH2D * h2 = (TH2D *) fs2->Get("EgThetaJpsi");
    h1->Add(h2);
    h1->Scale(47*24*0.8);

    TH1D * ha = h1->ProjectionY("7.2-7.4", h1->GetXaxis()->FindBin(7.2+eps), h1->GetXaxis()->FindBin(7.4-eps));
    ha->SetTitle("E_{#gamma}: [7.2,7.4]GeV;#theta_{J/#psi} (deg);Counts in 47 days");
    TH1D * hb = h1->ProjectionY("7.4-7.6", h1->GetXaxis()->FindBin(7.4+eps), h1->GetXaxis()->FindBin(7.6-eps));
    hb->SetTitle("E_{#gamma}: [7.4,7.6]GeV;#theta_{J/#psi} (deg);Counts in 47 days");
    TH1D * hc = h1->ProjectionY("7.6-7.8", h1->GetXaxis()->FindBin(7.6+eps), h1->GetXaxis()->FindBin(7.8-eps));
    hc->SetTitle("E_{#gamma}: [7.6,7.8]GeV;#theta_{J/#psi} (deg);Counts in 47 days");
    TH1D * hd = h1->ProjectionY("7.8-8.0", h1->GetXaxis()->FindBin(7.8+eps), h1->GetXaxis()->FindBin(8.0-eps));
    hd->SetTitle("E_{#gamma}: [7.8,8.0]GeV;#theta_{J/#psi} (deg);Counts in 47 days");
    TH1D * he = h1->ProjectionY("8.0-8.2", h1->GetXaxis()->FindBin(8.0+eps), h1->GetXaxis()->FindBin(8.2-eps));
    he->SetTitle("E_{#gamma}: [8.0,8.2]GeV;#theta_{J/#psi} (deg);Counts in 47 days");
    TH1D * hf = h1->ProjectionY("8.2-8.4", h1->GetXaxis()->FindBin(8.2+eps), h1->GetXaxis()->FindBin(8.4-eps));
    hf->SetTitle("E_{#gamma}: [8.2,8.4]GeV;#theta_{J/#psi} (deg);Counts in 47 days");
    TH1D * hg = h1->ProjectionY("8.4-8.6", h1->GetXaxis()->FindBin(8.4+eps), h1->GetXaxis()->FindBin(8.6-eps));
    hg->SetTitle("E_{#gamma} > 8.4GeV;#theta_{J/#psi} (deg);Counts in 47 days");
    
    
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);

    SetStyle(ha);
    SetStyle(hb);
    SetStyle(hc);
    SetStyle(hd);
    SetStyle(he);
    SetStyle(hf);
    SetStyle(hg);

 
    ha->DrawClone("");
    c0->Print("figures202009-E8.5/gamma-thetajpsi.pdf(", "pdf");

    hb->DrawClone("");
    c0->Print("figures202009-E8.5/gamma-thetajpsi.pdf", "pdf");

    hc->DrawClone("");
    c0->Print("figures202009-E8.5/gamma-thetajpsi.pdf", "pdf");

    hd->DrawClone("");
    c0->Print("figures202009-E8.5/gamma-thetajpsi.pdf", "pdf");

    he->DrawClone("");
    c0->Print("figures202009-E8.5/gamma-thetajpsi.pdf", "pdf");
    
    hf->DrawClone("");
    c0->Print("figures202009-E8.5/gamma-thetajpsi.pdf", "pdf");

    hg->DrawClone("");
    c0->Print("figures202009-E8.5/gamma-thetajpsi.pdf)", "pdf");
  }

  
  if (opt == 4){
    double eps = 1e-10;
    TH2D * h1 = (TH2D *) fs1->Get("EgPp");
    TH2D * h2 = (TH2D *) fs2->Get("EgPp");
    h1->Add(h2);
    h1->Scale(47*24*0.8);

    TH1D * ha = h1->ProjectionY("7.2-7.4", h1->GetXaxis()->FindBin(7.2+eps), h1->GetXaxis()->FindBin(7.4-eps));
    ha->SetTitle("E_{#gamma}: [7.2,7.4]GeV;P_{p*} (GeV);Counts in 47 days");
    TH1D * hb = h1->ProjectionY("7.4-7.6", h1->GetXaxis()->FindBin(7.4+eps), h1->GetXaxis()->FindBin(7.6-eps));
    hb->SetTitle("E_{#gamma}: [7.4,7.6]GeV;P_{p*} (GeV);Counts in 47 days");
    TH1D * hc = h1->ProjectionY("7.6-7.8", h1->GetXaxis()->FindBin(7.6+eps), h1->GetXaxis()->FindBin(7.8-eps));
    hc->SetTitle("E_{#gamma}: [7.6,7.8]GeV;P_{p*} (GeV);Counts in 47 days");
    TH1D * hd = h1->ProjectionY("7.8-8.0", h1->GetXaxis()->FindBin(7.8+eps), h1->GetXaxis()->FindBin(8.0-eps));
    hd->SetTitle("E_{#gamma}: [7.8,8.0]GeV;P_{p*} (GeV);Counts in 47 days");
    TH1D * he = h1->ProjectionY("8.0-8.2", h1->GetXaxis()->FindBin(8.0+eps), h1->GetXaxis()->FindBin(8.2-eps));
    he->SetTitle("E_{#gamma}: [8.0,8.2]GeV;P_{p*} (GeV);Counts in 47 days");
    TH1D * hf = h1->ProjectionY("8.2-8.4", h1->GetXaxis()->FindBin(8.2+eps), h1->GetXaxis()->FindBin(8.4-eps));
    hf->SetTitle("E_{#gamma}: [8.2,8.4]GeV;P_{p*} (GeV);Counts in 47 days");
    TH1D * hg = h1->ProjectionY("8.4-8.6", h1->GetXaxis()->FindBin(8.4+eps), h1->GetXaxis()->FindBin(8.6-eps));
    hg->SetTitle("E_{#gamma} > 8.4GeV;P_{p*} (GeV);Counts in 47 days");
    
    
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);

    SetStyle(ha);
    SetStyle(hb);
    SetStyle(hc);
    SetStyle(hd);
    SetStyle(he);
    SetStyle(hf);
    SetStyle(hg);

 
    ha->DrawClone("");
    c0->Print("figures202009-E8.5/gamma-pp.pdf(", "pdf");

    hb->DrawClone("");
    c0->Print("figures202009-E8.5/gamma-pp.pdf", "pdf");

    hc->DrawClone("");
    c0->Print("figures202009-E8.5/gamma-pp.pdf", "pdf");

    hd->DrawClone("");
    c0->Print("figures202009-E8.5/gamma-pp.pdf", "pdf");

    he->DrawClone("");
    c0->Print("figures202009-E8.5/gamma-pp.pdf", "pdf");
    
    hf->DrawClone("");
    c0->Print("figures202009-E8.5/gamma-pp.pdf", "pdf");

    hg->DrawClone("");
    c0->Print("figures202009-E8.5/gamma-pp.pdf)", "pdf");
  }

  if (opt == 5){
    double eps = 1e-10;
    TH2D * h1 = (TH2D *) fs1->Get("Egkappa");
    TH2D * h2 = (TH2D *) fs2->Get("Egkappa");
    h1->Add(h2);
    h1->Scale(47*24*0.8);

    TH1D * ha = h1->ProjectionY("7.2-7.4", h1->GetXaxis()->FindBin(7.2+eps), h1->GetXaxis()->FindBin(7.4-eps));
    ha->SetTitle("E_{#gamma}: [7.2,7.4]GeV;#kappa_{J/#psi} (GeV);Counts in 47 days");
    TH1D * hb = h1->ProjectionY("7.4-7.6", h1->GetXaxis()->FindBin(7.4+eps), h1->GetXaxis()->FindBin(7.6-eps));
    hb->SetTitle("E_{#gamma}: [7.4,7.6]GeV;#kappa_{J/#psi} (GeV);Counts in 47 days");
    TH1D * hc = h1->ProjectionY("7.6-7.8", h1->GetXaxis()->FindBin(7.6+eps), h1->GetXaxis()->FindBin(7.8-eps));
    hc->SetTitle("E_{#gamma}: [7.6,7.8]GeV;#kappa_{J/#psi} (GeV);Counts in 47 days");
    TH1D * hd = h1->ProjectionY("7.8-8.0", h1->GetXaxis()->FindBin(7.8+eps), h1->GetXaxis()->FindBin(8.0-eps));
    hd->SetTitle("E_{#gamma}: [7.8,8.0]GeV;#kappa_{J/#psi} (GeV);Counts in 47 days");
    TH1D * he = h1->ProjectionY("8.0-8.2", h1->GetXaxis()->FindBin(8.0+eps), h1->GetXaxis()->FindBin(8.2-eps));
    he->SetTitle("E_{#gamma}: [8.0,8.2]GeV;#kappa_{J/#psi} (GeV);Counts in 47 days");
    TH1D * hf = h1->ProjectionY("8.2-8.4", h1->GetXaxis()->FindBin(8.2+eps), h1->GetXaxis()->FindBin(8.4-eps));
    hf->SetTitle("E_{#gamma}: [8.2,8.4]GeV;#kappa_{J/#psi} (GeV);Counts in 47 days");
    TH1D * hg = h1->ProjectionY("8.4-8.6", h1->GetXaxis()->FindBin(8.4+eps), h1->GetXaxis()->FindBin(8.6-eps));
    hg->SetTitle("E_{#gamma} > 8.4GeV;#kappa_{J/#psi} (GeV);Counts in 47 days");
    
    
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);

    SetStyle(ha);
    SetStyle(hb);
    SetStyle(hc);
    SetStyle(hd);
    SetStyle(he);
    SetStyle(hf);
    SetStyle(hg);

 
    ha->DrawClone("");
    c0->Print("figures202009-E8.5/gamma-kappa.pdf(", "pdf");

    hb->DrawClone("");
    c0->Print("figures202009-E8.5/gamma-kappa.pdf", "pdf");

    hc->DrawClone("");
    c0->Print("figures202009-E8.5/gamma-kappa.pdf", "pdf");

    hd->DrawClone("");
    c0->Print("figures202009-E8.5/gamma-kappa.pdf", "pdf");

    he->DrawClone("");
    c0->Print("figures202009-E8.5/gamma-kappa.pdf", "pdf");
    
    hf->DrawClone("");
    c0->Print("figures202009-E8.5/gamma-kappa.pdf", "pdf");

    hg->DrawClone("");
    c0->Print("figures202009-E8.5/gamma-kappa.pdf)", "pdf");
  }

    
  return 0;
}

  
