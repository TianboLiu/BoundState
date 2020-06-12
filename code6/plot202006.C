/* plots for proposal 2020/05 */

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
    cout << "./plot202006 <opt>" << endl;
    cout << "1: mass" << endl;
    return 1;
  }

  const int opt = atoi(argv[1]);

  gStyle->SetPalette(55);
  gStyle->SetOptStat(0);

  if (opt == 0){// rates
    TFile * fs1 = new TFile("result202006/photo-default-smeared.root");
    TFile * fs2 = new TFile("result202006/electro-default-smeared.root");
  
    TH1D * f1 = (TH1D *) fs1->Get("Mass_e+e-");
    TH1D * f2 = (TH1D *) fs2->Get("Mass_e+e-");
    
    int binA = f1->FindBin(3.0969 - 1.06);
    int binB = f1->FindBin(3.0969 + 1.06);

    cout << f1->Integral(binA, binB) * 0.8 << endl;
    cout << f2->Integral(binA, binB) * 0.8 << endl;
  }
     
  
  if (opt == -1){// Mass spectrum of e+e-, photo-production
    TFile * fs1 = new TFile("result202006/photo-default-detected.root", "r");
    TFile * fs2 = new TFile("result202006/photo-default-smeared.root", "r");
    TH1D * h1 = (TH1D *) fs1->Get("Mass_e+e-");
    TH1D * h2 = (TH1D *) fs2->Get("Mass_e+e-");
  
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);

    SetStyle(h1);
    SetStyle(h2);

    h1->SetLineColor(1);
    h2->SetLineColor(2);

    h1->DrawClone("");
    h2->DrawClone("same");

    c0->Print("figures202006/mass-photon.pdf");
  }

  if (opt == -2){// Momentum distribution, photo-production
    TFile * fs = new TFile("result202006/photo-default-smeared.root", "r");
    TH2D * h1 = (TH2D *) fs->Get("ThetaP_e+");
    TH2D * h2 = (TH2D *) fs->Get("ThetaP_p");
  
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);
    c0->SetRightMargin(0.15);

    //c0->SetLogz();

    SetStyle(h1);
    SetStyle(h2);

    h1->DrawClone("colz");
    
    c0->Print("figures202006/thetap-photon.pdf(", "pdf");

    h2->DrawClone("colz");

    c0->Print("figures202006/thetap-photon.pdf)", "pdf");

  }

  if (opt == -3){// Jpsi momentum, photo-production
    TFile * fs = new TFile("result202006/photo-default-smeared.root", "r");
    TH2D * h1 = (TH2D *) fs->Get("ThetaP_Jpsi");
    TH2D * h2 = (TH2D *) fs->Get("Angle");
  
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);
    c0->SetRightMargin(0.15);

    //c0->SetLogz();

    SetStyle(h1);
    SetStyle(h2);

    h1->DrawClone("colz");

    c0->Print("figures202006/jpsi-photon.pdf(", "pdf");

    h2->DrawClone("colz");
    c0->Print("figures202006/jpsi-photon.pdf)", "pdf");

  }

  if (opt == -4){// kappa, photo-production
    TFile * fs = new TFile("result202006/photo-default-smeared.root", "r");
    TH1D * h1 = (TH1D *) fs->Get("kappa");
    TH2D * h2 = (TH2D *) fs->Get("kappaP");
    TH2D * h3 = (TH2D *) fs->Get("kappaTheta");
 
    TH1D * h1a = h2->ProjectionX("kappa_lowP", 1, h2->GetYaxis()->FindBin(1.0));
    TH1D * h1b = h2->ProjectionX("kappa_highP", h2->GetYaxis()->FindBin(1.0), -1);

    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);

    SetStyle(h1);
    SetStyle(h2);
    SetStyle(h3);

    h1->SetLineColor(1);
    h1a->SetLineColor(4);
    h1b->SetLineColor(2);
    h1->DrawClone();
    h1a->DrawClone("same");
    h1b->DrawClone("same");

    c0->Print("figures202006/kappa-photon.pdf(", "pdf");

    h2->DrawClone("colz");

    c0->SetRightMargin(0.15);
    //c0->SetLogz();

    c0->Print("figures202006/kappa-photon.pdf", "pdf");

    h3->DrawClone("colz");

    c0->Print("figures202006/kappa-photon.pdf)", "pdf");

  }

  if (opt == 1){// Mass spectrum of e+e-, photo-production
    TFile * fs1 = new TFile("result202006/electro-default-detected.root", "r");
    TFile * fs2 = new TFile("result202006/electro-default-smeared.root", "r");
    TH1D * h1 = (TH1D *) fs1->Get("Mass_e+e-");
    TH1D * h2 = (TH1D *) fs2->Get("Mass_e+e-");
   
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);

    SetStyle(h1);
    SetStyle(h2);

    h1->SetLineColor(1);
    h2->SetLineColor(2);

    h1->DrawClone("");
    h2->DrawClone("same");

    c0->Print("figures202006/mass-quasi.pdf");

  }

  if (opt == 2){// Momentum distribution, photo-production
    TFile * fs = new TFile("result202006/electro-default-smeared.root", "r");
    TH2D * h1 = (TH2D *) fs->Get("ThetaP_e+");
    TH2D * h2 = (TH2D *) fs->Get("ThetaP_p");

    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);
    c0->SetRightMargin(0.15);

    //c0->SetLogz();

    SetStyle(h1);
    SetStyle(h2);

    h1->DrawClone("colz");
    
    c0->Print("figures202006/thetap-quasi.pdf(", "pdf");

    h2->DrawClone("colz");

    c0->Print("figures202006/thetap-quasi.pdf)", "pdf");

  }

  if (opt == 3){// Jpsi momentum, photo-production
    TFile * fs = new TFile("result202006/electro-default-smeared.root", "r");
    TH2D * h1 = (TH2D *) fs->Get("ThetaP_Jpsi");
    TH2D * h2 = (TH2D *) fs->Get("Angle");

    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);
    c0->SetRightMargin(0.15);

    //c0->SetLogz();

    SetStyle(h1);
    SetStyle(h2);

    h1->DrawClone("colz");

    c0->Print("figures202006/jpsi-quasi.pdf(", "pdf");

    h2->DrawClone("colz");
    c0->Print("figures202006/jpsi-quasi.pdf)", "pdf");

  }

  if (opt == 4){// kappa, photo-production
    TFile * fs = new TFile("result202006/electro-default-smeared.root", "r");
    TH1D * h1 = (TH1D *) fs->Get("kappa");
    TH2D * h2 = (TH2D *) fs->Get("kappaP");
    TH2D * h3 = (TH2D *) fs->Get("kappaTheta");

    TH1D * h1a = h2->ProjectionX("kappa_lowP", 1, h2->GetYaxis()->FindBin(1.0));
    TH1D * h1b = h2->ProjectionX("kappa_highP", h2->GetYaxis()->FindBin(1.0), -1);

    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);

    SetStyle(h1);
    SetStyle(h2);
    SetStyle(h3);

    h1->SetLineColor(1);
    h1a->SetLineColor(4);
    h1b->SetLineColor(2);
    h1->DrawClone();
    h1a->DrawClone("same");
    h1b->DrawClone("same");

    c0->Print("figures202006/kappa-quasi.pdf(", "pdf");

    h2->DrawClone("colz");

    c0->SetRightMargin(0.15);
    //c0->SetLogz();

    c0->Print("figures202006/kappa-quasi.pdf", "pdf");

    h3->DrawClone("colz");

    c0->Print("figures202006/kappa-quasi.pdf)", "pdf");
  }

  if (opt == 5){
    TFile * fs1 = new TFile("result202006/photo-default-smeared.root", "r");
    TH2D * h1 = (TH2D *) fs1->Get("kappaP");
    TH1D * h1a = h1->ProjectionX("kappa_highP", h1->GetYaxis()->FindBin(1.0), -1);
       
    TFile * fs2 = new TFile("result202006/electro-default-smeared.root", "r");
    TH2D * h2 = (TH2D *) fs2->Get("kappaP");
    TH1D * h2a = h2->ProjectionX("kappa_highP", h2->GetYaxis()->FindBin(1.0), -1);

    h1a->Scale(50*24*0.8);
    h2a->Scale(50*24*0.8);

    h1a->Rebin(10);
    h2a->Rebin(10);
    h2a->SetMinimum(0);

    for (int i = 1; i <= 10; i++){
      h1a->SetBinError(i, sqrt(h1a->GetBinContent(i)));
      h2a->SetBinError(i, sqrt(h2a->GetBinContent(i)));
    }
      
    
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);

    SetStyle(h1a);
    SetStyle(h2a);

    h2a->GetYaxis()->SetTitle("Counts in 50 days");

    h1a->SetLineColor(2);
    h2a->SetLineColor(4);

    h1a->SetLineWidth(3);
    h2a->SetLineWidth(3);

    h1a->SetFillColorAlpha(2, 0.35);
    h2a->SetFillColorAlpha(4, 0.35);

    //h2a->SetLineStyle(7);

    h2a->DrawClone("pe2");
    h1a->DrawClone("pe2same");
    
    h2a->DrawClone("esame");
    h1a->DrawClone("esame");

    c0->Print("figures202006/kappaproj.pdf");

  }

  if (opt == 6){
    TFile * fs1 = new TFile("result202006/photo-default-smeared.root", "r");
    TH2D * h1 = (TH2D *) fs1->Get("kappaEgHigh");
    TH1D * h1a = h1->ProjectionX("kappaEg1", h1->GetYaxis()->FindBin(7.2), h1->GetYaxis()->FindBin(7.5));
    TH1D * h1b = h1->ProjectionX("kappaEg2", h1->GetYaxis()->FindBin(7.5), h1->GetYaxis()->FindBin(7.9));
    TH1D * h1c = h1->ProjectionX("kappaEg3", h1->GetYaxis()->FindBin(7.9), -1);
       
    TFile * fs2 = new TFile("result202006/electro-default-smeared.root", "r");
    TH2D * h2 = (TH2D *) fs2->Get("kappaEgHigh");
    TH1D * h2a = h2->ProjectionX("kappaEg1", h2->GetYaxis()->FindBin(7.2), h2->GetYaxis()->FindBin(7.5));
    TH1D * h2b = h2->ProjectionX("kappaEg2", h2->GetYaxis()->FindBin(7.5), h2->GetYaxis()->FindBin(7.9));
    TH1D * h2c = h2->ProjectionX("kappaEg3", h2->GetYaxis()->FindBin(7.9), -1);

    h1a->Scale(50*24*0.8);
    h1b->Scale(50*24*0.8);
    h1c->Scale(50*24*0.8);
    h2a->Scale(50*24*0.8);
    h2b->Scale(50*24*0.8);
    h2c->Scale(50*24*0.8);

    h1a->Rebin(20);
    h1b->Rebin(20);
    h1c->Rebin(20);
    h2a->Rebin(20);
    h2b->Rebin(20);
    h2c->Rebin(20);

    h1c->SetMinimum(0);
    h2c->SetMinimum(0);
    
    for (int i = 1; i <= 5; i++){
      h1a->SetBinError(i, sqrt(h1a->GetBinContent(i)));
      h1b->SetBinError(i, sqrt(h1b->GetBinContent(i)));
      h1c->SetBinError(i, sqrt(h1c->GetBinContent(i)));
      h2a->SetBinError(i, sqrt(h2a->GetBinContent(i)));
      h2b->SetBinError(i, sqrt(h2b->GetBinContent(i)));
      h2c->SetBinError(i, sqrt(h2c->GetBinContent(i)));
    }
      
    
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);
    
    SetStyle(h1a);
    SetStyle(h1b);
    SetStyle(h1c);
    SetStyle(h2a);
    SetStyle(h2b);
    SetStyle(h2c);

    h1a->GetYaxis()->SetTitle("Counts in 50 days");
    h2a->GetYaxis()->SetTitle("Counts in 50 days");
    h1b->GetYaxis()->SetTitle("Counts in 50 days");
    h2b->GetYaxis()->SetTitle("Counts in 50 days");
    h1c->GetYaxis()->SetTitle("Counts in 50 days");
    h2c->GetYaxis()->SetTitle("Counts in 50 days");

    h1a->SetLineColor(4);
    h1b->SetLineColor(kGreen+2);
    h1c->SetLineColor(2);
    h2a->SetLineColor(4);
    h2b->SetLineColor(kGreen+2);
    h2c->SetLineColor(2);

    h1a->SetLineWidth(3);
    h1b->SetLineWidth(3);
    h1c->SetLineWidth(3);
    h2a->SetLineWidth(3);
    h2b->SetLineWidth(3);
    h2c->SetLineWidth(3);

    h1a->SetFillColorAlpha(4, 0.35);
    h1b->SetFillColorAlpha(kGreen+2, 0.35);
    h1c->SetFillColorAlpha(2, 0.35);
    h2a->SetFillColorAlpha(4, 0.35);
    h2b->SetFillColorAlpha(kGreen+2, 0.35);
    h2c->SetFillColorAlpha(2, 0.35);

    h1c->DrawClone("pe2");
    h1b->DrawClone("pe2same");
    h1a->DrawClone("pe2same");
    h1a->DrawClone("esame");
    h1b->DrawClone("esame");
    h1c->DrawClone("esame");
    c0->Print("figures202006/kappagamma-photon.pdf");

    h2c->DrawClone("pe2");
    h2b->DrawClone("pe2same");
    h2a->DrawClone("pe2same");
    h2a->DrawClone("esame");
    h2b->DrawClone("esame");
    h2c->DrawClone("esame");
    c0->Print("figures202006/kappagamma-electro.pdf");
  
  }

  if (opt == 7){
    TFile * fs1a = new TFile("result202006/photo-default-smeared.root", "r");
    TFile * fs1b = new TFile("result202006/photo-nv2Ia-smeared.root", "r");
    TFile * fs1c = new TFile("result202006/photo-v18-smeared.root", "r");
    TH2D * h1A = (TH2D *) fs1a->Get("kappaEgHigh");
    TH2D * h1B = (TH2D *) fs1b->Get("kappaEgHigh");
    TH2D * h1C = (TH2D *) fs1c->Get("kappaEgHigh");
    TH1D * h1a = h1A->ProjectionX("kappadefault");
    TH1D * h1b = h1B->ProjectionX("kappanv2Ia");
    TH1D * h1c = h1C->ProjectionX("kappav18");

    TFile * fs2a = new TFile("result202006/electro-default-smeared.root", "r");
    TFile * fs2b = new TFile("result202006/electro-nv2Ia-smeared.root", "r");
    TFile * fs2c = new TFile("result202006/electro-v18-smeared.root", "r");
    TH2D * h2A = (TH2D *) fs2a->Get("kappaEgHigh");
    TH2D * h2B = (TH2D *) fs2b->Get("kappaEgHigh");
    TH2D * h2C = (TH2D *) fs2c->Get("kappaEgHigh");
    TH1D * h2a = h2A->ProjectionX("kappadefault");
    TH1D * h2b = h2B->ProjectionX("kappanv2Ia");
    TH1D * h2c = h2C->ProjectionX("kappav18");

    h1a->Scale(50*24*0.8);
    h1b->Scale(50*24*0.8);
    h1c->Scale(50*24*0.8);
    h2a->Scale(50*24*0.8);
    h2b->Scale(50*24*0.8);
    h2c->Scale(50*24*0.8);

    h1a->Rebin(10);
    h1b->Rebin(10);
    h1c->Rebin(10);
    h2a->Rebin(10);
    h2b->Rebin(10);
    h2c->Rebin(10);


    for (int i = 1; i <= 10; i++){
      h1a->SetBinError(i, sqrt(h1a->GetBinContent(i)));
      h1b->SetBinError(i, sqrt(h1b->GetBinContent(i)));
      h1c->SetBinError(i, sqrt(h1c->GetBinContent(i)));
      h2a->SetBinError(i, sqrt(h2a->GetBinContent(i)));
      h2b->SetBinError(i, sqrt(h2b->GetBinContent(i)));
      h2c->SetBinError(i, sqrt(h2c->GetBinContent(i)));
    }
      
    
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);

    SetStyle(h1a);
    SetStyle(h1b);
    SetStyle(h1c);
    SetStyle(h2a);
    SetStyle(h2b);
    SetStyle(h2c);

    h1a->GetYaxis()->SetTitle("Counts in 50 days");
    h2a->GetYaxis()->SetTitle("Counts in 50 days");
    h1b->GetYaxis()->SetTitle("Counts in 50 days");
    h2b->GetYaxis()->SetTitle("Counts in 50 days");
    h1c->GetYaxis()->SetTitle("Counts in 50 days");
    h2c->GetYaxis()->SetTitle("Counts in 50 days");    

    h1a->SetLineColor(1);
    h1b->SetLineColor(4);
    h1c->SetLineColor(2);
    h2a->SetLineColor(1);
    h2b->SetLineColor(4);
    h2c->SetLineColor(2);

    h1a->SetLineWidth(3);
    h1b->SetLineWidth(3);
    h1c->SetLineWidth(3);
    h2a->SetLineWidth(3);
    h2b->SetLineWidth(3);
    h2c->SetLineWidth(3);

    h1a->SetFillColorAlpha(1, 0.35);
    h1b->SetFillColorAlpha(4, 0.35);
    h1c->SetFillColorAlpha(2, 0.35);
    h2a->SetFillColorAlpha(1, 0.35);
    h2b->SetFillColorAlpha(4, 0.35);
    h2c->SetFillColorAlpha(2, 0.35);


    h1b->SetMaximum(130);
    h1b->DrawClone("pe2");
    h1c->DrawClone("pe2same");
    h1a->DrawClone("pe2same");
    h1a->DrawClone("esame");
    h1b->DrawClone("esame");
    h1c->DrawClone("esame");
    c0->Print("figures202006/models-photon.pdf");

    h2b->SetMaximum(130);
    h2b->DrawClone("pe2");
    h2c->DrawClone("pe2same");
    h2a->DrawClone("pe2same");
    h2a->DrawClone("esame");
    h2b->DrawClone("esame");
    h2c->DrawClone("esame");
    c0->Print("figures202006/models-electro.pdf");
  
  }

    
    
  return 0;
}

  
