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
    cout << "./plot202005 <opt>" << endl;
    cout << "1: mass" << endl;
    return 1;
  }

  const int opt = atoi(argv[1]);

  gStyle->SetPalette(55);
  gStyle->SetOptStat(0);

  if (opt == 1){// Mass spectrum of e+e- system
    TFile * fs1 = new TFile("result-photon/Dsolid-4dim-202005.root", "r");//photo-production
    TFile * fs2 = new TFile("result-electro/Dsolid-4dim-202005.root", "r");//elecro-production
    TFile * fs3 = new TFile("result-photon/Dsolid-4dim-smear-202005.root", "r");//photo-production smeared
    TFile * fs4 = new TFile("result-electro/Dsolid-4dim-smear-202005.root", "r");//elecro-production smeared
    TH1D * h1 = (TH1D *) fs1->Get("Mass_e+e-_Jpsi");
    TH1D * h2 = (TH1D *) fs2->Get("Mass_e+e-_Jpsi");
    TH1D * h3 = (TH1D *) fs3->Get("Mass_e+e-_Jpsi");
    TH1D * h4 = (TH1D *) fs4->Get("Mass_e+e-_Jpsi");
    
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);

    SetStyle(h1);
    SetStyle(h2);
    SetStyle(h3);
    SetStyle(h4);

    h1->SetLineColor(1);
    h2->SetLineColor(1);
    h3->SetLineColor(4);
    h4->SetLineColor(4);

    h1->SetMaximum(5.0);
    h3->SetMaximum(5.0);

    h1->DrawClone("");
    h3->DrawClone("same");

    c0->Print("figures202005/mass.pdf(", "pdf");

    h2->DrawClone("");
    h4->DrawClone("same");

    c0->Print("figures202005/mass.pdf)", "pdf");

    
  }

  if (opt == 2){//momentum-theta distribution
    TFile * fs1 = new TFile("result-photon/Dsolid-4dim-202005.root", "r");//photo-production
    TFile * fs2 = new TFile("result-electro/Dsolid-4dim-202005.root", "r");//elecro-production
    TFile * fs3 = new TFile("result-photon/Dsolid-4dim-smear-202005.root", "r");//photo-production smeared
    TFile * fs4 = new TFile("result-electro/Dsolid-4dim-smear-202005.root", "r");//elecro-production smeared
    TH2D * h1 = (TH2D *) fs1->Get("ThetaP_e-_Jpsi");
    TH2D * h2 = (TH2D *) fs2->Get("ThetaP_e-_Jpsi");
    TH2D * h3 = (TH2D *) fs3->Get("ThetaP_e-_Jpsi");
    TH2D * h4 = (TH2D *) fs4->Get("ThetaP_e-_Jpsi");

    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);
    c0->SetLogz();

    SetStyle(h1);
    SetStyle(h2);
    SetStyle(h3);
    SetStyle(h4);
    
    h1->GetXaxis()->SetRangeUser(0,50);
    h1->GetYaxis()->SetRangeUser(0,7);
    h1->SetMinimum(1e-4);
    h1->SetMaximum(5e-1);
    h1->DrawClone("colz");
    c0->Print("figures202005/momentum-theta.pdf(", "pdf");

    h2->GetXaxis()->SetRangeUser(0,50);
    h2->GetYaxis()->SetRangeUser(0,7);
    h2->SetMinimum(1e-4);
    h2->SetMaximum(5e-1);
    h2->DrawClone("colz");
    c0->Print("figures202005/momentum-theta.pdf", "pdf");

    h3->GetXaxis()->SetRangeUser(0,50);
    h3->GetYaxis()->SetRangeUser(0,7);
    h3->SetMinimum(1e-4);
    h3->SetMaximum(5e-1);
    h3->DrawClone("colz");
    c0->Print("figures202005/momentum-theta.pdf", "pdf");

    h4->GetXaxis()->SetRangeUser(0,50);
    h4->GetYaxis()->SetRangeUser(0,7);
    h4->SetMinimum(1e-4);
    h4->SetMaximum(5e-1);
    h4->DrawClone("colz");
    c0->Print("figures202005/momentum-theta.pdf)", "pdf");

    
  }

  if (opt == 3){//jpsi distribution
    TFile * fs1 = new TFile("result-photon/Dsolid-4dim-202005.root", "r");//photo-production
    TFile * fs2 = new TFile("result-electro/Dsolid-4dim-202005.root", "r");//elecro-production
    TFile * fs3 = new TFile("result-photon/Dsolid-4dim-smear-202005.root", "r");//photo-production smeared
    TFile * fs4 = new TFile("result-electro/Dsolid-4dim-smear-202005.root", "r");//elecro-production smeared
    TH2D * h1 = (TH2D *) fs1->Get("ThetaPJpsi");
    TH2D * h2 = (TH2D *) fs2->Get("ThetaPJpsi");
    TH2D * h3 = (TH2D *) fs3->Get("ThetaPJpsi");
    TH2D * h4 = (TH2D *) fs4->Get("ThetaPJpsi");

    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);
    c0->SetLogz();

    SetStyle(h1);
    SetStyle(h2);
    SetStyle(h3);
    SetStyle(h4);
    
    h1->GetXaxis()->SetRangeUser(0,12);
    h1->GetYaxis()->SetRangeUser(0,10);
    h1->SetMinimum(1e-4);
    h1->SetMaximum(1e0);
    h1->DrawClone("colz");
    c0->Print("figures202005/jpsi-distribution.pdf(", "pdf");

    h2->GetXaxis()->SetRangeUser(0,12);
    h2->GetYaxis()->SetRangeUser(0,10);
    h2->SetMinimum(1e-4);
    h2->SetMaximum(1e0);
    h2->DrawClone("colz");
    c0->Print("figures202005/jpsi-distribution.pdf", "pdf");

    h3->GetXaxis()->SetRangeUser(0,12);
    h3->GetYaxis()->SetRangeUser(0,10);
    h3->SetMinimum(1e-4);
    h3->SetMaximum(1e0);
    h3->DrawClone("colz");
    c0->Print("figures202005/jpsi-distribution.pdf", "pdf");

    h4->GetXaxis()->SetRangeUser(0,12);
    h4->GetYaxis()->SetRangeUser(0,10);
    h4->SetMinimum(1e-4);
    h4->SetMaximum(1e0);
    h4->DrawClone("colz");
    c0->Print("figures202005/jpsi-distribution.pdf)", "pdf");
  }

  if (opt == 4){//internal momentum kappa
    TFile * fs1 = new TFile("result-photon/Dsolid-4dim-202005.root", "r");//photo-production
    TFile * fs2 = new TFile("result-electro/Dsolid-4dim-202005.root", "r");//elecro-production
    TFile * fs3 = new TFile("result-photon/Dsolid-4dim-smear-202005.root", "r");//photo-production smeared
    TFile * fs4 = new TFile("result-electro/Dsolid-4dim-smear-202005.root", "r");//elecro-production smeared
    TH1D * h1 = (TH1D *) fs1->Get("kappa_Jpsi");
    TH1D * h2 = (TH1D *) fs2->Get("kappa_Jpsi");
    TH1D * h3 = (TH1D *) fs3->Get("kappa_Jpsi");
    TH1D * h4 = (TH1D *) fs4->Get("kappa_Jpsi");

    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);

    SetStyle(h1);
    SetStyle(h2);
    SetStyle(h3);
    SetStyle(h4);

    h1->SetLineColor(1);
    h2->SetLineColor(1);
    h3->SetLineColor(4);
    h4->SetLineColor(4);

    h1->SetMaximum(0.15);
    h3->SetMaximum(0.15);

    h1->DrawClone("");
    h3->DrawClone("same");

    c0->Print("figures202005/kappa.pdf(", "pdf");

    h2->DrawClone("");
    h4->DrawClone("same");

    c0->Print("figures202005/kappa.pdf)", "pdf");

    
  }
    
  return 0;
}

  
