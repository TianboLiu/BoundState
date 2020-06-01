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

  if (opt == 0){// rates
    TFile * fs1 = new TFile("result202005/photo-smeared.root");
    TFile * fs2a = new TFile("result202005/quasi-0-5-smeared.root");
    TFile * fs2b = new TFile("result202005/quasi-5-15-smeared.root");
    TFile * fs2c = new TFile("result202005/quasi-5-15-smeared.root");
    TFile * fs3a = new TFile("result202005/electro-0-5-smeared.root");
    TFile * fs3b = new TFile("result202005/electro-5-15-smeared.root");
    TFile * fs3c = new TFile("result202005/electro-5-15-smeared.root");

    TH1D * f1 = (TH1D *) fs1->Get("Mass_e+e-");
    TH1D * f2 = (TH1D *) fs2a->Get("Mass_e+e-");
    TH1D * f3 = (TH1D *) fs3a->Get("Mass_e+e-");

    f2->Add( (TH1D *) fs2b->Get("Mass_e+e-"));
    f2->Add( (TH1D *) fs2c->Get("Mass_e+e-"));
    f3->Add( (TH1D *) fs3b->Get("Mass_e+e-"));
    f3->Add( (TH1D *) fs3c->Get("Mass_e+e-"));

    f1->Scale(0.5);
    f2->Scale(0.5);
    f3->Scale(0.5);

    cout << f1->Integral() << endl;
    cout << f2->Integral() << endl;
    cout << f3->Integral() << endl;
  }
   
	     
    
    
  
  if (opt == -1){// Mass spectrum of e+e-, photo-production
    TFile * fs1 = new TFile("result202005/photo-detected.root", "r");
    TFile * fs2 = new TFile("result202005/photo-smeared.root", "r");
    TH1D * h1 = (TH1D *) fs1->Get("Mass_e+e-");
    TH1D * h2 = (TH1D *) fs2->Get("Mass_e+e-");
    h1->Scale(0.5);
    h2->Scale(0.5);

    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);

    SetStyle(h1);
    SetStyle(h2);

    h1->SetLineColor(1);
    h2->SetLineColor(2);

    h1->DrawClone("");
    h2->DrawClone("same");

    c0->Print("figures202005/mass-photon.pdf");

  }

  if (opt == -2){// Momentum distribution, photo-production
    TFile * fs = new TFile("result202005/photo-smeared.root", "r");
    TH2D * h1 = (TH2D *) fs->Get("ThetaP_e+");
    TH2D * h2 = (TH2D *) fs->Get("ThetaP_p");
    h1->Scale(0.5);
    h2->Scale(0.5);

    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);

    //c0->SetLogz();

    SetStyle(h1);
    SetStyle(h2);

    h1->DrawClone("colz");
    
    c0->Print("figures202005/thetap-photon.pdf(", "pdf");

    h2->DrawClone("colz");

    c0->Print("figures202005/thetap-photon.pdf)", "pdf");

  }

  if (opt == -3){// Jpsi momentum, photo-production
    TFile * fs = new TFile("result202005/photo-smeared.root", "r");
    TH2D * h1 = (TH2D *) fs->Get("ThetaP_Jpsi");
    TH2D * h2 = (TH2D *) fs->Get("Angle");
    h1->Scale(0.5);
    h2->Scale(0.5);

    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);

    //c0->SetLogz();

    SetStyle(h1);
    SetStyle(h2);

    h1->DrawClone("colz");

    c0->Print("figures202005/jpsi-photon.pdf(", "pdf");

    h2->DrawClone("colz");
    c0->Print("figures202005/jpsi-photon.pdf)", "pdf");

  }

  if (opt == -4){// kappa, photo-production
    TFile * fs = new TFile("result202005/photo-smeared.root", "r");
    TH1D * h1 = (TH1D *) fs->Get("kappa");
    TH2D * h2 = (TH2D *) fs->Get("kappaP");
    TH2D * h3 = (TH2D *) fs->Get("kappaTheta");
    h1->Scale(0.5);
    h2->Scale(0.5);
    h3->Scale(0.5);

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

    c0->Print("figures202005/kappa-photon.pdf(", "pdf");

    h2->DrawClone("colz");

    c0->SetRightMargin(0.15);
    //c0->SetLogz();

    c0->Print("figures202005/kappa-photon.pdf", "pdf");

    h3->DrawClone("colz");

    c0->Print("figures202005/kappa-photon.pdf)", "pdf");

  }

  if (opt == 1){// Mass spectrum of e+e-, photo-production
    TFile * fs1a = new TFile("result202005/quasi-0-5-detected.root", "r");
    TFile * fs2a = new TFile("result202005/quasi-0-5-smeared.root", "r");
    TFile * fs1b = new TFile("result202005/quasi-5-15-detected.root", "r");
    TFile * fs2b = new TFile("result202005/quasi-5-15-smeared.root", "r");
    TFile * fs1c = new TFile("result202005/quasi-15-40-detected.root", "r");
    TFile * fs2c = new TFile("result202005/quasi-15-40-smeared.root", "r");
    TH1D * h1 = (TH1D *) fs1a->Get("Mass_e+e-");
    TH1D * h2 = (TH1D *) fs2a->Get("Mass_e+e-");
    h1->Add( (TH1D *) fs1b->Get("Mass_e+e-"));
    h1->Add( (TH1D *) fs1c->Get("Mass_e+e-"));
    h2->Add( (TH1D *) fs2b->Get("Mass_e+e-"));
    h2->Add( (TH1D *) fs2c->Get("Mass_e+e-"));
   
    h1->Scale(0.5);
    h2->Scale(0.5);

    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);

    SetStyle(h1);
    SetStyle(h2);

    h1->SetLineColor(1);
    h2->SetLineColor(2);

    h1->DrawClone("");
    h2->DrawClone("same");

    c0->Print("figures202005/mass-quasi.pdf");

  }

  if (opt == 2){// Momentum distribution, photo-production
    TFile * fsa = new TFile("result202005/quasi-0-5-smeared.root", "r");
    TFile * fsb = new TFile("result202005/quasi-0-5-smeared.root", "r");
    TFile * fsc = new TFile("result202005/quasi-0-5-smeared.root", "r");
    TH2D * h1 = (TH2D *) fsa->Get("ThetaP_e+");
    TH2D * h2 = (TH2D *) fsa->Get("ThetaP_p");

    h1->Add( (TH2D *) fsb->Get("ThetaP_e+"));
    h1->Add( (TH2D *) fsc->Get("ThetaP_e+"));
    h2->Add( (TH2D *) fsb->Get("ThetaP_p"));
    h2->Add( (TH2D *) fsc->Get("ThetaP_p"));	     
    
    h1->Scale(0.5);
    h2->Scale(0.5);

    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);

    //c0->SetLogz();

    SetStyle(h1);
    SetStyle(h2);

    h1->DrawClone("colz");
    
    c0->Print("figures202005/thetap-quasi.pdf(", "pdf");

    h2->DrawClone("colz");

    c0->Print("figures202005/thetap-quasi.pdf)", "pdf");

  }

  if (opt == 3){// Jpsi momentum, photo-production
    TFile * fsa = new TFile("result202005/quasi-0-5-smeared.root", "r");
    TFile * fsb = new TFile("result202005/quasi-0-5-smeared.root", "r");
    TFile * fsc = new TFile("result202005/quasi-0-5-smeared.root", "r");
    TH2D * h1 = (TH2D *) fsa->Get("ThetaP_Jpsi");
    TH2D * h2 = (TH2D *) fsa->Get("Angle");

    h1->Add( (TH2D *) fsb->Get("ThetaP_Jpsi"));
    h1->Add( (TH2D *) fsc->Get("ThetaP_Jpsi"));
    h2->Add( (TH2D *) fsb->Get("Angle"));
    h2->Add( (TH2D *) fsc->Get("Angle"));
    
    h1->Scale(0.5);
    h2->Scale(0.5);

    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);

    //c0->SetLogz();

    SetStyle(h1);
    SetStyle(h2);

    h1->DrawClone("colz");

    c0->Print("figures202005/jpsi-quasi.pdf(", "pdf");

    h2->DrawClone("colz");
    c0->Print("figures202005/jpsi-quasi.pdf)", "pdf");

  }

  if (opt == 4){// kappa, photo-production
    TFile * fsa = new TFile("result202005/quasi-0-5-smeared.root", "r");
    TFile * fsb = new TFile("result202005/quasi-0-5-smeared.root", "r");
    TFile * fsc = new TFile("result202005/quasi-0-5-smeared.root", "r");
    TH1D * h1 = (TH1D *) fsa->Get("kappa");
    TH2D * h2 = (TH2D *) fsa->Get("kappaP");
    TH2D * h3 = (TH2D *) fsa->Get("kappaTheta");

    h1->Add( (TH1D *) fsb->Get("kappa"));
    h1->Add( (TH1D *) fsc->Get("kappa"));
    h2->Add( (TH2D *) fsb->Get("kappaP"));
    h2->Add( (TH2D *) fsc->Get("kappaP"));
    h3->Add( (TH2D *) fsb->Get("kappaTheta"));
    h3->Add( (TH2D *) fsc->Get("kappaTheta"));
	     
    
    h1->Scale(0.5);
    h2->Scale(0.5);
    h3->Scale(0.5);

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

    c0->Print("figures202005/kappa-quasi.pdf(", "pdf");

    h2->DrawClone("colz");

    c0->SetRightMargin(0.15);
    //c0->SetLogz();

    c0->Print("figures202005/kappa-quasi.pdf", "pdf");

    h3->DrawClone("colz");

    c0->Print("figures202005/kappa-quasi.pdf)", "pdf");

  }

    
    
  return 0;
}

  
