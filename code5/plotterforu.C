#include "Lcore.h"

int main(const int argc, const char * argv[]){

  if (argc < 4) {
    cout << "./plotterforu <opt> <loadpath> <savepath>" << endl;
    return 1;
  }

  const int opt = atoi(argv[1]);
  gStyle->SetPalette(53);

  double lumi = 3.0e34 * 1.0e-26 * pow(Phys::hbar, 2) / GOLD::NA;//eA GeV^2 s^-1
  double time = 3600.0;

  TString loadpath = argv[2];
  TString savepath = argv[3];

  TH1D * hB = new TH1D("hB", "", 2000, 0.0, 200.0);
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

  
  if(opt == -311){//missing mass
    TFile * fs = new TFile(loadpath + "breakcuttest.root", "r");
    TFile * fs1 = new TFile(loadpath + "modifiedtwopinew.root", "r");
    TH1D * h0 = (TH1D *) fs->Get("Missing_BoundStateAll");
    TH1D * h5 = (TH1D *) fs1->Get("Missing_PiPi");
    
    double mispid = 2500.0;//10000.0;//2500.0;
    //
    h0->Scale(lumi*time/h0->GetXaxis()->GetBinWidth(1)/1000);
    h5->Scale(lumi*time/1000);
    
    h0->SetLineColor(1);h5->SetLineColor(3);
    
    h5->Scale(1.0/mispid);
    h5->Scale(0.01);
    
    TLegend * leg = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg->AddEntry(h0, "pKK from N-#phi", "l");
    leg->AddEntry(h5, "misidentified #pi#pi", "l");
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);
    hB->SetMaximum(h0->GetMaximum() * 1.1);
    hB->SetAxisRange(181.71, 184.61, "X");
    hB->GetXaxis()->SetTitle("Missing Mass(pKK) (GeV)");
    hB->GetYaxis()->SetTitle("counts / hour / MeV");
    hB->DrawClone("axis");
    h0->DrawClone("same");
    h5->DrawClone("same");
    leg->DrawClone("same");
    c0->Print(savepath + "Massbreakfinalcut.pdf","pdf");//
  }
  
  
  
if (opt == -401){//invariant mass cut
    TFile * fs = new TFile(loadpath + "breakcuttest.root", "r");
    TFile * fs1 = new TFile(loadpath + "modifiedtwopinew.root", "r");
    TH1D * h0 = (TH1D *) fs->Get("MpKK_BoundStateAll");
    TH1D * h1 = (TH1D *) fs->Get("MpKK_BoundStateKK");
    TH1D * h2 = (TH1D *) fs->Get("MpKK_phi");
    TH1D * h3 = (TH1D *) fs->Get("MpKK_Lambda1520");
    TH1D * h4 = (TH1D *) fs->Get("MpKK_directKK");
    TH1D * h5 = (TH1D *) fs1->Get("MpKK_PiPi");
    
    TH1D * h0a = (TH1D *) fs->Get("MpKp_BoundStateAll");
    TH1D * h1a = (TH1D *) fs->Get("MpKp_BoundStateKK");
    TH1D * h2a = (TH1D *) fs->Get("MpKp_phi");
    TH1D * h3a = (TH1D *) fs->Get("MpKp_Lambda1520");
    TH1D * h4a = (TH1D *) fs->Get("MpKp_directKK");
    TH1D * h5a = (TH1D *) fs1->Get("MpKp_PiPi");
 
    TH1D * h0b = (TH1D *) fs->Get("MpKm_BoundStateAll");
    TH1D * h1b = (TH1D *) fs->Get("MpKm_BoundStateKK");
    TH1D * h2b = (TH1D *) fs->Get("MpKm_phi");
    TH1D * h3b = (TH1D *) fs->Get("MpKm_Lambda1520");
    TH1D * h4b = (TH1D *) fs->Get("MpKm_directKK");
    TH1D * h5b = (TH1D *) fs1->Get("MpKm_PiPi");
    
    TH1D * h0c = (TH1D *) fs->Get("MKK_BoundStateAll");
    TH1D * h1c = (TH1D *) fs->Get("MKK_BoundStateKK");
    TH1D * h2c = (TH1D *) fs->Get("MKK_phi");
    TH1D * h3c = (TH1D *) fs->Get("MKK_Lambda1520");
    TH1D * h4c = (TH1D *) fs->Get("MKK_directKK");
    TH1D * h5c = (TH1D *) fs1->Get("MKK_PiPi");
    
    TH1D * h0d = (TH1D *) fs->Get("Theta_pKp_BoundStateAll");
    TH1D * h0dz = (TH1D *) fs->Get("Theta_z_pKp_BoundStateAll");
    //TH1D * h1d = (TH1D *) fs->Get();
    //TH1D * h2 = (TH1D *) fs->Get();
    //TH1D * h3 = (TH1D *) fs->Get();
    //TH1D * h4 = (TH1D *) fs->Get();
    TH1D * h5d = (TH1D *) fs1->Get("Theta_pKp_PiPi");
    TH1D * h5dz = (TH1D *) fs1->Get("Theta_z_pKp_PiPi");
    
    TH1D * h0e = (TH1D *) fs->Get("Theta_pKm_BoundStateAll");
    TH1D * h0ez = (TH1D *) fs->Get("Theta_z_pKm_BoundStateAll");
    TH1D * h5e = (TH1D *) fs1->Get("Theta_pKm_PiPi");
    TH1D * h5ez = (TH1D *) fs1->Get("Theta_z_pKm_PiPi");
    
    TH1D * h0f = (TH1D *) fs->Get("Theta_KpKm_BoundStateAll");
    TH1D * h0fz = (TH1D *) fs->Get("Theta_z_KpKm_BoundStateAll");
    TH1D * h5f = (TH1D *) fs1->Get("Theta_KpKm_PiPi");
    TH1D * h5fz = (TH1D *) fs1->Get("Theta_z_KpKm_PiPi");
    
    double mispid = 2500.0;//10000.0;//2500.0;
    //
    h0->Scale(lumi*time/h0->GetXaxis()->GetBinWidth(1)/1000);
    h1->Scale(lumi*time/h1->GetXaxis()->GetBinWidth(1)/1000);
    h2->Scale(lumi*time/h2->GetXaxis()->GetBinWidth(1)/1000);
    h3->Scale(lumi*time/h3->GetXaxis()->GetBinWidth(1)/1000);
    h4->Scale(lumi*time/h4->GetXaxis()->GetBinWidth(1)/1000);
    h5->Scale(lumi*time/1000);
    
    h0a->Scale(lumi*time/h0a->GetXaxis()->GetBinWidth(1)/1000);
    h1a->Scale(lumi*time/h1a->GetXaxis()->GetBinWidth(1)/1000);
    h2a->Scale(lumi*time/h2a->GetXaxis()->GetBinWidth(1)/1000);
    h3a->Scale(lumi*time/h3a->GetXaxis()->GetBinWidth(1)/1000);
    h4a->Scale(lumi*time/h4a->GetXaxis()->GetBinWidth(1)/1000);
    h5a->Scale(lumi*time/1000);
    
    h0b->Scale(lumi*time/h0b->GetXaxis()->GetBinWidth(1)/1000);
    h1b->Scale(lumi*time/h1b->GetXaxis()->GetBinWidth(1)/1000);
    h2b->Scale(lumi*time/h2b->GetXaxis()->GetBinWidth(1)/1000);
    h3b->Scale(lumi*time/h3b->GetXaxis()->GetBinWidth(1)/1000);
    h4b->Scale(lumi*time/h4b->GetXaxis()->GetBinWidth(1)/1000);
    h5b->Scale(lumi*time/1000);
    
    h0c->Scale(lumi*time/h0c->GetXaxis()->GetBinWidth(1)/1000);
    h1c->Scale(lumi*time/h1c->GetXaxis()->GetBinWidth(1)/1000);
    h2c->Scale(lumi*time/h2c->GetXaxis()->GetBinWidth(1)/1000);
    h3c->Scale(lumi*time/h3c->GetXaxis()->GetBinWidth(1)/1000);
    h4c->Scale(lumi*time/h4c->GetXaxis()->GetBinWidth(1)/1000);
    h5c->Scale(lumi*time/1000);
    //

    
    h0d->Scale(lumi*time/h0d->GetXaxis()->GetBinWidth(1));
    h5d->Scale(lumi*time);
    
    h0dz->Scale(lumi*time/h0d->GetXaxis()->GetBinWidth(1));
    h5dz->Scale(lumi*time);
    
    h0e->Scale(lumi*time/h0e->GetXaxis()->GetBinWidth(1));
    h5e->Scale(lumi*time);
    
    h0ez->Scale(lumi*time/h0e->GetXaxis()->GetBinWidth(1));
    h5ez->Scale(lumi*time);
    
    h0f->Scale(lumi*time/h0f->GetXaxis()->GetBinWidth(1));
    h5f->Scale(lumi*time);
    
    h0fz->Scale(lumi*time/h0f->GetXaxis()->GetBinWidth(1));
    h5fz->Scale(lumi*time);
    //h5d->Scale(0.01);
    //h5e->Scale(0.01);
    //h5f->Scale(0.01);
    
    //
    //h5->Scale(0.01);
    //h5a->Scale(0.01);
    //h5b->Scale(0.01);
    //h5c->Scale(0.01);
    //
    h0->SetLineColor(1); h1->SetLineColor(4); h2->SetLineColor(2); h3->SetLineColor(6); h4->SetLineColor(7);h5->SetLineColor(3);
    h0a->SetLineColor(1); h1a->SetLineColor(4); h2a->SetLineColor(2); h3a->SetLineColor(6); h4a->SetLineColor(7);h5a->SetLineColor(3);
    h0b->SetLineColor(1); h1b->SetLineColor(4); h2b->SetLineColor(2); h3b->SetLineColor(6); h4b->SetLineColor(7);h5b->SetLineColor(3);
    h0c->SetLineColor(1); h1c->SetLineColor(4); h2c->SetLineColor(2); h3c->SetLineColor(6); h4c->SetLineColor(7);h5c->SetLineColor(3);
    h0d->SetLineColor(1);h5d->SetLineColor(3);
    h0e->SetLineColor(1);h5e->SetLineColor(3);
    h0f->SetLineColor(1);h5f->SetLineColor(3);
    
    h0dz->SetLineColor(1);h5dz->SetLineColor(3);
    h0ez->SetLineColor(1);h5ez->SetLineColor(3);
    h0fz->SetLineColor(1);h5fz->SetLineColor(3);
    //
    h5->Scale(1.0/mispid);
    h5a->Scale(1.0/mispid);
    h5b->Scale(1.0/mispid);
    h5c->Scale(1.0/mispid);
    h5d->Scale(1.0/mispid);
    h5e->Scale(1.0/mispid);
    h5f->Scale(1.0/mispid);
    h5dz->Scale(1.0/mispid);
    h5ez->Scale(1.0/mispid);
    h5fz->Scale(1.0/mispid);
    
    TLegend * leg = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg->AddEntry(h0, "pKK from N-#phi", "l");
    leg->AddEntry(h1, "only KK from N-#phi", "l");
    leg->AddEntry(h2, "#phi production", "l");
    leg->AddEntry(h3, "#Lambda(1520) production", "l");
    leg->AddEntry(h4, "direct KK production", "l");
    leg->AddEntry(h5, "misidentified #pi#pi", "l");
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);
    hB->SetMaximum(h0->GetMaximum() * 1.1);
    //hB->SetAxisRange(1.88, 2.31, "X");
    hB->GetXaxis()->SetRangeUser(1.88,2.31);
    hB->GetXaxis()->SetTitle("M(pKK) (GeV)");
    hB->GetYaxis()->SetTitle("counts / hour / MeV");
    hB->DrawClone("axis");
    h0->DrawClone("same");
    h1->DrawClone("same");
    h2->DrawClone("same");
    h3->DrawClone("same");
    h4->DrawClone("same");
    h5->DrawClone("same");
    leg->DrawClone("same");
    c0->Print(savepath + "Massbreakfinalcut.pdf(","pdf");//
    
    hB->SetMaximum(h5a->GetMaximum() * 1.1);
    hB->SetAxisRange(1.38, 1.81, "X");
    hB->GetXaxis()->SetTitle("M(pK^{+}) (GeV)");
    hB->DrawClone("axis");
    h0a->DrawClone("same");
    h1a->DrawClone("same");
    h2a->DrawClone("same");
    h3a->DrawClone("same");
    h4a->DrawClone("same");
    h5a->DrawClone("same");
    leg->DrawClone("same");
    c0->Print(savepath + "Massbreakfinalcut.pdf","pdf");//
    
    hB->SetMaximum(h5b->GetMaximum() * 1.1);
    hB->SetAxisRange(1.38, 1.81, "X");
    hB->GetXaxis()->SetTitle("M(pK^{-}) (GeV)");
    hB->DrawClone("axis");
    h0b->DrawClone("same");
    h1b->DrawClone("same");
    h2b->DrawClone("same");
    h3b->Scale(1/20.);h3b->DrawClone("same");
    h4b->DrawClone("same");
    h5b->DrawClone("same");
    leg->Clear();
    leg->AddEntry(h0, "pKK from N-#phi", "l");
    leg->AddEntry(h1, "only KK from N-#phi", "l");
    leg->AddEntry(h2, "#phi production", "l");
    leg->AddEntry(h3, "#Lambda(1520) production / 20", "l");
    leg->AddEntry(h4, "direct KK production", "l");
    leg->AddEntry(h5, "misidentified #pi#pi", "l");
    leg->DrawClone("same");
    c0->Print(savepath + "Massbreakfinalcut.pdf","pdf");//
    
    hB->SetMaximum(h5c->GetMaximum() * 1.1);
    hB->SetAxisRange(0.98, 1.41, "X");
    hB->GetXaxis()->SetTitle("M(K^{+}K^{-}) (GeV)");
    hB->DrawClone("axis");
    h0c->DrawClone("same");
    h1c->DrawClone("same");
    h2c->Scale(1/50.);h2c->DrawClone("same");
    h3c->DrawClone("same");
    h4c->DrawClone("same");
    h5c->DrawClone("same");
    leg->Clear();
    leg->AddEntry(h0, "pKK from N-#phi", "l");
    leg->AddEntry(h1, "only KK from N-#phi", "l");
    leg->AddEntry(h2, "#phi production / 50", "l");
    leg->AddEntry(h3, "#Lambda(1520) production", "l");
    leg->AddEntry(h4, "direct KK production", "l");
    leg->AddEntry(h5, "misidentified #pi#pi", "l");
    leg->DrawClone("same");
    c0->Print(savepath + "Massbreakfinalcut.pdf)","pdf");//
    
/*
    hB->SetMaximum(20.0);
    hB->SetAxisRange(0.0,  M_PI, "X");
    hB->GetXaxis()->SetTitle("Theta(pK^{+}) (rad)");
    hB->GetYaxis()->SetTitle("counts / hour / rad");
    hB->DrawClone("axis");
    h0d->DrawClone("same");
    h5d->DrawClone("same");
    leg->Clear();
    leg->AddEntry(h0, "pKK from N-#phi", "l");
    leg->AddEntry(h5, "misidentified #pi#pi", "l");
    leg->DrawClone("same");
    c0->Print(savepath + "Massbreakfinalcut.pdf","pdf");//
     
    
    hB->SetMaximum(20.0);
    hB->SetAxisRange(0.0,  M_PI, "X");
    hB->GetXaxis()->SetTitle("Theta(pK^{-}) (rad)");
    hB->DrawClone("axis");
    h0e->DrawClone("same");
    h5e->DrawClone("same");
    leg->Clear();
    leg->AddEntry(h0, "pKK from N-#phi", "l");
    leg->AddEntry(h5, "misidentified #pi#pi", "l");
    leg->DrawClone("same");
    c0->Print(savepath + "Massbreakfinalcut.pdf","pdf");//
    
    hB->SetMaximum(20.0);
    hB->SetAxisRange(0.0,  M_PI, "X");
    hB->GetXaxis()->SetTitle("Theta(K^{+}K^{-}) (rad)");
    hB->DrawClone("axis");
    h0f->DrawClone("same");
    h5f->DrawClone("same");
    leg->Clear();
    leg->AddEntry(h0, "pKK from N-#phi", "l");
    leg->AddEntry(h5, "misidentified #pi#pi", "l");
    leg->DrawClone("same");
    c0->Print(savepath + "Massbreakfinalcut.pdf","pdf");//
    
    
    
    
    
    
    hB->SetMaximum(20.0);
    hB->SetTitle("c.m. angle betweem p and K^{+}");
    hB->SetAxisRange(0.0,  M_PI, "X");
    hB->GetXaxis()->SetTitle("Theta(pK^{+}) (rad)");
    hB->GetYaxis()->SetTitle("counts / hour / rad");
    hB->DrawClone("axis");
    h0dz->DrawClone("same");
    h5dz->DrawClone("same");
    leg->Clear();
    leg->AddEntry(h0, "pKK from N-#phi", "l");
    leg->AddEntry(h5, "misidentified #pi#pi", "l");
    leg->DrawClone("same");
    c0->Print(savepath + "Massbreakfinalcut.pdf","pdf");//
     
    
    hB->SetMaximum(20.0);
    hB->SetTitle("c.m. angle betweem p and K^{-}");
    hB->SetAxisRange(0.0,  M_PI, "X");
    hB->GetXaxis()->SetTitle("Theta(pK^{-}) (rad)");
    hB->DrawClone("axis");
    h0ez->DrawClone("same");
    h5ez->DrawClone("same");
    leg->Clear();
    leg->AddEntry(h0, "pKK from N-#phi", "l");
    leg->AddEntry(h5, "misidentified #pi#pi", "l");
    leg->DrawClone("same");
    c0->Print(savepath + "Massbreakfinalcut.pdf","pdf");//
    
    hB->SetMaximum(20.0);
    hB->SetTitle("c.m. angle betweem K^{+} and K^{-}");
    hB->SetAxisRange(0.0,  M_PI, "X");
    hB->GetXaxis()->SetTitle("Theta(K^{+}K^{-}) (rad)");
    hB->DrawClone("axis");
    h0fz->DrawClone("same");
    h5fz->DrawClone("same");
    leg->Clear();
    leg->AddEntry(h0, "pKK from N-#phi", "l");
    leg->AddEntry(h5, "misidentified #pi#pi", "l");
    leg->DrawClone("same");
    c0->Print(savepath + "Massbreakfinalcut.pdf)","pdf");//
    
  */ 
    
    
    int i1, i2; 
    double count;
    double count1,count2,count3,count4;
    i1 = h0->GetXaxis()->FindBin(1.94);
    i2 = h0->GetXaxis()->FindBin(1.96);
    count = h0->Integral(i1,i2,"width");
    cout<< "signal rate is " << count * 1000.0 << "/h" <<endl;
    
    i1 = h5->GetXaxis()->FindBin(1.94);
    i2 = h5->GetXaxis()->FindBin(1.96);
    count = h5->Integral(i1,i2,"width");
    cout<< "misidentified twopion rate is " << count * 1000.0 << "/h" <<endl;
    
    //i1 = h5a->GetXaxis()->FindBin(1.43);
    //i2 = h5a->GetXaxis()->FindBin(2.2);
    //count = h5a->Integral(i1,i2,"width");
    //cout<< "the area of pKp is " << count <<endl;
    
    i1 = h1->GetXaxis()->FindBin(1.94);
    i2 = h1->GetXaxis()->FindBin(1.96);
    count1 = h1->Integral(i1,i2,"width");
    
    i1 = h2->GetXaxis()->FindBin(1.94);
    i2 = h2->GetXaxis()->FindBin(1.96);
    count2 = h2->Integral(i1,i2,"width");
    
    i1 = h3->GetXaxis()->FindBin(1.94);
    i2 = h3->GetXaxis()->FindBin(1.96);
    count3 = h3->Integral(i1,i2,"width");
    
    i1 = h4->GetXaxis()->FindBin(1.94);
    i2 = h4->GetXaxis()->FindBin(1.96);
    count4 = h4->Integral(i1,i2,"width");
    
    cout<< "uncorrelated rate is " << (count1 + count2 + count3 + count4) * 1000.0 << "/h" <<endl; //
  }
 
 
  if (opt == -511){// momentum cut    center of mass frame
    TFile * fs = new TFile(loadpath + "breakcuttest.root", "r");
    TFile * fs1 = new TFile(loadpath + "modifiedtwopinew.root", "r");
    TH2D * d0a = (TH2D *) fs->Get("Momentum_z_p_Kp_BoundStateAll");
    TH2D * d1a = (TH2D *) fs->Get("Momentum_z_p_Kp_BoundStateKK");
    TH2D * d2a = (TH2D *) fs->Get("Momentum_z_p_Kp_phi");
    TH2D * d3a = (TH2D *) fs->Get("Momentum_z_p_Kp_Lambda1520");
    TH2D * d4a = (TH2D *) fs->Get("Momentum_z_p_Kp_KK");
    TH2D * d5a = (TH2D *) fs1->Get("Momentum_z_p_Pip_PiPi");
    
    TH2D * d0b = (TH2D *) fs->Get("Momentum_z_p_Km_BoundStateAll");
    TH2D * d1b = (TH2D *) fs->Get("Momentum_z_p_Km_BoundStateKK");
    TH2D * d2b = (TH2D *) fs->Get("Momentum_z_p_Km_phi");
    TH2D * d3b = (TH2D *) fs->Get("Momentum_z_p_Km_Lambda1520");
    TH2D * d4b = (TH2D *) fs->Get("Momentum_z_p_Km_KK");
    TH2D * d5b = (TH2D *) fs1->Get("Momentum_z_p_Pim_PiPi");
    
    TH2D * d0c = (TH2D *) fs->Get("Momentum_z_Kp_Km_BoundStateAll");
    TH2D * d1c = (TH2D *) fs->Get("Momentum_z_Kp_Km_BoundStateKK");
    TH2D * d2c = (TH2D *) fs->Get("Momentum_z_Kp_Km_phi");
    TH2D * d3c = (TH2D *) fs->Get("Momentum_z_Kp_Km_Lambda1520");
    TH2D * d4c = (TH2D *) fs->Get("Momentum_z_Kp_Km_KK");
    TH2D * d5c = (TH2D *) fs1->Get("Momentum_z_Pip_Pim_PiPi");
    
    
    
    d0a->GetXaxis()->SetTitle("P(K^{+}) (GeV)");
    d0a->GetXaxis()->CenterTitle(true);
    d0a->GetXaxis()->SetTitleSize(0.06);
    d0a->GetXaxis()->SetTitleOffset(1.15);
    d0a->GetXaxis()->SetLabelSize(0.06);
    d0a->GetXaxis()->SetRangeUser(0.0, 2.0);
    d0a->GetXaxis()->SetNdivisions(6, 5, 0);
    
    d0a->GetYaxis()->SetTitle("P(p) (GeV)");
    d0a->GetYaxis()->CenterTitle(true);
    d0a->GetYaxis()->SetTitleSize(0.06);
    d0a->GetYaxis()->SetTitleOffset(1.15);
    d0a->GetYaxis()->SetLabelSize(0.06);
    d0a->GetYaxis()->SetRangeUser(0.0, 2.0);
    d0a->GetYaxis()->SetNdivisions(5, 5, 0);
    
    d0b->GetXaxis()->SetTitle("P(K^{-}) (GeV)");
    d0b->GetXaxis()->CenterTitle(true);
    d0b->GetXaxis()->SetTitleSize(0.06);
    d0b->GetXaxis()->SetTitleOffset(1.15);
    d0b->GetXaxis()->SetLabelSize(0.06);
    d0b->GetXaxis()->SetRangeUser(0.0, 2.0);
    d0b->GetXaxis()->SetNdivisions(6, 5, 0);
    d0b->GetYaxis()->SetTitle("P(p) (GeV)");
    d0b->GetYaxis()->CenterTitle(true);
    d0b->GetYaxis()->SetTitleSize(0.06);
    d0b->GetYaxis()->SetTitleOffset(1.15);
    d0b->GetYaxis()->SetLabelSize(0.06);
    d0b->GetYaxis()->SetRangeUser(0.0, 2.0);
    d0b->GetYaxis()->SetNdivisions(5, 5, 0);
    d0c->GetXaxis()->SetTitle("P(K^{+}) (GeV)");
    d0c->GetXaxis()->CenterTitle(true);
    d0c->GetXaxis()->SetTitleSize(0.06);
    d0c->GetXaxis()->SetTitleOffset(1.15);
    d0c->GetXaxis()->SetLabelSize(0.06);
    d0c->GetXaxis()->SetRangeUser(0.0, 2.0);
    d0c->GetXaxis()->SetNdivisions(6, 5, 0);
    d0c->GetYaxis()->SetTitle("P(K^{-}) (GeV)");
    d0c->GetYaxis()->CenterTitle(true);
    d0c->GetYaxis()->SetTitleSize(0.06);
    d0c->GetYaxis()->SetTitleOffset(1.15);
    d0c->GetYaxis()->SetLabelSize(0.06);
    d0c->GetYaxis()->SetRangeUser(0.0, 2.0);
    d0c->GetYaxis()->SetNdivisions(5, 5, 0);
    //
    d0a->SetStats(0);
    d0b->SetStats(0);
    d0c->SetStats(0);
   
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetRightMargin(0.12);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);
    //
    d0a->SetTitle("pKK from N-#phi");
    d0a->DrawClone("axis");
    d0a->DrawClone("colsame");
    c0->Print(savepath + "cmMomentumbreakfinalcut.pdf(","pdf");//
    
    d0a->SetTitle("only KK from N-#phi");
    d0a->DrawClone("axis");
    d1a->DrawClone("colsame");
    c0->Print(savepath + "cmMomentumbreakfinalcut.pdf","pdf");//
    
    d0a->SetTitle("#phi production");
    d0a->DrawClone("axis");
    d2a->DrawClone("colsame");
    c0->Print(savepath + "cmMomentumbreakfinalcut.pdf","pdf");//
    
    d0a->SetTitle("#Lambda(1520) production");
    d0a->DrawClone("axis");
    d3a->DrawClone("colsame");
    c0->Print(savepath + "cmMomentumbreakfinalcut.pdf","pdf");//    
    
    d0a->SetTitle("direct KK production");
    d0a->DrawClone("axis");
    d4a->DrawClone("colsame");
    c0->Print(savepath + "cmMomentumbreakfinalcut.pdf","pdf");//
    
    d0a->SetTitle("misidentified #pi#pi");
    d0a->DrawClone("axis");
    d5a->DrawClone("colsame");
    c0->Print(savepath + "cmMomentumbreakfinalcut.pdf","pdf");//
    //
    d0b->SetTitle("pKK from N-#phi");
    d0b->DrawClone("axis");
    d0b->DrawClone("colsame");
    c0->Print(savepath + "cmMomentumbreakfinalcut.pdf(","pdf");//
    
    d0b->SetTitle("only KK from N-#phi");
    d0b->DrawClone("axis");
    d1b->DrawClone("colsame");
    c0->Print(savepath + "cmMomentumbreakfinalcut.pdf","pdf");//
    
    d0b->SetTitle("#phi production");
    d0b->DrawClone("axis");
    d2b->DrawClone("colsame");
    c0->Print(savepath + "cmMomentumbreakfinalcut.pdf","pdf");//
    
    d0b->SetTitle("#Lambda(1520) production");
    d0b->DrawClone("axis");
    d3b->DrawClone("colsame");
    c0->Print(savepath + "cmMomentumbreakfinalcut.pdf","pdf");//    
    
    d0b->SetTitle("direct KK production");
    d0b->DrawClone("axis");
    d4b->DrawClone("colsame");
    c0->Print(savepath + "cmMomentumbreakfinalcut.pdf","pdf");//
    
    d0b->SetTitle("misidentified #pi#pi");
    d0b->DrawClone("axis");
    d5b->DrawClone("colsame");
    c0->Print(savepath + "cmMomentumbreakfinalcut.pdf","pdf");//
    //
    d0c->SetTitle("pKK from N-#phi");
    d0c->DrawClone("axis");
    d0c->DrawClone("colsame");
    c0->Print(savepath + "cmMomentumbreakfinalcut.pdf(","pdf");//
    
    d0c->SetTitle("only KK from N-#phi");
    d0c->DrawClone("axis");
    d1c->DrawClone("colsame");
    c0->Print(savepath + "cmMomentumbreakfinalcut.pdf","pdf");//
    
    d0c->SetTitle("#phi production");
    d0c->DrawClone("axis");
    d2c->DrawClone("colsame");
    c0->Print(savepath + "cmMomentumbreakfinalcut.pdf","pdf");//
    
    d0c->SetTitle("#Lambda(1520) production");
    d0c->DrawClone("axis");
    d3c->DrawClone("colsame");
    c0->Print(savepath + "cmMomentumbreakfinalcut.pdf","pdf");//    
    
    d0c->SetTitle("direct KK production");
    d0c->DrawClone("axis");
    d4c->DrawClone("colsame");
    c0->Print(savepath + "cmMomentumbreakfinalcut.pdf","pdf");//
    
    d0c->SetTitle("misidentified #pi#pi");
    d0c->DrawClone("axis");
    d5c->DrawClone("colsame");
    c0->Print(savepath + "cmMomentumbreakfinalcut.pdf)","pdf");//
    
  }

  if (opt == -411){//Momentum cut
    TFile * fs = new TFile(loadpath + "breakcuttest.root", "r");
    TFile * fs1 = new TFile(loadpath + "modifiedtwopinew.root", "r");
    TH2D * d0a = (TH2D *) fs->Get("Momentum_p_Kp_BoundStateAll");
    TH2D * d1a = (TH2D *) fs->Get("Momentum_p_Kp_BoundStateKK");
    TH2D * d2a = (TH2D *) fs->Get("Momentum_p_Kp_phi");
    TH2D * d3a = (TH2D *) fs->Get("Momentum_p_Kp_Lambda1520");
    TH2D * d4a = (TH2D *) fs->Get("Momentum_p_Kp_KK");
    TH2D * d5a = (TH2D *) fs1->Get("Momentum_p_Pip_PiPi");
    
    TH2D * d0b = (TH2D *) fs->Get("Momentum_p_Km_BoundStateAll");
    TH2D * d1b = (TH2D *) fs->Get("Momentum_p_Km_BoundStateKK");
    TH2D * d2b = (TH2D *) fs->Get("Momentum_p_Km_phi");
    TH2D * d3b = (TH2D *) fs->Get("Momentum_p_Km_Lambda1520");
    TH2D * d4b = (TH2D *) fs->Get("Momentum_p_Km_KK");
    TH2D * d5b = (TH2D *) fs1->Get("Momentum_p_Pim_PiPi");
    
    TH2D * d0c = (TH2D *) fs->Get("Momentum_Kp_Km_BoundStateAll");
    TH2D * d1c = (TH2D *) fs->Get("Momentum_Kp_Km_BoundStateKK");
    TH2D * d2c = (TH2D *) fs->Get("Momentum_Kp_Km_phi");
    TH2D * d3c = (TH2D *) fs->Get("Momentum_Kp_Km_Lambda1520");
    TH2D * d4c = (TH2D *) fs->Get("Momentum_Kp_Km_KK");
    TH2D * d5c = (TH2D *) fs1->Get("Momentum_Pip_Pim_PiPi");
    
    TGraph * g0a = new TGraph(1000);
    TGraph * g1a = new TGraph(1000);
    TGraph * g2a = new TGraph(1000);
    TGraph * g3a = new TGraph(1000);
    TGraph * g4a = new TGraph(1000);
    TGraph * g5a = new TGraph(1000);
    TGraph * g0b = new TGraph(1000);
    TGraph * g1b = new TGraph(1000);
    TGraph * g2b = new TGraph(1000);
    TGraph * g3b = new TGraph(1000);
    TGraph * g4b = new TGraph(1000);
    TGraph * g5b = new TGraph(1000);
    TGraph * g0c = new TGraph(1000);
    TGraph * g1c = new TGraph(1000);
    TGraph * g2c = new TGraph(1000);
    TGraph * g3c = new TGraph(1000);
    TGraph * g4c = new TGraph(1000);
    TGraph * g5c = new TGraph(1000);
    double x, y;
    for (int i = 0; i < 1000; i++){
      d0a->GetRandom2(x, y);
      g0a->SetPoint(i, x, y);
      d1a->GetRandom2(x, y);
      g1a->SetPoint(i, x, y);
      d2a->GetRandom2(x, y);
      g2a->SetPoint(i, x, y);
      d3a->GetRandom2(x, y);
      g3a->SetPoint(i, x, y);
      d4a->GetRandom2(x, y);
      g4a->SetPoint(i, x, y);
      d5a->GetRandom2(x, y);
      g5a->SetPoint(i, x, y);
      
      d0b->GetRandom2(x, y);
      g0b->SetPoint(i, x, y);
      d1b->GetRandom2(x, y);
      g1b->SetPoint(i, x, y);
      d2b->GetRandom2(x, y);
      g2b->SetPoint(i, x, y);
      d3b->GetRandom2(x, y);
      g3b->SetPoint(i, x, y);
      d4b->GetRandom2(x, y);
      g4b->SetPoint(i, x, y);
      d5b->GetRandom2(x, y);
      g5b->SetPoint(i, x, y);
      
      d0c->GetRandom2(x, y);
      g0c->SetPoint(i, x, y);
      d1c->GetRandom2(x, y);
      g1c->SetPoint(i, x, y);
      d2c->GetRandom2(x, y);
      g2c->SetPoint(i, x, y);
      d3c->GetRandom2(x, y);
      g3c->SetPoint(i, x, y);
      d4c->GetRandom2(x, y);
      g4c->SetPoint(i, x, y);
      d5c->GetRandom2(x, y);
      g5c->SetPoint(i, x, y);
    }
    //
    d0a->GetXaxis()->SetTitle("P(K^{+}) (GeV)");
    d0a->GetXaxis()->CenterTitle(true);
    d0a->GetXaxis()->SetTitleSize(0.06);
    d0a->GetXaxis()->SetTitleOffset(1.15);
    d0a->GetXaxis()->SetLabelSize(0.06);
    d0a->GetXaxis()->SetRangeUser(0.0, 2.0);
    d0a->GetXaxis()->SetNdivisions(6, 5, 0);
    
    d0a->GetYaxis()->SetTitle("P(p) (GeV)");
    d0a->GetYaxis()->CenterTitle(true);
    d0a->GetYaxis()->SetTitleSize(0.06);
    d0a->GetYaxis()->SetTitleOffset(1.15);
    d0a->GetYaxis()->SetLabelSize(0.06);
    d0a->GetYaxis()->SetRangeUser(0.0, 2.0);
    d0a->GetYaxis()->SetNdivisions(5, 5, 0);
    
    d0b->GetXaxis()->SetTitle("P(K^{-}) (GeV)");
    d0b->GetXaxis()->CenterTitle(true);
    d0b->GetXaxis()->SetTitleSize(0.06);
    d0b->GetXaxis()->SetTitleOffset(1.15);
    d0b->GetXaxis()->SetLabelSize(0.06);
    d0b->GetXaxis()->SetRangeUser(0.0, 2.0);
    d0b->GetXaxis()->SetNdivisions(6, 5, 0);
    d0b->GetYaxis()->SetTitle("P(p) (GeV)");
    d0b->GetYaxis()->CenterTitle(true);
    d0b->GetYaxis()->SetTitleSize(0.06);
    d0b->GetYaxis()->SetTitleOffset(1.15);
    d0b->GetYaxis()->SetLabelSize(0.06);
    d0b->GetYaxis()->SetRangeUser(0.0, 2.0);
    d0b->GetYaxis()->SetNdivisions(5, 5, 0);
    d0c->GetXaxis()->SetTitle("P(K^{+}) (GeV)");
    d0c->GetXaxis()->CenterTitle(true);
    d0c->GetXaxis()->SetTitleSize(0.06);
    d0c->GetXaxis()->SetTitleOffset(1.15);
    d0c->GetXaxis()->SetLabelSize(0.06);
    d0c->GetXaxis()->SetRangeUser(0.0, 2.0);
    d0c->GetXaxis()->SetNdivisions(6, 5, 0);
    d0c->GetYaxis()->SetTitle("P(K^{-}) (GeV)");
    d0c->GetYaxis()->CenterTitle(true);
    d0c->GetYaxis()->SetTitleSize(0.06);
    d0c->GetYaxis()->SetTitleOffset(1.15);
    d0c->GetYaxis()->SetLabelSize(0.06);
    d0c->GetYaxis()->SetRangeUser(0.0, 2.0);
    d0c->GetYaxis()->SetNdivisions(5, 5, 0);
    //
    d0a->SetStats(0);
    d0b->SetStats(0);
    d0c->SetStats(0);
    //
    g0a->SetMarkerStyle(20);
    g0a->SetMarkerColor(1);
    g0a->SetMarkerSize(0.6);
    g1a->SetMarkerStyle(21);
    g1a->SetMarkerColor(4);
    g1a->SetMarkerSize(0.6);
    g2a->SetMarkerStyle(22);
    g2a->SetMarkerColor(2);
    g2a->SetMarkerSize(0.6);
    g3a->SetMarkerStyle(23);
    g3a->SetMarkerColor(6);
    g3a->SetMarkerSize(0.6);
    g4a->SetMarkerStyle(25);
    g4a->SetMarkerColor(7);
    g4a->SetMarkerSize(0.6);
    g5a->SetMarkerStyle(33);
    g5a->SetMarkerColor(3);
    g5a->SetMarkerSize(0.6);
    
    g0b->SetMarkerStyle(20);
    g0b->SetMarkerColor(1);
    g0b->SetMarkerSize(0.6);
    g1b->SetMarkerStyle(21);
    g1b->SetMarkerColor(4);
    g1b->SetMarkerSize(0.6);
    g2b->SetMarkerStyle(22);
    g2b->SetMarkerColor(2);
    g2b->SetMarkerSize(0.6);
    g3b->SetMarkerStyle(23);
    g3b->SetMarkerColor(6);
    g3b->SetMarkerSize(0.6);
    g4b->SetMarkerStyle(25);
    g4b->SetMarkerColor(7);
    g4b->SetMarkerSize(0.6);
    g5b->SetMarkerStyle(33);
    g5b->SetMarkerColor(3);
    g5b->SetMarkerSize(0.6);
    
    g0c->SetMarkerStyle(20);
    g0c->SetMarkerColor(1);
    g0c->SetMarkerSize(0.6);
    g1c->SetMarkerStyle(21);
    g1c->SetMarkerColor(4);
    g1c->SetMarkerSize(0.6);
    g2c->SetMarkerStyle(22);
    g2c->SetMarkerColor(2);
    g2c->SetMarkerSize(0.6);
    g3c->SetMarkerStyle(23);
    g3c->SetMarkerColor(6);
    g3c->SetMarkerSize(0.6);
    g4c->SetMarkerStyle(25);
    g4c->SetMarkerColor(7);
    g4c->SetMarkerSize(0.6);
    g5c->SetMarkerStyle(33);
    g5c->SetMarkerColor(3);
    g5c->SetMarkerSize(0.6);
    
    TLegend * leg = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg->AddEntry(g0a, "pKK from N-#phi", "p");
    leg->AddEntry(g1a, "only KK from N-#phi", "p");
    leg->AddEntry(g2a, "#phi production", "p");
    leg->AddEntry(g3a, "#Lambda(1520) production", "p");
    leg->AddEntry(g4a, "direct KK production", "p");
    leg->AddEntry(g5a, "misidentified PiPi", "p");
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetRightMargin(0.12);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);
    //
    d0a->DrawClone("axis");
    g0a->DrawClone("psame");
    g1a->DrawClone("psame");
    g2a->DrawClone("psame");
    g3a->DrawClone("psame");
    g4a->DrawClone("psame");
    g5a->DrawClone("psame");
    leg->DrawClone("same");
    c0->Print(savepath + "Momentumbreakfinalcut.pdf(","pdf");//
    
    d0a->SetTitle("pKK from N-#phi");
    d0a->DrawClone("axis");
    d0a->DrawClone("colsame");
    c0->Print(savepath + "Momentumbreakfinalcut.pdf","pdf");//
    
    d0a->SetTitle("only KK from N-#phi");
    d0a->DrawClone("axis");
    d1a->DrawClone("colsame");
    c0->Print(savepath + "Momentumbreakfinalcut.pdf","pdf");//
    
    d0a->SetTitle("#phi production");
    d0a->DrawClone("axis");
    d2a->DrawClone("colsame");
    c0->Print(savepath + "Momentumbreakfinalcut.pdf","pdf");//
    
    d0a->SetTitle("#Lambda(1520) production");
    d0a->DrawClone("axis");
    d3a->DrawClone("colsame");
    c0->Print(savepath + "Momentumbreakfinalcut.pdf","pdf");//    
    
    d0a->SetTitle("direct KK production");
    d0a->DrawClone("axis");
    d4a->DrawClone("colsame");
    c0->Print(savepath + "Momentumbreakfinalcut.pdf","pdf");//
    
    d0a->SetTitle("misidentified #pi#pi");
    d0a->DrawClone("axis");
    d5a->DrawClone("colsame");
    c0->Print(savepath + "Momentumbreakfinalcut.pdf","pdf");//
    //
    d0b->DrawClone("axis");
    g0b->DrawClone("psame");
    g1b->DrawClone("psame");
    g2b->DrawClone("psame");
    g3b->DrawClone("psame");
    g4b->DrawClone("psame");
    g5b->DrawClone("psame");
    leg->DrawClone("same");
    c0->Print(savepath + "Momentumbreakfinalcut.pdf(","pdf");//
    
    d0b->SetTitle("pKK from N-#phi");
    d0b->DrawClone("axis");
    d0b->DrawClone("colsame");
    c0->Print(savepath + "Momentumbreakfinalcut.pdf","pdf");//
    
    d0b->SetTitle("only KK from N-#phi");
    d0b->DrawClone("axis");
    d1b->DrawClone("colsame");
    c0->Print(savepath + "Momentumbreakfinalcut.pdf","pdf");//
    
    d0b->SetTitle("#phi production");
    d0b->DrawClone("axis");
    d2b->DrawClone("colsame");
    c0->Print(savepath + "Momentumbreakfinalcut.pdf","pdf");//
    
    d0b->SetTitle("#Lambda(1520) production");
    d0b->DrawClone("axis");
    d3b->DrawClone("colsame");
    c0->Print(savepath + "Momentumbreakfinalcut.pdf","pdf");//    
    
    d0b->SetTitle("direct KK production");
    d0b->DrawClone("axis");
    d4b->DrawClone("colsame");
    c0->Print(savepath + "Momentumbreakfinalcut.pdf","pdf");//
    
    d0b->SetTitle("misidentified #pi#pi");
    d0b->DrawClone("axis");
    d5b->DrawClone("colsame");
    c0->Print(savepath + "Momentumbreakfinalcut.pdf","pdf");//
    //
    d0c->DrawClone("axis");
    g0c->DrawClone("psame");
    g1c->DrawClone("psame");
    g2c->DrawClone("psame");
    g3c->DrawClone("psame");
    g4c->DrawClone("psame");
    g5c->DrawClone("psame");
    leg->DrawClone("same");
    c0->Print(savepath + "Momentumbreakfinalcut.pdf(","pdf");//
    
    d0c->SetTitle("pKK from N-#phi");
    d0c->DrawClone("axis");
    d0c->DrawClone("colsame");
    c0->Print(savepath + "Momentumbreakfinalcut.pdf","pdf");//
    
    d0c->SetTitle("only KK from N-#phi");
    d0c->DrawClone("axis");
    d1c->DrawClone("colsame");
    c0->Print(savepath + "Momentumbreakfinalcut.pdf","pdf");//
    
    d0c->SetTitle("#phi production");
    d0c->DrawClone("axis");
    d2c->DrawClone("colsame");
    c0->Print(savepath + "Momentumbreakfinalcut.pdf","pdf");//
    
    d0c->SetTitle("#Lambda(1520) production");
    d0c->DrawClone("axis");
    d3c->DrawClone("colsame");
    c0->Print(savepath + "Momentumbreakfinalcut.pdf","pdf");//    
    
    d0c->SetTitle("direct KK production");
    d0c->DrawClone("axis");
    d4c->DrawClone("colsame");
    c0->Print(savepath + "Momentumbreakfinalcut.pdf","pdf");//
    
    d0c->SetTitle("misidentified #pi#pi");
    d0c->DrawClone("axis");
    d5c->DrawClone("colsame");
    c0->Print(savepath + "Momentumbreakfinalcut.pdf)","pdf");//
  }
  if (opt == -521){//P-theta cut    center of mass frame
    TFile * fs = new TFile(loadpath + "breakcuttest.root", "r");
    TFile * fs1 = new TFile(loadpath + "modifiedtwopinew.root", "r");
    TH2D * r0a = (TH2D *) fs->Get("PTheta_z_p_BoundStateAll");
    TH2D * r0b = (TH2D *) fs->Get("PTheta_z_Kp_BoundStateAll");
    TH2D * r0c = (TH2D *) fs->Get("PTheta_z_Km_BoundStateAll");
    
    TH2D * r1a = (TH2D *) fs->Get("PTheta_z_p_BoundStateKK");
    TH2D * r1b = (TH2D *) fs->Get("PTheta_z_Kp_BoundStateKK");
    TH2D * r1c = (TH2D *) fs->Get("PTheta_z_Km_BoundStateKK");
    
    TH2D * r2a = (TH2D *) fs->Get("PTheta_z_p_phi");
    TH2D * r2b = (TH2D *) fs->Get("PTheta_z_Kp_phi");
    TH2D * r2c = (TH2D *) fs->Get("PTheta_z_Km_phi");
    TH2D * r3a = (TH2D *) fs->Get("PTheta_z_p_Lambda1520");
    TH2D * r3b = (TH2D *) fs->Get("PTheta_z_Kp_Lambda1520");
    TH2D * r3c = (TH2D *) fs->Get("PTheta_z_Km_Lambda1520");
    
    TH2D * r4a = (TH2D *) fs->Get("PTheta_z_p_KK");
    TH2D * r4b = (TH2D *) fs->Get("PTheta_z_Kp_KK");
    TH2D * r4c = (TH2D *) fs->Get("PTheta_z_Km_KK");
    
    TH2D * r5a = (TH2D *) fs1->Get("PTheta_z_p_PiPi");
    TH2D * r5b = (TH2D *) fs1->Get("PTheta_z_Kp_PiPi");
    TH2D * r5c = (TH2D *) fs1->Get("PTheta_z_Km_PiPi");
    
    TH2D * rB = new TH2D("rB", "", 1, 0.0, 90.0, 1, 0.0, 2.0);
  
    rB->SetStats(0);
    rB->GetXaxis()->SetTitle("#theta (deg)");
    rB->GetXaxis()->CenterTitle(true);
    rB->GetXaxis()->SetTitleSize(0.06);
    rB->GetXaxis()->SetTitleOffset(1.15);
    rB->GetXaxis()->SetLabelSize(0.06);
    rB->GetXaxis()->SetRangeUser(0.0, 90.0);
    rB->GetXaxis()->SetNdivisions(6, 5, 0);
    rB->GetYaxis()->SetTitle("P (GeV)");
    rB->GetYaxis()->CenterTitle(true);
    rB->GetYaxis()->SetTitleSize(0.06);
    rB->GetYaxis()->SetTitleOffset(1.15);
    rB->GetYaxis()->SetLabelSize(0.06);
    rB->GetYaxis()->SetRangeUser(0.0, 2.0);
    rB->GetYaxis()->SetNdivisions(5, 5, 0);
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetRightMargin(0.12);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);
    
    rB->SetTitle("proton from N-#phi");
    rB->DrawClone("axis");
    r0a->Draw("colsame");
    c0->Print(savepath + "cmPThetabreakfinalcut.pdf(", "pdf");
   
    rB->SetTitle("K^{+} from N-#phi");
    rB->DrawClone("axis");
    r0b->Draw("colsame");
    c0->Print(savepath + "cmPThetabreakfinalcut.pdf", "pdf");
    
    rB->SetTitle("K^{-} from N-#phi");
    rB->DrawClone("axis");
    r0c->Draw("colsame");
    c0->Print(savepath + "cmPThetabreakfinalcut.pdf", "pdf");//
    
    rB->SetTitle("proton with N-#phi production");
    rB->DrawClone("axis");
    r1a->Draw("colsame");
    c0->Print(savepath + "cmPThetabreakfinalcut.pdf", "pdf");
    
    rB->SetTitle("K^{+} from N-#phi");
    rB->DrawClone("axis");
    r1b->Draw("colsame");
    c0->Print(savepath + "cmPThetabreakfinalcut.pdf", "pdf");
    
    rB->SetTitle("K^{-} from N-#phi");
    rB->DrawClone("axis");
    r1c->Draw("colsame");
    c0->Print(savepath + "cmPThetabreakfinalcut.pdf", "pdf");//
    
    rB->SetTitle("proton with #phi production");
    rB->DrawClone("axis");
    r2a->Draw("colsame");
    c0->Print(savepath + "cmPThetabreakfinalcut.pdf", "pdf");
    
    rB->SetTitle("K^{+} from #phi");
    rB->DrawClone("axis");
    r2b->Draw("colsame");
    c0->Print(savepath + "cmPThetabreakfinalcut.pdf", "pdf");
    
    rB->SetTitle("K^{-} from #phi");
    rB->DrawClone("axis");
    r2c->Draw("colsame");
    c0->Print(savepath + "cmPThetabreakfinalcut.pdf", "pdf");//
    
    rB->SetTitle("proton from #Lambda(1520)");
    rB->DrawClone("axis");
    r3a->Draw("colsame");
    c0->Print(savepath + "cmPThetabreakfinalcut.pdf", "pdf");
    
    rB->SetTitle("K^{+} with #Lambda(1520) production");
    rB->DrawClone("axis");
    r3b->Draw("colsame");
    c0->Print(savepath + "cmPThetabreakfinalcut.pdf", "pdf");
    
    rB->SetTitle("K^{-} from #Lambda(1520)");
    rB->DrawClone("axis");
    r3c->Draw("colsame");
    c0->Print(savepath + "cmPThetabreakfinalcut.pdf", "pdf");//
    
    rB->SetTitle("proton with direct KK production");
    rB->DrawClone("axis");
    r4a->Draw("colsame");
    c0->Print(savepath + "cmPThetabreakfinalcut.pdf", "pdf");
    
    rB->SetTitle("K^{+} from direct KK production");
    rB->DrawClone("axis");
    r4b->Draw("colsame");
    c0->Print(savepath + "cmPThetabreakfinalcut.pdf", "pdf");
    
    rB->SetTitle("K^{-} from direct KK production");
    rB->DrawClone("axis");
    r4c->Draw("colsame");
    c0->Print(savepath + "cmPThetabreakfinalcut.pdf", "pdf");//
    
    rB->SetTitle("proton with misidentified #pi#pi");
    rB->DrawClone("axis");
    r5a->Draw("colsame");
    c0->Print(savepath + "cmPThetabreakfinalcut.pdf", "pdf");
    
    rB->SetTitle("K^{+} from misidentified #pi#pi");
    rB->DrawClone("axis");
    r5b->Draw("colsame");
    c0->Print(savepath + "cmPThetabreakfinalcut.pdf", "pdf");
    
    rB->SetTitle("K^{-} from misidentified #pi#pi");
    rB->DrawClone("axis");
    r5c->Draw("colsame");
    c0->Print(savepath + "cmPThetabreakfinalcut.pdf)", "pdf");//
    
    
  }
  if (opt == -421){//P-theta cut
    TFile * fs = new TFile(loadpath + "breakcuttest.root", "r");
    TFile * fs1 = new TFile(loadpath + "modifiedtwopinew.root", "r");
    TH2D * r0a = (TH2D *) fs->Get("PTheta_p_BoundStateAll");
    TH2D * r0b = (TH2D *) fs->Get("PTheta_Kp_BoundStateAll");
    TH2D * r0c = (TH2D *) fs->Get("PTheta_Km_BoundStateAll");
    
    TH2D * r1a = (TH2D *) fs->Get("PTheta_p_BoundStateKK");
    TH2D * r1b = (TH2D *) fs->Get("PTheta_Kp_BoundStateKK");
    TH2D * r1c = (TH2D *) fs->Get("PTheta_Km_BoundStateKK");
    
    TH2D * r2a = (TH2D *) fs->Get("PTheta_p_phi");
    TH2D * r2b = (TH2D *) fs->Get("PTheta_Kp_phi");
    TH2D * r2c = (TH2D *) fs->Get("PTheta_Km_phi");
    TH2D * r3a = (TH2D *) fs->Get("PTheta_p_Lambda1520");
    TH2D * r3b = (TH2D *) fs->Get("PTheta_Kp_Lambda1520");
    TH2D * r3c = (TH2D *) fs->Get("PTheta_Km_Lambda1520");
    
    TH2D * r4a = (TH2D *) fs->Get("PTheta_p_KK");
    TH2D * r4b = (TH2D *) fs->Get("PTheta_Kp_KK");
    TH2D * r4c = (TH2D *) fs->Get("PTheta_Km_KK");
    
    TH2D * r5a = (TH2D *) fs1->Get("PTheta_p_PiPi");
    TH2D * r5b = (TH2D *) fs1->Get("PTheta_Kp_PiPi");
    TH2D * r5c = (TH2D *) fs1->Get("PTheta_Km_PiPi");
    
    TH2D * rB = new TH2D("rB", "", 1, 0.0, 90.0, 1, 0.0, 2.0);
    rB->SetStats(0);
    rB->GetXaxis()->SetTitle("#theta (deg)");
    rB->GetXaxis()->CenterTitle(true);
    rB->GetXaxis()->SetTitleSize(0.06);
    rB->GetXaxis()->SetTitleOffset(1.15);
    rB->GetXaxis()->SetLabelSize(0.06);
    rB->GetXaxis()->SetRangeUser(0.0, 90.0);
    rB->GetXaxis()->SetNdivisions(6, 5, 0);
    rB->GetYaxis()->SetTitle("P (GeV)");
    rB->GetYaxis()->CenterTitle(true);
    rB->GetYaxis()->SetTitleSize(0.06);
    rB->GetYaxis()->SetTitleOffset(1.15);
    rB->GetYaxis()->SetLabelSize(0.06);
    rB->GetYaxis()->SetRangeUser(0.0, 2.0);
    rB->GetYaxis()->SetNdivisions(5, 5, 0);
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
      c0->SetRightMargin(0.12);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);
    
    rB->SetTitle("proton from N-#phi");
    rB->DrawClone("axis");
    r0a->Draw("colsame");
    c0->Print(savepath + "PThetabreakfinalcut.pdf(", "pdf");
   
    rB->SetTitle("K^{+} from N-#phi");
    rB->DrawClone("axis");
    r0b->Draw("colsame");
    c0->Print(savepath + "PThetabreakfinalcut.pdf", "pdf");
    
    rB->SetTitle("K^{-} from N-#phi");
    rB->DrawClone("axis");
    r0c->Draw("colsame");
    c0->Print(savepath + "PThetabreakfinalcut.pdf", "pdf");//
    
    rB->SetTitle("proton with N-#phi production");
    rB->DrawClone("axis");
    r1a->Draw("colsame");
    c0->Print(savepath + "PThetabreakfinalcut.pdf", "pdf");
    
    rB->SetTitle("K^{+} from N-#phi");
    rB->DrawClone("axis");
    r1b->Draw("colsame");
    c0->Print(savepath + "PThetabreakfinalcut.pdf", "pdf");
    
    rB->SetTitle("K^{-} from N-#phi");
    rB->DrawClone("axis");
    r1c->Draw("colsame");
    c0->Print(savepath + "PThetabreakfinalcut.pdf", "pdf");//
    
    rB->SetTitle("proton with #phi production");
    rB->DrawClone("axis");
    r2a->Draw("colsame");
    c0->Print(savepath + "PThetabreakfinalcut.pdf", "pdf");
    
    rB->SetTitle("K^{+} from #phi");
    rB->DrawClone("axis");
    r2b->Draw("colsame");
    c0->Print(savepath + "PThetabreakfinalcut.pdf", "pdf");
    
    rB->SetTitle("K^{-} from #phi");
    rB->DrawClone("axis");
    r2c->Draw("colsame");
    c0->Print(savepath + "PThetabreakfinalcut.pdf", "pdf");//
    
    rB->SetTitle("proton from #Lambda(1520)");
    rB->DrawClone("axis");
    r3a->Draw("colsame");
    c0->Print(savepath + "PThetabreakfinalcut.pdf", "pdf");
    
    rB->SetTitle("K^{+} with #Lambda(1520) production");
    rB->DrawClone("axis");
    r3b->Draw("colsame");
    c0->Print(savepath + "PThetabreakfinalcut.pdf", "pdf");
    
    rB->SetTitle("K^{-} from #Lambda(1520)");
    rB->DrawClone("axis");
    r3c->Draw("colsame");
    c0->Print(savepath + "PThetabreakfinalcut.pdf", "pdf");//
    
    rB->SetTitle("proton with direct KK production");
    rB->DrawClone("axis");
    r4a->Draw("colsame");
    c0->Print(savepath + "PThetabreakfinalcut.pdf", "pdf");
    
    rB->SetTitle("K^{+} from direct KK production");
    rB->DrawClone("axis");
    r4b->Draw("colsame");
    c0->Print(savepath + "PThetabreakfinalcut.pdf", "pdf");
    
    rB->SetTitle("K^{-} from direct KK production");
    rB->DrawClone("axis");
    r4c->Draw("colsame");
    c0->Print(savepath + "PThetabreakfinalcut.pdf", "pdf");//
    
    rB->SetTitle("proton with misidentified #pi#pi");
    rB->DrawClone("axis");
    r5a->Draw("colsame");
    c0->Print(savepath + "PThetabreakfinalcut.pdf", "pdf");
    
    rB->SetTitle("K^{+} from misidentified #pi#pi");
    rB->DrawClone("axis");
    r5b->Draw("colsame");
    c0->Print(savepath + "PThetabreakfinalcut.pdf", "pdf");
    
    rB->SetTitle("K^{-} from misidentified #pi#pi");
    rB->DrawClone("axis");
    r5c->Draw("colsame");
    c0->Print(savepath + "PThetabreakfinalcut.pdf)", "pdf");//
  }
  
  if(opt == -431){//phi-theta cut
   TFile * fs = new TFile(loadpath + "breakcuttest.root", "r");
    TFile * fs1 = new TFile(loadpath + "modifiedtwopinew.root", "r");
    
    TH2D * p0a = (TH2D *) fs->Get("phitheta_p_BoundStateAll");
    TH2D * p0b = (TH2D *) fs->Get("phitheta_Kp_BoundStateAll");
    TH2D * p0c = (TH2D *) fs->Get("phitheta_Km_BoundStateAll");
    
    TH2D * p5a = (TH2D *) fs1->Get("phitheta_p_PiPi");
    TH2D * p5b = (TH2D *) fs1->Get("phitheta_Kp_PiPi");
    TH2D * p5c = (TH2D *) fs1->Get("phitheta_Km_PiPi");
    
    TH2D * pB = new TH2D("pB", "", 1, 0.0, 150.0, 1, 0.0, 200.0);
    pB->SetStats(0);
    pB->GetXaxis()->SetTitle("#theta (rad)");
    pB->GetXaxis()->CenterTitle(true);
    pB->GetXaxis()->SetTitleSize(0.06);
    pB->GetXaxis()->SetTitleOffset(1.15);
    pB->GetXaxis()->SetLabelSize(0.06);
    pB->GetXaxis()->SetRangeUser(0.0, M_PI);
    pB->GetXaxis()->SetNdivisions(6, 5, 0);
    pB->GetYaxis()->SetTitle("#phi (rad)");
    pB->GetYaxis()->CenterTitle(true);
    pB->GetYaxis()->SetTitleSize(0.06);
    pB->GetYaxis()->SetTitleOffset(1.15);
    pB->GetYaxis()->SetLabelSize(0.06);
    pB->GetYaxis()->SetRangeUser(0.0, 2 * M_PI);
    pB->GetYaxis()->SetNdivisions(5, 5, 0);
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
    c0->SetRightMargin(0.12);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);
    
    pB->SetTitle("proton from N-#phi");
    pB->DrawClone("axis");
    p0a->Draw("colsame");
    c0->Print(savepath + "PhiThetabreakfinalcut.pdf(", "pdf");
   
    pB->SetTitle("K^{+} from N-#phi");
    pB->DrawClone("axis");
    p0b->Draw("colsame");
    c0->Print(savepath + "PhiThetabreakfinalcut.pdf", "pdf");
    
    pB->SetTitle("K^{-} from N-#phi");
    pB->DrawClone("axis");
    p0c->Draw("colsame");
    c0->Print(savepath + "PhiThetabreakfinalcut.pdf", "pdf");//
  
    pB->SetTitle("proton with misidentified #pi#pi");
    pB->DrawClone("axis");
    p5a->Draw("colsame");
    c0->Print(savepath + "PhiThetabreakfinalcut.pdf", "pdf");
    
    pB->SetTitle("K^{+} from misidentified #pi#pi");
    pB->DrawClone("axis");
    p5b->Draw("colsame");
    c0->Print(savepath + "PhiThetabreakfinalcut.pdf", "pdf");
    
    pB->SetTitle("K^{-} from misidentified #pi#pi");
    pB->DrawClone("axis");
    p5c->Draw("colsame");
    c0->Print(savepath + "PhiThetabreakfinalcut.pdf)", "pdf");//
  
  }
  if(opt == -441){//p-phi cut
    TFile * fs = new TFile(loadpath + "breakcuttest.root", "r");
    TFile * fs1 = new TFile(loadpath + "modifiedtwopinew.root", "r");
    TH2D * s0a = (TH2D *) fs->Get("Pphi_p_BoundStateAll");
    TH2D * s0b = (TH2D *) fs->Get("Pphi_Kp_BoundStateAll");
    TH2D * s0c = (TH2D *) fs->Get("Pphi_Km_BoundStateAll");
    
    TH2D * s5a = (TH2D *) fs1->Get("Pphi_p_PiPi");
    TH2D * s5b = (TH2D *) fs1->Get("Pphi_Kp_PiPi");
    TH2D * s5c = (TH2D *) fs1->Get("Pphi_Km_PiPi");
    
    TH2D * sB = new TH2D("rB", "", 1, 0.0, 180.0, 1, 0.0, 2.0);
    sB->SetStats(0);
    sB->GetXaxis()->SetTitle("#phi (deg)");
    sB->GetXaxis()->CenterTitle(true);
    sB->GetXaxis()->SetTitleSize(0.06);
    sB->GetXaxis()->SetTitleOffset(1.15);
    sB->GetXaxis()->SetLabelSize(0.06);
    sB->GetXaxis()->SetRangeUser(0.0, 90.0);
    sB->GetXaxis()->SetNdivisions(6, 5, 0);
    sB->GetYaxis()->SetTitle("P (GeV)");
    sB->GetYaxis()->CenterTitle(true);
    sB->GetYaxis()->SetTitleSize(0.06);
    sB->GetYaxis()->SetTitleOffset(1.15);
    sB->GetYaxis()->SetLabelSize(0.06);
    sB->GetYaxis()->SetRangeUser(0.0, 2.0);
    sB->GetYaxis()->SetNdivisions(5, 5, 0);
    TCanvas * c0 = new TCanvas("c0", "", 800, 600);
      c0->SetRightMargin(0.12);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);
    
    sB->SetTitle("proton from N-#phi");
    sB->DrawClone("axis");
    s0a->Draw("colsame");
    c0->Print(savepath + "Pphibreakfinalcut.pdf(", "pdf");
   
    sB->SetTitle("K^{+} from N-#phi");
    sB->DrawClone("axis");
    s0b->Draw("colsame");
    c0->Print(savepath + "Pphibreakfinalcut.pdf", "pdf");
    
    sB->SetTitle("K^{-} from N-#phi");
    sB->DrawClone("axis");
    s0c->Draw("colsame");
    c0->Print(savepath + "Pphibreakfinalcut.pdf", "pdf");//
    
    sB->SetTitle("proton with misidentified #pi#pi");
    sB->DrawClone("axis");
    s5a->Draw("colsame");
    c0->Print(savepath + "Pphibreakfinalcut.pdf", "pdf");
    
    sB->SetTitle("K^{+} from misidentified #pi#pi");
    sB->DrawClone("axis");
    s5b->Draw("colsame");
    c0->Print(savepath + "Pphibreakfinalcut.pdf", "pdf");
    
    sB->SetTitle("K^{-} from misidentified #pi#pi");
    sB->DrawClone("axis");
    s5c->Draw("colsame");
    c0->Print(savepath + "Pphibreakfinalcut.pdf)", "pdf");//
  }
 return 0;
}
