#include "Lsidis0.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"

using namespace std;

TFile * file_negative = new TFile("/var/phy/project/mepg/tl190/SoLIDacceptance/acceptance_solid_CLEO_SIDIS_3he_negative_output.root","r");
TFile * file_positive = new TFile("/var/phy/project/mepg/tl190/SoLIDacceptance/acceptance_solid_CLEO_SIDIS_3he_positive_output.root","r");
TH2F * acc_FA_negative = (TH2F *) file_negative->Get("acceptance_forwardangle");
TH2F * acc_LA_negative = (TH2F *) file_negative->Get("acceptance_largeangle");
TH2F * acc_FA_positive = (TH2F *) file_positive->Get("acceptance_forwardangle");
TH2F * acc_LA_positive = (TH2F *) file_positive->Get("acceptance_largeangle");

double GetAcceptance(const TLorentzVector p, const char * dfile){//
  double theta = p.Theta() / M_PI * 180.0;
  double mom = p.P();
  double acc = 0;
  if (theta > 8.0 && mom > 0.8){
    if (strcmp(dfile, "e-") == 0){
      acc += acc_FA_negative->GetBinContent(acc_FA_negative->GetXaxis()->FindBin(theta), acc_FA_negative->GetYaxis()->FindBin(mom));
      if (mom > 3.5){
	acc += acc_LA_negative->GetBinContent(acc_LA_negative->GetXaxis()->FindBin(theta), acc_LA_negative->GetYaxis()->FindBin(mom));
      }
    }
    else if (strcmp(dfile, "h+") == 0){
      acc += acc_FA_positive->GetBinContent(acc_FA_positive->GetXaxis()->FindBin(theta), acc_FA_positive->GetYaxis()->FindBin(mom));
    }    
    else if (strcmp(dfile, "h-") == 0){
      acc += acc_FA_negative->GetBinContent(acc_FA_negative->GetXaxis()->FindBin(theta), acc_FA_negative->GetYaxis()->FindBin(mom));
    }
  }
  return acc;
}

int main(int argc, char * argv[]){

  if (argc < 6) {
    cout << "missing inputs" << endl;
    return -1;
  }

  gRandom->SetSeed(1);
  gStyle->SetOptStat(0);

  Lsidis mysidis;
  TLorentzVector l(0, 0, 11.0, 11.0);
  TLorentzVector P(0, 0, 0, 0.938272);

  mysidis.SetNucleus(2, 1);//helium-3
  mysidis.SetHadron(argv[1]);

  int ic = mysidis.GetHadronCharge();

  mysidis.SetInitialState(l, P);
  mysidis.SetPDFset("CT14lo");

  double Q2min = atof(argv[2]);
  double Q2max = atof(argv[3]);
  double zmin = atof(argv[4]);
  double zmax = atof(argv[5]);
  //double Ptmin = atof(argv[6]);
  //double Ptmax = atof(argv[7]);
  
  double Xmin[6] = {Q2min/24.0, Q2min, zmin, 0.0, -M_PI, -M_PI};//x, Q2, z, Pt, phih, phiS
  double Xmax[6] = {0.7, Q2max, zmax, 1.6, M_PI, M_PI};
  mysidis.SetRange(Xmin, Xmax);

  cout << "Q2: " << Q2min << " -- " << Q2max << "     z: " << zmin << " -- " << zmax << endl; 

  double weight = 0.0;
  TLorentzVector lp(0, 0, 0, 0);
  TLorentzVector Ph(0, 0, 0, 0);

  double lumi = 1.0e+10 * pow(0.197327, 2);
  double vz;
  
  Long64_t Nsim = 0;
  Long64_t Nrec = 0;
  double acc, accr;

  TH1D * h0 = new TH1D("h0", "", 35, 0.0, 0.7);
  TH1D * h1 = new TH1D("h1", "", 35, 0.0, 0.7);
  TH1D * h2 = new TH1D("h2", "", 35, 0.0, 0.7);
  TH1D * h3 = new TH1D("h3", "", 35, 0.0, 0.7);
  TH1D * h4 = new TH1D("h4", "", 35, 0.0, 0.7);
  TH1D * h5 = new TH1D("h5", "", 35, 0.0, 0.7);

  TH1D * r0 = new TH1D("r0", "", 35, 0.0, 0.7);
  TH1D * r1 = new TH1D("r1", "", 35, 0.0, 0.7);
  TH1D * r2 = new TH1D("r2", "", 35, 0.0, 0.7);
  TH1D * r3 = new TH1D("r3", "", 35, 0.0, 0.7);
  TH1D * r4 = new TH1D("r4", "", 35, 0.0, 0.7);
  TH1D * r5 = new TH1D("r5", "", 35, 0.0, 0.7);

  TH1D * rh0 = new TH1D("rh0", "", 35, 0.0, 0.7);
  TH1D * rh1 = new TH1D("rh1", "", 35, 0.0, 0.7);
  TH1D * rh2 = new TH1D("rh2", "", 35, 0.0, 0.7);
  TH1D * rh3 = new TH1D("rh3", "", 35, 0.0, 0.7);
  TH1D * rh4 = new TH1D("rh4", "", 35, 0.0, 0.7);
  TH1D * rh5 = new TH1D("rh5", "", 35, 0.0, 0.7);

  h0->GetXaxis()->SetTitle("x");
  h0->GetXaxis()->SetLabelSize(0.04);
  h0->GetYaxis()->SetTitle("rate/Hz");
  h0->SetLineColor(4);
  h0->SetLineWidth(1.5);
  h1->GetXaxis()->SetTitle("x");
  h1->GetYaxis()->SetTitle("rate/Hz");
  h1->SetLineColor(4);
  h1->SetLineWidth(1.5);
  h2->GetXaxis()->SetTitle("x");
  h2->GetYaxis()->SetTitle("rate/Hz");
  h2->SetLineColor(4);
  h2->SetLineWidth(1.5);
  h3->GetXaxis()->SetTitle("x");
  h3->GetYaxis()->SetTitle("rate/Hz");
  h3->SetLineColor(4);
  h3->SetLineWidth(1.5);
  h4->GetXaxis()->SetTitle("x");
  h4->GetYaxis()->SetTitle("rate/Hz");
  h4->SetLineColor(4);
  h4->SetLineWidth(1.5);
  h5->GetXaxis()->SetTitle("x");
  h5->GetYaxis()->SetTitle("rate/Hz");
  h5->SetLineColor(4);
  h5->SetLineWidth(1.5);

  r0->SetFillColor(2);
  r1->SetFillColor(2);
  r2->SetFillColor(2);
  r3->SetFillColor(2);
  r4->SetFillColor(2);
  r5->SetFillColor(2);

  r0->SetLineColor(2);
  r1->SetLineColor(2);
  r2->SetLineColor(2);
  r3->SetLineColor(2);
  r4->SetLineColor(2);
  r5->SetLineColor(2);

  rh0->GetXaxis()->SetTitle("x");
  rh0->SetLineColor(4);
  rh1->GetXaxis()->SetTitle("x");
  rh1->SetLineColor(4);
  rh2->GetXaxis()->SetTitle("x");
  rh2->SetLineColor(4);
  rh3->GetXaxis()->SetTitle("x");
  rh3->SetLineColor(4);
  rh4->GetXaxis()->SetTitle("x");
  rh4->SetLineColor(4);
  rh5->GetXaxis()->SetTitle("x");
  rh5->SetLineColor(4);

  double x, Q2, z, W, Wp, Pt;

  while (Nrec < 1000000){
    Nsim++;
    if (Nsim%1000000 == 0) cout << Nsim << endl;
    vz = gRandom->Uniform(-370, -330);
    weight = mysidis.GenerateEvent(0, 1);
    if (weight > 0){
      acc = 0;
      accr = 0;
      lp = mysidis.GetLorentzVector("lp");
      Ph = mysidis.GetLorentzVector("Ph");

      x = mysidis.GetVariable("x");
      Q2 = mysidis.GetVariable("Q2");
      z = mysidis.GetVariable("z");
      Pt = mysidis.GetVariable("Pt");
      W = mysidis.GetVariable("W");
      Wp = mysidis.GetVariable("Wp");

      if (Q2 < 1.0) continue;
      if (W < 2.3) continue;
      if (Wp < 1.6) continue;
      
      if (ic == 1)
	acc = GetAcceptance(lp, "e-") * GetAcceptance(Ph, "h+");
      else if (ic == -1)
	acc = GetAcceptance(lp, "e-") * GetAcceptance(Ph, "h-");

      if (acc > 0){
	Nrec++;
	accr = acc;
	if(lp.Angle(Ph.Vect()) < 2.5 / 180.0 * M_PI)
	  accr = acc / 10.0;
	
	//cout << acc << endl;
	if (Pt < 0.2){
	  h0->Fill(x, weight * acc);
	  r0->Fill(x, weight * accr);
	}
	else if (Pt < 0.4){
	  h1->Fill(x, weight * acc);
	  r1->Fill(x, weight * accr);
	}
	else if (Pt < 0.6){
	  h2->Fill(x, weight * acc);
	  r2->Fill(x, weight * accr);
	}
	else if (Pt < 0.8){
	  h3->Fill(x, weight * acc);
	  r3->Fill(x, weight * accr);
	}
	else if (Pt < 1.0){
	  h4->Fill(x, weight * acc);
	  r4->Fill(x, weight * accr);
	}
	else {
	  h5->Fill(x, weight * acc);
	  r5->Fill(x, weight * accr);
	}
      }
    }
  }

  for (int i = 1; i < 36; i++){
    if (h0->GetBinContent(i) > 0){
      rh0->SetBinContent(i, r0->GetBinContent(i)/h0->GetBinContent(i));
    }
    if (h1->GetBinContent(i) > 0){
      rh1->SetBinContent(i, r1->GetBinContent(i)/h1->GetBinContent(i));
    }
    if (h2->GetBinContent(i) > 0){
      rh2->SetBinContent(i, r2->GetBinContent(i)/h2->GetBinContent(i));
    }
    if (h3->GetBinContent(i) > 0){
      rh3->SetBinContent(i, r3->GetBinContent(i)/h3->GetBinContent(i));
    }
    if (h4->GetBinContent(i) > 0){
      rh4->SetBinContent(i, r4->GetBinContent(i)/h4->GetBinContent(i));
    }
    if (h5->GetBinContent(i) > 0){
      rh5->SetBinContent(i, r5->GetBinContent(i)/h5->GetBinContent(i));
    }
  }
  rh0->GetXaxis()->SetTitle("x");
  rh0->GetYaxis()->SetRangeUser(0.1, 1.1);
  rh1->GetXaxis()->SetTitle("x");
  rh1->GetYaxis()->SetRangeUser(0.1, 1.1);
  rh2->GetXaxis()->SetTitle("x");
  rh2->GetYaxis()->SetRangeUser(0.1, 1.1);
  rh3->GetXaxis()->SetTitle("x");
  rh3->GetYaxis()->SetRangeUser(0.1, 1.1);
  rh4->GetXaxis()->SetTitle("x");
  rh4->GetYaxis()->SetRangeUser(0.1, 1.1);
  rh5->GetXaxis()->SetTitle("x");
  rh5->GetYaxis()->SetRangeUser(0.1, 1.1);
  


  h0->Scale(lumi/Nsim);
  h1->Scale(lumi/Nsim);
  h2->Scale(lumi/Nsim);
  h3->Scale(lumi/Nsim);
  h4->Scale(lumi/Nsim);
  h5->Scale(lumi/Nsim);

  r0->Scale(lumi/Nsim);
  r1->Scale(lumi/Nsim);
  r2->Scale(lumi/Nsim);
  r3->Scale(lumi/Nsim);
  r4->Scale(lumi/Nsim);
  r5->Scale(lumi/Nsim);

  TLegend * leg = new TLegend(0.72, 0.75, 0.9, 0.9);
  leg->AddEntry(h0, "before cut", "l");
  leg->AddEntry(r0, "after cut", "f");

  TCanvas * c0 = new TCanvas("c0", "", 800, 900);
  c0->SetTitle(Form("Rate before and after 2.5deg cut [%.1f<Q2<%.1f, %.2f<z<%.2f]", Q2min, Q2max, zmin, zmax));
  c0->SetTitle("title");

  c0->Divide(2,3);

  c0->cd(1);
  leg->SetHeader("0.0<P_{hT}<0.2");
  h0->Draw();
  r0->Draw("same");
  leg->DrawClone("same");

  c0->cd(2);
  leg->SetHeader("0.2<P_{hT}<0.4");
  h1->Draw();
  r1->Draw("bsame");
  leg->DrawClone("same");
    
  c0->cd(3);
  leg->SetHeader("0.4<P_{hT}<0.6");
  h2->Draw();
  r2->Draw("bsame");
  leg->DrawClone("same");

  c0->cd(4);
  leg->SetHeader("0.6<P_{hT}<0.8");
  h3->Draw();
  r3->Draw("bsame");
  leg->DrawClone("same");

  c0->cd(5);
  leg->SetHeader("0.8<P_{hT}<1.0");
  h4->Draw();
  r4->Draw("bsame");
  leg->DrawClone("same");

  c0->cd(6);
  leg->SetHeader("P_{hT}>1.0");
  h5->Draw();
  r5->Draw("bsame");
  leg->DrawClone("same");

  if (strcmp(argv[1], "pi+") == 0){
    c0->Print(Form("gallery/sidis_3he_pip_11GeV_Q2(%.1f,%.1f)_z(%.2f,%.2f)_2.5degcut.pdf", Q2min, Q2max, zmin, zmax));
  }
  else if (strcmp(argv[1], "pi-") == 0){
    c0->Print(Form("gallery/sidis_3he_pim_11GeV_Q2(%.1f,%.1f)_z(%.2f,%.2f)_2.5degcut.pdf", Q2min, Q2max, zmin, zmax));
  }

  TLegend * leg1 = new TLegend(0.72, 0.75, 0.9, 0.9);
  TCanvas * c1 = new TCanvas("c1", "", 800, 900);
  c1->Divide(2,3);

  c1->cd(1);
  rh0->Draw();
  leg1->SetHeader("0.0<P_{hT}<0.2");
  leg1->DrawClone("same");
  c1->cd(2);
  rh1->Draw();
  leg1->SetHeader("0.2<P_{hT}<0.4");
  leg1->DrawClone("same");
  c1->cd(3);
  rh2->Draw();
  leg1->SetHeader("0.4<P_{hT}<0.6");
  leg1->DrawClone("same");
  c1->cd(4);
  rh3->Draw();
  leg1->SetHeader("0.6<P_{hT}<0.8");
  leg1->DrawClone("same");
  c1->cd(5);
  rh4->Draw();
  leg1->SetHeader("0.8<P_{hT}<1.0");
  leg1->DrawClone("same");
  c1->cd(6);
  rh5->Draw();
  leg1->SetHeader("P_{hT}>1.0");
  leg1->DrawClone("same");

  if (strcmp(argv[1], "pi+") == 0){
    c1->Print(Form("gallery/sidis_3he_pip_11GeV_Q2(%.1f,%.1f)_z(%.2f,%.2f)_2.5degcut_ratio.pdf", Q2min, Q2max, zmin, zmax));
  }
  else if (strcmp(argv[1], "pi-") == 0){
    c1->Print(Form("gallery/sidis_3he_pim_11GeV_Q2(%.1f,%.1f)_z(%.2f,%.2f)_2.5degcut_ratio.pdf", Q2min, Q2max, zmin, zmax));
  }

      
  return 0;
}
