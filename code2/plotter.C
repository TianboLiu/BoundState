#include "bound.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLatex.h"

using namespace std;

int main(){
  
  if (true){//wave function and potential
    readgrid("wf.dat", "Veff.dat");
    urone();
    TGraph * gwf = new TGraph(100);
    double r, u, V;
    for (int i = 0; i < 100; i++){
      r = i * 0.05 * GeVfm;
      u = ur(r);
      gwf->SetPoint(i, 0.05*i, u);
    }
    TCanvas * cwf = new TCanvas("cwf", "u(r)", 800, 600);
    cwf->SetLeftMargin(0.15);
    cwf->SetBottomMargin(0.15);
    
    gwf->SetTitle("wave function");
    gwf->SetLineColor(4);
    gwf->SetLineWidth(2.5);
    gwf->GetXaxis()->SetTitle("r / fm");
    gwf->GetXaxis()->SetTitleSize(0.06);
    gwf->GetXaxis()->SetLabelSize(0.055);
    gwf->GetXaxis()->CenterTitle(true);
    gwf->GetXaxis()->SetRangeUser(0.0, 5.0);
    gwf->GetYaxis()->SetTitle("u(r) / GeV^{1/2}");
    gwf->GetYaxis()->SetTitleSize(0.06);
    gwf->GetYaxis()->SetTitleOffset(1.15);
    gwf->GetYaxis()->SetLabelSize(0.055);
    gwf->GetYaxis()->CenterTitle(true);
    gwf->GetYaxis()->SetRangeUser(0.0, 0.5);
    gwf->GetYaxis()->SetNdivisions(5, 5, 0, kTRUE);
    gwf->Draw("AC");
    cwf->Print("wf.pdf");
    cwf->Close();
    gwf->Delete();

    TGraph * gV = new TGraph(100);
    for (int i = 0; i < 100; i++){
      r = i * 0.05 * GeVfm;
      V = Veff(r);
      gV->SetPoint(i, 0.05*i, V);
    }
    TCanvas * cV = new TCanvas("V", "V(r)", 800, 600);
    gV->SetTitle("effective potential");
    
    gV->SetLineColor(4);
    gV->SetLineWidth(2.5);
    gV->GetXaxis()->SetTitle("r / fm");
    gV->GetXaxis()->SetTitleSize(0.06);
    gV->GetXaxis()->SetLabelSize(0.055);
    gV->GetXaxis()->CenterTitle();
    gV->GetXaxis()->SetTitleOffset(1.15);
    gV->GetXaxis()->SetRangeUser(0.0, 5.0);
    gV->GetYaxis()->SetTitle("V_{eff }/GeV");
    gV->GetYaxis()->SetTitleSize(0.06);
    gV->GetYaxis()->SetTitleOffset(1.15);
    gV->GetYaxis()->SetLabelSize(0.055);
    gV->GetYaxis()->CenterTitle();
    gV->GetYaxis()->SetRangeUser(-0.6, 0.2);
    gV->GetYaxis()->SetNdivisions(4, 5, 0, kTRUE);
    gV->Draw("AC");
    cV->Print("Veff.pdf");
    cV->Close();
    gV->Delete();
  }

  return 0;
}
