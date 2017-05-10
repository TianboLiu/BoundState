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

  return 0;
}
