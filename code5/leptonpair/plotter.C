#include "../Lcore.h"

int main(const int argc, const char * argv[]){

  if (argc < 4){
    cout << "./plotter <loadfile> <savefile> <opt>" << endl;
    return 0;
  }

  const int opt = atoi(argv[3]);
  gStyle->SetPalette(55);

  TH2D * hB = new TH2D("hB", "", 200, 0.0, 200.0, 200, 0.0, 200.0);
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

  if (opt == 1){
    TString savefile = argv[2];
    TFile * fs = new TFile(argv[1], "r");
    TH2D * PTheta_phi = (TH2D *) fs->Get("PTheta_phi");
    TH2D * PTheta_electron = (TH2D *) fs->Get("PTheta_electron");
    TH2D * PTheta_positron = (TH2D *) fs->Get("PTheta_positron");
    TH2D * ThetaTheta = (TH2D *) fs->Get("ThetaTheta");
    TH2D * PP = (TH2D *) fs->Get("PP");
    TH2D * ThetaAngle = (TH2D *) fs->Get("ThetaAngle_pair");
    TH2D * PAngle = (TH2D *) fs->Get("PAngle_pair");
    TCanvas * c0 =  new TCanvas("c0", "", 800, 600);
    c0->SetLeftMargin(0.15);
    c0->SetBottomMargin(0.15);
    hB->SetTitle("");
    hB->GetXaxis()->SetTitle("#theta_{#phi} (deg)");
    hB->GetXaxis()->SetRangeUser(0.0, 180.0);
    hB->GetYaxis()->SetTitle("P_{#phi} (GeV)");
    hB->GetYaxis()->SetRangeUser(0.0, 3.0);
    hB->SetMaximum(PTheta_phi->GetMaximum());
    hB->DrawClone("axis");
    PTheta_phi->DrawClone("colsame");
    c0->Print(savefile + "(", "pdf");
    hB->GetXaxis()->SetTitle("#theta_{e^{-}} (deg)");
    hB->GetYaxis()->SetTitle("P_{e^{-}} (GeV)");
    hB->SetMaximum(PTheta_electron->GetMaximum());
    hB->DrawClone("axis");
    PTheta_electron->DrawClone("colsame");
    c0->Print(savefile, "pdf");
    hB->GetXaxis()->SetTitle("#theta_{e^{+}} (deg)");
    hB->GetYaxis()->SetTitle("P_{e^{+}} (GeV)");
    hB->SetMaximum(PTheta_positron->GetMaximum());
    hB->DrawClone("axis");
    PTheta_positron->DrawClone("colsame");
    c0->Print(savefile, "pdf");
    hB->GetXaxis()->SetTitle("#theta_{e^{-}} (deg)");
    hB->GetYaxis()->SetTitle("#theta_{e^{+}} (deg)");
    hB->GetYaxis()->SetRangeUser(0.0, 180.0);
    hB->SetMaximum(ThetaTheta->GetMaximum());
    hB->DrawClone("axis");
    ThetaTheta->DrawClone("colsame");
    c0->Print(savefile, "pdf");
    hB->GetXaxis()->SetTitle("#theta_{e^{-}} (deg)");
    hB->GetYaxis()->SetTitle("#alpha_{<e^{-}e^{+}>} (deg)");
    hB->SetMaximum(ThetaAngle->GetMaximum());
    hB->DrawClone("axis");
    ThetaAngle->DrawClone("colsame");
    c0->Print(savefile, "pdf");
    hB->GetXaxis()->SetTitle("P_{e^{-}} (GeV)");
    hB->GetXaxis()->SetRangeUser(0.0, 3.0);
    hB->GetYaxis()->SetTitle("#alpha_{<e^{-}e^{+}>} (deg)");
    hB->SetMaximum(PAngle->GetMaximum());
    hB->DrawClone("axis");
    PAngle->DrawClone("colsame");
    c0->Print(savefile, "pdf");
    hB->GetXaxis()->SetTitle("P_{e^{-}} (GeV)");
    hB->GetYaxis()->SetRangeUser(0.0, 3.0);
    hB->GetYaxis()->SetTitle("P_{e^{+}} (GeV)");
    hB->SetMaximum(PAngle->GetMaximum());
    hB->DrawClone("axis");
    PP->DrawClone("colsame");
    c0->Print(savefile + ")", "pdf");
  }

  return 0;
}
    
    
    
