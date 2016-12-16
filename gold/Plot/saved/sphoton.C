{
//=========Macro generated from canvas: cwf/
//=========  (Fri Dec 16 03:27:29 2016) by ROOT version5.34/21
   TCanvas *cwf = new TCanvas("cwf", "",0,0,800,600);
   cwf->SetHighLightColor(2);
   cwf->Range(-0.4004392,-0.032,2.269155,0.1813333);
   cwf->SetFillColor(0);
   cwf->SetBorderMode(0);
   cwf->SetBorderSize(2);
   cwf->SetLeftMargin(0.15);
   cwf->SetBottomMargin(0.15);
   cwf->SetFrameBorderMode(0);
   cwf->SetFrameBorderMode(0);
   
   cwf->SetLogy();

   TGraph *graph = new TGraph(11);
   graph->SetName("Graph0");
   graph->SetTitle("");
   graph->SetFillColor(1);
   graph->SetLineColor(4);
   graph->SetLineWidth(2);

   graph->SetPoint(0,1.0,0.676939/2);
   graph->SetPoint(1,1.1,64.144/2);
   graph->SetPoint(2,1.2,795.547/2);
   graph->SetPoint(3,1.3,3424.02/2);
   graph->SetPoint(4,1.4,6171.96/2);
   graph->SetPoint(5,1.5,7131.48/2);
   graph->SetPoint(6,1.6,4981.91/2);  
   graph->SetPoint(7,1.7,3805.06/2);
   graph->SetPoint(8,1.8,2255.33/2);
   graph->SetPoint(9,1.9,1783.44/2);
   graph->SetPoint(10,2,755.647/2);

   
   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","",1000,0,2.1978);
   Graph_Graph1->SetMinimum(0);
   Graph_Graph1->SetMaximum(0.16);
   Graph_Graph1->SetDirectory(0);
   Graph_Graph1->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   Graph_Graph1->SetLineColor(ci);
   Graph_Graph1->GetXaxis()->SetTitle("E_{#gamma} / GeV");
   Graph_Graph1->GetXaxis()->SetRangeUser(1,2);
   Graph_Graph1->GetXaxis()->CenterTitle(true);
   Graph_Graph1->GetXaxis()->SetLabelFont(42);
   Graph_Graph1->GetXaxis()->SetLabelSize(0.055);
   Graph_Graph1->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph1->GetXaxis()->SetTitleFont(42);
   Graph_Graph1->GetYaxis()->SetTitle("#sigma / nb");
   Graph_Graph1->GetYaxis()->CenterTitle(true);
   //Graph_Graph1->GetYaxis()->SetNdivisions(505);
   Graph_Graph1->GetYaxis()->SetLabelFont(42);
   Graph_Graph1->GetYaxis()->SetLabelSize(0.055);
   Graph_Graph1->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph1->GetYaxis()->SetTitleOffset(1.15);
   Graph_Graph1->GetYaxis()->SetTitleFont(42);
   Graph_Graph1->GetYaxis()->SetRangeUser(1,10000);
   Graph_Graph1->GetZaxis()->SetLabelFont(42);
   Graph_Graph1->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph1->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph1);
   
   graph->Draw("ac");
   cwf->Modified();
   cwf->cd();
   cwf->SetSelected(cwf);
   
   cwf->Print("sphoton.pdf");

}
