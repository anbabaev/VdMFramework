{
//=========Macro generated from canvas: c1_n2/c1_n2
//=========  (Sun Feb  7 16:33:53 2016) by ROOT version5.34/32
   TCanvas *c1_n2 = new TCanvas("c1_n2", "c1_n2",0,0,700,500);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   c1_n2->SetHighLightColor(2);
   c1_n2->Range(-1.25,26.63501,12.25,34.46125);
   c1_n2->SetFillColor(0);
   c1_n2->SetBorderMode(0);
   c1_n2->SetBorderSize(2);
   c1_n2->SetFrameBorderMode(0);
   c1_n2->SetFrameBorderMode(0);
   
   TGraphErrors *gre = new TGraphErrors(50);
   gre->SetName("");
   gre->SetTitle("4266 PCC DGConst area");
   gre->SetFillColor(1);
   gre->SetMarkerStyle(8);
   gre->SetMarkerSize(0.4);
   gre->SetPoint(0,10,28.83778);
   gre->SetPointError(0,0,0);
   gre->SetPoint(1,10,28.52146);
   gre->SetPointError(1,0,0);
   gre->SetPoint(2,10,28.35742);
   gre->SetPointError(2,0,0);
   gre->SetPoint(3,10,28.20292);
   gre->SetPointError(3,0,0);
   gre->SetPoint(4,10,29.48801);
   gre->SetPointError(4,0,0);
   gre->SetPoint(5,4,31.70236);
   gre->SetPointError(5,0,0);
   gre->SetPoint(6,4,31.46482);
   gre->SetPointError(6,0,0);
   gre->SetPoint(7,4,31.11724);
   gre->SetPointError(7,0,0);
   gre->SetPoint(8,4,30.85766);
   gre->SetPointError(8,0,0);
   gre->SetPoint(9,4,32.52116);
   gre->SetPointError(9,0,0);
   gre->SetPoint(10,5,31.37033);
   gre->SetPointError(10,0,0);
   gre->SetPoint(11,5,31.25076);
   gre->SetPointError(11,0,0);
   gre->SetPoint(12,5,30.74589);
   gre->SetPointError(12,0,0);
   gre->SetPoint(13,5,31.07207);
   gre->SetPointError(13,0,0);
   gre->SetPoint(14,5,31.73379);
   gre->SetPointError(14,0,0);
   gre->SetPoint(15,6,28.491);
   gre->SetPointError(15,0,0);
   gre->SetPoint(16,6,28.57869);
   gre->SetPointError(16,0,0);
   gre->SetPoint(17,6,28.04914);
   gre->SetPointError(17,0,0);
   gre->SetPoint(18,6,27.97835);
   gre->SetPointError(18,0,0);
   gre->SetPoint(19,6,29.37269);
   gre->SetPointError(19,0,0);
   gre->SetPoint(20,7,32.16322);
   gre->SetPointError(20,0,0);
   gre->SetPoint(21,7,31.78767);
   gre->SetPointError(21,0,0);
   gre->SetPoint(22,7,31.29637);
   gre->SetPointError(22,0,0);
   gre->SetPoint(23,7,31.49588);
   gre->SetPointError(23,0,0);
   gre->SetPoint(24,7,32.97174);
   gre->SetPointError(24,0,0);
   gre->SetPoint(25,1,31.26983);
   gre->SetPointError(25,0,0);
   gre->SetPoint(26,1,30.93735);
   gre->SetPointError(26,0,0);
   gre->SetPoint(27,1,30.63303);
   gre->SetPointError(27,0,0);
   gre->SetPoint(28,1,30.68255);
   gre->SetPointError(28,0,0);
   gre->SetPoint(29,1,31.96409);
   gre->SetPointError(29,0,0);
   gre->SetPoint(30,2,28.88354);
   gre->SetPointError(30,0,0);
   gre->SetPoint(31,2,28.77894);
   gre->SetPointError(31,0,0);
   gre->SetPoint(32,2,28.5452);
   gre->SetPointError(32,0,0);
   gre->SetPoint(33,2,28.38309);
   gre->SetPointError(33,0,0);
   gre->SetPoint(34,2,29.7435);
   gre->SetPointError(34,0,0);
   gre->SetPoint(35,3,28.83448);
   gre->SetPointError(35,0,0);
   gre->SetPoint(36,3,28.70569);
   gre->SetPointError(36,0,0);
   gre->SetPoint(37,3,28.22357);
   gre->SetPointError(37,0,0);
   gre->SetPoint(38,3,28.15763);
   gre->SetPointError(38,0,0);
   gre->SetPoint(39,3,29.55198);
   gre->SetPointError(39,0,0);
   gre->SetPoint(40,8,28.67321);
   gre->SetPointError(40,0,0);
   gre->SetPoint(41,8,28.45883);
   gre->SetPointError(41,0,0);
   gre->SetPoint(42,8,28.05395);
   gre->SetPointError(42,0,0);
   gre->SetPoint(43,8,27.93939);
   gre->SetPointError(43,0,0);
   gre->SetPoint(44,8,29.06939);
   gre->SetPointError(44,0,0);
   gre->SetPoint(45,9,32.2132);
   gre->SetPointError(45,0,0);
   gre->SetPoint(46,9,31.90743);
   gre->SetPointError(46,0,0);
   gre->SetPoint(47,9,31.91434);
   gre->SetPointError(47,0,0);
   gre->SetPoint(48,9,31.72686);
   gre->SetPointError(48,0,0);
   gre->SetPoint(49,9,33.15687);
   gre->SetPointError(49,0,0);
   
   TH1F *Graph_Graph95 = new TH1F("Graph_Graph95","4266 PCC DGConst area",100,0.1,10.9);
   Graph_Graph95->SetMinimum(27.41764);
   Graph_Graph95->SetMaximum(33.67862);
   Graph_Graph95->SetDirectory(0);
   Graph_Graph95->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   Graph_Graph95->SetLineColor(ci);
   Graph_Graph95->GetXaxis()->SetTitle("Scan");
   Graph_Graph95->GetXaxis()->SetLabelFont(42);
   Graph_Graph95->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph95->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph95->GetXaxis()->SetTitleFont(42);
   Graph_Graph95->GetYaxis()->SetTitle("area");
   Graph_Graph95->GetYaxis()->SetLabelFont(42);
   Graph_Graph95->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph95->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph95->GetYaxis()->SetTitleFont(42);
   Graph_Graph95->GetZaxis()->SetLabelFont(42);
   Graph_Graph95->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph95->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph95->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph95);
   
   gre->Draw("ap");
   
   TPaveText *pt = new TPaveText(0.2910057,0.94,0.7089943,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   TText *text = pt->AddText("4266 PCC DGConst area");
   pt->Draw();
   c1_n2->Modified();
   c1_n2->cd();
   c1_n2->SetSelected(c1_n2);
}
