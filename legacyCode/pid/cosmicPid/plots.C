void plots(){
  gStyle->SetOptStat(0);

  // TO USE DO THE FOLLOWING 3 THINGS...

  // 1. Specify the two input files, the nue_all and numu_all output files from the PID...
  const char* signalFile = "/unix/chips/jtingey/CHIPS/code/WCSimAnalysis/pid/cosmicPid/numu_cc_combined.root";
  const char* bkgFile = "/unix/chips/jtingey/CHIPS/code/WCSimAnalysis/pid/cosmicPid/numu_cr_combined.root";


  TChain * PIDTree_ann_signal = new TChain("PIDTree_ann", "PIDTree_ann");
  TChain * PIDTree_ann_bkg = new TChain("PIDTree_ann", "PIDTree_ann");

  PIDTree_ann_signal->Add(signalFile);
  PIDTree_ann_bkg->Add(bkgFile);

  //Define all the histograms we want...

  //First we want All ANN variables plots without the ANN cuts and just preselection...

  /*
  // All events
  TH1F* hNumuCC= new TH1F("hNumuCC",";vtxZ;Fraction of Events", 100, -1000, 1000);
  hNumuCC->SetFillColor(kRed);
  hNumuCC->SetFillStyle(3005);
  hNumuCC->SetLineColor(kRed);
  hNumuCC->GetYaxis()->CenterTitle();
  hNumuCC->GetYaxis()->SetTitleOffset(1.3);
  hNumuCC->GetXaxis()->CenterTitle();
  TH1F* hNumuCR = new TH1F("hNumuCR",";vtxZ;Fraction of Events", 100, -1000, 1000);
  hNumuCR->SetFillColor(kBlue);
  hNumuCR->SetFillStyle(3004);
  hNumuCR->SetLineColor(kBlue);
  hNumuCR->GetYaxis()->CenterTitle();
  hNumuCR->GetXaxis()->CenterTitle();

  // Define selection and plot from PIDTree the All ANN variables plots without the ANN cuts but with preselection...
  PIDTree_ann_signal->Draw("fVtxZ>>hNumuCC");
  PIDTree_ann_bkg->Draw("fVtxZ>>hNumuCR");


  // All events
  TH1F* hNumuCC= new TH1F("hNumuCC",";dirX;Fraction of Events", 100, -1, 1);
  hNumuCC->SetFillColor(kRed);
  hNumuCC->SetFillStyle(3005);
  hNumuCC->SetLineColor(kRed);
  hNumuCC->GetYaxis()->CenterTitle();
  hNumuCC->GetYaxis()->SetTitleOffset(1.3);
  hNumuCC->GetXaxis()->CenterTitle();
  TH1F* hNumuCR = new TH1F("hNumuCR",";dirX;Fraction of Events", 100, -1, 1);
  hNumuCR->SetFillColor(kBlue);
  hNumuCR->SetFillStyle(3004);
  hNumuCR->SetLineColor(kBlue);
  hNumuCR->GetYaxis()->CenterTitle();
  hNumuCR->GetXaxis()->CenterTitle();

  // Define selection and plot from PIDTree the All ANN variables plots without the ANN cuts but with preselection...
  PIDTree_ann_signal->Draw("fDirX>>hNumuCC");
  PIDTree_ann_bkg->Draw("fDirX>>hNumuCR");


  // All events
  TH1F* hNumuCC= new TH1F("hNumuCC",";dirZ;Fraction of Events", 100, -1, 1);
  hNumuCC->SetFillColor(kRed);
  hNumuCC->SetFillStyle(3005);
  hNumuCC->SetLineColor(kRed);
  hNumuCC->GetYaxis()->CenterTitle();
  hNumuCC->GetYaxis()->SetTitleOffset(1.3);
  hNumuCC->GetXaxis()->CenterTitle();
  TH1F* hNumuCR = new TH1F("hNumuCR",";dirZ;Fraction of Events", 100, -1, 1);
  hNumuCR->SetFillColor(kBlue);
  hNumuCR->SetFillStyle(3004);
  hNumuCR->SetLineColor(kBlue);
  hNumuCR->GetYaxis()->CenterTitle();
  hNumuCR->GetXaxis()->CenterTitle();

  // Define selection and plot from PIDTree the All ANN variables plots without the ANN cuts but with preselection...
  PIDTree_ann_signal->Draw("fDirZ>>hNumuCC");
  PIDTree_ann_bkg->Draw("fDirZ>>hNumuCR");
  */

  // All events
  TH1F* hNumuCC= new TH1F("hNumuCC",";vtxR;Fraction of Events", 100, 0, 1500);
  hNumuCC->SetFillColor(kRed);
  hNumuCC->SetFillStyle(3005);
  hNumuCC->SetLineColor(kRed);
  hNumuCC->GetYaxis()->CenterTitle();
  hNumuCC->GetYaxis()->SetTitleOffset(1.3);
  hNumuCC->GetXaxis()->CenterTitle();
  TH1F* hNumuCR = new TH1F("hNumuCR",";vtxR;Fraction of Events", 100, 0, 1500);
  hNumuCR->SetFillColor(kBlue);
  hNumuCR->SetFillStyle(3004);
  hNumuCR->SetLineColor(kBlue);
  hNumuCR->GetYaxis()->CenterTitle();
  hNumuCR->GetXaxis()->CenterTitle();

  // Define selection and plot from PIDTree the All ANN variables plots without the ANN cuts but with preselection...
  PIDTree_ann_signal->Draw("fvtxR>>hNumuCC");
  PIDTree_ann_bkg->Draw("fvtxR>>hNumuCR");

  // Draw The All ANN ElMu Variable Plot without cuts but with preselection...
  TCanvas *c1 = new TCanvas("c1", "", 800, 600);
  c1->cd();

  hNumuCC->Scale(1/hNumuCC->GetEntries());
  hNumuCR->Scale(1/hNumuCR->GetEntries());

  hNumuCC->Draw();
  hNumuCR->Draw("SAME");

  TLegend * leg = new TLegend(0.35, 0.83, 0.65, 0.58);
  leg->AddEntry(hNumuCC, "Numu Beam Events", "FL");
  leg->AddEntry(hNumuCR, "Numu Cosmic Events", "FL");
  leg->SetTextFont(42);
  leg->Draw("SAME");

  c1->Update();
  //c1->SaveAs("vtxZ.png");
  //c1->SaveAs("plot.root");
  //c1->SaveAs("plot.pdf");
  //c1->SaveAs("plot.C");
}
