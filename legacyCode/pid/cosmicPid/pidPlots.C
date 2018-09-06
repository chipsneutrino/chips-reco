void pidPlots(){
  gStyle->SetOptStat(0);

  // TO USE DO THE FOLLOWING 3 THINGS...

  // 1. Specify the two input files, the nue_all and numu_all output files from the PID...
  const char* tmvaFile = "/unix/chips/jtingey/CHIPS/code/WCSimAnalysis/pid/cosmicPid/weights/tmva_cosmic.root";

  TChain * PIDTree_ann = new TChain("TestTree", "TestTree");
  PIDTree_ann->Add(tmvaFile);


  //Define all the histograms we want...

  //First we want All ANN variables plots without the ANN cuts and just preselection...

  // All events
  TH1F* hNumuCC= new TH1F("hNumuCC",";numuBeam vs. numuCosmic ANN value;Fraction of Events", 100, -0.5, 1.5);
  hNumuCC->SetFillColor(kRed);
  hNumuCC->SetFillStyle(3001);
  hNumuCC->SetLineColor(kRed);
  hNumuCC->GetYaxis()->CenterTitle();
  hNumuCC->GetYaxis()->SetTitleOffset(1.3);
  hNumuCC->GetXaxis()->CenterTitle();
  TH1F* hNumuCR = new TH1F("hNumuCR",";numuBeam vs. numuCosmic ANN value;Fraction of Events", 100, -0.5, 1.5);
  hNumuCR->SetFillColor(kBlue);
  hNumuCR->SetFillStyle(3001);
  hNumuCR->SetLineColor(kBlue);
  hNumuCR->GetYaxis()->CenterTitle();
  hNumuCR->GetXaxis()->CenterTitle();

  // Define selection and plot from PIDTree the All ANN variables plots without the ANN cuts but with preselection...
  TString string_hNumuCC = Form("classID==0");
  TString string_hNumuCR = Form("classID==1");
  PIDTree_ann->Draw("MLP>>hNumuCC", string_hNumuCC.Data());
  PIDTree_ann->Draw("MLP>>hNumuCR", string_hNumuCR.Data());

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
  c1->SaveAs("plot.png");
  c1->SaveAs("plot.root");
  c1->SaveAs("plot.pdf");
  c1->SaveAs("plot.C");
}
