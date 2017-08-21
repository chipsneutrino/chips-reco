void pidPlots(){
  gStyle->SetOptStat(0);

  // TO USE DO THE FOLLOWING 3 THINGS...

  // 1. Specify the two input files, the nue_all and numu_all output files from the PID...
  const char* nueFileName = "/unix/chips/jtingey/CHIPS/data/geometries/CHIPS_10kton_3inch_HQE_6perCent/PID/output/CHIPS_10kton_3inch_HQE_6perCent_nuelAll_readOutput.root";
  const char* numuFileName = "/unix/chips/jtingey/CHIPS/data/geometries/CHIPS_10kton_3inch_HQE_6perCent/PID/output/CHIPS_10kton_3inch_HQE_6perCent_numuAll_readOutput.root";

  // 2. Give the cut values on both PID variables, use wc_scanANNCuts.py to find these...
  float ElMuCut = 0.96;
  float ElNcCut = 0.81;

  // 3. Change the cuts in these two TPaveText's so they are printed nicely on the cut plots...
  TPaveText * pt1 = new TPaveText(0.18, 0.89, 0.845, 0.84, "blNDC");
  pt1->AddText("Preselection & e vs. NC ANN > 0.81 cuts applied to initial sample of 95% #nu_{#mu}, 5% #nu_{e}");

  TPaveText * pt2 = new TPaveText(0.18, 0.89, 0.845, 0.84, "blNDC");
  pt2->AddText("Preselection & e vs. #mu ANN > 0.96 cuts applied to initial sample of 95% #nu_{#mu}, 5% #nu_{e}");


  TChain * PIDTree_ann = new TChain("PIDTree_ann", "PIDTree_ann");
  PIDTree_ann->Add(nueFileName);
  PIDTree_ann->Add(numuFileName);

  //Define all the histograms we want...

  //First we want All ANN variables plots without the ANN cuts and just preselection...

  // All events
  TH1F* hAll_ElMu= new TH1F("hAll_ElMu",";numuCCQE vs. nueCCQE ANN value;Events (arb. units)", 100, -0.5, 1.5);
  hAll_ElMu->SetFillColor(kGray);
  hAll_ElMu->SetFillStyle(3001);
  hAll_ElMu->SetLineColor(kGray);
  hAll_ElMu->GetYaxis()->CenterTitle();
  hAll_ElMu->GetXaxis()->CenterTitle();
  TH1F* hAll_ElNc = new TH1F("hAll_ElNc",";NC vs. nueCCQE ANN value;Events (arb. units)", 100, -0.5, 1.5);
  hAll_ElNc->SetFillColor(kGray);
  hAll_ElNc->SetFillStyle(3001);
  hAll_ElNc->SetLineColor(kGray);
  hAll_ElNc->GetYaxis()->CenterTitle();
  hAll_ElNc->GetXaxis()->CenterTitle();
  // All NC events
  TH1F* hNC_ElMu = new TH1F("hNC_ElMu",";numuCCQE vs. nueCCQE ANN value;Events (arb. units)", 100, -0.5, 1.5);
  hNC_ElMu->SetLineColor(kRed-2);
  hNC_ElMu->GetYaxis()->CenterTitle();
  hNC_ElMu->GetXaxis()->CenterTitle();
  TH1F* hNC_ElNc = new TH1F("hNC_ElNc",";NC vs. nueCCQE ANN value;Events (arb. units)", 100, -0.5, 1.5);
  hNC_ElNc->SetLineColor(kRed-2);
  hNC_ElNc->GetYaxis()->CenterTitle();
  hNC_ElNc->GetXaxis()->CenterTitle();
  // Muon CC events
  TH1F* hNumuCC_ElMu = new TH1F("hNumuCC_ElMu",";numuCCQE vs. nueCCQE ANN value;Events (arb. units)", 100, -0.5, 1.5);
  hNumuCC_ElMu->SetLineColor(kRed);
  hNumuCC_ElMu->GetYaxis()->CenterTitle();
  hNumuCC_ElMu->GetXaxis()->CenterTitle();
  TH1F* hNumuCC_ElNc = new TH1F("hNumuCC_ElNc",";NC vs. nueCCQE ANN value;Events (arb. units)", 100, -0.5, 1.5);
  hNumuCC_ElNc->SetLineColor(kRed);
  hNumuCC_ElNc->GetYaxis()->CenterTitle();
  hNumuCC_ElNc->GetXaxis()->CenterTitle();
  // Electron CC events
  TH1F* hNueCC_ElMu = new TH1F("hNueCC_ElMu",";numuCCQE vs. nueCCQE ANN value;Events (arb. units)", 100, -0.5, 1.5);
  hNueCC_ElMu->SetLineColor(kSpring+2);
  hNueCC_ElMu->SetFillColor(kSpring+2);
  hNueCC_ElMu->GetYaxis()->CenterTitle();
  hNueCC_ElMu->GetXaxis()->CenterTitle();
  TH1F* hNueCC_ElNc = new TH1F("hNueCC_ElNc",";NC vs. nueCCQE ANN value;Events (arb. units)", 100, -0.5, 1.5);
  hNueCC_ElNc->SetLineColor(kSpring+2);
  hNueCC_ElNc->SetFillColor(kSpring+2);
  hNueCC_ElNc->GetYaxis()->CenterTitle();
  hNueCC_ElNc->GetXaxis()->CenterTitle();
  // Electron CCQE events
  TH1F* hNueCCQE_ElMu = new TH1F("hNueCCQE_ElMu",";numuCCQE vs. nueCCQE ANN value;Events (arb. units)", 100, -0.5, 1.5);
  hNueCCQE_ElMu->SetLineColor(kGreen+2);
  hNueCCQE_ElMu->SetFillColor(kGreen+2);
  hNueCCQE_ElMu->GetYaxis()->CenterTitle();
  hNueCCQE_ElMu->GetXaxis()->CenterTitle();
  TH1F* hNueCCQE_ElNc = new TH1F("hNueCCQE_ElNc",";NC vs. nueCCQE ANN value;Events (arb. units)", 100, -0.5, 1.5);
  hNueCCQE_ElNc->SetLineColor(kGreen+2);
  hNueCCQE_ElNc->SetFillColor(kGreen+2);
  hNueCCQE_ElNc->GetYaxis()->CenterTitle();
  hNueCCQE_ElNc->GetXaxis()->CenterTitle();

  //Now we want All ANN variables plots with the ANN cuts and preselection...

  // All events
  TH1F* hAll_ElMu_cuts= new TH1F("hAll_ElMu_cuts",";numuCCQE vs. nueCCQE ANN value;Events (arb. units)", 100, -0.5, 1.5);
  hAll_ElMu_cuts->SetFillColor(kGray);
  hAll_ElMu_cuts->SetFillStyle(3001);
  hAll_ElMu_cuts->SetLineColor(kGray);
  hAll_ElMu_cuts->GetYaxis()->CenterTitle();
  hAll_ElMu_cuts->GetXaxis()->CenterTitle();
  TH1F* hAll_ElNc_cuts = new TH1F("hAll_ElNc_cuts",";NC vs. nueCCQE ANN value;Events (arb. units)", 100, -0.5, 1.5);
  hAll_ElNc_cuts->SetFillColor(kGray);
  hAll_ElNc_cuts->SetFillStyle(3001);
  hAll_ElNc_cuts->SetLineColor(kGray);
  hAll_ElNc_cuts->GetYaxis()->CenterTitle();
  hAll_ElNc_cuts->GetXaxis()->CenterTitle();
  // All NC events
  TH1F* hNC_ElMu_cuts = new TH1F("hNC_ElMu_cuts",";numuCCQE vs. nueCCQE ANN value;Events (arb. units)", 100, -0.5, 1.5);
  hNC_ElMu_cuts->SetLineColor(kRed-2);
  hNC_ElMu_cuts->GetYaxis()->CenterTitle();
  hNC_ElMu_cuts->GetXaxis()->CenterTitle();
  TH1F* hNC_ElNc_cuts = new TH1F("hNC_ElNc_cuts",";NC vs. nueCCQE ANN value;Events (arb. units)", 100, -0.5, 1.5);
  hNC_ElNc_cuts->SetLineColor(kRed-2);
  hNC_ElNc_cuts->GetYaxis()->CenterTitle();
  hNC_ElNc_cuts->GetXaxis()->CenterTitle();
  // Muon CC events
  TH1F* hNumuCC_ElMu_cuts = new TH1F("hNumuCC_ElMu_cuts",";numuCCQE vs. nueCCQE ANN value;Events (arb. units)", 100, -0.5, 1.5);
  hNumuCC_ElMu_cuts->SetLineColor(kRed);
  hNumuCC_ElMu_cuts->GetYaxis()->CenterTitle();
  hNumuCC_ElMu_cuts->GetXaxis()->CenterTitle();
  TH1F* hNumuCC_ElNc_cuts = new TH1F("hNumuCC_ElNc_cuts",";NC vs. nueCCQE ANN value;Events (arb. units)", 100, -0.5, 1.5);
  hNumuCC_ElNc_cuts->SetLineColor(kRed);
  hNumuCC_ElNc_cuts->GetYaxis()->CenterTitle();
  hNumuCC_ElNc_cuts->GetXaxis()->CenterTitle();
  // Electron CC events
  TH1F* hNueCC_ElMu_cuts = new TH1F("hNueCC_ElMu_cuts",";numuCCQE vs. nueCCQE ANN value;Events (arb. units)", 100, -0.5, 1.5);
  hNueCC_ElMu_cuts->SetLineColor(kSpring+2);
  hNueCC_ElMu_cuts->SetFillColor(kSpring+2);
  hNueCC_ElMu_cuts->GetYaxis()->CenterTitle();
  hNueCC_ElMu_cuts->GetXaxis()->CenterTitle();
  TH1F* hNueCC_ElNc_cuts = new TH1F("hNueCC_ElNc_cuts",";NC vs. nueCCQE ANN value;Events (arb. units)", 100, -0.5, 1.5);
  hNueCC_ElNc_cuts->SetLineColor(kSpring+2);
  hNueCC_ElNc_cuts->SetFillColor(kSpring+2);
  hNueCC_ElNc_cuts->GetYaxis()->CenterTitle();
  hNueCC_ElNc_cuts->GetXaxis()->CenterTitle();
  // Electron CCQE events
  TH1F* hNueCCQE_ElMu_cuts = new TH1F("hNueCCQE_ElMu_cuts",";numuCCQE vs. nueCCQE ANN value;Events (arb. units)", 100, -0.5, 1.5);
  hNueCCQE_ElMu_cuts->SetLineColor(kGreen+2);
  hNueCCQE_ElMu_cuts->SetFillColor(kGreen+2);
  hNueCCQE_ElMu_cuts->GetYaxis()->CenterTitle();
  hNueCCQE_ElMu_cuts->GetXaxis()->CenterTitle();
  TH1F* hNueCCQE_ElNc_cuts = new TH1F("hNueCCQE_ElNc_cuts",";NC vs. nueCCQE ANN value;Events (arb. units)", 100, -0.5, 1.5);
  hNueCCQE_ElNc_cuts->SetLineColor(kGreen+2);
  hNueCCQE_ElNc_cuts->SetFillColor(kGreen+2);
  hNueCCQE_ElNc_cuts->GetYaxis()->CenterTitle();
  hNueCCQE_ElNc_cuts->GetXaxis()->CenterTitle();

  //Define Scaling Variables
  double normNumEvents = 1000;
  double scaleNue = 0.05;
  double scaleNumu = 0.95;

  //Calculate the WEIGHTS!!!

  Long64_t numNueEvents = PIDTree_ann->GetEntries("trueBeamPDG==12");
  Long64_t numNumuEvents = PIDTree_ann->GetEntries("trueBeamPDG==14");

  Long64_t numNueCCEvents = PIDTree_ann->GetEntries("trueBeamPDG==12 && trueCCEvent");
  Long64_t numNueNCEvents = PIDTree_ann->GetEntries("trueBeamPDG==12 && trueNCEvent");
  Long64_t numNumuCCEvents = PIDTree_ann->GetEntries("trueBeamPDG==14 && trueCCEvent");
  Long64_t numNumuNCEvents = PIDTree_ann->GetEntries("trueBeamPDG==14 && trueNCEvent");

  Long64_t numNueCCQEEvents = PIDTree_ann->GetEntries("trueBeamPDG==12 && trueCCEvent && trueQEEvent");
  Long64_t numNumuCCQEEvents = PIDTree_ann->GetEntries("trueBeamPDG==14 && trueCCEvent && trueQEEvent");
  Long64_t numNueCCnonQEEvents = PIDTree_ann->GetEntries("trueBeamPDG==12 && trueCCEvent && !trueQEEvent");
  Long64_t numNumuCCnonQEEvents = PIDTree_ann->GetEntries("trueBeamPDG==14 && trueCCEvent && !trueQEEvent");

  float weightNueEvents = normNumEvents * 0.05/numNueEvents;
  float weightNumuEvents = normNumEvents * 0.95/numNumuEvents;

  float weightNueCCEvents  = normNumEvents * (0.05 * 0.7)/numNueCCEvents;
  float weightNueNCEvents  = normNumEvents * (0.05 * 0.3)/numNueNCEvents;
  float weightNumuCCEvents = normNumEvents * (0.95 * 0.7)/numNumuCCEvents;
  float weightNumuNCEvents = normNumEvents * (0.95 * 0.3)/numNumuNCEvents;

  float weightNueCCQEEvents  = normNumEvents * (0.05 * 0.7 * 0.2)/numNueCCQEEvents;
  float weightNumuCCQEEvents = normNumEvents * (0.95 * 0.7 * 0.2)/numNumuCCQEEvents;
  float weightNueCCnonQEEvents  = normNumEvents * (0.05 * 0.7 * 0.8)/numNueCCnonQEEvents;
  float weightNumuCCnonQEEvents = normNumEvents * (0.95 * 0.7 * 0.8)/numNumuCCnonQEEvents;

  // Define selection and plot from PIDTree the All ANN variables plots without the ANN cuts but with preselection...
  TString nuAllString_ElMu = Form("((trueBeamPDG==14)*(trueCCEvent*(%f*trueQEEvent + %f*!trueQEEvent) + %f*trueNCEvent) + (trueBeamPDG == 12)*(trueCCEvent * (%f*trueQEEvent + %f * !trueQEEvent) + %f*trueNCEvent)) * (preselected)", weightNumuCCQEEvents, weightNumuCCnonQEEvents, weightNumuNCEvents, weightNueCCQEEvents, weightNueCCnonQEEvents, weightNueNCEvents);
  TString nuAllString_ElNc = Form("((trueBeamPDG==14)*(trueCCEvent*(%f*trueQEEvent + %f*!trueQEEvent) + %f*trueNCEvent) + (trueBeamPDG == 12)*(trueCCEvent * (%f*trueQEEvent + %f * !trueQEEvent) + %f*trueNCEvent)) * (preselected)", weightNumuCCQEEvents, weightNumuCCnonQEEvents, weightNumuNCEvents, weightNueCCQEEvents, weightNueCCnonQEEvents, weightNueNCEvents);
  PIDTree_ann->Draw("annNueCCQEvsNumuCCQE>>hAll_ElMu", nuAllString_ElMu.Data());
  PIDTree_ann->Draw("annNueCCQEvsNC>>hAll_ElNc", nuAllString_ElNc.Data());
  TString nuNCString_ElMu = Form("(%f*(trueBeamPDG == 14) + %f*(trueBeamPDG == 12)) * (preselected && trueNCEvent)", weightNumuNCEvents, weightNueNCEvents);
  TString nuNCString_ElNc = Form("(%f*(trueBeamPDG == 14) + %f*(trueBeamPDG == 12)) * (preselected && trueNCEvent)", weightNumuNCEvents, weightNueNCEvents);
  PIDTree_ann->Draw("annNueCCQEvsNumuCCQE>>hNC_ElMu", nuNCString_ElMu.Data());
  PIDTree_ann->Draw("annNueCCQEvsNC>>hNC_ElNc", nuNCString_ElNc.Data());
  TString numuCCString_ElMu  = Form("(trueBeamPDG==14 * trueCCEvent) * (%f*trueQEEvent + %f*(!trueQEEvent)) * (preselected)", weightNumuCCQEEvents, weightNumuCCnonQEEvents);
  TString numuCCString_ElNc  = Form("(trueBeamPDG==14 * trueCCEvent) * (%f*trueQEEvent + %f*(!trueQEEvent)) * (preselected)", weightNumuCCQEEvents, weightNumuCCnonQEEvents);
  PIDTree_ann->Draw("annNueCCQEvsNumuCCQE>>hNumuCC_ElMu", numuCCString_ElMu.Data());
  PIDTree_ann->Draw("annNueCCQEvsNC>>hNumuCC_ElNc", numuCCString_ElNc.Data());
  TString nueCCString_ElMu   = Form("(trueBeamPDG==12 * trueCCEvent) * (%f*trueQEEvent + %f*(!trueQEEvent)) * (preselected)", weightNueCCQEEvents, weightNueCCnonQEEvents);
  TString nueCCString_ElNc   = Form("(trueBeamPDG==12 * trueCCEvent) * (%f*trueQEEvent + %f*(!trueQEEvent)) * (preselected)", weightNueCCQEEvents, weightNueCCnonQEEvents);
  PIDTree_ann->Draw("annNueCCQEvsNumuCCQE>>hNueCC_ElMu", nueCCString_ElMu.Data());
  PIDTree_ann->Draw("annNueCCQEvsNC>>hNueCC_ElNc", nueCCString_ElNc.Data());
  TString nueCCQEString_ElMu = Form("(trueBeamPDG==12 * trueCCEvent) * (%f*trueQEEvent) * (preselected)", weightNueCCQEEvents);
  TString nueCCQEString_ElNc = Form("(trueBeamPDG==12 * trueCCEvent) * (%f*trueQEEvent) * (preselected)", weightNueCCQEEvents);
  PIDTree_ann->Draw("annNueCCQEvsNumuCCQE>>hNueCCQE_ElMu", nueCCQEString_ElMu.Data());
  PIDTree_ann->Draw("annNueCCQEvsNC>>hNueCCQE_ElNc", nueCCQEString_ElNc.Data());

  // Define selection and plot from PIDTree the All ANN variables plots with the ANN cuts and preselection...
  TString nuAllString_ElMu_cuts = Form("((trueBeamPDG==14)*(trueCCEvent*(%f*trueQEEvent + %f*!trueQEEvent) + %f*trueNCEvent) + (trueBeamPDG == 12)*(trueCCEvent * (%f*trueQEEvent + %f * !trueQEEvent) + %f*trueNCEvent)) * (preselected && (annNueCCQEvsNC > %f))", weightNumuCCQEEvents, weightNumuCCnonQEEvents, weightNumuNCEvents, weightNueCCQEEvents, weightNueCCnonQEEvents, weightNueNCEvents, ElNcCut);
  TString nuAllString_ElNc_cuts = Form("((trueBeamPDG==14)*(trueCCEvent*(%f*trueQEEvent + %f*!trueQEEvent) + %f*trueNCEvent) + (trueBeamPDG == 12)*(trueCCEvent * (%f*trueQEEvent + %f * !trueQEEvent) + %f*trueNCEvent)) * (preselected && (annNueCCQEvsNumuCCQE > %f))", weightNumuCCQEEvents, weightNumuCCnonQEEvents, weightNumuNCEvents, weightNueCCQEEvents, weightNueCCnonQEEvents, weightNueNCEvents, ElMuCut);
  PIDTree_ann->Draw("annNueCCQEvsNumuCCQE>>hAll_ElMu_cuts", nuAllString_ElMu_cuts.Data());
  PIDTree_ann->Draw("annNueCCQEvsNC>>hAll_ElNc_cuts", nuAllString_ElNc_cuts.Data());
  TString nuNCString_ElMu_cuts = Form("(%f*(trueBeamPDG == 14) + %f*(trueBeamPDG == 12)) * (preselected && trueNCEvent && (annNueCCQEvsNC > %f))", weightNumuNCEvents, weightNueNCEvents, ElNcCut);
  TString nuNCString_ElNc_cuts = Form("(%f*(trueBeamPDG == 14) + %f*(trueBeamPDG == 12)) * (preselected && trueNCEvent && (annNueCCQEvsNumuCCQE > %f))", weightNumuNCEvents, weightNueNCEvents, ElMuCut);
  PIDTree_ann->Draw("annNueCCQEvsNumuCCQE>>hNC_ElMu_cuts", nuNCString_ElMu_cuts.Data());
  PIDTree_ann->Draw("annNueCCQEvsNC>>hNC_ElNc_cuts", nuNCString_ElNc_cuts.Data());
  TString numuCCString_ElMu_cuts  = Form("(trueBeamPDG==14 * trueCCEvent) * (%f*trueQEEvent + %f*(!trueQEEvent)) * (preselected && (annNueCCQEvsNC > %f))", weightNumuCCQEEvents, weightNumuCCnonQEEvents, ElNcCut);
  TString numuCCString_ElNc_cuts  = Form("(trueBeamPDG==14 * trueCCEvent) * (%f*trueQEEvent + %f*(!trueQEEvent)) * (preselected && (annNueCCQEvsNumuCCQE > %f))", weightNumuCCQEEvents, weightNumuCCnonQEEvents, ElMuCut);
  PIDTree_ann->Draw("annNueCCQEvsNumuCCQE>>hNumuCC_ElMu_cuts", numuCCString_ElMu_cuts.Data());
  PIDTree_ann->Draw("annNueCCQEvsNC>>hNumuCC_ElNc_cuts", numuCCString_ElNc_cuts.Data());
  TString nueCCString_ElMu_cuts   = Form("(trueBeamPDG==12 * trueCCEvent) * (%f*trueQEEvent + %f*(!trueQEEvent)) * (preselected && (annNueCCQEvsNC > %f))", weightNueCCQEEvents, weightNueCCnonQEEvents, ElNcCut);
  TString nueCCString_ElNc_cuts   = Form("(trueBeamPDG==12 * trueCCEvent) * (%f*trueQEEvent + %f*(!trueQEEvent)) * (preselected && (annNueCCQEvsNumuCCQE > %f))", weightNueCCQEEvents, weightNueCCnonQEEvents, ElMuCut);
  PIDTree_ann->Draw("annNueCCQEvsNumuCCQE>>hNueCC_ElMu_cuts", nueCCString_ElMu_cuts.Data());
  PIDTree_ann->Draw("annNueCCQEvsNC>>hNueCC_ElNc_cuts", nueCCString_ElNc_cuts.Data());
  TString nueCCQEString_ElMu_cuts = Form("(trueBeamPDG==12 * trueCCEvent) * (%f*trueQEEvent) * (preselected && (annNueCCQEvsNC > %f))", weightNueCCQEEvents, ElNcCut);
  TString nueCCQEString_ElNc_cuts = Form("(trueBeamPDG==12 * trueCCEvent) * (%f*trueQEEvent) * (preselected && (annNueCCQEvsNumuCCQE > %f))", weightNueCCQEEvents, ElMuCut);
  PIDTree_ann->Draw("annNueCCQEvsNumuCCQE>>hNueCCQE_ElMu_cuts", nueCCQEString_ElMu_cuts.Data());
  PIDTree_ann->Draw("annNueCCQEvsNC>>hNueCCQE_ElNc_cuts", nueCCQEString_ElNc_cuts.Data());

  double maximum;
  double scaleAll;

  // Draw The All ANN ElMu Variable Plot without cuts but with preselection...
  TCanvas *c1 = new TCanvas("c1", "", 800, 600);
  c1->cd();
  maximum = hAll_ElMu->GetMaximum();
  scaleAll = 1.0/maximum;

  hAll_ElMu->Draw();
  hNueCC_ElMu->Draw("SAME");
  hNueCCQE_ElMu->Draw("SAME");
  hNC_ElMu->Draw("SAME");
  hNumuCC_ElMu->Draw("SAME");

  hAll_ElMu->Scale(scaleAll);
  hNC_ElMu->Scale(scaleAll);
  hNumuCC_ElMu->Scale(scaleAll);
  hNueCC_ElMu->Scale(scaleAll);
  hNueCCQE_ElMu->Scale(scaleAll);
  hAll_ElMu->SetMaximum(1.19);

  TPaveText * pt = new TPaveText(0.18, 0.89, 0.845, 0.84, "blNDC");
  pt->AddText("Preselection applied to initial sample of 95% #nu_{#mu}, 5% #nu_{e}");
  pt->Draw("SAME");
  pt->SetBorderSize(0.0);

  TLegend * leg = new TLegend(0.35, 0.83, 0.65, 0.58);
  leg->AddEntry(hAll_ElMu, "All events", "F");
  leg->AddEntry(hNC_ElMu, "All NC events", "L");
  leg->AddEntry(hNumuCC_ElMu, "#nu_{#mu} CC events", "L");
  leg->AddEntry(hNueCC_ElMu, "#nu_{e} CC events", "F");
  leg->AddEntry(hNueCCQE_ElMu, "#nu_{e} CCQE events only", "F");
  leg->SetTextFont(42);
  leg->Draw("SAME");

  c1->Update();
  c1->SaveAs("annNueCCQEvsNumuCCQE.png");
  //c1->SaveAs("annNueCCQEvsNumuCCQE.root");
  //c1->SaveAs("annNueCCQEvsNumuCCQE.pdf");
  c1->SaveAs("annNueCCQEvsNumuCCQE.C");

  // Draw The All ANN ElNc Variable Plot without cuts but with preselection...
  TCanvas *c2 = new TCanvas("c2", "", 800, 600);
  c2->cd();
  maximum = hAll_ElNc->GetMaximum();
  scaleAll = 1.0/maximum;

  hAll_ElNc->Draw();
  hNueCC_ElNc->Draw("SAME");
  hNueCCQE_ElNc->Draw("SAME");
  hNC_ElNc->Draw("SAME");
  hNumuCC_ElNc->Draw("SAME");

  hAll_ElNc->Scale(scaleAll);
  hNC_ElNc->Scale(scaleAll);
  hNumuCC_ElNc->Scale(scaleAll);
  hNueCC_ElNc->Scale(scaleAll);
  hNueCCQE_ElNc->Scale(scaleAll);
  hAll_ElNc->SetMaximum(1.19);

  TPaveText * pt = new TPaveText(0.18, 0.89, 0.845, 0.84, "blNDC");
  pt->AddText("Preselection applied to initial sample of 95% #nu_{#mu}, 5% #nu_{e}");
  pt->Draw("SAME");
  pt->SetBorderSize(0.0);

  TLegend * leg = new TLegend(0.2, 0.83, 0.43, 0.58);
  leg->AddEntry(hAll_ElNc, "All events", "F");
  leg->AddEntry(hNC_ElNc, "All NC events", "L");
  leg->AddEntry(hNumuCC_ElNc, "#nu_{#mu} CC events", "L");
  leg->AddEntry(hNueCC_ElNc, "#nu_{e} CC events", "F");
  leg->AddEntry(hNueCCQE_ElNc, "#nu_{e} CCQE events only", "F");
  leg->SetTextFont(42);
  leg->Draw("SAME");

  c2->Update();
  c2->SaveAs("annNueCCQEvsNC.png");
  //c2->SaveAs("annNueCCQEvsNC.root");
  //c2->SaveAs("annNueCCQEvsNC.pdf");
  c2->SaveAs("annNueCCQEvsNC.C");

  // Draw The All ANN ElMu Variable Plot With cuts and preselection...
  TCanvas *c3 = new TCanvas("c3", "", 800, 600);
  c3->cd();
  maximum = hAll_ElMu_cuts->GetMaximum();
  scaleAll = 1.0/maximum;

  hAll_ElMu_cuts->Draw();
  hNueCC_ElMu_cuts->Draw("SAME");
  hNueCCQE_ElMu_cuts->Draw("SAME");
  hNC_ElMu_cuts->Draw("SAME");
  hNumuCC_ElMu_cuts->Draw("SAME");

  hAll_ElMu_cuts->Scale(scaleAll);
  hNC_ElMu_cuts->Scale(scaleAll);
  hNumuCC_ElMu_cuts->Scale(scaleAll);
  hNueCC_ElMu_cuts->Scale(scaleAll);
  hNueCCQE_ElMu_cuts->Scale(scaleAll);
  hAll_ElMu_cuts->SetMaximum(1.19);

  TLine * line = new TLine(ElMuCut, 0.0, ElMuCut, 0.9);
  line->SetLineWidth(2);
  line->SetLineColor(kBlack);
  line->Draw("SAME");

  //TArrow * arr = new TArrow(0.9, 0.9, 1.05, 0.9, 0.05, "|>");
  //arr->SetLineWidth(2);
  //arr->SetLineColor(kBlack);
  //arr->Draw("SAME");

  pt1->Draw("SAME");
  pt1->SetBorderSize(0.0);

  TLegend * leg = new TLegend(0.35, 0.83, 0.65, 0.58);
  leg->AddEntry(hAll_ElMu_cuts, "All events", "F");
  leg->AddEntry(hNC_ElMu_cuts, "All NC events", "L");
  leg->AddEntry(hNumuCC_ElMu_cuts, "#nu_{#mu} CC events", "L");
  leg->AddEntry(hNueCC_ElMu_cuts, "#nu_{e} CC events", "F");
  leg->AddEntry(hNueCCQE_ElMu_cuts, "#nu_{e} CCQE events only", "F");
  leg->SetTextFont(42);
  leg->Draw("SAME");

  c3->Update();
  c3->SaveAs("annNueCCQEvsNumuCCQE_cuts.png");
  //c3->SaveAs("annNueCCQEvsNumuCCQE_cuts.root");
  //c3->SaveAs("annNueCCQEvsNumuCCQE_cuts.pdf");
  c3->SaveAs("annNueCCQEvsNumuCCQE_cuts.C");

  // Draw The All ANN ElNc Variable Plot With cuts and preselection...
  TCanvas *c4 = new TCanvas("c4", "", 800, 600);
  c4->cd();
  maximum = hAll_ElNc_cuts->GetMaximum();
  scaleAll = 1.0/maximum;

  hAll_ElNc_cuts->Draw();
  hNueCC_ElNc_cuts->Draw("SAME");
  hNueCCQE_ElNc_cuts->Draw("SAME");
  hNC_ElNc_cuts->Draw("SAME");
  hNumuCC_ElNc_cuts->Draw("SAME");

  hAll_ElNc_cuts->Scale(scaleAll);
  hNC_ElNc_cuts->Scale(scaleAll);
  hNumuCC_ElNc_cuts->Scale(scaleAll);
  hNueCC_ElNc_cuts->Scale(scaleAll);
  hNueCCQE_ElNc_cuts->Scale(scaleAll);
  hAll_ElNc_cuts->SetMaximum(1.19);

  TLine * line = new TLine(ElNcCut, 0.0, ElNcCut, 0.9);
  line->SetLineWidth(2);
  line->SetLineColor(kBlack);
  line->Draw("SAME");

  //TArrow * arr = new TArrow(0.9, 0.9, 1.05, 0.9, 0.05, "|>");
  //arr->SetLineWidth(2);
  //arr->SetLineColor(kBlack);
  //arr->Draw("SAME");

  pt2->Draw("SAME");
  pt2->SetBorderSize(0.0);

  TLegend * leg = new TLegend(0.2, 0.83, 0.43, 0.58);
  leg->AddEntry(hAll_ElMu_cuts, "All events", "F");
  leg->AddEntry(hNC_ElMu_cuts, "All NC events", "L");
  leg->AddEntry(hNumuCC_ElMu_cuts, "#nu_{#mu} CC events", "L");
  leg->AddEntry(hNueCC_ElMu_cuts, "#nu_{e} CC events", "F");
  leg->AddEntry(hNueCCQE_ElMu_cuts, "#nu_{e} CCQE events only", "F");
  leg->SetTextFont(42);
  leg->Draw("SAME");

  c4->Update();
  c4->SaveAs("annNueCCQEvsNC_cuts.png");
  //c4->SaveAs("annNueCCQEvsNC_cuts.root");
  //c4->SaveAs("annNueCCQEvsNC_cuts.pdf");
  c4->SaveAs("annNueCCQEvsNC_cuts.C");


  // Now I want to make the Efficiency and Purity Plot
  // First find all the efficiencies we want...

  TH1F* nueCCQEAll = new TH1F("nueCCQEAll", "", 10, 50, 5050);
  TH1F* nueCCQESel = new TH1F("nueCCQESel", "", 10, 50, 5050);

  TH1F* nueCCAll = new TH1F("nueCCAll", "", 10, 50, 5050);
  TH1F* nueCCSel = new TH1F("nueCCSel", "", 10, 50, 5050);

  TH1F* numuCCAll = new TH1F("numuCCAll", "", 10, 50, 5050);
  TH1F* numuCCSel = new TH1F("numuCCSel", "", 10, 50, 5050);

  TH1F* NCAll = new TH1F("NCAll", "", 10, 50, 5050);
  TH1F* NCSel = new TH1F("NCSel", "", 10, 50, 5050);

  TString nueCCQEAll_string = Form("1 && trueBeamPDG==12 && trueCCEvent && trueQEEvent");
  TString nueCCQESel_string = Form("1 && trueBeamPDG==12 && trueCCEvent && trueQEEvent && preselected && (annNueCCQEvsNumuCCQE > %f) && (annNueCCQEvsNC > %f)", ElMuCut, ElNcCut);
  PIDTree_ann->Draw("recoE_el>>nueCCQEAll", nueCCQEAll_string.Data());
  PIDTree_ann->Draw("recoE_el>>nueCCQESel", nueCCQESel_string.Data());
  TEfficiency* nueCCQEEff = new TEfficiency(*nueCCQESel,*nueCCQEAll);

  TString nueCCAll_string = Form("1 && trueBeamPDG==12 && trueCCEvent");
  TString nueCCSel_string = Form("1 && trueBeamPDG==12 && trueCCEvent && preselected && (annNueCCQEvsNumuCCQE > %f) && (annNueCCQEvsNC > %f)", ElMuCut, ElNcCut);
  PIDTree_ann->Draw("recoE_el>>nueCCAll", nueCCAll_string.Data());
  PIDTree_ann->Draw("recoE_el>>nueCCSel", nueCCSel_string.Data());
  TEfficiency* nueCCEff = new TEfficiency(*nueCCSel,*nueCCAll);

  TString numuCCAll_string = Form("1 && trueBeamPDG==14 && trueCCEvent");
  TString numuCCSel_string = Form("1 && trueNCEvent && preselected && (annNueCCQEvsNumuCCQE > %f) && (annNueCCQEvsNC > %f)", ElMuCut, ElNcCut);
  PIDTree_ann->Draw("recoE_el>>numuCCAll", numuCCAll_string.Data());
  PIDTree_ann->Draw("recoE_el>>numuCCSel", numuCCSel_string.Data());
  TEfficiency* numuCCEff = new TEfficiency(*numuCCSel,*numuCCAll);

  TString NCAll_string = Form("1 && trueNCEvent");
  TString NCSel_string = Form("1 && trueNCEvent && preselected && (annNueCCQEvsNumuCCQE > %f) && (annNueCCQEvsNC > %f)", ElMuCut, ElNcCut);
  PIDTree_ann->Draw("recoE_el>>NCAll", NCAll_string.Data());
  PIDTree_ann->Draw("recoE_el>>NCSel", NCSel_string.Data());
  TEfficiency* NCEff = new TEfficiency(*NCSel,*NCAll);

  // Now find the signal purity
  TH1F* nueCCSignal = new TH1F("nueCCSignal", "", 10, 50, 5050);
  TH1F* nueCCTotal = new TH1F("nueCCTotal", "", 10, 50, 5050);

  TString signal_string = Form("((trueBeamPDG==12)*(trueCCEvent * (%f*trueQEEvent + %f * !trueQEEvent))) * (preselected && (annNueCCQEvsNumuCCQE > %f) && (annNueCCQEvsNC > %f))", weightNueCCQEEvents, weightNueCCnonQEEvents, ElMuCut, ElNcCut);
  TString total_string  = Form("((trueBeamPDG==14)*(trueCCEvent*(%f*trueQEEvent + %f*!trueQEEvent) + %f*trueNCEvent) + (trueBeamPDG == 12)*(trueCCEvent * (%f*trueQEEvent + %f * !trueQEEvent) + %f*trueNCEvent)) * (preselected && (annNueCCQEvsNumuCCQE > %f) && (annNueCCQEvsNC > %f))", weightNumuCCQEEvents, weightNumuCCnonQEEvents, weightNumuNCEvents, weightNueCCQEEvents, weightNueCCnonQEEvents, weightNueNCEvents, ElMuCut, ElNcCut);
  PIDTree_ann->Draw("recoE_el>>nueCCSignal", signal_string.Data());
  PIDTree_ann->Draw("recoE_el>>nueCCTotal", total_string.Data());

  nueCCSignal->Divide(nueCCTotal);

  TCanvas *c5 = new TCanvas("c5", "", 800, 600);
  c5->cd();

  TH2F *hempty = new TH2F("hempty", ";recoE [MeV]; Efficiency or Purity", 10, 50, 5050, 10,  0, 1 );
  hempty->GetXaxis()->SetTitleSize(0.06);    hempty->GetYaxis()->SetTitleSize(0.06);
  hempty->GetXaxis()->SetTitleOffset(0.8);   hempty->GetYaxis()->SetTitleOffset(0.8);
  hempty->GetXaxis()->SetLabelSize(0.05);    hempty->GetYaxis()->SetLabelSize(0.05);
  hempty->Draw();

  nueCCQEEff->SetLineColor(kBlack);       nueCCQEEff->SetLineWidth(2);   nueCCQEEff->SetMarkerSize(1.2);    nueCCQEEff->SetMarkerStyle(20);    nueCCQEEff->SetMarkerColor(kBlack);    nueCCQEEff->Draw("sameP");
  nueCCEff->SetLineColor(kBlue);          nueCCEff->SetLineWidth(2);     nueCCEff->SetMarkerSize(1.2);      nueCCEff->SetMarkerStyle(20);      nueCCEff->SetMarkerColor(kBlue);       nueCCEff->Draw("sameP");
  numuCCEff->SetLineColor(kGreen+2);      numuCCEff->SetLineWidth(2);    numuCCEff->SetMarkerSize(1.2);     numuCCEff->SetMarkerStyle(20);     numuCCEff->SetMarkerColor(kGreen+2);   numuCCEff->Draw("sameP");
  NCEff->SetLineColor(kMagenta);          NCEff->SetLineWidth(2);        NCEff->SetMarkerSize(1.2);         NCEff->SetMarkerStyle(20);         NCEff->SetMarkerColor(kMagenta);       NCEff->Draw("sameP");
  nueCCSignal->SetLineColor(kRed);        nueCCSignal->SetLineWidth(0);  nueCCSignal->SetMarkerSize(1.2);   nueCCSignal->SetMarkerStyle(20);   nueCCSignal->SetMarkerColor(kRed);     nueCCSignal->Draw("sameP");

  TLegend *leg = new TLegend(0.12, 0.70, 0.32, 0.9, "Efficiency");
  leg->AddEntry(nueCCEff, "#nu_{e} CC all", "P");
  leg->AddEntry(nueCCQEEff, "#nu_{e} CC QE", "P");
  leg->AddEntry(numuCCEff, "#nu_{#mu} CC", "P");
  leg->AddEntry(NCEff, "NC", "P");
  leg->SetTextSize(0.03);
  leg->SetTextFont(42);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->Draw();

  TLegend *leg2 = new TLegend(0.12, 0.62, 0.32, 0.7, "Purity");
  leg2->AddEntry(nueCCSignal, "#nu_{e} CC", "P");
  leg2->SetTextFont(42);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(0);
  leg2->SetTextSize(0.03);
  leg2->Draw();

  c5->Update();
  c5->SaveAs("effPur.png");
  //c5->SaveAs("effPur.root");
  //c5->SaveAs("effPur.pdf");
  c5->SaveAs("effPur.C");

  /*
  TFile * mainOutput = new TFile("test.root","RECREATE");

  hAllEvents_ElMu_cuts->Write();
  hAllEvents_ElNc_cuts->Write();
  hNCEvents_ElMu_cuts->Write();
  hNCEvents_ElNc_cuts->Write();
  hMuCCEvents_ElMu_cuts->Write();
  hMuCCEvents_ElNc_cuts->Write();
  hElCCEvents_ElMu_cuts->Write();
  hElCCEvents_ElNc_cuts->Write();
  hElCCQEEvents_ElMu_cuts->Write();
  hElCCQEEvents_ElNc_cuts->Write();

  mainOutput->Close();
  */


} // eof
