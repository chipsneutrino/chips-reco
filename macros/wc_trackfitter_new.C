void wc_trackfitter_new(const char * infile = "", int start=0, int fit=100){
  // Path to WCSim ROOT file
  // =======================
  TString filename(infile);
  if(filename.CompareTo("") == 0 )
  {
    filename = TString("localfile.root");
  }
  gApplication->ProcessLine(".except");

  // Load libraries
  // ==============
  gSystem->Load("libGeom");
  gSystem->Load("libEve");
  gSystem->Load("libMinuit");

  TString libWCSimRoot = TString::Format("%s%s",gSystem->Getenv("WCSIMHOME"), "/libWCSimRoot.so");
  TString libWCSimAnalysis = TString::Format("%s%s",gSystem->Getenv("WCSIMANAHOME"), "/lib/libWCSimAnalysis.so");
  gSystem->Load(libWCSimRoot.Data());
  gSystem->Load(libWCSimAnalysis.Data());


  // Load Data
  // =========
  WCSimInterface::LoadData(filename.Data());
  
  WCSimFitterInterface::Instance()->SetInputFileName(filename.Data()); // For inclusion in the fitterPlots file
  WCSimFitterInterface::Instance()->SetNumTracks(1);
  WCSimFitterInterface::Instance()->SetTrackType(0, "MuonLike");

  // Set parameter(track number, "name", minimum, maximum, start, is fixed?)
  // Names are: kVtxX, kVtxY, kVtxZ, kDirTh, kDirPhi, kEnergy
  WCSimFitterInterface::Instance()->SetParameter(0, "kVtxX", -1200, 1200, 0, false);
  WCSimFitterInterface::Instance()->SetParameter(0, "kVtxY", -1200, 1200, 0, false);
  WCSimFitterInterface::Instance()->SetParameter(0, "kVtxZ", -900, 900, 0, false);
  WCSimFitterInterface::Instance()->SetParameter(0, "kVtxT", -900, 900, 0, true);
  WCSimFitterInterface::Instance()->SetParameter(0, "kDirTh", 0, TMath::Pi(), 0.5*TMath::Pi(), false);
  WCSimFitterInterface::Instance()->SetParameter(0, "kDirPhi", -TMath::Pi(), TMath::Pi(), 0.0*TMath::Pi(), false);
  WCSimFitterInterface::Instance()->SetParameter(0, "kEnergy", 500, 2400, 1500, false);


  // Plot best-fit results
  WCSimFitterInterface::Instance()->PlotForEachEvent("kVtxX");
  WCSimFitterInterface::Instance()->PlotForEachEvent("kVtxY");
  WCSimFitterInterface::Instance()->PlotForEachEvent("kVtxZ");
  WCSimFitterInterface::Instance()->PlotForEachEvent("kVtxT");
  WCSimFitterInterface::Instance()->PlotForEachEvent("kDirTh");
  WCSimFitterInterface::Instance()->PlotForEachEvent("kDirPhi");
  WCSimFitterInterface::Instance()->PlotForEachEvent("kEnergy");

  // Plot reco - true values
  WCSimFitterInterface::Instance()->PlotRecoMinusTrue("kVtxX");
  WCSimFitterInterface::Instance()->PlotRecoMinusTrue("kVtxY");
  WCSimFitterInterface::Instance()->PlotRecoMinusTrue("kVtxZ");
  WCSimFitterInterface::Instance()->PlotRecoMinusTrue("kVtxT");
  WCSimFitterInterface::Instance()->PlotRecoMinusTrue("kDirTh");
  WCSimFitterInterface::Instance()->PlotRecoMinusTrue("kDirPhi");
  WCSimFitterInterface::Instance()->PlotRecoMinusTrue("kEnergy");

  // Sweep out one variable
  WCSimFitterInterface::Instance()->SetNumSurfaceBins(50);
  WCSimFitterInterface::Instance()->Make1DSurface(0, "kEnergy");
  WCSimFitterInterface::Instance()->SetMakeFits(kTRUE);
  WCSimFitterInterface::Instance()->SetMakeSurfaces(kFALSE);



  WCSimFitterInterface::Instance()->Print();
  WCSimFitterInterface::Instance()->SetNumEventsToFit(fit);
  WCSimFitterInterface::Instance()->SetFirstEventToFit(start);
  WCSimFitterInterface::Instance()->Run();
  
  std::cout << "Done!" << std::endl;

}
