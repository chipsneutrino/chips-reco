void wc_trackfitter_twoTracks(const char * infile = "", int start=0, int fit=100){
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
  WCSimFitterInterface::Instance()->SetNumTracks(2);
  WCSimFitterInterface::Instance()->SetTrackType(0, "MuonLike");
  WCSimFitterInterface::Instance()->SetTrackType(1, "MuonLike");

  // Set parameter(track number, "name", minimum, maximum, start, is fixed?)
  // Names are: kVtxX, kVtxY, kVtxZ, kDirTh, kDirPhi, kEnergy
  WCSimFitterInterface::Instance()->SetParameter(0, "kVtxX", -1200, 1200, 226.27, true);
  WCSimFitterInterface::Instance()->SetParameter(0, "kVtxY", -1200, 1200, -575.35, true);
  WCSimFitterInterface::Instance()->SetParameter(0, "kVtxZ", -900, 900, 762.81, true);
  WCSimFitterInterface::Instance()->SetParameter(0, "kVtxT", -900, 900, 0, true);
  WCSimFitterInterface::Instance()->SetParameter(0, "kDirTh", 0, TMath::Pi(), 0.0, false);
  WCSimFitterInterface::Instance()->SetParameter(0, "kDirPhi", -TMath::Pi(), TMath::Pi(), 0.0, false);
  WCSimFitterInterface::Instance()->SetParameter(0, "kEnergy", 900, 3100, 1500, false);

  WCSimFitterInterface::Instance()->SetParameter(1, "kVtxX", -1200, 1200, 226.27, true);
  WCSimFitterInterface::Instance()->SetParameter(1, "kVtxY", -1200, 1200, -575.35, true);
  WCSimFitterInterface::Instance()->SetParameter(1, "kVtxZ", -900, 900, 762.81, true);
  WCSimFitterInterface::Instance()->SetParameter(1, "kVtxT", -900, 900, 0, true);
  WCSimFitterInterface::Instance()->SetParameter(1, "kDirTh", 0, TMath::Pi(), 0.0, false);
  WCSimFitterInterface::Instance()->SetParameter(1, "kDirPhi", -TMath::Pi(), TMath::Pi(), 0.0, false);
  WCSimFitterInterface::Instance()->SetParameter(1, "kEnergy", 900, 3100, 850, false);

  WCSimFitterInterface::Instance()->JoinParametersTogether(0,1,"kVtxX");
  WCSimFitterInterface::Instance()->JoinParametersTogether(0,1,"kVtxY");
  WCSimFitterInterface::Instance()->JoinParametersTogether(0,1,"kVtxZ");
  WCSimFitterInterface::Instance()->JoinParametersTogether(0,1,"kVtxT");

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
  WCSimFitterInterface::Instance()->SetNumSurfaceBins(22);
  WCSimFitterInterface::Instance()->Make1DSurface(0, "kEnergy");
  //WCSimFitterInterface::Instance()->Make2DSurface(0, "kVtxX", 0, "kVtxY");
  WCSimFitterInterface::Instance()->SetMakeFits(kTRUE);
  WCSimFitterInterface::Instance()->SetMakeSurfaces(kFALSE);



  WCSimFitterInterface::Instance()->Print();
  WCSimFitterInterface::Instance()->SetNumEventsToFit(fit);
  WCSimFitterInterface::Instance()->SetFirstEventToFit(start);
  WCSimFitterInterface::Instance()->Run();
  
  std::cout << "Done!" << std::endl;

}
