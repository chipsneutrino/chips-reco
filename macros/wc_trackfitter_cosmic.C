void wc_trackfitter_cosmic(const char * infile = "", int start=0, int fit=100){
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

  WCSimFitterInterface myInterface;  
  myInterface.SetInputFileName(filename.Data()); // For inclusion in the fitterPlots file
  myInterface.SetNumTracks(2);
  myInterface.SetTrackType(0, "ElectronLike");
  myInterface.SetTrackType(1, "MuonLike");

  // Set parameter(track number, "name", minimum, maximum, start, is fixed?)
  // Names are: kVtxX, kVtxY, kVtxZ, kDirTh, kDirPhi, kEnergy
  
  // Convenience settings
  bool fFixVtx = false;
  bool fFixDir = false;
  double fVtxX = 0.0;
  double fVtxY = 0.0;
  double fVtxZ = 0.0;
  double fDirX1 = 0.0;
  double fDirY1 = 0.0;
  double fDirZ1 = 0.0;
  double fDirX2 = 0.0;
  double fDirY2 = 0.0;
  double fDirZ2 = 0.0;
  double fTheta1 = TMath::ACos(fDirZ1);
  double fPhi1 = TMath::ATan2(fDirY1,fDirX1);
  double fTheta2 = TMath::ACos(fDirZ2);
  double fPhi2 = TMath::ATan2(fDirY2,fDirX2);
  //
  myInterface.SetParameter(0, "kVtxX", -1267, 1267, fVtxX, fFixVtx, 10.0);
  myInterface.SetParameter(0, "kVtxY", -1267, 1267, fVtxY, fFixVtx, 10.0);
  myInterface.SetParameter(0, "kVtxZ", -1020, 1020, fVtxZ, fFixVtx, 10.0);
  myInterface.SetParameter(0, "kVtxT", 0, 10000, 0, false, 1.0);
  myInterface.SetParameter(0, "kDirTh", 0, TMath::Pi(), fTheta1, fFixDir, 0.01);
  myInterface.SetParameter(0, "kDirPhi", -TMath::Pi(), TMath::Pi(), fPhi1, fFixDir, 0.02);
  myInterface.SetParameter(0, "kEnergy", 500, 5000, 1500, false, 250.0);

  myInterface.SetParameter(1, "kVtxX", -1267, 1267, fVtxX, fFixVtx, 10.0);
  myInterface.SetParameter(1, "kVtxY", -1267, 1267, fVtxY, fFixVtx, 10.0);
  myInterface.SetParameter(1, "kVtxZ", -1020, 1020, fVtxZ, fFixVtx, 10.0);
  myInterface.SetParameter(1, "kVtxT", 0, 10000, 0, false, 1.0);
  myInterface.SetParameter(1, "kDirTh", 0, TMath::Pi(), fTheta2, fFixDir, 0.01);
  myInterface.SetParameter(1, "kDirPhi", -TMath::Pi(), TMath::Pi(), fPhi2, fFixDir, 0.02);
  myInterface.SetParameter(1, "kEnergy", 500, 5000, 850, false, 250.0);

  // Standard clustering
  WCSimParameters::Instance()->SetSlicerClusterDistance(250);
  WCSimParameters::Instance()->SetSlicerMinSize(25);
  WCSimParameters::Instance()->SetSlicerChargeCut(0.0);
  WCSimParameters::Instance()->SetSlicerTimeCut(100);
  WCSimParameters::Instance()->SetIterateSlicing(true);

  // Veto clustering
  WCSimParameters::Instance()->SetVetoClusterDistance(500);
  WCSimParameters::Instance()->SetVetoMinSize(5);
  WCSimParameters::Instance()->SetVetoMinChargeCut(2.0);
  WCSimParameters::Instance()->SetVetoMaxChargeCut(30.0);
  WCSimParameters::Instance()->SetVetoTimeCut(200);

  // Make sure we know this is a cosmic event
  myInterface.SetIsCosmicFit(true);

  // Plot best-fit results
  myInterface.PlotForEachEvent("kVtxX");
  myInterface.PlotForEachEvent("kVtxY");
  myInterface.PlotForEachEvent("kVtxZ");
  myInterface.PlotForEachEvent("kVtxT");
  myInterface.PlotForEachEvent("kDirTh");
  myInterface.PlotForEachEvent("kDirPhi");
  myInterface.PlotForEachEvent("kEnergy");

  // Plot reco - true values
  myInterface.PlotRecoMinusTrue("kVtxX");
  myInterface.PlotRecoMinusTrue("kVtxY");
  myInterface.PlotRecoMinusTrue("kVtxZ");
  myInterface.PlotRecoMinusTrue("kVtxT");
  myInterface.PlotRecoMinusTrue("kDirTh");
  myInterface.PlotRecoMinusTrue("kDirPhi");
  myInterface.PlotRecoMinusTrue("kEnergy");

  // Sweep out one variable
  myInterface.SetNumSurfaceBins(22);
  myInterface.Make1DSurface(0, "kEnergy");
  //WCSimFitterInterface::Instance()->Make2DSurface(0, "kVtxX", 0, "kVtxY");
  myInterface.SetMakeFits(kTRUE);
  myInterface.SetMakeSurfaces(kFALSE);

  myInterface.Print();
  myInterface.SetNumEventsToFit(fit);
  myInterface.SetFirstEventToFit(start);
  myInterface.Run();
  
  std::cout << "Done!" << std::endl;

}
