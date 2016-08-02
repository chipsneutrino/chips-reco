void wc_trackfitter_piZero(const char * infile = "", int start=0, int fit=100){
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
  myInterface.SetTrackType(0, "PhotonLike");
  myInterface.SetTrackType(1, "PhotonLike");
  myInterface.SetIsPiZeroFit(true);
  myInterface.SetForcePiZeroMass(false);

//  myInterface.SetTrackType(0, "MuonLike");
//  myInterface.SetTrackType(1, "MuonLike");

  // Set parameter(track number, "name", minimum, maximum, start, is fixed?)
  // Names are: kVtxX, kVtxY, kVtxZ, kDirTh, kDirPhi, kEnergy
  
  // Convenience settings
  bool fFixVtx = false;
  bool fFixDir = false;
  double fVtxX = 0.0;
  double fVtxY = 0.0;
  double fVtxZ = 0.0;
  TVector3 dir1(1,1,0);
  TVector3 dir2(1,-1,0);
  double fTheta1 = dir1.Theta();
  double fPhi1 = dir1.Phi();
  double fTheta2 = dir2.Theta();
  double fPhi2 = dir2.Phi();


  myInterface.SetParameter(0, "kVtxX", -1200, 1200, fVtxX, fFixVtx, 10.0);
  myInterface.SetParameter(0, "kVtxY", -1200, 1200, fVtxY, fFixVtx, 10.0);
  myInterface.SetParameter(0, "kVtxZ", -900, 900, fVtxZ, fFixVtx, 10.0);
  myInterface.SetParameter(0, "kVtxT", -100, 101000, 5000, false, 1.0);
  myInterface.SetParameter(0, "kDirTh", 0, TMath::Pi(), fTheta1, fFixDir, 0.01);
  myInterface.SetParameter(0, "kDirPhi", -TMath::Pi(), TMath::Pi(), fPhi1, fFixDir, 0.02);
  myInterface.SetParameter(0, "kEnergy", 100, 5000, 600, false, 250.);
  myInterface.SetParameter(0, "kConversionDistance", 0, 400, 50, false, 10.0);

  myInterface.SetParameter(1, "kVtxX", -1200, 1200, fVtxX, fFixVtx, 10.0);
  myInterface.SetParameter(1, "kVtxY", -1200, 1200, fVtxY, fFixVtx, 10.0);
  myInterface.SetParameter(1, "kVtxZ", -900, 900, fVtxZ, fFixVtx, 10.0);
  myInterface.SetParameter(1, "kVtxT", -100, 101000, 5000, false, 1.0);
  myInterface.SetParameter(1, "kDirTh", 0, TMath::Pi(), fTheta2, fFixDir, 0.01);
  myInterface.SetParameter(1, "kDirPhi", -TMath::Pi(), TMath::Pi(), fPhi2, fFixDir, 0.02);
  myInterface.SetParameter(1, "kEnergy", 100, 5000, 600, false, 250.);
  myInterface.SetParameter(1, "kConversionDistance", 0, 400, 50, false, 10.0);

  myInterface.JoinParametersTogether(0,1,"kVtxX");
  myInterface.JoinParametersTogether(0,1,"kVtxY");
  myInterface.JoinParametersTogether(0,1,"kVtxZ");
  myInterface.JoinParametersTogether(0,1,"kVtxT");

  // Set up some initial parameters
  /*if(myInterface.GetTrackType(0) == TrackType::MuonLike){
    WCSimParameters::Instance()->SetSlicerClusterDistance(250);
  }
  else{
    WCSimParameters::Instance()->SetSlicerClusterDistance(100);
  }
  WCSimParameters::Instance()->SetSlicerMinSize(10);
  WCSimParameters::Instance()->SetSlicerChargeCut(0.0);
  */

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
  //myInterface.Make2DSurface(0, "kVtxX", 0, "kVtxY");
  myInterface.SetMakeFits(kTRUE);
  myInterface.SetMakeSurfaces(kFALSE);



  myInterface.Print();
  myInterface.SetNumEventsToFit(fit);
  myInterface.SetFirstEventToFit(start);
  myInterface.Run();
  
  std::cout << "Done!" << std::endl;

}
