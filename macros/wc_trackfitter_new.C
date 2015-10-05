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

  WCSimFitterInterface myFitter;
  
  myFitter.SetInputFileName(filename.Data()); // For inclusion in the fitterPlots file
  myFitter.SetNumTracks(1);
  myFitter.SetTrackType(0, "MuonLike");

  // Set parameter(track number, "name", minimum, maximum, start, is fixed?)
  // Names are: kVtxX, kVtxY, kVtxZ, kDirTh, kDirPhi, kEnergy
  myFitter.SetParameter(0, "kVtxZ", -900, 900, 0, false);
  myFitter.SetParameter(0, "kVtxT", 900, 1000, 935, true);
  myFitter.SetParameter(0, "kDirTh", 0, TMath::Pi(), 0.0, false);
  myFitter.SetParameter(0, "kDirPhi", -TMath::Pi(), TMath::Pi(), 0.0, false);
  myFitter.SetParameter(0, "kEnergy", 100, 4500, 1000, false);
  myFitter.SetParameter(0, "kConversionDistance", 0, 250, 50, true);


  // Plot best-fit results
  myFitter.PlotForEachEvent("kVtxX");
  myFitter.PlotForEachEvent("kVtxY");
  myFitter.PlotForEachEvent("kVtxZ");
  myFitter.PlotForEachEvent("kVtxT");
  myFitter.PlotForEachEvent("kDirTh");
  myFitter.PlotForEachEvent("kDirPhi");
  myFitter.PlotForEachEvent("kEnergy");

  // Plot reco - true values
  myFitter.PlotRecoMinusTrue("kVtxX");
  myFitter.PlotRecoMinusTrue("kVtxY");
  myFitter.PlotRecoMinusTrue("kVtxZ");
  myFitter.PlotRecoMinusTrue("kVtxT");
  myFitter.PlotRecoMinusTrue("kDirTh");
  myFitter.PlotRecoMinusTrue("kDirPhi");
  myFitter.PlotRecoMinusTrue("kEnergy");

  // Sweep out one variable
  myFitter.SetNumSurfaceBins(22);
  myFitter.Make1DSurface(0, "kEnergy");
  //myFitter.Make2DSurface(0, "kVtxX", 0, "kVtxY");
  myFitter.SetMakeFits(kTRUE);
  myFitter.SetMakeSurfaces(kFALSE);



  myFitter.Print();
  myFitter.SetNumEventsToFit(fit);
  myFitter.SetFirstEventToFit(start);
  myFitter.Run();
  
  std::cout << "Done!" << std::endl;

}
