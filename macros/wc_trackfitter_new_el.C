void wc_trackfitter_new_el(const char * infile = "", int start=0, int fit=100){
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
  myFitter.SetTrackType(0, "ElectronLike");

  TVector3 dir(1,0,0);
  double theta = dir.Theta();
  double phi = dir.Phi();
 
  // Set parameter(track number, "name", minimum, maximum, start, is fixed?)
  // Names are: kVtxX, kVtxY, kVtxZ, kDirTh, kDirPhi, kEnergy
  myFitter.SetParameter(0, "kVtxX", -1200, 1200, 500, false, 10.0);
  myFitter.SetParameter(0, "kVtxY", -1200, 1200, 500, false, 10.0);
  myFitter.SetParameter(0, "kVtxZ",  -950, 950, 500, false, 10.0);
  myFitter.SetParameter(0, "kVtxT",  -100, 10100, 5000, false, 1.0);
  myFitter.SetParameter(0, "kDirTh",  0.0, TMath::Pi(), 0.25*TMath::Pi(), false, 0.01);
  myFitter.SetParameter(0, "kDirPhi", -1.0*TMath::Pi(), 1.0*TMath::Pi(), 0, false, 0.02);
  myFitter.SetParameter(0, "kEnergy", 500, 5000, 1000, false, 250.0);
  myFitter.SetParameter(0, "kConversionDistance", -50, 50, 0, true);


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
  myFitter.SetNumSurfaceBins(20);
  myFitter.Make1DSurface(0, "kEnergy");
  myFitter.Make1DSurface(0, "kVtxX");
  myFitter.Make1DSurface(0, "kVtxY");
  myFitter.Make1DSurface(0, "kVtxZ");
  myFitter.Make1DSurface(0, "kVtxT");
  myFitter.Make1DSurface(0, "kDirTh");
  myFitter.Make1DSurface(0, "kDirPhi");
  myFitter.Make1DSurface(0, "kEnergy");
  myFitter.Make2DSurface(0, "kVtxX", 0, "kVtxY");
  myFitter.Make2DSurface(0, "kVtxX", 0, "kEnergy");
  myFitter.Make2DSurface(0, "kVtxT", 0, "kEnergy");


  myFitter.SetMakeFits(kTRUE);
  myFitter.SetMakeSurfaces(kFALSE);

  myFitter.Print();
  myFitter.SetNumEventsToFit(fit);
  myFitter.SetFirstEventToFit(start);
  myFitter.Run();
  
  std::cout << "Done!" << std::endl;

}
