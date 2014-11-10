void wc_trackfitter_new(){
  // Path to WCSim ROOT file
  // =======================
  TString filename("/home/ajperch/CHIPS/WCSim/muons_1357_testEInterp.root");
  
  gApplication->ProcessLine(".except");

  // Load libraries
  // ==============
  gSystem->Load("libGeom");
  gSystem->Load("libEve");
  gSystem->Load("libMinuit");

  gSystem->Load("../WCSim/libWCSimRoot.so");
  gSystem->Load("./lib/libWCSimAnalysis.so");


  // Load Data
  // =========
  WCSimInterface::LoadData(filename.Data());
  
  WCSimFitterInterface::Instance()->SetNumTracks(1);

  // Set parameter(track number, "name", minimum, maximum, start, is fixed?)
  // Names are: kVtxX, kVtxY, kVtxZ, kDirTh, kDirPhi, kEnergy
  WCSimFitterInterface::Instance()->SetParameter(0, "kVtxX", -1900, 1900, 0, true);
  WCSimFitterInterface::Instance()->SetParameter(0, "kVtxY", -1900, 1900, 0, true);
  WCSimFitterInterface::Instance()->SetParameter(0, "kVtxZ", -900, 900, 0, true);
  WCSimFitterInterface::Instance()->SetParameter(0, "kVtxT", -900, 900, 0, true);
  WCSimFitterInterface::Instance()->SetParameter(0, "kDirTh", 0, TMath::Pi(), 0, true);
  WCSimFitterInterface::Instance()->SetParameter(0, "kDirPhi", -TMath::Pi(), TMath::Pi(), 0, true);
  WCSimFitterInterface::Instance()->SetParameter(0, "kEnergy", 1000, 2000, 1500, false);


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
  WCSimFitterInterface::Instance()->Set
  WCSimFitterInterface::Instance()->Make1DSurface(0, "kEnergy");



  WCSimFitterInterface::Instance()->Print();
  WCSimFitterInterface::Instance()->SetNumEventsToFit(1);
  WCSimFitterInterface::Instance()->Run();
  
  std::cout << "Done!" << std::endl;

}
