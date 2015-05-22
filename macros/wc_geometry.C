{
 
  gSystem->Load("libGeom");
  gSystem->Load("libEve");
  gSystem->Load("libMinuit");

  gSystem->Load("../WCSim_fixGeneration/libWCSimRoot.so");
  gSystem->Load("./lib/libWCSimAnalysis.so");


  WCSimInterface::LoadData("/unix/fnu/ajperch/chips_10kton_20perCent_muons_10in.root");

  WCSimGeometry::WriteGeometry("wcsim.geometry.test.root");




}
