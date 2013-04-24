{
 
  gSystem->Load("libGeom");
  gSystem->Load("libEve");
  gSystem->Load("libMinuit");

  gSystem->Load("../WCSim/libWCSimRoot.so");
  gSystem->Load("./lib/libWCSimAnalysis.so");


  WCSimInterface::LoadData("../mcSamples/numi_numu_NC_100kT_10pC_Leigh.root");

  WCSimGeometry::WriteGeometry("wcsim.geometry.test.root");




}
