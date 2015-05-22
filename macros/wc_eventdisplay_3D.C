{
  gSystem->Load("libGeom");
  gSystem->Load("libEve");
  gSystem->Load("libMinuit.so");
  gSystem->Load("~/CHIPS/WCSim_fixGeneration/libWCSimRoot.so");
  gSystem->Load("lib/libWCSimAnalysis.so");


  // Load Data
  // =========
  //  WCSimInterface::LoadData("/lbne/software2010/WCSim/wcsim.test.beam.root");
  WCSimInterface::LoadData("/home/ajperch/CHIPS/WCSim_fixGeneration/localfile_muons.root");///unix/fnu/ajperch/sampleEvents/WCSimOutput/nue_chips_me_1000_rightVtx.root");//../mcSamples/muonZ.root");//numi_numu_NC_100kT_10pC_Leigh.root");


  // create viewer: must be called 'viewer'
  // =====================================
  WCSimDisplayViewer* viewer = new WCSimDisplayViewer();
 
  // for EVE display
  // ===============
  viewer->UseDisplay("EVE"); 

  // configuration
  // =============
  // viewer->DisplayRecoClusters();  // display clusters

  // not implemented:
  // ================
  // viewer->DisplayRecoVertex(); // display vertex
  // viewer->DisplayTrueEvent();  // display truth

}
