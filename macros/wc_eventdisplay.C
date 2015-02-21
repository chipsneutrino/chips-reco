 {
  gSystem->Load("libGeom");
  gSystem->Load("libEve");
  gSystem->Load("libMinuit.so");
  gSystem->Load("$FNU/software/WCSim_github/libWCSimRoot.so");
  gSystem->Load("lib/libWCSimAnalysis.so");

  // Load Data
  // =========
  //  WCSimInterface::LoadData("../mcSamples/numi_numu_1_100kT_10pC.root");
  WCSimInterface::LoadData("/unix/fnu/ajperch/WCSim/pions.root");
  WCSimInterface::LoadData("fitterPlots_2015216_1627.root");

  // create viewer: must be called 'viewer'
  // =====================================
  WCSimDisplayViewer* viewer = new WCSimDisplayViewer();
 
  // type of display
  // ===============
  //viewer->UseDisplay("AB"); 
  //viewer->UseDisplay("Vtx"); 
  //viewer->UseDisplay("Eve"); 

  // configuration
  // =============
  // viewer->DisplayRecoClusters();  // display clusters
  //viewer->DisplayRecoEvent();        // display vertex
  viewer->DisplayTrueEvent();        // display truth

  // display configuration
  // =====================
  viewer->SetPulseHeightCut(2.0);
  viewer->PrintGIF();
  viewer->PrintEPS();

  // configuration
  // =============
 //  WCSimRingFinder::Instance()->SetUsingTSpectrum2();
   WCSimVertexFinder::UseTrueVertex();	 // suppress vertex finder
   WCSimRingFinder::UseRecoVertex();   // suppress ring finder
  
}
