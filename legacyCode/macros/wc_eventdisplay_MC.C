 {
  gSystem->Load("libGeom");
  gSystem->Load("libEve");
  gSystem->Load("libMinuit.so");
  gSystem->Load("../WCSim/libWCSimRoot.so");
  gSystem->Load("./lib/libWCSimAnalysis.so");

  // Load Data
  // =========
  // WCSimInterface::LoadData("/lbne/software2010/WCSim/wcsim.test.pizero.root");
  WCSimInterface::LoadData("../mcSamples/piZero.root");
  
  // create viewer: must be called 'viewer'
  // =====================================
  WCSimDisplayViewer* viewer = new WCSimDisplayViewer();
 
  // type of display
  // ===============
  viewer->UseDisplay("AB"); 
  // viewer->UseDisplay("Vtx"); 
  // viewer->UseDisplay("Eve"); 

  // configuration
  // =============
  // viewer->DisplayRecoClusters();  // display clusters
  viewer->DisplayRecoEvent();        // display vertex
  viewer->DisplayTrueEvent();        // display truth

  // display configuration
  // =====================
  // viewer->SetPulseHeightCut(1.0);
  viewer->PrintGIF();
  // viewer->PrintEPS();

  // configuration
  // =============
  WCSimVertexFinder::UseTrueVertex(); // suppress vertex finder
  WCSimRingFinder::UseRecoVertex();   // suppress ring finder
  
}
