{

  gSystem->Load("libGeom");
  gSystem->Load("libEve");
  gSystem->Load("libMinuit");
  gSystem->Load("../WCSim/libWCSimRoot.so");
  gSystem->Load("lib/libWCSimAnalysis.so");

  // Configuration
  // =============
   WCSimVertexFinder::UseTrueVertex(); // suppress vertex finder
   WCSimRingFinder::UseRecoVertex();      // suppress ring finder

  // Print Parameters
  // ================
  //WCSimParameters::PrintParameters();
  //WCSimDataCleaner::PrintParameters();
  //WCSimVertexFinder::PrintParameters();
  WCSimRingFinder::PrintParameters();

  // Load Data
  // =========
  WCSimInterface::LoadData("/home/ajperch/MINOS/waterCherenkov/mcSamples/electronGun_3GeV.root");
  //WCSimInterface::LoadData("/r01/lbne/water/wcsim_root_files/DUSEL_100kton_10inch_15perCent/muon_vtx_12001.root");

  // Create Ntuple Writer
  // ====================
  WCSimNtupleWriter* ntpwriter = new WCSimNtupleWriter();

  // type of ntuple
  // ==============
//  ntpwriter->UseNtuple("Reco");
   ntpwriter->UseNtuple("Vertex");
  // ntpwriter->UseNtuple("VertexSeed");

  // set file name
  // =============
  ntpwriter->SetFileName("/home/ajperch/MINOS/waterCherenkov/recoNtp/electronGun_3GeV.root");

  // loop over events
  // =================
  ntpwriter->Run(1000);

  WCSimRecoObjectTable::Instance()->Print();

}
