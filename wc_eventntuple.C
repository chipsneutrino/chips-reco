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
  WCSimParameters::PrintParameters();
  //WCSimDataCleaner::PrintParameters();
  WCSimVertexFinder::PrintParameters();
  WCSimRingFinder::PrintParameters();

  // Load Data
  // =========
  //WCSimInterface::LoadData("/home/ajperch/MINOS/waterCherenkov/mcSamples/numi_nue_100x20x30_10pC_1000evts_fullyRandom.root");
  WCSimInterface::LoadData("/home/ajperch/MINOS/waterCherenkov/mcSamples/cyl_tuning_nue_5mrad_1000_CC.root");
	//WCSimInterface::LoadData("/r01/lbne/water/wcsim_root_files/DUSEL_100kton_10inch_15perCent/muon_vtx_12001.root");

  // Create Ntuple Writer
  // ====================
  WCSimNtupleWriter* ntpwriter = new WCSimNtupleWriter();

  // type of ntuple
  // ==============
  //ntpwriter->UseNtuple("Reco");
   ntpwriter->UseNtuple("Vertex");
  // ntpwriter->UseNtuple("VertexSeed");

  // set file name
  // =============

 // ntpwriter->SetFileName("/home/ajperch/MINOS/waterCherenkov/recoNtp/split/mc_5mrad_numu_100evts_RECO.root");
  ntpwriter->SetFileName("/home/ajperch/MINOS/waterCherenkov/recoNtp/cyl_tuning_nue_5mrad_1000_CC.root");


  // loop over events
  // =================
  ntpwriter->Run(1000);

  WCSimRecoObjectTable::Instance()->Print();

}
