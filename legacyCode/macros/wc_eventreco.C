{
  
  
  // Path to WCSim ROOT file
  // =======================
  TString filename("../mcSamples/numi_numu_NC_100kT_10pC_Leigh.root");
  

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


  // Initialize Reconstruction
  // =========================
  WCSimReco* myReco = WCSimRecoFactory::Instance()->MakeReco();
  

  // Reconstruct Event
  // =================

  // get first entry
  WCSimInterface::LoadEvent(0);

  WCSimRecoEvent* recoEvent = WCSimInterface::RecoEvent();
  WCSimTrueEvent* trueEvent = WCSimInterface::TrueEvent();

  // run reconstruction
  myReco->Run(recoEvent);

  // get reconstructed quantities
  Double_t recoVtxX = recoEvent->GetVtxX();
  Double_t recoVtxY = recoEvent->GetVtxY();
  Double_t recoVtxZ = recoEvent->GetVtxZ();
  Double_t recoDirX = recoEvent->GetDirX();
  Double_t recoDirY = recoEvent->GetDirY();
  Double_t recoDirZ = recoEvent->GetDirZ();

  std::cout << "  recoVtx=(" << recoVtxX << ", " << recoVtxY << ", " << recoVtxZ << std::endl
            << "           " << recoDirX << ", " << recoDirY << ", " << recoDirZ << ") " << std::endl;

  // Interface to MC Truth
  //======================

  // get true quantities
  Double_t trueVtxX = trueEvent->GetVtxX();
  Double_t trueVtxY = trueEvent->GetVtxY();
  Double_t trueVtxZ = trueEvent->GetVtxZ();  
  Double_t trueDirX = trueEvent->GetDirX();
  Double_t trueDirY = trueEvent->GetDirY();
  Double_t trueDirZ = trueEvent->GetDirZ();

  std::cout << "  trueVtx=(" << trueVtxX << ", " << trueVtxY << ", " << trueVtxZ << std::endl
            << "           " << trueDirX << ", " << trueDirY << ", " << trueDirZ << ") " << std::endl;
}
