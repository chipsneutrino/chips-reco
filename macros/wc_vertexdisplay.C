void wc_vertexdisplay( Int_t ievent )
{
  
  gSystem->Load("libGeom");
  gSystem->Load("libEve");
  gSystem->Load("libMinuit");

  gSystem->Load("../WCSim/libWCSimRoot.so");
  gSystem->Load("./lib/libWCSimAnalysis.so");


  // Load Data
  // =========
  WCSimInterface::LoadData("../mcSamples/numi_numu_NC_100kT_10pC_Leigh.root");

  // Load Event
  // ==========
  WCSimInterface::LoadEvent(ievent);

  // ring finder configuration
  // =========================
  WCSimRingFinder::UseRecoVertex();  // suppress ring finder

  // Initialize Viewer
  // =================
  WCSimVertexViewer* myVertexViewer = new WCSimVertexViewer();

  // Initialize Reconstruction
  // =========================
  WCSimReco* myReco = WCSimRecoFactory::Instance()->MakeReco();

  WCSimRecoEvent* myRecoEvent = WCSimInterface::RecoEvent();
  WCSimTrueEvent* myTrueEvent = WCSimInterface::TrueEvent();

  // run reconstruction
  myReco->Run(myRecoEvent);

  // truth
  Double_t trueVtxX = myTrueEvent->GetVtxX();
  Double_t trueVtxY = myTrueEvent->GetVtxY();
  Double_t trueVtxZ = myTrueEvent->GetVtxZ();  
  Double_t trueDirX = myTrueEvent->GetDirX();
  Double_t trueDirY = myTrueEvent->GetDirY();
  Double_t trueDirZ = myTrueEvent->GetDirZ();

  cout << " Monte Carlo Truth: " << endl;
  cout << "  TrueVtx=(" << trueVtxX << "," << trueVtxY << "," << trueVtxZ << ")" << endl;

  // Reconstruct Vertex
  // ==================
  WCSimVertexFinder* myVertexFinder = WCSimVertexFinder::Instance();
  
  
 

  // Extended Fit
  // ============
  WCSimRecoVertex* myExtendedVertex = myVertexFinder->RunExtendedFit(myRecoEvent);
  //WCSimRecoVertex* myExtendedVertex = myVertexFinder->RunExtendedFit(myRecoEvent,trueVtxX,trueVtxY,trueVtxZ,trueDirX,trueDirY,trueDirZ); // USING TRUTH
  myVertexViewer->DrawNewRecoEvent(myRecoEvent,myExtendedVertex);

  cout << " Result of Extended Fit: " << endl;
  cout << "  ExtendedVtx=(" << myExtendedVertex->GetX() << "," <<  myExtendedVertex->GetY() << "," << myExtendedVertex->GetZ() << "," << myExtendedVertex->GetTime() << "," << myExtendedVertex->GetFOM() << ")" << endl;

 

}
