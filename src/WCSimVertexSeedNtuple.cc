#include "WCSimVertexSeedNtuple.hh"
#include "WCSimVertexGeometry.hh"

#include "TDirectory.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

#include <iostream>
#include <cmath>

ClassImp(WCSimVertexSeedNtuple)

WCSimVertexSeedNtuple::WCSimVertexSeedNtuple() : WCSimNtuple(),
  fRecoFile(0),
  fRecoTree(0)
{

}

WCSimVertexSeedNtuple::~WCSimVertexSeedNtuple()
{

}

void WCSimVertexSeedNtuple::Reset()
{
  std::cout << " *** WCSimVertexSeedNtuple::Reset() *** " << std::endl;

  fRecoFile = 0;
  fRecoTree = 0;

  return;
}

void WCSimVertexSeedNtuple::Fill( WCSimTrueEvent* trueEvent, WCSimRecoEvent* recoEvent )
{
  this->ResetVariables();
  this->WriteVariables( trueEvent, recoEvent );
}

void WCSimVertexSeedNtuple::ResetVariables()
{
 // reset variables
  // ===============
  fRunNum = -1;
  fEventNum = -1;
  fTriggerNum = -1;

  index = -1;
  last = -1;

  fTrueVtxX = 0.0;
  fTrueVtxY = 0.0;
  fTrueVtxZ = 0.0;
  fTrueDirX = 0.0;
  fTrueDirY = 0.0;
  fTrueDirZ = 0.0;

  fRecoVtxX = 0.0;
  fRecoVtxY = 0.0;
  fRecoVtxZ = 0.0;
  fRecoVtxTime = 0.0;

  fRecoDelta = -99999.9;
  fRecoDeltaTrans = -99999.9;
  fRecoDeltaLong = -99999.9;

  return;
}

void WCSimVertexSeedNtuple::WriteVariables( WCSimTrueEvent* fTrueEvent, WCSimRecoEvent* fRecoEvent )
{
  std::cout << " *** WCSimVertexSeedNtuple::WriteVariables() *** " << std::endl;

  // number of seeds (try 1000)
  // ==========================
  Int_t fNumSeeds = 1000;

  // calculate vertex seeds
  // ======================
  WCSimVertexGeometry::Instance()->CalcVertexSeeds(fRecoEvent,fNumSeeds);

  // header information
  // ==================
  fRunNum = fRecoEvent->GetRun();
  fEventNum = fRecoEvent->GetEvent();
  fTriggerNum = fRecoEvent->GetTrigger();

  // true vertex position
  // ====================
  fTrueVtxX  = fTrueEvent->GetVtxX();
  fTrueVtxY  = fTrueEvent->GetVtxY();
  fTrueVtxZ  = fTrueEvent->GetVtxZ();
  fTrueDirX  = fTrueEvent->GetDirX();
  fTrueDirY  = fTrueEvent->GetDirY();
  fTrueDirZ  = fTrueEvent->GetDirZ();

  // loop over vertex seeds
  // ======================
  UInt_t nlast = WCSimVertexGeometry::Instance()->GetNSeeds();

  for( UInt_t n=0; n<nlast; n++ ){
    index = (Int_t)n;
    last = (Int_t)nlast;

    fRecoVtxX = WCSimVertexGeometry::Instance()->GetSeedVtxX(n);
    fRecoVtxY = WCSimVertexGeometry::Instance()->GetSeedVtxY(n);
    fRecoVtxZ = WCSimVertexGeometry::Instance()->GetSeedVtxZ(n);
    fRecoVtxTime = WCSimVertexGeometry::Instance()->GetSeedVtxTime(n);

    Double_t dx = fRecoVtxX - fTrueVtxX;
    Double_t dy = fRecoVtxY - fTrueVtxY;
    Double_t dz = fRecoVtxZ - fTrueVtxZ;
    Double_t ds = sqrt( dx*dx + dy*dy + dz*dz );
    Double_t px = fTrueDirX;
    Double_t py = fTrueDirY;
    Double_t pz = fTrueDirZ;

    Double_t epsilon = 1.0e-7;

    if( ds>epsilon ){
      px = dx/ds;
      py = dy/ds;
      pz = dz/ds;
    }

    Double_t costh = px*fTrueDirX+py*fTrueDirY+pz*fTrueDirZ;
    Double_t sinth = 0.0;
    if( costh<1.0-epsilon ) sinth = sqrt(1.0-costh*costh);
    else costh = 1.0;

    fRecoDelta = ds;
    fRecoDeltaTrans = ds*sinth;
    fRecoDeltaLong = ds*costh;

    this->WriteToFile();
  }

  return;
}

void WCSimVertexSeedNtuple::WriteToFile()
{
  TDirectory* tmpd = 0;

  if( fRecoFile==0 ){
    tmpd = gDirectory;
    std::cout << " *** WCSimVertexSeedNtuple::OpenFile(...) *** " << std::endl;
    std::cout << "  opening file: " << fNtpFileName.Data() << std::endl;
    fRecoFile = new TFile(fNtpFileName.Data(),"recreate");
    fRecoTree = new TTree("ntuple","my analysis ntuple");
    fRecoTree->Branch("run",&fRunNum,"run/I");
    fRecoTree->Branch("event",&fEventNum,"event/I");
    fRecoTree->Branch("trigger",&fTriggerNum,"trigger/I");    
    fRecoTree->Branch("index",&index,"index/I");
    fRecoTree->Branch("last",&last,"last/I");
    fRecoTree->Branch("trueVtxX",&fTrueVtxX,"trueVtxX/D");
    fRecoTree->Branch("trueVtxY",&fTrueVtxY,"trueVtxY/D");
    fRecoTree->Branch("trueVtxZ",&fTrueVtxZ,"trueVtxZ/D");   
    fRecoTree->Branch("trueDirX",&fTrueDirX,"trueDirX/D");
    fRecoTree->Branch("trueDirY",&fTrueDirY,"trueDirY/D");
    fRecoTree->Branch("trueDirZ",&fTrueDirZ,"trueDirZ/D");
    fRecoTree->Branch("recoVtxX",&fRecoVtxX,"recoVtxX/D");
    fRecoTree->Branch("recoVtxY",&fRecoVtxY,"recoVtxY/D");
    fRecoTree->Branch("recoVtxZ",&fRecoVtxZ,"recoVtxZ/D"); 
    fRecoTree->Branch("recoVtxTime",&fRecoVtxTime,"recoVtxTime/D");  
    fRecoTree->Branch("recoDelta",&fRecoDelta,"recoDelta/D");
    fRecoTree->Branch("recoDeltaTrans",&fRecoDeltaTrans,"recoDeltaTrans/D");
    fRecoTree->Branch("recoDeltaLong",&fRecoDeltaLong,"recoDeltaLong/D");
    gDirectory = tmpd;
  }

  if( fRecoFile ){
    tmpd = gDirectory;
    fRecoFile->cd();
    fRecoTree->Fill();
    gDirectory = tmpd;
  }

  return;
}

void WCSimVertexSeedNtuple::CloseFile()
{
  TDirectory* tmpd = 0;

  if( fRecoFile ){
    tmpd = gDirectory;
    std::cout << " *** WCSimVertexSeedNtuple::CloseFile() *** " << std::endl;
    std::cout << "  closing file: " << fNtpFileName.Data() << std::endl;
    fRecoFile->cd();
    fRecoTree->Write();
    fRecoFile->Close();
    gDirectory = tmpd;
  }

  return;
}

