#ifndef WCSIMVERTEXSEEDNTUPLE_HH
#define WCSIMVERTEXSEEDNTUPLE_HH

#include "TObject.h"
#include "TString.h"

#include "TFile.h"
#include "TTree.h"

#include "WCSimNtuple.hh"

class WCSimVertexSeedNtuple : public WCSimNtuple
{
 public:

  WCSimVertexSeedNtuple();
  ~WCSimVertexSeedNtuple();

  void Fill( WCSimTrueEvent* trueEvent,
             WCSimRecoEvent* recoEvent );

  void CloseFile();
  void Reset();

 private:

  void ResetVariables();
  void WriteVariables( WCSimTrueEvent* trueEvent,
                       WCSimRecoEvent* recoEvent);
  void WriteToFile();

  // TFile/TTree
  TFile* fRecoFile;
  TTree* fRecoTree;
  TString fRecoFileName;

  // True Quantities
  Int_t fRunNum;
  Int_t fEventNum;
  Int_t fTriggerNum;

  // Entries
  Int_t index;
  Int_t last;

  // Vertex
  Double_t fTrueVtxX;
  Double_t fTrueVtxY;
  Double_t fTrueVtxZ;
  Double_t fTrueDirX;
  Double_t fTrueDirY;
  Double_t fTrueDirZ;

  // Reco Vertex
  Double_t fRecoVtxX;
  Double_t fRecoVtxY;
  Double_t fRecoVtxZ;
  Double_t fRecoVtxTime;

  // Reco vs True
  Double_t fRecoDelta;
  Double_t fRecoDeltaTrans;
  Double_t fRecoDeltaLong;

  ClassDef(WCSimVertexSeedNtuple,0)
};

#endif
