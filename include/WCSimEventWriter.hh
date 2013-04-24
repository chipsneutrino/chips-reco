#ifndef WCSIMEVENTWRITER_HH
#define WCSIMEVENTWRITER_HH

#include "TObject.h"

class TFile;
class TTree;
class TChain;

class WCSimRootEvent;
class WCSimRootGeom;

class WCSimEventWriter : public TObject {

 public:
  WCSimEventWriter();
  ~WCSimEventWriter();

  void AddFile(const char* file);

  void WriteOutGeometry(Bool_t yesno = 1){
    fWriteOutGeometry = yesno;
  }

  void WriteOutTracks(Bool_t yesno = 1){
    fWriteOutTracks = yesno;
  }
  
  void WriteOutDigits(Bool_t yesno = 1){
    fWriteOutDigits = yesno;
  }

  void Run(Int_t entries = -1);

 private:

  void Initialize();

  void WriteThisDigit();
  void WriteThisTrack();

  void CloseFiles();

  Double_t fScale;

  Int_t fPMTs;
  Double_t* fPmtX;
  Double_t* fPmtY;
  Double_t* fPmtZ;
  Int_t* fPmtRegion;

  Int_t fRun;
  Int_t fEventNumber;
  Int_t fSubEventNumber;

  Int_t fTube;
  Int_t fRegion;
  Double_t fXpos;
  Double_t fYpos;
  Double_t fZpos;
  Double_t fQPEs;
  Double_t fTime;

  Int_t fTrackNumber;
  Int_t fPdg;
  Int_t fFlag;
  Double_t fM;
  Double_t fP;
  Double_t fE;
  Double_t fTrkDirX;
  Double_t fTrkDirY;
  Double_t fTrkDirZ;
  Double_t fVtxX;
  Double_t fVtxY;
  Double_t fVtxZ;
  Double_t fVtxTime;
  Double_t fEndX;
  Double_t fEndY;
  Double_t fEndZ;
  Int_t fPdgParent;

  TFile* fFileDigits;
  TTree* fTreeDigits;

  TFile* fFileTracks;
  TTree* fTreeTracks;

  TFile* fFileGeometry;
  TTree* fTreeGeometry;  

  Bool_t fWriteOutGeometry;
  Bool_t fWriteOutTracks;
  Bool_t fWriteOutDigits;

  TChain* fChain;
  TChain* fChainGeom;

  WCSimRootEvent* fEvent;
  WCSimRootGeom* fGeometry;

  ClassDef(WCSimEventWriter,0)

};

#endif
