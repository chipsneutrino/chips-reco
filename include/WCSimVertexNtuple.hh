#ifndef WCSIMVERTEXNTUPLE_HH
#define WCSIMVERTEXNTUPLE_HH

#include "TObject.h"
#include "TString.h"

#include "TFile.h"
#include "TTree.h"

#include "WCSimNtuple.hh"

class WCSimVertexNtuple : public WCSimNtuple
{
 public:

  WCSimVertexNtuple();
  ~WCSimVertexNtuple();

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

  // Vertex
  Double_t fVtxX;
  Double_t fVtxY;
  Double_t fVtxZ;
  Double_t fVtxTime;
  Double_t fDirX;
  Double_t fDirY;
  Double_t fDirZ;

  Double_t fTrueVtxX;
  Double_t fTrueVtxY;
  Double_t fTrueVtxZ;
  Double_t fTrueDirX;
  Double_t fTrueDirY;
  Double_t fTrueDirZ;
  Double_t fTrueLength;

  // Digits
  Int_t index;
  Int_t last;

  Int_t fIsFiltered;

  Double_t fDigitX;
  Double_t fDigitY;
  Double_t fDigitZ;
  Double_t fDigitT;
  Double_t fDigitQ;

  Double_t fConeAngle;
  Double_t fAngle;
  Double_t fZenith;
  Double_t fAzimuth;
  Double_t fSolidAngle;

  Double_t fDistPoint;
  Double_t fDistTrack;
  Double_t fDistPhoton;
  Double_t fDistScatter;

  Double_t fDeltaTime;
  Double_t fDeltaSigma;

  Double_t fDeltaAngle;
  Double_t fDeltaPoint;
  Double_t fDeltaTrack;
  Double_t fDeltaPhoton;
  Double_t fDeltaScatter;

  Double_t fPointPath;
  Double_t fExtendedPath;
  Double_t fPointResidual;
  Double_t fExtendedResidual;

  // --- debug ---
  Double_t fExtendedResidualCorrected;

  ClassDef(WCSimVertexNtuple,0)
};

#endif 

