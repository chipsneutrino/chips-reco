#ifndef WCSIMRECONTUPLE_HH
#define WCSIMRECONTUPLE_HH

#include "TObject.h"
#include "TString.h"

#include "TFile.h"
#include "TTree.h"

#include "WCSimNtuple.hh"

class WCSimRecoNtuple : public WCSimNtuple
{
 public:

  WCSimRecoNtuple();
  ~WCSimRecoNtuple();

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

  // True Quantities
  Int_t fRunNum;
  Int_t fEventNum;
  Int_t fTriggerNum;

  Int_t fTruePDG;
  Double_t fTrueP;
  Double_t fTrueE; 
  Double_t fTrueG4VtxX;
  Double_t fTrueG4VtxY;
  Double_t fTrueG4VtxZ;
  Double_t fTrueG4EndX;
  Double_t fTrueG4EndY;
  Double_t fTrueG4EndZ;
  Double_t fTrueVtxX;
  Double_t fTrueVtxY;
  Double_t fTrueVtxZ;
  Double_t fTrueEndX;
  Double_t fTrueEndY;
  Double_t fTrueEndZ;
  Double_t fTrueDirX;
  Double_t fTrueDirY;
  Double_t fTrueDirZ;
  Double_t fTrueLength;
  
  Double_t fTrueVtxTraceFwd;
  Double_t fTrueVtxTraceBwd;
  Double_t fTrueVtxDistToEdge;
  Double_t fTrueEndDistToEdge;

	Int_t fTrueTrackPDG[1000];				// Note that if you change these sizes
	Int_t fTrueTrackParentPDG[1000];	// you also have to change WCSimRecoNtuple::ResetVariables()
	Double_t fTrueTrackVtxX[1000];			// to zero the right number of entries
	Double_t fTrueTrackVtxY[1000];
	Double_t fTrueTrackVtxZ[1000];
	Double_t fTrueTrackDirX[1000];
	Double_t fTrueTrackDirY[1000];
	Double_t fTrueTrackDirZ[1000];

  // Reconstructed Quantities
  Int_t fNDigits;
  Int_t fNFilterDigits;
  Int_t fNRings;

  Int_t fFoundVertex;
  Int_t fFoundDirection;
  Int_t fFoundRings;

  Double_t fRecoVtxX;
  Double_t fRecoVtxY;
  Double_t fRecoVtxZ;
  Double_t fRecoVtxTime;
  Double_t fRecoDirX;
  Double_t fRecoDirY;
  Double_t fRecoDirZ;
  Double_t fRecoConeAngle;
  Double_t fRecoTrackLength;

  Double_t fRecoVtxFOM;
  Int_t fRecoVtxSteps;
  Int_t fRecoVtxPass;
  Int_t fRecoVtxStatus;

  Double_t fRecoVtxTraceFwd;
  Double_t fRecoVtxTraceBwd;
  Double_t fRecoVtxDistToEdge;

  Double_t fRecoRingVtxX;
  Double_t fRecoRingVtxY;
  Double_t fRecoRingVtxZ;  
  Double_t fRecoRingDirX;
  Double_t fRecoRingDirY;
  Double_t fRecoRingDirZ;
  Double_t fRecoRingAngle;
  Double_t fRecoRingHeight;
  Double_t fRecoRingHeightFrac;

  // Vertex: Reco minus True
  Double_t fRecoDelta;
  Double_t fRecoDeltaTrans;
  Double_t fRecoDeltaLong;
  Double_t fRecoDeltaAngle;
  Double_t fRecoRingDeltaAngle;

  // Vertex: Reconstruction Chain
  Double_t fSimplePositionX;
  Double_t fSimplePositionY;
  Double_t fSimplePositionZ;
  Double_t fSimplePositionDelta;
  Double_t fSimplePositionDeltaTrans;
  Double_t fSimplePositionDeltaLong;

  Double_t fPointPositionX;
  Double_t fPointPositionY;
  Double_t fPointPositionZ;
  Double_t fPointPositionFOM;
  Int_t fPointPositionSteps;
  Int_t fPointPositionPass;
  Double_t fPointPositionDelta;
  Double_t fPointPositionDeltaTrans;
  Double_t fPointPositionDeltaLong;

  Double_t fSimpleDirectionX;
  Double_t fSimpleDirectionY;
  Double_t fSimpleDirectionZ;
  Double_t fSimpleDirectionDeltaAngle;

  Double_t fPointDirectionX;
  Double_t fPointDirectionY;
  Double_t fPointDirectionZ;
  Double_t fPointDirectionFOM;
  Int_t fPointDirectionSteps;
  Int_t fPointDirectionPass;  
  Double_t fPointDirectionDeltaAngle;

  Double_t fPointVtxX;
  Double_t fPointVtxY;
  Double_t fPointVtxZ;
  Double_t fPointDirX;
  Double_t fPointDirY;
  Double_t fPointDirZ; 
  Double_t fPointConeAngle;
  Double_t fPointTrackLength;
  Double_t fPointVtxFOM;
  Int_t fPointVtxSteps;
  Int_t fPointVtxPass;  
  Double_t fPointVtxDelta;
  Double_t fPointVtxDeltaTrans;
  Double_t fPointVtxDeltaLong;
  Double_t fPointVtxDeltaAngle;

  Double_t fExtendedVtxX;
  Double_t fExtendedVtxY;
  Double_t fExtendedVtxZ;
  Double_t fExtendedDirX;
  Double_t fExtendedDirY;
  Double_t fExtendedDirZ;  
  Double_t fExtendedConeAngle;
  Double_t fExtendedTrackLength;
  Double_t fExtendedVtxFOM;
  Int_t fExtendedVtxSteps;
  Int_t fExtendedVtxPass;  
  Double_t fExtendedVtxDelta;
  Double_t fExtendedVtxDeltaTrans;
  Double_t fExtendedVtxDeltaLong;
  Double_t fExtendedVtxDeltaAngle;

  // --- debug ---
  Double_t fExtendedVtxTimeParam0;
  Double_t fExtendedVtxConeParam0;
  Double_t fExtendedVtxConeParam1;
  Double_t fExtendedVtxConeParam2;

  ClassDef(WCSimRecoNtuple,0)
};

#endif 

