#ifndef WCSIMTRUETRACK_HH
#define WCSIMTRUETRACK_HH

#include "TObject.h"

class WCSimTrueTrack : public TObject {

 public:
  WCSimTrueTrack(Int_t ipdg, Int_t ipdgParent,
                 Double_t g4vx, Double_t g4vy, Double_t g4vz,
                 Double_t g4ex, Double_t g4ey, Double_t g4ez,
                 Double_t vx, Double_t vy, Double_t vz,
                 Double_t ex, Double_t ey, Double_t ez,
                 Double_t px, Double_t py, Double_t pz,
                 Double_t trkE, Double_t trkP);
  ~WCSimTrueTrack();

  Double_t GetG4VtxX() { return fG4VtxX; }
  Double_t GetG4VtxY() { return fG4VtxY; }
  Double_t GetG4VtxZ() { return fG4VtxZ; }

  Double_t GetG4EndX() { return fG4EndX; }
  Double_t GetG4EndY() { return fG4EndY; }
  Double_t GetG4EndZ() { return fG4EndZ; }

  Double_t GetVtxX() { return fVtxX; }
  Double_t GetVtxY() { return fVtxY; }
  Double_t GetVtxZ() { return fVtxZ; } 

  Double_t GetEndX() { return fEndX; }
  Double_t GetEndY() { return fEndY; }
  Double_t GetEndZ() { return fEndZ; }

  Double_t GetDirX() { return fDirX; }
  Double_t GetDirY() { return fDirY; }
  Double_t GetDirZ() { return fDirZ; }

  Double_t GetMomentum() { return fTrkP; }
  Double_t GetEnergy()   { return fTrkE; }

  Int_t GetTrackPDG()  { return fIpdg; }
  Int_t GetParentPDG() { return fIpdgParent; }
  
  void PrintTrack();

 private:
  
  Int_t fIpdg;
  Int_t fIpdgParent;

  Double_t fTrkP;
  Double_t fTrkE;

  Double_t fG4VtxX;
  Double_t fG4VtxY;
  Double_t fG4VtxZ; 

  Double_t fG4EndX;
  Double_t fG4EndY;
  Double_t fG4EndZ; 

  Double_t fVtxX;
  Double_t fVtxY;
  Double_t fVtxZ;
  
  Double_t fEndX;
  Double_t fEndY;
  Double_t fEndZ;

  Double_t fDirX;
  Double_t fDirY;
  Double_t fDirZ;

  ClassDef(WCSimTrueTrack,0)

};

#endif
