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

  Double_t fG4VtxX;  //< Where the Geant4 track object actually started (x)
  Double_t fG4VtxY;  //< Where the Geant4 track object actually started (y) 
  Double_t fG4VtxZ;  //< Where the Geant4 track object actually started (z) 
  
  Double_t fG4EndX;  //< Where the Geant4 track object actually ended (x)
  Double_t fG4EndY;  //< Where the Geant4 track object actually ended (y)
  Double_t fG4EndZ;  //< Where the Geant4 track object actually ended (z)

  Double_t fVtxX;  //< Where the track first entered the detector (x)
  Double_t fVtxY;  //< Where the track first entered the detector (y)
  Double_t fVtxZ;  //< Where the track first entered the detector (z) 
  
  Double_t fEndX; //< Where the track exited the detector (x)
  Double_t fEndY; //< Where the track exited the detector (y)
  Double_t fEndZ; //< Where the track exited the detector (z)

  Double_t fDirX;
  Double_t fDirY;
  Double_t fDirZ;

  ClassDef(WCSimTrueTrack,0)

};

#endif
