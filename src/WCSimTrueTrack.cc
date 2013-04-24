#include "WCSimTrueTrack.hh"

#include <iostream>

ClassImp(WCSimTrueTrack)

WCSimTrueTrack::WCSimTrueTrack( Int_t ipdg, Int_t ipdgParent,
                                Double_t g4vx, Double_t g4vy, Double_t g4vz,
                                Double_t g4ex, Double_t g4ey, Double_t g4ez,
                                Double_t vx, Double_t vy, Double_t vz,
                                Double_t ex, Double_t ey, Double_t ez,
                                Double_t px, Double_t py, Double_t pz,
                                Double_t trkE, Double_t trkP )
{
  fIpdg = ipdg;
  fIpdgParent = ipdgParent;

  fG4VtxX = g4vx;
  fG4VtxY = g4vy;
  fG4VtxZ = g4vz;
 
  fG4EndX = g4ex;
  fG4EndY = g4ey;
  fG4EndZ = g4ez; 

  fVtxX = vx;
  fVtxY = vy;
  fVtxZ = vz;
 
  fEndX = ex;
  fEndY = ey;
  fEndZ = ez; 

  fDirX = px;
  fDirY = py;
  fDirZ = pz;

  fTrkP = trkP;
  fTrkE = trkE;
}

WCSimTrueTrack::~WCSimTrueTrack()
{

}
 
void WCSimTrueTrack::PrintTrack()
{
  std::cout << " *** WCSimTrueTrack::PrintTrack() *** " << std::endl
            << " * fIpdg = " << fIpdg << std::endl
            << " * fIpdgParent=" << fIpdgParent << std::endl         
            << " * fG4VtxX = " << fG4VtxX << std::endl
	    << " * fG4VtxY = " << fG4VtxY << std::endl
	    << " * fG4VtxZ = " << fG4VtxZ << std::endl
            << " * fG4EndX = " << fG4EndX << std::endl
	    << " * fG4EndY = " << fG4EndY << std::endl
	    << " * fG4EndZ = " << fG4EndZ << std::endl
	    << " * fVtxX = " << fVtxX << std::endl
	    << " * fVtxY = " << fVtxY << std::endl
	    << " * fVtxZ = " << fVtxZ << std::endl
	    << " * fEndX = " << fEndX << std::endl
	    << " * fEndY = " << fEndY << std::endl
	    << " * fEndZ = " << fEndZ << std::endl
	    << " * fDirX = " << fDirX << std::endl
	    << " * fDirY = " << fDirY << std::endl
	    << " * fDirZ = " << fDirZ << std::endl
	    << " * fTrkP = " << fTrkP << std::endl
	    << " * fTrkE = " << fTrkE << std::endl
            << " *********************************** " << std::endl;

  return;
}

