#include "WCSimLikelihoodTrack.hh"
#include "TMath.h"
#include "TVector3.h"
#include <string>
#include <iostream>

#ifndef REFLEX_DICTIONARY
ClassImp(WCSimLikelihoodTrack)
#endif

/*
 * Constructors
 */
WCSimLikelihoodTrack::WCSimLikelihoodTrack()
{   
	fVtx[0] = 0;
	fVtx[1] = 0;
	fVtx[2] = 0;
	fT0 = 0;
	fTheta0 = 0;
	fPhi0 = 0;
	fE0 = 0;	
    fType = WCSimLikelihoodTrack::Unknown;
	return;
}

WCSimLikelihoodTrack::WCSimLikelihoodTrack( double x, double y, double z, double t, double theta, double phi, double E, WCSimLikelihoodTrack::TrackType myType )
{   
	fVtx[0] = x;
	fVtx[1] = y;
	fVtx[2] = z;
	fT0     = t;
	fTheta0 = theta;
	fPhi0   = phi;
	fE0     = E;
    fType   = myType;
	return;
}

void WCSimLikelihoodTrack::Print()
{
  std::cout << "WCSimLikelihoodTrack::Print():" << std::endl
            << "fVtx[0] = " << fVtx[0] << std::endl
            << "fVtx[1] = " << fVtx[1] << std::endl
            << "fVtx[2] = " << fVtx[2] << std::endl
            << "fT0     = " << fT0     << std::endl
            << "fTheta0 = " << fTheta0 << std::endl
            << "fPhi0   = " << fPhi0   << std::endl
            << "fE0     = " << fE0     << std::endl
            << "fType   = " << fType   << std::endl << std::endl;
}

///////////////////////////////////////////////////////////////////////////
// Setters
///////////////////////////////////////////////////////////////////////////
void WCSimLikelihoodTrack::SetX(double x){ fVtx[0] = x; }
void WCSimLikelihoodTrack::SetY(double y){ fVtx[1] = y; }
void WCSimLikelihoodTrack::SetZ(double z){ fVtx[2] = z; }
void WCSimLikelihoodTrack::SetT(double t){ fT0     = t; }
void WCSimLikelihoodTrack::SetTheta(double th){ fTheta0 = th; }
void WCSimLikelihoodTrack::SetPhi(double phi){ fPhi0 = phi; }
void WCSimLikelihoodTrack::SetE(double E){ fE0 = E; }
void WCSimLikelihoodTrack::SetType(WCSimLikelihoodTrack::TrackType type){ fType = type; }

///////////////////////////////////////////////////////////////////////////
// Getters
///////////////////////////////////////////////////////////////////////////
const double WCSimLikelihoodTrack::GetX(){return fVtx[0]; }
const double WCSimLikelihoodTrack::GetY(){return fVtx[1]; }
const double WCSimLikelihoodTrack::GetZ(){return fVtx[2]; }
const TVector3 WCSimLikelihoodTrack::GetVtx(){ return TVector3(fVtx[0], fVtx[1], fVtx[2]); }

const double WCSimLikelihoodTrack::GetT(){return fT0; }
const double WCSimLikelihoodTrack::GetTheta(){return fTheta0; }
const double WCSimLikelihoodTrack::GetPhi(){return fPhi0; }
const double WCSimLikelihoodTrack::GetE(){return fE0; }

const double WCSimLikelihoodTrack::GetDirX(){ return TMath::Sin(fTheta0) * TMath::Cos(fPhi0); }
const double WCSimLikelihoodTrack::GetDirY(){ return TMath::Sin(fTheta0) * TMath::Sin(fPhi0); }
const double WCSimLikelihoodTrack::GetDirZ(){ return TMath::Cos(fTheta0); }
const TVector3 WCSimLikelihoodTrack::GetDir(){ return TVector3( this->GetDirX(), this->GetDirY(), this->GetDirZ()); }

const WCSimLikelihoodTrack::TrackType WCSimLikelihoodTrack::GetType(){ return fType;}

///////////////////////////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////////////////////////
WCSimLikelihoodTrack::~WCSimLikelihoodTrack()
{
}

/**
 * Output strings of our tracktype enums
 * @param myType The particle type of the track
 * @return A string corresponding the the particle type of the track
 */
std::string WCSimLikelihoodTrack::TrackTypeToString( WCSimLikelihoodTrack::TrackType myType )
{
  std::string type;

  switch(myType)
  {
    case WCSimLikelihoodTrack::ElectronLike:
      type = "electron";
      break;
    case WCSimLikelihoodTrack::MuonLike:
      type = "muon";
      break;
    default:
      type = "UNKOWN";
      break;
  }
  return type;
}


const TVector3 WCSimLikelihoodTrack::GetPropagatedPos(Double_t s)
{
    return (this->GetVtx() + s * this->GetDir());
}
