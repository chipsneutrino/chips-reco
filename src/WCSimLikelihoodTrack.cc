#include "WCSimLikelihoodTrack.hh"
#include "WCSimTrueTrack.hh"
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

WCSimLikelihoodTrack::WCSimLikelihoodTrack(WCSimTrueTrack* trueTrack) {
	fVtx[0] = trueTrack->GetG4VtxX();
	fVtx[1] = trueTrack->GetG4VtxY();
	fVtx[2] = trueTrack->GetG4VtxZ();
	fT0 = 0;
	fTheta0 = TMath::ACos( trueTrack->GetDirZ() );
	fPhi0 = TMath::ATan2(trueTrack->GetVtxY(), trueTrack->GetVtxZ());
	fE0 = trueTrack->GetEnergy();
	fType = GetTypeFromPDG(trueTrack->GetTrackPDG());

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
double WCSimLikelihoodTrack::GetX() const {return fVtx[0]; }
double WCSimLikelihoodTrack::GetY() const {return fVtx[1]; }
double WCSimLikelihoodTrack::GetZ() const {return fVtx[2]; }
TVector3 WCSimLikelihoodTrack::GetVtx() const { return TVector3(fVtx[0], fVtx[1], fVtx[2]); }

double WCSimLikelihoodTrack::GetT() const {return fT0; }
double WCSimLikelihoodTrack::GetTheta() const {return fTheta0; }
double WCSimLikelihoodTrack::GetPhi() const {return fPhi0; }
double WCSimLikelihoodTrack::GetE() const {return fE0; }

double WCSimLikelihoodTrack::GetDirX() const { return TMath::Sin(fTheta0) * TMath::Cos(fPhi0); }
double WCSimLikelihoodTrack::GetDirY() const { return TMath::Sin(fTheta0) * TMath::Sin(fPhi0); }
double WCSimLikelihoodTrack::GetDirZ() const { return TMath::Cos(fTheta0); }
TVector3 WCSimLikelihoodTrack::GetDir() const { return TVector3( this->GetDirX(), this->GetDirY(), this->GetDirZ()); }

WCSimLikelihoodTrack::TrackType WCSimLikelihoodTrack::GetType() const { return fType;}

double WCSimLikelihoodTrack::GetTrackParameter(
		FitterParameterType::Type type) const {
	if( type == FitterParameterType::kVtxX ){ return GetX(); }
	if( type == FitterParameterType::kVtxY ){ return GetY(); }
	if( type == FitterParameterType::kVtxZ ){ return GetZ(); }
	if( type == FitterParameterType::kVtxT ){ return GetT(); }
	if( type == FitterParameterType::kDirTh ){ return GetPhi(); }
	if( type == FitterParameterType::kDirPhi ){ return GetTheta(); }
	if( type == FitterParameterType::kEnergy ){ return GetE(); }
	assert(0);
	return -99999;
}

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


TVector3 WCSimLikelihoodTrack::GetPropagatedPos(Double_t s) const
{
    return (this->GetVtx() + s * this->GetDir());
}

Int_t WCSimLikelihoodTrack::GetPDG() const {
	if( fType == MuonLike ){ return 13; }
	else if( fType == ElectronLike ){ return 11; }
	else
	{
		assert(    (std::cerr << "Could not get track PDG type for " << TrackTypeToString(fType) << std::endl)
				&& false);
	}
	return -999;
}
