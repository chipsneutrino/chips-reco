#include "WCSimLikelihoodTrack.hh"
#include "WCSimTrackParameterEnums.hh"
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
WCSimLikelihoodTrack::WCSimLikelihoodTrack() : WCSimLikelihoodTrackBase()
{
}

WCSimLikelihoodTrack::WCSimLikelihoodTrack( double x, double y, double z, double t, double theta, double phi, double E, TrackType::Type myType )
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
	  std::cout << "WCSimLikelihoodTrack::Print():" << std::endl;
	  printf("Vertex = (%.02fcm,%.02fcm,%.02fcm,%.02fns)   Dir = (%.04f,%.04f)   E = %.03f  Type = %s",
           fVtx[0], fVtx[1], fVtx[2], fT0, fTheta0, fPhi0, fE0, TrackType::AsString(fType).c_str());

}

double WCSimLikelihoodTrack::GetTrackParameter(
		const FitterParameterType::Type &type) const {
	if( type == FitterParameterType::kVtxX ){ return GetX(); }
	if( type == FitterParameterType::kVtxY ){ return GetY(); }
	if( type == FitterParameterType::kVtxZ ){ return GetZ(); }
	if( type == FitterParameterType::kVtxT ){ return GetT(); }
	if( type == FitterParameterType::kDirTh ){ return GetTheta(); }
	if( type == FitterParameterType::kDirPhi ){ return GetPhi(); }
	if( type == FitterParameterType::kEnergy ){ return GetE(); }
	assert(0);
	return -99999;
}

TVector3 WCSimLikelihoodTrack::GetPropagatedPos(const Double_t &s) const
{
	return (this->GetVtx() + s * this->GetDir());
}

void WCSimLikelihoodTrack::SetType(const TrackType::Type &type)
{
	if(		type == TrackType::ElectronLike
		||	type == TrackType::MuonLike )
	{
		fType = type;
	}
	else
	{
		std::cerr << "Error in WCSimLikelihoodTrack::SetType: track type must be ElectronLike or MuonLike" << std::endl;
		assert(		type == TrackType::ElectronLike
				||	type == TrackType::MuonLike );
	}
}

