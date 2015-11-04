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
	fDir[0] = (TMath::Sin(fTheta0) * TMath::Cos(fPhi0));
	fDir[1] = (TMath::Sin(fTheta0) * TMath::Sin(fPhi0));
	fDir[2] = (TMath::Cos(fTheta0));
	fE0     = E;
	fType   = myType;
	fConversionDistance = 0.0;
	return;
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
  if( type == FitterParameterType::kConversionDistance ) { return GetConversionDistance(); }
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
		assert(		type == TrackType::ElectronLike	||	type == TrackType::MuonLike );
	}
}

TVector3 WCSimLikelihoodTrack::GetFirstEmissionVtx() const { return GetVtx(); } // Conversion distance is 0 so it emits straight away
