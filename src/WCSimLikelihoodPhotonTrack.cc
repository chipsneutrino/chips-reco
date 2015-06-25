/*
 * WCSimLikelihoodPhotonTrack.cc
 *
 *  Created on: 9 Jun 2015
 *      Author: andy
 */

#include "WCSimLikelihoodPhotonTrack.hh"

WCSimLikelihoodPhotonTrack::WCSimLikelihoodPhotonTrack() : WCSimLikelihoodTrackBase(){
	// TODO Auto-generated constructor stub
	fConversionDistance = 0.0;
	SetType(TrackType::PhotonLike);
}

WCSimLikelihoodPhotonTrack::WCSimLikelihoodPhotonTrack(double x, double y,
		double z, double t, double theta, double phi, double E,
		double convDistance) {
	fVtx[0] = x;
	fVtx[1] = y;
	fVtx[2] = z;
	fT0     = t;
	fTheta0 = theta;
	fPhi0   = phi;
	fE0     = E;
    fType   = TrackType::PhotonLike;
	fConversionDistance = convDistance;
}

WCSimLikelihoodPhotonTrack::~WCSimLikelihoodPhotonTrack() {
	// TODO Auto-generated destructor stub
}

double WCSimLikelihoodPhotonTrack::GetConversionDistance() const {
	return fConversionDistance;
}

void WCSimLikelihoodPhotonTrack::SetConversionDistance(const double& convDist) {
	fConversionDistance = convDist;
}

void WCSimLikelihoodPhotonTrack::Print() {
	  std::cout << "WCSimLikelihoodTrack::Print():" << std::endl
	            << "fVtx[0] = " << fVtx[0] << std::endl
	            << "fVtx[1] = " << fVtx[1] << std::endl
	            << "fVtx[2] = " << fVtx[2] << std::endl
	            << "fT0     = " << fT0     << std::endl
	            << "fTheta0 = " << fTheta0 << std::endl
	            << "fPhi0   = " << fPhi0   << std::endl
	            << "fE0     = " << fE0     << std::endl
	            << "fType   = " << fType   << " (" << TrackType::AsString(fType) << ")" << std::endl << std::endl
	  	  	  	<< "fConversionDistance = " << fConversionDistance << std::endl;
}

double WCSimLikelihoodPhotonTrack::GetTrackParameter(
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

TVector3 WCSimLikelihoodPhotonTrack::GetPropagatedPos(const Double_t &s) const
{
	return (this->GetVtx() + s * this->GetDir());
}

void WCSimLikelihoodPhotonTrack::SetType(
		const TrackType::Type& type) {
	if(	type == TrackType::PhotonLike )
	{
		fType = type;
	}
	else
	{
		std::cerr << "Error in WCSimLikelihoodPhotonTrack::SetType: track type must be PhotonLike" << std::endl;
		assert(	type == TrackType::PhotonLike );
	}
}
