/*
 * WCSimLikelihoodPhotonTrack.cc
 *
 *  Created on: 9 Jun 2015
 *      Author: andy
 */

#include "WCSimLikelihoodPhotonTrack.hh"

WCSimLikelihoodPhotonTrack::WCSimLikelihoodPhotonTrack() : WCSimLikelihoodTrackBase(){
	// TODO Auto-generated constructor stub
	SetType(TrackType::PhotonLike);
	SetConversionDistance(0.0);
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
  SetConversionDistance(convDistance);

}

WCSimLikelihoodPhotonTrack::~WCSimLikelihoodPhotonTrack() {
	// TODO Auto-generated destructor stub
}

void WCSimLikelihoodPhotonTrack::SetConversionDistance(const double& convDist) {
	fConversionDistance = convDist;
  fFirstEmissionVtx = GetPropagatedPos(convDist); // Where are Cherenkov photons first emitted?
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

// Get the point at which the track first releases Cherenkov photons
// For a photon, this is found by travelling away from the vertex in the distance
// of track propagation by the conversion distance.  Work this out when we set the
// conversion distance so that this is just a lookup
TVector3 WCSimLikelihoodPhotonTrack::GetFirstEmissionVtx() const
{
  return fFirstEmissionVtx;
}
