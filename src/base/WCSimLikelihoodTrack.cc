#include "WCSimLikelihoodTrack.hh"
#include "WCSimTrackParameterEnums.hh"
#include "WCSimTrueTrack.hh"
#include "TMath.h"
#include "TVector3.h"
#include <string>
#include <iostream>

#ifndef REFLEX_DICTIONARY
ClassImp (WCSimLikelihoodTrack)
ClassImp (WCSimLikelihoodPhotonTrack)
ClassImp (WCSimLikelihoodUnknownTrack)
#endif

///////////////////////////////////////////////////////////////
//  METHODS FOR WCSIMLIKELIHOODTRACK CLASS
//////////////////////////////////////////////////////////////

WCSimLikelihoodTrack::WCSimLikelihoodTrack() :
		WCSimLikelihoodTrackBase() {
	fVtx[0] = -999;
	fVtx[1] = -999;
	fVtx[2] = -999;
	fT0 = -999;
	fTheta0 = -999;
	fPhi0 = -999;
	fDir[0] = (TMath::Sin(fTheta0) * TMath::Cos(fPhi0));
	fDir[1] = (TMath::Sin(fTheta0) * TMath::Sin(fPhi0));
	fDir[2] = (TMath::Cos(fTheta0));
	fE0 = -999;
	fType = TrackType::Unknown;
	fConversionDistance = 0.0;
	return;
}

WCSimLikelihoodTrack::WCSimLikelihoodTrack(double x, double y, double z, double t, double theta, double phi, double E,
		TrackType::Type myType) {
	fVtx[0] = x;
	fVtx[1] = y;
	fVtx[2] = z;
	fT0 = t;
	fTheta0 = theta;
	fPhi0 = phi;
	fDir[0] = (TMath::Sin(fTheta0) * TMath::Cos(fPhi0));
	fDir[1] = (TMath::Sin(fTheta0) * TMath::Sin(fPhi0));
	fDir[2] = (TMath::Cos(fTheta0));
	fE0 = E;
	fType = myType;
	fConversionDistance = 0.0;
	return;
}

WCSimLikelihoodTrack::~WCSimLikelihoodTrack() {
	// TODO Auto-generated destructor stub
}

double WCSimLikelihoodTrack::GetTrackParameter(const FitterParameterType::Type &type) const {
	if (type == FitterParameterType::kVtxX) {
		return GetX();
	}
	if (type == FitterParameterType::kVtxY) {
		return GetY();
	}
	if (type == FitterParameterType::kVtxZ) {
		return GetZ();
	}
	if (type == FitterParameterType::kVtxT) {
		return GetT();
	}
	if (type == FitterParameterType::kDirTh) {
		return GetTheta();
	}
	if (type == FitterParameterType::kDirPhi) {
		return GetPhi();
	}
	if (type == FitterParameterType::kEnergy) {
		return GetE();
	}
	if (type == FitterParameterType::kConversionDistance) {
		return GetConversionDistance();
	}
	assert(0);
	return -99999;
}

TVector3 WCSimLikelihoodTrack::GetPropagatedPos(const Double_t &s) const {
	return (this->GetVtx() + s * this->GetDir());
}

void WCSimLikelihoodTrack::SetType(const TrackType::Type &type) {
	if (type == TrackType::ElectronLike || type == TrackType::MuonLike) {
		fType = type;
	} else {
		std::cerr << "Error in WCSimLikelihoodTrack::SetType: track type must be ElectronLike or MuonLike" << std::endl;
		assert(type == TrackType::ElectronLike || type == TrackType::MuonLike);
	}
}

TVector3 WCSimLikelihoodTrack::GetFirstEmissionVtx() const {
	return GetVtx();
} // Conversion distance is 0 so it emits straight away

///////////////////////////////////////////////////////////////
//  METHODS FOR WCSIMLIKELIHOODPHOTONTRACK CLASS
//////////////////////////////////////////////////////////////

WCSimLikelihoodPhotonTrack::WCSimLikelihoodPhotonTrack() :
		WCSimLikelihoodTrackBase() {
	fVtx[0] = -999;
	fVtx[1] = -999;
	fVtx[2] = -999;
	fT0 = -999;
	fTheta0 = -999;
	fPhi0 = -999;
	fDir[0] = (TMath::Sin(fTheta0) * TMath::Cos(fPhi0));
	fDir[1] = (TMath::Sin(fTheta0) * TMath::Sin(fPhi0));
	fDir[2] = (TMath::Cos(fTheta0));
	fE0 = -999;
	fType = TrackType::Unknown;
	fConversionDistance = 0.0;
	return;
}

WCSimLikelihoodPhotonTrack::WCSimLikelihoodPhotonTrack(double x, double y, double z, double t, double theta, double phi,
		double E, double convDistance) {
	fVtx[0] = x;
	fVtx[1] = y;
	fVtx[2] = z;
	fTheta0 = theta;
	fPhi0 = phi;
	fDir[0] = (TMath::Sin(fTheta0) * TMath::Cos(fPhi0));
	fDir[1] = (TMath::Sin(fTheta0) * TMath::Sin(fPhi0));
	fDir[2] = (TMath::Cos(fTheta0));
	fT0 = t;
	fE0 = E;
	fType = TrackType::PhotonLike;
	SetConversionDistance(convDistance);
}

WCSimLikelihoodPhotonTrack::~WCSimLikelihoodPhotonTrack() {
	// TODO Auto-generated destructor stub
}

void WCSimLikelihoodPhotonTrack::SetConversionDistance(const double& convDist) {
	fConversionDistance = convDist;
	fFirstEmissionVtx = GetPropagatedPos(convDist); // Where are Cherenkov photons first emitted?
}

double WCSimLikelihoodPhotonTrack::GetTrackParameter(const FitterParameterType::Type &type) const {
	if (type == FitterParameterType::kVtxX) {
		return GetX();
	}
	if (type == FitterParameterType::kVtxY) {
		return GetY();
	}
	if (type == FitterParameterType::kVtxZ) {
		return GetZ();
	}
	if (type == FitterParameterType::kVtxT) {
		return GetT();
	}
	if (type == FitterParameterType::kDirTh) {
		return GetTheta();
	}
	if (type == FitterParameterType::kDirPhi) {
		return GetPhi();
	}
	if (type == FitterParameterType::kEnergy) {
		return GetE();
	}
	if (type == FitterParameterType::kConversionDistance) {
		return GetConversionDistance();
	}
	assert(0);
	return -99999;
}

TVector3 WCSimLikelihoodPhotonTrack::GetPropagatedPos(const Double_t &s) const {
	return (this->GetVtx() + s * this->GetDir());
}

void WCSimLikelihoodPhotonTrack::SetType(const TrackType::Type& type) {
	if (type == TrackType::PhotonLike) {
		fType = type;
	} else {
		std::cerr << "Error in WCSimLikelihoodPhotonTrack::SetType: track type must be PhotonLike" << std::endl;
		assert(type == TrackType::PhotonLike);
	}
}

// Get the point at which the track first releases Cherenkov photons
// For a photon, this is found by travelling away from the vertex in the distance
// of track propagation by the conversion distance.  Work this out when we set the
// conversion distance so that this is just a lookup
TVector3 WCSimLikelihoodPhotonTrack::GetFirstEmissionVtx() const {
	return fFirstEmissionVtx;
}

///////////////////////////////////////////////////////////////
//  METHODS FOR WCSIMLIKELIHOODUNKNOWNTRACK CLASS
//////////////////////////////////////////////////////////////

WCSimLikelihoodUnknownTrack::WCSimLikelihoodUnknownTrack() :
		WCSimLikelihoodTrackBase() {
	fVtx[0] = -999;
	fVtx[1] = -999;
	fVtx[2] = -999;
	fT0 = -999;
	fTheta0 = -999;
	fPhi0 = -999;
	fDir[0] = (TMath::Sin(fTheta0) * TMath::Cos(fPhi0));
	fDir[1] = (TMath::Sin(fTheta0) * TMath::Sin(fPhi0));
	fDir[2] = (TMath::Cos(fTheta0));
	fE0 = -999;
	fType = TrackType::Unknown;
	fConversionDistance = 0.0;
	return;
}

WCSimLikelihoodUnknownTrack::WCSimLikelihoodUnknownTrack(double x, double y, double z, double t, double theta,
		double phi, double E, double convDistance) {
	fVtx[0] = x;
	fVtx[1] = y;
	fVtx[2] = z;
	fTheta0 = theta;
	fPhi0 = phi;
	fDir[0] = (TMath::Sin(fTheta0) * TMath::Cos(fPhi0));
	fDir[1] = (TMath::Sin(fTheta0) * TMath::Sin(fPhi0));
	fDir[2] = (TMath::Cos(fTheta0));
	fT0 = t;
	fE0 = E;
	fType = TrackType::PhotonLike;
	SetConversionDistance(convDistance);
}

WCSimLikelihoodUnknownTrack::~WCSimLikelihoodUnknownTrack() {
	// TODO Auto-generated destructor stub
}

void WCSimLikelihoodUnknownTrack::SetConversionDistance(const double& convDist) {
	fConversionDistance = convDist;
	fFirstEmissionVtx = GetPropagatedPos(convDist); // Where are Cherenkov photons first emitted?
}

double WCSimLikelihoodUnknownTrack::GetTrackParameter(const FitterParameterType::Type &type) const {
	if (type == FitterParameterType::kVtxX) {
		return GetX();
	}
	if (type == FitterParameterType::kVtxY) {
		return GetY();
	}
	if (type == FitterParameterType::kVtxZ) {
		return GetZ();
	}
	if (type == FitterParameterType::kVtxT) {
		return GetT();
	}
	if (type == FitterParameterType::kDirTh) {
		return GetTheta();
	}
	if (type == FitterParameterType::kDirPhi) {
		return GetPhi();
	}
	if (type == FitterParameterType::kEnergy) {
		return GetE();
	}
	if (type == FitterParameterType::kConversionDistance) {
		return GetConversionDistance();
	}
	assert(0);
	return -99999;
}

TVector3 WCSimLikelihoodUnknownTrack::GetPropagatedPos(const Double_t &s) const {
	return (this->GetVtx() + s * this->GetDir());
}

void WCSimLikelihoodUnknownTrack::SetType(const TrackType::Type& type) {
	if (type == TrackType::Unknown) {
		fType = type;
	} else {
		std::cerr << "Error in WCSimLikelihoodPhotonTrack::SetType: track type must be Unknown" << std::endl;
		assert(type == TrackType::Unknown);
	}
}

// Get the point at which the track first releases Cherenkov photons
// For a photon, this is found by travelling away from the vertex in the distance
// of track propagation by the conversion distance.  Work this out when we set the
// conversion distance so that this is just a lookup
TVector3 WCSimLikelihoodUnknownTrack::GetFirstEmissionVtx() const {
	return fFirstEmissionVtx;
}
