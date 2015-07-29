/*
 * WCSimLikelihoodTrackBase.cc
 *
 *  Created on: 10 Jun 2015
 *      Author: andy
 */

#include "WCSimLikelihoodTrackBase.hh"
#include "WCSimTrackParameterEnums.hh"

WCSimLikelihoodTrackBase::WCSimLikelihoodTrackBase()
{
	fVtx[0] = 0;
	fVtx[1] = 0;
	fVtx[2] = 0;
	fT0 = 0;
	fTheta0 = 0;
	fPhi0 = 0;
	fE0 = 0;
    fType = TrackType::Unknown;
	return;
}



bool WCSimLikelihoodTrackBase::operator == (const WCSimLikelihoodTrackBase &b) const
{
	return (   fVtx[0] == b.fVtx[0]
              && fVtx[1] == b.fVtx[1]
              && fVtx[2] == b.fVtx[2]
              && fT0 == b.fT0
              && fTheta0 == b.fTheta0
              && fPhi0 == b.fPhi0
              && fE0 == b.fE0
              && fType == b.fType );
}

///////////////////////////////////////////////////////////////////////////
// Setters
///////////////////////////////////////////////////////////////////////////
void WCSimLikelihoodTrackBase::SetX(double x){ fVtx[0] = x; }
void WCSimLikelihoodTrackBase::SetY(double y){ fVtx[1] = y; }
void WCSimLikelihoodTrackBase::SetZ(double z){ fVtx[2] = z; }
void WCSimLikelihoodTrackBase::SetT(double t){ fT0     = t; }
void WCSimLikelihoodTrackBase::SetTheta(double th){ fTheta0 = th; }
void WCSimLikelihoodTrackBase::SetPhi(double phi){ fPhi0 = phi; }
void WCSimLikelihoodTrackBase::SetE(double E){ fE0 = E; }

///////////////////////////////////////////////////////////////////////////
// Getters
///////////////////////////////////////////////////////////////////////////
double WCSimLikelihoodTrackBase::GetX() const {return fVtx[0]; }
double WCSimLikelihoodTrackBase::GetY() const {return fVtx[1]; }
double WCSimLikelihoodTrackBase::GetZ() const {return fVtx[2]; }
TVector3 WCSimLikelihoodTrackBase::GetVtx() const { return TVector3(fVtx[0], fVtx[1], fVtx[2]); }

double WCSimLikelihoodTrackBase::GetT() const {return fT0; }
double WCSimLikelihoodTrackBase::GetTheta() const {return fTheta0; }
double WCSimLikelihoodTrackBase::GetPhi() const {return fPhi0; }
double WCSimLikelihoodTrackBase::GetE() const {return fE0; }

double WCSimLikelihoodTrackBase::GetDirX() const { return (TMath::Sin(fTheta0) * TMath::Cos(fPhi0)); }
double WCSimLikelihoodTrackBase::GetDirY() const { return (TMath::Sin(fTheta0) * TMath::Sin(fPhi0)); }
double WCSimLikelihoodTrackBase::GetDirZ() const { return (TMath::Cos(fTheta0)); }
TVector3 WCSimLikelihoodTrackBase::GetDir() const { return TVector3( this->GetDirX(), this->GetDirY(), this->GetDirZ()); }

TrackType WCSimLikelihoodTrackBase::GetType() const{ return fType; }

///////////////////////////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////////////////////////
WCSimLikelihoodTrackBase::~WCSimLikelihoodTrackBase()
{
}


Bool_t WCSimLikelihoodTrackBase::EnergyGreaterThanOrEqual(const WCSimLikelihoodTrackBase &a, const WCSimLikelihoodTrackBase &b)
{
  return (a.GetE() >= b.GetE());
}

Bool_t WCSimLikelihoodTrackBase::EnergyGreaterThanOrEqualPtrs(WCSimLikelihoodTrackBase *a, WCSimLikelihoodTrackBase *b)
{
  return EnergyGreaterThanOrEqual(*a, *b);
}

bool WCSimLikelihoodTrackBase::IsSameTrack(WCSimLikelihoodTrackBase * b) const
{
	return (   fVtx[0] == b->fVtx[0]
              && fVtx[1] == b->fVtx[1]
              && fVtx[2] == b->fVtx[2]
              && fT0 == b->fT0
              && fTheta0 == b->fTheta0
              && fPhi0 == b->fPhi0
              && fE0 == b->fE0
              && fType == b->fType );
}

int WCSimLikelihoodTrackBase::GetPDG() const
{
  return TrackType::GetPDGFromType(fType);
}