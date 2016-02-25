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
	fDir[0] = 0;
	fDir[1] = 0;
	fDir[2] = 0;
	fT0 = 0;
	fTheta0 = 0;
	fPhi0 = 0;
	fE0 = 0;
	fType = TrackType::Unknown;
	fConversionDistance = 0.0;
	return;
}



bool WCSimLikelihoodTrackBase::operator == (const WCSimLikelihoodTrackBase &b) const
{
	return (   fVtx[0] == b.fVtx[0]
              && fVtx[1] == b.fVtx[1]
              && fVtx[2] == b.fVtx[2]
              && fDir[0] == b.fDir[0]
              && fDir[1] == b.fDir[1]
              && fDir[2] == b.fDir[2]
              && fT0 == b.fT0
              && fTheta0 == b.fTheta0
              && fPhi0 == b.fPhi0
              && fE0 == b.fE0
              && fType == b.fType 
              && fConversionDistance == b.fConversionDistance );
}

///////////////////////////////////////////////////////////////////////////
// Setters
///////////////////////////////////////////////////////////////////////////
void WCSimLikelihoodTrackBase::SetX(double x){ fVtx[0] = x; }
void WCSimLikelihoodTrackBase::SetY(double y){ fVtx[1] = y; }
void WCSimLikelihoodTrackBase::SetZ(double z){ fVtx[2] = z; }
void WCSimLikelihoodTrackBase::SetT(double t){ fT0     = t; }
void WCSimLikelihoodTrackBase::SetTheta(double th)
{
	fTheta0 = th;
	fDir[0] = (TMath::Sin(fTheta0) * TMath::Cos(fPhi0));
	fDir[1] = (TMath::Sin(fTheta0) * TMath::Sin(fPhi0));
	fDir[2] = (TMath::Cos(fTheta0));
}
void WCSimLikelihoodTrackBase::SetPhi(double phi)
{
	fPhi0 = phi;
	fDir[0] = (TMath::Sin(fTheta0) * TMath::Cos(fPhi0));
	fDir[1] = (TMath::Sin(fTheta0) * TMath::Sin(fPhi0));
	fDir[2] = (TMath::Cos(fTheta0));
}

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
double WCSimLikelihoodTrackBase::GetConversionDistance() const { return fConversionDistance; }

double WCSimLikelihoodTrackBase::GetDirX() const { return fDir[0]; }
double WCSimLikelihoodTrackBase::GetDirY() const { return fDir[1]; }
double WCSimLikelihoodTrackBase::GetDirZ() const { return fDir[2]; }
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
              && fDir[0] == b->fDir[0]
              && fDir[1] == b->fDir[1]
              && fDir[2] == b->fDir[2]
              && fT0 == b->fT0
              && fTheta0 == b->fTheta0
              && fPhi0 == b->fPhi0
              && fE0 == b->fE0
              && fType == b->fType 
              && fConversionDistance == b->fConversionDistance);
}

int WCSimLikelihoodTrackBase::GetPDG() const
{
  return TrackType::GetPDGFromType(fType);
}

void WCSimLikelihoodTrackBase::Print()
{
	  printf("Vertex = (%.02fcm,%.02fcm,%.02fcm,%.02fns)   Dir = (%.04f,%.04f)   E = %.03f  Type = %s  Conv = %.02fcm\n",
           fVtx[0], fVtx[1], fVtx[2], fT0, fTheta0, fPhi0, fE0, TrackType::AsString(fType).c_str(), fConversionDistance);

}

double WCSimLikelihoodTrackBase::GetPropagationSpeedFrac() const{
    return GetPropagationSpeedFrac(fType);
}

double WCSimLikelihoodTrackBase::GetPropagationSpeedFrac(const TrackType &type)
{
  double speed = 1.0;
  switch(type){
    case TrackType::PhotonLike:
      // Fall through
    case TrackType::ElectronLike:
      speed = 0.8127;
      break;
    case TrackType::MuonLike:
      speed =  0.9395;
      break;
    default:
      break;
  }
  return speed;
}
