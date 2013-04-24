#include "WCSimLikelihoodTrack.hh"

#ifndef REFLEX_DICTIONARY
ClassImp(WCSimLikelihoodTrack)
#endif

///////////////////////////////////////////////////////////////////////////
// Constructors
///////////////////////////////////////////////////////////////////////////
WCSimLikelihoodTrack::WCSimLikelihoodTrack()
{   
	fVtx[0] = 0;
	fVtx[1] = 0;
	fVtx[2] = 0;
	fT0 = 0;
	fTheta0 = 0;
	fPhi0 = 0;
	fE0 = 0;	
	return;
}

WCSimLikelihoodTrack::WCSimLikelihoodTrack( double x, double y, double z, double t, double theta, double phi, double E )
{   
	fVtx[0] = x;
	fVtx[1] = y;
	fVtx[2] = z;
	fT0 = t;
	fTheta0 = theta;
	fPhi0 = phi;
	fE0 = E;	
	return;
}

///////////////////////////////////////////////////////////////////////////
// Setters
///////////////////////////////////////////////////////////////////////////
void WCSimLikelihoodTrack::SetX(double x){ fVtx[0] = x; }
void WCSimLikelihoodTrack::SetY(double y){ fVtx[1] = y; }
void WCSimLikelihoodTrack::SetZ(double z){ fVtx[2] = z; }
void WCSimLikelihoodTrack::SetT(double t){ fT0 = t; }
void WCSimLikelihoodTrack::SetTheta(double th){ fTheta0 = th; }
void WCSimLikelihoodTrack::SetPhi(double phi){ fPhi0 = phi; }
void WCSimLikelihoodTrack::SetE(double E){ fE0 = E; }

///////////////////////////////////////////////////////////////////////////
// Getters
///////////////////////////////////////////////////////////////////////////
const double WCSimLikelihoodTrack::GetX(){return fVtx[0]; }
const double WCSimLikelihoodTrack::GetY(){return fVtx[1]; }
const double WCSimLikelihoodTrack::GetZ(){return fVtx[2]; }
const double WCSimLikelihoodTrack::GetT(){return fT0; }
const double WCSimLikelihoodTrack::GetTheta(){return fTheta0; }
const double WCSimLikelihoodTrack::GetPhi(){return fPhi0; }
const double WCSimLikelihoodTrack::GetE(){return fE0; }


///////////////////////////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////////////////////////
WCSimLikelihoodTrack::~WCSimLikelihoodTrack()
{
}

