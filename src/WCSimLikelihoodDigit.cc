#include "WCSimLikelihoodDigit.hh"
#include "WCSimGeometry.hh"
#include "WCSimRootGeom.hh"
#include <iostream>
#include <cmath>

#ifndef REFLEX_DICTIONARY
ClassImp(WCSimLikelihoodDigit)
#endif

///////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////
WCSimLikelihoodDigit::WCSimLikelihoodDigit( Double_t x, Double_t y, Double_t z, Double_t t, Double_t Q, Int_t tubeId, Double_t faceX, Double_t faceY, Double_t faceZ )
{ 
  if(tubeId == 0) std::cerr << "Error, trying to create a PMT with tubeID = 0.  WCSim numbers from 1." << std::endl;
	fTubeId = tubeId;
    fPos[0] = x;
    fPos[1] = y;
    fPos[2] = z;
    fT = t;
    fQ = Q;
	fFace[0] = faceX;
	fFace[1] = faceY;
	fFace[2] = faceZ;
}

WCSimLikelihoodDigit::WCSimLikelihoodDigit( WCSimRootCherenkovDigiHit * myDigiHit )
{
	fTubeId = myDigiHit->GetTubeId();
  fQ = myDigiHit->GetQ();
  fT = myDigiHit->GetT();	
	WCSimRootGeom * myGeom = (WCSimRootGeom*)(WCSimGeometry::Instance())->GetWCSimGeometry();
    WCSimRootPMT myPMT = myGeom->GetPMTFromTubeID(fTubeId);
    fPos[0] = myPMT.GetPosition(0);
    fPos[1] = myPMT.GetPosition(1);
    fPos[2] = myPMT.GetPosition(2);
    fFace[0] = myPMT.GetOrientation(0);
    fFace[1] = myPMT.GetOrientation(1);
    fFace[2] = myPMT.GetOrientation(2);
}

///////////////////////////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////////////////////////
WCSimLikelihoodDigit::~WCSimLikelihoodDigit()
{
}

///////////////////////////////////////////////////////////////////////////
// Getters
///////////////////////////////////////////////////////////////////////////
int WCSimLikelihoodDigit::GetTubeId(){ return fTubeId;}
double WCSimLikelihoodDigit::GetQ(){ return fQ;}
double WCSimLikelihoodDigit::GetT(){ return fT;}
double WCSimLikelihoodDigit::GetX(){ return fPos[0]; }
double WCSimLikelihoodDigit::GetY(){ return fPos[1]; }
double WCSimLikelihoodDigit::GetZ(){ return fPos[2]; }
double WCSimLikelihoodDigit::GetFaceX(){ return fFace[0]; }
double WCSimLikelihoodDigit::GetFaceY(){ return fFace[1]; }
double WCSimLikelihoodDigit::GetFaceZ(){ return fFace[2]; }

void WCSimLikelihoodDigit::Print()
{
	Double_t fQ;
	Double_t fT;
	Double_t fPos[3];
	Double_t fFace[3];
	std::cout << "WCSimLikelihoodDigit::Print()" << std::endl
			  << "     fTubeId = " << fTubeId << std::endl
			  << "     fQ      = " << fQ << std::endl
			  << "     fT      = " << fT << std::endl
			  << "     fPos    = " << "(" << fPos[0] << ", " << fPos[1] << ", " << fPos[2] << ")" << std::endl
			  << "     fFace   = " << "(" << fFace[0] << ", " << fFace[1] << ", " << fFace[2] << ")" << std::endl << std::endl;
			
}
