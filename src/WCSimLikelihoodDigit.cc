#include "WCSimLikelihoodDigit.hh"
#include "WCSimGeometry.hh"
#include "WCSimRootGeom.hh"
#include "TGraph.h"
#include <iostream>
#include <cmath>
#include <cstdlib>

#ifndef REFLEX_DICTIONARY
ClassImp(WCSimLikelihoodDigit)
#endif

WCSimLikelihoodDigit::WCSimLikelihoodDigit(Double_t x, Double_t y, Double_t z,
        Double_t t, Double_t Q, Int_t tubeId, Double_t faceX, Double_t faceY,
        Double_t faceZ, TString pmtName, TGraph * wlWeightedQE, double wlWeightedRefIndex, double exposeHeight)
{
    if (tubeId == 0)
    {
        std::cerr << "Error, trying to create a PMT with tubeID = 0.  WCSim numbers from 1."
                  << std::endl;
        exit(EXIT_FAILURE);
    }
    fTubeId  = tubeId;
    fPos[0]  = x;
    fPos[1]  = y;
    fPos[2]  = z;
    fT       = t;
    fQ       = Q;
    fFace[0] = faceX;
    fFace[1] = faceY;
    fFace[2] = faceZ;
    fPMTName = pmtName;
    fAverageQE = wlWeightedQE;
    fAverageRefIndex = wlWeightedRefIndex;
    fExposeHeight = exposeHeight;
}


WCSimLikelihoodDigit::WCSimLikelihoodDigit( WCSimRootCherenkovDigiHit * myDigiHit)
{
    fTubeId = myDigiHit->GetTubeId();
    fQ      = myDigiHit->GetQ();
    fT      = myDigiHit->GetT();

    WCSimRootGeom * myGeom = (WCSimRootGeom*) (WCSimGeometry::Instance())->GetWCSimGeometry();
    WCSimRootPMT myPMT     = myGeom->GetPMTFromTubeID(fTubeId);
    fPos[0]  = myPMT.GetPosition(0);
    fPos[1]  = myPMT.GetPosition(1);
    fPos[2]  = myPMT.GetPosition(2);
    fFace[0] = myPMT.GetOrientation(0);
    fFace[1] = myPMT.GetOrientation(1);
    fFace[2] = myPMT.GetOrientation(2);
    fPMTName = myPMT.GetPMTName();
    fAverageQE = WCSimDetectorParameters::WavelengthAveragedQE(fPMTName.Data());
    fAverageRefIndex = WCSimDetectorParameters::QEAveragedRefIndex(fPMTName.Data());
    fExposeHeight = WCSimDetectorParameters::PMTExposeHeight(fPMTName.Data());
}

WCSimLikelihoodDigit::WCSimLikelihoodDigit(const WCSimLikelihoodDigit &otherLikelihoodDigit)
{
  fTubeId = otherLikelihoodDigit.fTubeId;
  fQ = otherLikelihoodDigit.fQ;
  fT = otherLikelihoodDigit.fT;

  fPos[0] = otherLikelihoodDigit.fPos[0];
  fPos[1] = otherLikelihoodDigit.fPos[1];
  fPos[2] = otherLikelihoodDigit.fPos[2];

  fFace[0] = otherLikelihoodDigit.fFace[0];
  fFace[1] = otherLikelihoodDigit.fFace[1];
  fFace[2] = otherLikelihoodDigit.fFace[2];

  fPMTName = otherLikelihoodDigit.fPMTName;
  fAverageQE = otherLikelihoodDigit.fAverageQE;
  fAverageRefIndex = otherLikelihoodDigit.fAverageRefIndex;
  fExposeHeight = otherLikelihoodDigit.fExposeHeight;
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
int WCSimLikelihoodDigit::GetTubeId() const
{
    return fTubeId;
}
double WCSimLikelihoodDigit::GetQ() const
{
    return fQ;
}
double WCSimLikelihoodDigit::GetT() const
{
    return fT;
}

TVector3 WCSimLikelihoodDigit::GetPos() const
{
    return TVector3(fPos[0], fPos[1], fPos[2]);
}
double WCSimLikelihoodDigit::GetX() const
{
    return fPos[0];
}
double WCSimLikelihoodDigit::GetY() const
{
    return fPos[1];
}
double WCSimLikelihoodDigit::GetZ() const
{
    return fPos[2];
}

TVector3 WCSimLikelihoodDigit::GetFace() const
{
    return TVector3(fFace[0], fFace[1], fFace[2]);
}
double WCSimLikelihoodDigit::GetFaceX() const
{
    return fFace[0];
}
double WCSimLikelihoodDigit::GetFaceY() const
{
    return fFace[1];
}
double WCSimLikelihoodDigit::GetFaceZ() const
{
    return fFace[2];
}

double WCSimLikelihoodDigit::GetAverageQE(const double &distanceToPMT) const
{
    double xLo, yLo;
    double xHi, yHi;
    fAverageQE->GetPoint(0, xLo, yLo);
    fAverageQE->GetPoint(fAverageQE->GetN()-1, xHi, yHi);
    if( distanceToPMT < xLo )
    {
      std::cerr << "WCSimLikelihoodDigit::GetAverageQE - Warning: distance to PMT of " << distanceToPMT << " is less than the minimum x-value in the graph, of " << xLo << std::endl;
    }
    if( distanceToPMT > xHi )
    {
      std::cerr << "WCSimLikelihoodDigit::GetAverageQE - Warning: distance to PMT of " << distanceToPMT << " is greater than the maximum x-value in the graph, of " << xHi << std::endl;
    }
    return fAverageQE->Eval(distanceToPMT);
}

double WCSimLikelihoodDigit::GetAverageRefIndex() const
{
    return fAverageRefIndex;
}
double WCSimLikelihoodDigit::GetExposeHeight() const
{
  return fExposeHeight;
}

void WCSimLikelihoodDigit::Print() const
{
    std::cout << "WCSimLikelihoodDigit::Print()" << std::endl
            << "     fTubeId = " << fTubeId << std::endl
            << "     fQ      = " << fQ << std::endl
            << "     fT      = " << fT << std::endl
            << "     fPos    = " << "(" << fPos[0] 
            << ", " << fPos[1] << ", " << fPos[2] << ")" << std::endl
            << "     fFace   = " << "(" << fFace[0] << ", " << fFace[1] << ", " << fFace[2] << ")" << std::endl
            << "     fExposeHeight = " << fExposeHeight << std::endl;
    return;
}

TString WCSimLikelihoodDigit::GetPMTName() const {
	return fPMTName;
}

bool WCSimLikelihoodDigit::operator == (const WCSimLikelihoodDigit &other) const
{
	return (
		       fTubeId == other.fTubeId     ///< Unique PMT ID number from WCSim
		        && fQ == other.fQ       ///< Digitized charge (P.E.)
		        && fT == other.fT       ///< Time of hit
		        && fPos[0] == other.fPos[0] ///< (x,y,z) co-ordinates of the PMT location
            && fPos[1] == other.fPos[1]
            && fPos[2] == other.fPos[2]
		        && fFace[0] == other.fFace[0] ///< (x,y,z) components of the direction normal to the PMT
            && fFace[1] == other.fFace[1]
            && fFace[2] == other.fFace[2]
		        && fPMTName == other.fPMTName ///< Name of PMT type, e.g. 3_inch_HQE
            && fAverageQE == other.fAverageQE ///< Average QE of the PMT, from weighting QE(wavelength) by the average Cherenkov spectrum - as a function of photon travel distance
		        && fAverageRefIndex == other.fAverageRefIndex ///< Weight WCSim's refractive index by (wavelength * PMT QE(wavelength))
		        && fExposeHeight == other.fExposeHeight ///< Height of PMT dome expose through the detector liner
			  );
}
