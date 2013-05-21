#include "WCSimChargeLikelihood.hh"
#include "WCSimRecoDigit.hh"
#include "WCSimLikelihoodDigit.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimLikelihoodTrack.hh"

#include "TCollection.h"
#include "TMath.h"

#include <iostream>
#include <cmath>

///////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////
WCSimChargeLikelihood::WCSimChargeLikelihood(WCSimLikelihoodDigitArray * myDigitArray, WCSimLikelihoodTrack * myTrack )
{   
    std::cout << "*** WCSimChargeLikelihood::WCSimChargeLikelihood() *** Created charge likelihood object" << std::endl;
    this->Initialize( myDigitArray, myTrack );
}

///////////////////////////////////////////////////////////////////////////
// Set starting values. TODO: Eventually WCSimLikelihoodTuner will be used
// to work these out, but for now they're hardcoded
///////////////////////////////////////////////////////////////////////////
void WCSimChargeLikelihood::Initialize( WCSimLikelihoodDigitArray * myDigitArray, WCSimLikelihoodTrack * myTrack )
{
    std::cout << "*** WCSimChargeLikelihood::Initialize() *** Initializing charge likelihood with tuned values" << std::endl; 
    fCoeffsCh.push_back(1);
    fCoeffsCh.push_back(2);
    fCoeffsCh.push_back(3);

    fCoeffsInd.push_back(4);
    fCoeffsInd.push_back(5);
    fCoeffsInd.push_back(6);
 
    fTrack = myTrack;
    fGotTrackParameters = false;
    fDigitArray = myDigitArray;
    fEnergyInterval = 0.1;    // steps by which to increment variables 
    fR0Interval = 10;        // in the lookup table from one bin
    fCosTheta0Interval = 0.1; // to the next

    for(int i = 0; i < 100; ++i)
      { for(int j = 0; j < 1000; ++j)
        { for(int k = 0; k < 100; ++k)
          { fChIntegral[i][j][k] = 0.05; }
        }
        fIndIntegral[1] = 0.06;
        fIndIntegral[2] = 0.07;
      }

    fEnergy = fTrack->GetE();
    return;
}   

///////////////////////////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////////////////////////
WCSimChargeLikelihood::~WCSimChargeLikelihood()
{
}

///////////////////////////////////////////////////////////////////////////
// Calculate the (log) likelihood of the measured PMT hit
// for a given set of track parameters.  First we have to calculate the
// expected charge given the parameters, then it's just Poisson statistics
///////////////////////////////////////////////////////////////////////////
double WCSimChargeLikelihood::Calc2LnL()
{
  std::cout << "*** WCSimChargeLikelihood::Calc2LnL() *** Calculating the charge log likelihood" << std::endl; 
  Double_t Charge2LnL;
  // std::cout << fDigitArray->GetNDigits() << std::endl;
  for(int iDigit = 0; iDigit < fDigitArray->GetNDigits(); ++iDigit)
    {	
      WCSimLikelihoodDigit * myDigit = (WCSimLikelihoodDigit *)fDigitArray->GetDigit(iDigit);
      fGotTrackParameters = false;
      Double_t chargeObserved = myDigit->GetQ();
      // This is where we need to work out the digit to vertex parameters
      this->GetTrackParameters( myDigit );



  /*    std::cout << "Source position:   (" << sourcePos[0] << "," << sourcePos[1] << "," << sourcePos[2] << ")" << std::endl
                << "Source direction:  (" << TMath::Sin(fTrack->GetTheta()) * TMath::Cos(fTrack->GetPhi()) << "," << TMath::Sin(fTrack->GetTheta()) * TMath::Sin(fTrack->GetPhi()) << "," << TMath::Cos(fTrack->GetTheta()) << ")" << std::endl
                << "Digit position:    (" << myDigit->GetX() << "," << myDigit->GetY() << "," << myDigit->GetZ() << ")" << std::endl
                << "R0 = " << fR0 << std::endl << "CosTheta0 = " << fCosTheta0 << std::endl;
  */
      

      Charge2LnL += -2.0 * TMath::Log(TMath::Poisson( chargeObserved, this->ChargeExpectation() ));
  //    std::cout << iDigit << "   " << Charge2LnL << std::endl;
	 }
  return Charge2LnL;
}


///////////////////////////////////////////////////////////////////////////
// Charge at a given PMT is Poisson distributed.  This calculates the 
// expected mean by summing direct and indirect light contributions.
///////////////////////////////////////////////////////////////////////////
double WCSimChargeLikelihood::ChargeExpectation( WCSimLikelihoodDigit * myDigit)
{
    std::cout << "*** WCSimChargeLikelihood::ChargeExpectation() *** Calculating the total expected mean charge at the PMT" << std::endl; 
    if(!fGotTrackParameters) this->GetTrackParameters(myDigit);
    return this->ChargeExpectation();
}

double WCSimChargeLikelihood::ChargeExpectation()
{
    std::cout << "*** WCSimChargeLikelihood::ChargeExpectation() *** Calculating the total expected mean charge at the PMT" << std::endl; 
    double muDir=0, muIndir=0;
    if(fGotTrackParameters)
    {
      muDir = this->GetMuDirect();
      muIndir = this->GetMuIndirect(); 
    }
    else std::cout << "Error: did not get track parameters first" << std::endl; 
  return 0; //TODO: return muDir + muIndir; ?? /mp
}

///////////////////////////////////////////////////////////////////////////
// Calculate the expected contribution to the Poisson mean from direct
// Cherenkov light
///////////////////////////////////////////////////////////////////////////
double WCSimChargeLikelihood::GetMuDirect( WCSimLikelihoodDigit * myDigit )
{
    std::cout << "*** WCSimChargeLikelihood::GetMuDirect() *** Calculating the direct light contribution to the expected charge" << std::endl; 
    if(!fGotTrackParameters) this->GetTrackParameters( myDigit );
    return this->GetMuDirect(); 
}

double WCSimChargeLikelihood::GetMuDirect( )
{
  std::cout << "*** WCSimChargeLikelihood::GetMuDirect() *** Calculating the direct light contribution to the expected charge" << std::endl; 
    double muDirect = 0.0 ;
  if(fGotTrackParameters)
  {
    double lightFlux = this->GetLightFlux();
    double i0 = this->LookupChIntegrals(0);
    double i1 = this->LookupChIntegrals(1);
    double i2 = this->LookupChIntegrals(2);
    muDirect = lightFlux * ( i0 * fCoeffsCh[0] + i1 * fCoeffsCh[1] + i2 * fCoeffsCh[2] );
  }
  else std::cout << "Error: did not get track parameters first, returning 0" << std::endl; 
  return muDirect;

}

///////////////////////////////////////////////////////////////////////////
// Calculate the expected contribution to the Poisson mean from indirect
// (eg. scatter, reflected) Cherenkov light
///////////////////////////////////////////////////////////////////////////
double WCSimChargeLikelihood::GetMuIndirect( WCSimLikelihoodDigit * myDigit)
{
  if(!fGotTrackParameters) this->GetTrackParameters( myDigit );
  return this->GetMuIndirect();
}


double WCSimChargeLikelihood::GetMuIndirect()
{
  std::cout << "*** WCSimChargeLikelihood::GetMuIndirect() *** Calculating the indirect light contribution to the expected charge" << std::endl; 
  double muIndirect = 0.0;
  if(fGotTrackParameters)
  { 
    double lightFlux = this->GetLightFlux();
    double i1 = this->LookupIndIntegrals(1);
    double i2 = this->LookupIndIntegrals(2);
    muIndirect = lightFlux * ( 1.0 * fCoeffsInd[0] + i1 * fCoeffsInd[1] + i2 * fCoeffsInd[2] );
    std::cout << "Mu indirect = " << muIndirect << std::endl;
  }
  else std::cout << "Error: did not get track parameters first, returning 0" << std::endl;
  return muIndirect;
}

///////////////////////////////////////////////////////////////////////////
// Cherenkov flux is a function of energy and is prefactor to the expected
// number of PE - this calculates it
///////////////////////////////////////////////////////////////////////////
double WCSimChargeLikelihood::GetLightFlux()
{
  std::cout << "*** WCSimChargeLikelihood::GetLightFlux() *** Getting the light flux as a function of the track KE" << std::endl; 
  return 3.0*fEnergy;   // Return 3 PE per MeV. TODO: This isn't right, sort it later.
}

///////////////////////////////////////////////////////////////////////////
// We decouple the direct part into a source factor (flux x profile) and 
// an acceptance factor (solid angle x transmittance x PMT acceptance)
// This tabulates integrals for the acceptance part
///////////////////////////////////////////////////////////////////////////
double WCSimChargeLikelihood::LookupChIntegrals(int i)
{
  std::cout << "*** WCSimChargeLikelihood::LookupChIntegrals() *** Looking up the tabulated integrals for direct Cherenkov light" << std::endl; 
 
  if( i < 0 || i > 2 )
   {
      std::cout << "Error: integral to look up is out of range,   i = " << i << " (should be 0,1,2)" << std::endl;
      return 0; 
   }
   
    //TODO: maybe (trilinear) interpolation? /mp
    int energyBin;
    if( fmod(fEnergy/fEnergyInterval,1) > 0.5 ) energyBin = (int)ceil(fEnergy/fEnergyInterval);
    else energyBin = (int)floor(fEnergy/fEnergyInterval);

    int r0Bin;
    if( fmod(fR0/fR0Interval,1) > 0.5 ) r0Bin = (int)ceil(fR0/fR0Interval);
    else r0Bin = (int)floor(fR0/fR0Interval);

    int cosTheta0Bin;
    if( fmod(fCosTheta0/fCosTheta0Interval,1) > 0.5 ) cosTheta0Bin = (int)ceil(fCosTheta0/fCosTheta0Interval);
    else cosTheta0Bin = (int)floor(fCosTheta0/fCosTheta0Interval);
    std::cout << "Energy bin = " << energyBin << "    CosThetaBin = " << cosTheta0Bin << "     r0Bin = " << r0Bin << std::endl;

    return fChIntegral[energyBin][r0Bin][cosTheta0Bin];
}

///////////////////////////////////////////////////////////////////////////
// We decouple the indirect part into a source factor (flux x profile) and 
// an acceptance factor.  This tabulates integrals for the acceptance part
///////////////////////////////////////////////////////////////////////////
double WCSimChargeLikelihood::LookupIndIntegrals(int i)
{
   std::cout << "*** WCSimChargeLikelihood::LookupIndIntegrals() *** Looking up the tabulated integrals for indirect light" << std::endl; 
   if( i < 1 || i > 2 )
   {
      std::cout << "Error: integral to look up is out of range,   i = " << i << " (should be 1,2)" << std::endl;
      return 0; 
   }
   
    int energyBin;
    //TODO: maybe interpolation? /mp
    if( fmod(fEnergy/fEnergyInterval,1) > 0.5 ) energyBin = (int)ceil(fEnergy/fEnergyInterval);
    else energyBin = (int)floor(fEnergy/fEnergyInterval);
    std::cout << "Energy bin = " << energyBin << std::endl;
    return fIndIntegral[energyBin];
}

///////////////////////////////////////////////////////////////////////////
// Instead of working out another set of integrals for scattered light, we
// note that it must be proportional to the amount of direct light, and
// depend on the position in the tank of its emission, and on its source
///////////////////////////////////////////////////////////////////////////
Double_t WCSimChargeLikelihood::ScatteringTable()
{
  return 1;
}


void WCSimChargeLikelihood::GetTrackParameters( WCSimLikelihoodDigit * myDigit )
{
   Double_t s = 10.0;
  Double_t vtxToSource[3] = { s*TMath::Sin(fTrack->GetTheta()) * TMath::Cos(fTrack->GetPhi()),
                              s*TMath::Sin(fTrack->GetTheta()) * TMath::Sin(fTrack->GetPhi()),
                              s*TMath::Cos(fTrack->GetTheta()) };
  Double_t sourcePos[3]   = { fTrack->GetX() + vtxToSource[0],  // = vertex + distance travelled in particle's direction
                              fTrack->GetY() + vtxToSource[1],
                              fTrack->GetZ() + vtxToSource[2] };

  Double_t r[3]           = { (myDigit->GetX() - sourcePos[0]), // Vector from source to PMT
                              (myDigit->GetY() - sourcePos[1]),
                              (myDigit->GetZ() - sourcePos[2]) };

 
  fR0 = TMath::Sqrt( r[0] * r[0] + r[1] * r[1] + r[2] * r[2] );  // distance from source to PMT
  fCosTheta0 = (r[0]*vtxToSource[0] + r[1]*vtxToSource[1] + r[2]*vtxToSource[2]) / (fR0 * s); // angle between track and vector to PMT
  fEta = TMath::ACos( (myDigit->GetFaceX() * r[0] + myDigit->GetFaceY() * r[1] + myDigit->GetFaceZ() * r[2])/fR0 ); // angle of incidence wrt PMT normal


  fRadius = TMath::Sqrt( sourcePos[0] * sourcePos[0] + sourcePos[1] * sourcePos[1] + sourcePos[2] * sourcePos[2] ); // distance from source to centre of detector
  fTheta = TMath::ACos(fCosTheta0); // as above
  fAngle = TMath::ACos( (sourcePos[0] * myDigit->GetX() + sourcePos[1] * myDigit->GetY() + sourcePos[2] * myDigit->GetZ()) / 
                        (fRadius * TMath::Sqrt(myDigit->GetX() * myDigit->GetX() + myDigit->GetY() * myDigit->GetY() + myDigit->GetZ() * myDigit->GetZ()) )
                      ); // angle at the origin between vector to source and vector to PMT


  Double_t RxPMT[4];  // (vector from origin to source) x (vector from origin to PMT) -- fourth component is modulus
  RxPMT[0] = sourcePos[2]*myDigit->GetZ() - sourcePos[3]*myDigit->GetY();
  RxPMT[1] = sourcePos[2]*myDigit->GetX() - sourcePos[1]*myDigit->GetZ();
  RxPMT[2] = sourcePos[0]*myDigit->GetZ() - sourcePos[1]*myDigit->GetX();
  RxPMT[4] = TMath::Sqrt( RxPMT[0]*RxPMT[0] + RxPMT[1]*RxPMT[1] + RxPMT[2]*RxPMT[2] );

  Double_t RxS[3]; // (vector from origin to source) x (vector source's travel direction) -- fourth component is modulus
  RxS[0] = sourcePos[2]*vtxToSource[2] - sourcePos[3]*vtxToSource[1];
  RxS[1] = sourcePos[2]*vtxToSource[0] - sourcePos[1]*vtxToSource[2];
  RxS[2] = sourcePos[0]*vtxToSource[2] - sourcePos[1]*vtxToSource[0];
  RxS[4] = TMath::Sqrt( RxS[0]*RxS[0] + RxS[1]*RxS[1] + RxS[2]*RxS[2] );
  // These two vectors are normal to the plane contianing the origin, particle and PMT, and to the plane containing the origin and track direction vector
  // (the second should be sourcePos x (sourcePos+S) but sourcePos x sourcePos = 0 )
  // fPhi is the angle between them:  
  fPhi = TMath::ACos( (RxPMT[0]*RxS[0] + RxPMT[1]*RxS[1] + RxPMT[2]*RxS[2]) / (RxS[3]*RxPMT[3]) );

  fGotTrackParameters = true;
  return;

}
