#include "WCSimChargeLikelihood.hh"
#include "WCSimRecoDigit.hh"
#include "WCSimLikelihoodDigit.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimLikelihoodTrack.hh"
#include "WCSimLikelihoodTuner.hh"

#include "TCollection.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"

#include <iostream>
#include <cmath>

///////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////
WCSimChargeLikelihood::WCSimChargeLikelihood(WCSimLikelihoodDigitArray * myDigitArray, WCSimLikelihoodTrack * myTrack )
{   
    //std::cout << "*** WCSimChargeLikelihood::WCSimChargeLikelihood() *** Created charge likelihood object" << std::endl;
    this->Initialize( myDigitArray, myTrack );
}

///////////////////////////////////////////////////////////////////////////
// Set starting values.  Eventually WCSimLikelihoodTuner will be used
// to work these out, but for now they're hardcoded
///////////////////////////////////////////////////////////////////////////
void WCSimChargeLikelihood::Initialize( WCSimLikelihoodDigitArray * myDigitArray, WCSimLikelihoodTrack * myTrack )
{
    //std::cout << "*** WCSimChargeLikelihood::Initialize() *** Initializing charge likelihood with tuned values" << std::endl; 
    fCalculateIntegrals = true;    // Setting true will calculate integrals numerically instead of looking them up

    fCoeffsCh.push_back(1);
    fCoeffsCh.push_back(2);
    fCoeffsCh.push_back(3);

    fCoeffsInd.push_back(4);
    fCoeffsInd.push_back(5);
    fCoeffsInd.push_back(6);
 
    fTrack = myTrack;
    fGotTrackParameters = false;
    fDigitArray = myDigitArray;
    fEnergyInterval = 100;     // MeV steps by which to increment variables 
    fR0Interval = 1;           // cm in the lookup table from one bin
    fCosTheta0Interval = 0.02; // to the next

    fEnergy = fTrack->GetE();
    fTuner = new WCSimLikelihoodTuner();
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
  //std::cout << "*** WCSimChargeLikelihood::Calc2LnL() *** Calculating the charge log likelihood" << std::endl; 
  Double_t Charge2LnL=0.0;
  Double_t likelihood;

  Double_t X, Y, Z, Qobs, Qpred;
  TString * str = new TString("likelihood.root");
//  str->Form("likelihood_%f_%f.root", fEnergy, fTrack->GetZ());
  TFile * f = new TFile(str->Data(),"RECREATE");
  TTree * t = new TTree("likelihood","likelihood");
  t->Branch("X",&X,"X/D");
  t->Branch("Y",&Y,"Y/D");
  t->Branch("Z",&Z,"Z/D");
  t->Branch("Qobs",&Qobs,"Qobs/D");
  t->Branch("Qpred",&Qpred,"Qpred/D");


  for(int iDigit = 0; iDigit < fDigitArray->GetNDigits(); ++iDigit)
    {	
      fDigit = (WCSimLikelihoodDigit *)fDigitArray->GetDigit(iDigit);
      fGotTrackParameters = false;
      Double_t chargeObserved = fDigit->GetQ();
      //std::cout << "chargeObserved = " << chargeObserved << std::endl;
      // This is where we need to work out the digit to vertex parameters
      this->GetTrackParameters();



  /*    //std::cout << "Source position:   (" << sourcePos[0] << "," << sourcePos[1] << "," << sourcePos[2] << ")" << std::endl
                << "Source direction:  (" << TMath::Sin(fTrack->GetTheta()) * TMath::Cos(fTrack->GetPhi()) << "," << TMath::Sin(fTrack->GetTheta()) * TMath::Sin(fTrack->GetPhi()) << "," << TMath::Cos(fTrack->GetTheta()) << ")" << std::endl
                << "Digit position:    (" << myDigit->GetX() << "," << myDigit->GetY() << "," << myDigit->GetZ() << ")" << std::endl
                << "R0 = " << fR0 << std::endl << "CosTheta0 = " << fCosTheta0 << std::endl;
  */
      
  //    likelihood *= TMath::Poisson( chargeObserved, this->ChargeExpectation());
      Double_t chargeExpected = this->ChargeExpectation();
      
      Double_t likelihood = TMath::Poisson( chargeObserved, chargeExpected );
      if(likelihood == 0)
      {
        Double_t lnL = (chargeObserved * TMath::Log(chargeExpected) - chargeExpected) - TMath::LnGamma(chargeObserved+1);
                      // This is just log(Poisson) taking n! ~ sqrt((2n + 1/3)pi) * n^n * exp(-n) which is better than Stirling at low n
        Charge2LnL += -2.0 * lnL;
        std::cout << "Adding -2.0 * approx(lnL) = " << -2.0 * lnL << std::endl;
      }
      else
      {
        Charge2LnL += -2.0 * TMath::Log(likelihood);
      }

      if((iDigit % 100) == 0){ std::cout << "iDigit  = " << iDigit << "/" << fDigitArray->GetNDigits() << "  Observed charge = " << chargeObserved << "   Expected charge = " << chargeExpected << "    ... So -2lnL now = " << Charge2LnL << std::endl;  }
      X = fDigit->GetX();
      Y = fDigit->GetY();
      Z = fDigit->GetZ();
      Qobs = chargeObserved;
      Qpred = chargeExpected;
      t->Fill();
	 }
     t->Write();
     f->Close();
  return Charge2LnL;
}


///////////////////////////////////////////////////////////////////////////
// Charge at a given PMT is Poisson distributed.  This calculates the 
// expected mean by summing direct and indirect light contributions.
///////////////////////////////////////////////////////////////////////////
double WCSimChargeLikelihood::ChargeExpectation()
{

    //std::cout << "*** WCSimChargeLikelihood::ChargeExpectation() *** Calculating the total expected mean charge at the PMT" << std::endl; 
    double muDir=0, muIndir=0;
    if(fGotTrackParameters)
    {
      muDir = this->GetMuDirect();
      muIndir = this->GetMuIndirect(); 
      //std::cout << " Direct charge expectation = " << muDir << " and indirect expectation = " << muIndir << "    total = " << muDir+muIndir << std::endl;
    }
    else{std::cout << "Error: did not get track parameters first" << std::endl; }

    // Sometimes the integral is very slightly negative, we'll return zero in these cases:
    if( muDir < 0 )
    {
      std::cerr << "Warning: predicted direct charge was less than zero (" << muDir <<") - returning zero" << std::endl;
      muDir = 0;
    }
    if( muIndir < 0 )
    {
      std::cerr << "Warning: predicted indirect charge was less than zero (" << muIndir <<") - returning zero" << std::endl;
      muIndir = 0;
    }

  return (muDir+muIndir);
}

///////////////////////////////////////////////////////////////////////////
// Calculate the expected contribution to the Poisson mean from direct
// Cherenkov light
///////////////////////////////////////////////////////////////////////////

double WCSimChargeLikelihood::GetMuDirect( )
{
  //std::cout << "*** WCSimChargeLikelihood::GetMuDirect() *** Calculating the direct light contribution to the expected charge" << std::endl; 
    double muDirect = 0.0 ;
  if(!fGotTrackParameters) this->GetTrackParameters();
  if(fGotTrackParameters)
  {
    double lightFlux = this->GetLightFlux();
//    std::cout << "Energy = " << fEnergy << "  and light flux = " << lightFlux <<std::endl;
    //std::cout << " **** CALCULATING COEFFICIENTS ***** " << std::endl;
    fTuner->CalculateCoefficients(fTrack, fDigit);

    double i0, i1, i2;
    if(!fCalculateIntegrals)
    {
        i0 = this->LookupChIntegrals(0);
        i1 = this->LookupChIntegrals(1);
        i2 = this->LookupChIntegrals(2);
    }
    else
    {
        i0 = fTuner->CalculateChIntegrals(0, fTrack, fDigit);
        i1 = fTuner->CalculateChIntegrals(1, fTrack, fDigit);
        i2 = fTuner->CalculateChIntegrals(2, fTrack, fDigit);
    }

    //std::cout << "i0 = " << i0 << "  i1 = " << i1 << "   i2 = " << i2 << std::endl
    //          << "fCoeffsCh = " << fCoeffsCh[0] << "," << fCoeffsCh[1] << "," << fCoeffsCh[2] << std::endl;
    
    muDirect = lightFlux * ( i0 * fCoeffsCh[0] + i1 * fCoeffsCh[1] + i2 * fCoeffsCh[2] );
  }
  else std::cout << "Error: did not get track parameters first, returning 0" << std::endl; 
    //std::cout << "muDirect = " << muDirect << std::endl;
  return muDirect;

}

///////////////////////////////////////////////////////////////////////////
// Calculate the expected contribution to the Poisson mean from indirect
// (eg. scatter, reflected) Cherenkov light
///////////////////////////////////////////////////////////////////////////
double WCSimChargeLikelihood::GetMuIndirect()
{
  //std::cout << "*** WCSimChargeLikelihood::GetMuIndirect() *** Calculating the indirect light contribution to the expected charge" << std::endl; 
  double muIndirect = 0.0;
  if(!fGotTrackParameters) this->GetTrackParameters(  );
  if(fGotTrackParameters)
  { 
    double lightFlux = this->GetLightFlux();
    double i1;
    double i2;
    if(!fCalculateIntegrals)
    {
        i1 = this->LookupIndIntegrals(1);
        i2 = this->LookupIndIntegrals(2);
    }
    else
    {   
        i1 = fTuner->CalculateIndIntegrals(1, fTrack);
        i2 = fTuner->CalculateIndIntegrals(2, fTrack);
    }
    muIndirect = lightFlux * ( 1.0 * fCoeffsInd[0] + i1 * fCoeffsInd[1] + i2 * fCoeffsInd[2] );
    if(muIndirect < 0)
    {
        std::cout << "IT'S NEGATIVE! " << "  i1 = " << i1 << "   i2 = " << i2 << "    fCoeffsInd = " << fCoeffsInd[0] << "," << fCoeffsInd[1] << "," << fCoeffsInd[2] << std::endl;
        muIndirect = 0;        
    }
//    std::cout << "Mu indirect = " << muIndirect << std::endl;
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
  //std::cout << "*** WCSimChargeLikelihood::GetLightFlux() *** Getting the light flux as a function of the track KE" << std::endl; 
  
  
   Double_t factorFromG = 1/(4.0*TMath::Pi());
  // We normalize g to 1, but it's angular so is missing 1/(4pi)
  // I've absorbed this into g now
  
   Double_t factorFromSolidAngle = 2*TMath::Pi(); 
  // There's an overall normalization of 2pi neglected because the solid angle subtended by a cone is 2pi(1 - cosTheta)
  // I've absorbed this into SolidAngle now

   Double_t nPhotons = 341.726*fEnergy;    
  // From a linear fit to nPhotons as a function of energy (it's very linear)
  // 
  // Calculating these is a waste of time because they're so simple, but this is where the factors come from
  
  return factorFromG * factorFromSolidAngle * nPhotons;
  //return 341.726*fEnergy;///(2 * TMath::Pi());
}

///////////////////////////////////////////////////////////////////////////
// We decouple the direct part into a source factor (flux x profile) and 
// an acceptance factor (solid angle x transmittance x PMT acceptance)
// This tabulates integrals for the acceptance part
///////////////////////////////////////////////////////////////////////////
double WCSimChargeLikelihood::LookupChIntegrals(int i)
{
  //std::cout << "*** WCSimChargeLikelihood::LookupChIntegrals() *** Looking up the tabulated integrals for direct Cherenkov light" << std::endl; 
 
  if( i < 0 || i > 2 )
   {
      //std::cout << "Error: integral to look up is out of range,   i = " << i << " (should be 0,1,2)" << std::endl;
      return 0; 
   }
   
    int energyBin;
    if( fmod(fEnergy/fEnergyInterval,1) > 0.5 ) energyBin = (int)ceil(fEnergy/fEnergyInterval);
    else energyBin = (int)floor(fEnergy/fEnergyInterval);

    int r0Bin;
    if( fmod(fR0/fR0Interval,1) > 0.5 ) r0Bin = (int)ceil(fR0/fR0Interval);
    else r0Bin = (int)floor(fR0/fR0Interval);

    int cosTheta0Bin;
    if( fmod(fCosTheta0/fCosTheta0Interval,1) > 0.5 ) cosTheta0Bin = (int)ceil(fCosTheta0/fCosTheta0Interval);
    else cosTheta0Bin = (int)floor(fCosTheta0/fCosTheta0Interval);
    //std::cout << "Energy bin = " << energyBin << "    CosThetaBin = " << cosTheta0Bin << "     r0Bin = " << r0Bin << std::endl;
    return 1;
    return fChIntegral[energyBin][r0Bin][cosTheta0Bin];
}

///////////////////////////////////////////////////////////////////////////
// We decouple the indirect part into a source factor (flux x profile) and 
// an acceptance factor.  This tabulates integrals for the acceptance part
///////////////////////////////////////////////////////////////////////////
double WCSimChargeLikelihood::LookupIndIntegrals(int i)
{
   //std::cout << "*** WCSimChargeLikelihood::LookupIndIntegrals() *** Looking up the tabulated integrals for indirect light" << std::endl; 
   if( i < 1 || i > 2 )
   {
      //std::cout << "Error: integral to look up is out of range,   i = " << i << " (should be 1,2)" << std::endl;
      return 0; 
   }
   
    int energyBin;
    if( fmod(fEnergy/fEnergyInterval,1) > 0.5 ) energyBin = (int)ceil(fEnergy/fEnergyInterval);
    else energyBin = (int)floor(fEnergy/fEnergyInterval);
    //std::cout << "Energy bin = " << energyBin << std::endl;
    
    return fIndIntegral[energyBin];

}

///////////////////////////////////////////////////////////////////////////
// Instead of working out another set of integrals for scattered light, we
// note that it must be proportional to the amount of direct light, and
// depend on the position in the tank of its emission, and on its source
///////////////////////////////////////////////////////////////////////////
Double_t WCSimChargeLikelihood::ScatteringTable()
{
  // THIS ISN'T USED YET - THERE'S AN EQUIVALENT IN WCSIMLIKELIHOODTUNER
  // We'll turn this off for the time being...
  return 0.05;
}


void WCSimChargeLikelihood::GetTrackParameters( )
{
  // Get the coefficients to multiply the integrals by
  std::vector<Double_t> coeffs = fTuner->CalculateCoefficientsVector(fTrack, fDigit);
  
  fCoeffsCh[0] = coeffs[0];
  fCoeffsCh[1] = coeffs[1];
  fCoeffsCh[2] = coeffs[2];

  fCoeffsInd[0] = coeffs[3];
  fCoeffsInd[1] = coeffs[4];
  fCoeffsInd[2] = coeffs[5];
//  std::cout << "*** WCSimShcargeLikelihood::GetTrackParameters() *** " << std::endl
//            << "fCoeffsCh = " << fCoeffsCh[0] << "," << fCoeffsCh[1] << "," << fCoeffsCh[2] << std::endl
//            << "fCoeffsInd = " << fCoeffsInd[0] << "," << fCoeffsInd[1] << "," << fCoeffsInd[2] << std::endl;

  fGotTrackParameters = true;
  return;

}

//////////////////////////////////////////////////////////////////////////////
// Calculate the likelihood for a given PMT hit Q and a known track, using the 
// values of the individual functions rather than tabulated integrals.  
// We won't use this in the fit, but it's useful for validation
//////////////////////////////////////////////////////////////////////////////
Double_t WCSimChargeLikelihood::CalculateExactLikelihood()
{
    if(!fGotTrackParameters){ this->GetTrackParameters(); }
    fCalculateIntegrals = true;
    this->Calc2LnL();
    // Charge expectation = integral (flux * rho * solid angle * transmission * acceptance * angular emission * ds)
    return 0; 

}

Double_t WCSimChargeLikelihood::GetR0Interval()
{
  return fR0Interval;
}

Double_t WCSimChargeLikelihood::GetCosTheta0Interval()
{
  return fCosTheta0Interval;
}


