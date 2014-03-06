#include <iostream>
#include <cmath>
#include <sstream>
#include <string>
#include <exception>

#include "WCSimChargeLikelihood.hh"
#include "WCSimDigitizerLikelihood.hh"
#include "WCSimRecoDigit.hh"
#include "WCSimLikelihoodDigit.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimLikelihoodTrack.hh"
#include "WCSimLikelihoodTuner.hh"

#include "TAxis.h"
#include "TCollection.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TROOT.h"
#include "TTree.h"


///////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////
WCSimChargeLikelihood::WCSimChargeLikelihood(WCSimLikelihoodDigitArray * myDigitArray, Bool_t calculateIntegrals)
{
    // std::cout << "*** WCSimChargeLikelihood::WCSimChargeLikelihood() *** Created charge likelihood object" << std::endl;
    fTuner = NULL;
    fDigitArray = NULL;
    fDigit = NULL;
    fDigitizer = NULL;
    this->Initialize( myDigitArray, calculateIntegrals );
}

///////////////////////////////////////////////////////////////////////////
// Set starting values.  Eventually WCSimLikelihoodTuner will be used
// to work these out, but for now they're hardcoded
///////////////////////////////////////////////////////////////////////////
void WCSimChargeLikelihood::Initialize( WCSimLikelihoodDigitArray * myDigitArray, Bool_t calculateIntegrals)
{
    // std::cout << "*** WCSimChargeLikelihood::Initialize() *** Initializing charge likelihood with tuned values" << std::endl;
    fCalculateIntegrals = calculateIntegrals;    // Setting true will calculate integrals numerically instead of looking them up

    fCoeffsCh.push_back(1);
    fCoeffsCh.push_back(2);
    fCoeffsCh.push_back(3);

    fCoeffsInd.push_back(4);
    fCoeffsInd.push_back(5);
    fCoeffsInd.push_back(6);

    fGotTrackParameters = false;
    fDigitArray = myDigitArray;

    fTuner = new WCSimLikelihoodTuner(fDigitArray->GetExtent(0),
                                      fDigitArray->GetExtent(1),
                                      fDigitArray->GetExtent(2),
                                      fCalculateIntegrals);
    
    fDigitizer = new WCSimDigitizerLikelihood(WCSimDigitizerLikelihood::kSimple);

    return;
}

///////////////////////////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////////////////////////
WCSimChargeLikelihood::~WCSimChargeLikelihood()
{
    if(fDigit != NULL) delete fDigit;
    if(fDigitArray != NULL) delete fDigitArray;
    if(fTuner != NULL) delete fTuner;
    if(fDigitizer != NULL) delete fDigitizer;
}

void WCSimChargeLikelihood::ClearTracks()
{
  // std::cout << " *** WCSimChargeLikelihood::ClearTracks() *** " << std::endl;
  // This erases the vector of track pointers.  It doesn't delete the tracks that are pointed to.
  // This should be handled by the process that creates the tracks.
  fTracks.clear();
}

void WCSimChargeLikelihood::SetTracks( std::vector<WCSimLikelihoodTrack*> myTracks )
{
  // std::cout << " *** WCSimChargeLikelihood::SetTracks() *** " << std::endl;
  fTracks = myTracks;
  fGotTrackParameters = false;
}

void WCSimChargeLikelihood::AddTrack( WCSimLikelihoodTrack * myTrack )
{
  // std::cout << " *** WCSimChargeLikelihood::AddTrack() *** " << std::endl;
  fTracks.push_back(myTrack);
  fGotTrackParameters = false;
  return;
}

void WCSimChargeLikelihood::UpdateDigitArray( WCSimLikelihoodDigitArray * myDigitArray)
{
  // std::cout << " *** WCSimChargeLikelihood::UpdateDigitArray() *** " << std::endl;
  fDigitArray = myDigitArray;
  fTuner->UpdateDigitArray(fDigitArray);
  return;
}



///////////////////////////////////////////////////////////////////////////
// Calculate the (log) likelihood of the measured PMT hit
// for a given set of track parameters.  First we have to calculate the
// expected charge given the parameters, then it's just Poisson statistics
///////////////////////////////////////////////////////////////////////////
double WCSimChargeLikelihood::Calc2LnL()
{

//  std::cout << "fDigiPDF = " << fDigiPDF << "   " << fCoeffsCh.at(2) << std::endl;
//  fDigiPDF->Draw("COLZ");
  // std::cout << "*** WCSimChargeLikelihood::Calc2LnL() *** Calculating the charge log likelihood" << std::endl;
  Double_t Charge2LnL=0.0;
  Int_t trackNum = 0;

  Double_t X, Y, Z, Qobs, Qpred;
  std::string str("likelihoodDoubleBins.root");
//  str->Form("likelihood_%f_%f.root", fEnergy, fTrack->GetZ());
  TFile * file = new TFile(str.c_str(),"RECREATE");
  TTree * t = new TTree("likelihood","likelihood");
  t->Branch("X",&X,"X/D");
  t->Branch("Y",&Y,"Y/D");
  t->Branch("Z",&Z,"Z/D");
  t->Branch("Qobs",&Qobs,"Qobs/D");
  t->Branch("Qpred",&Qpred,"Qpred/D");
  t->Branch("trackNum",&trackNum,"TrackNum/I");


  Double_t * predictedCharges = new Double_t[fDigitArray->GetNDigits()];
  for(Int_t i = 0; i < fDigitArray->GetNDigits(); ++ i)
  {
    predictedCharges[i] = 0.0;
  }

  for( std::vector<WCSimLikelihoodTrack *>::iterator trackItr = fTracks.begin(); trackItr != fTracks.end(); ++trackItr)
  {
    trackNum = std::distance(fTracks.begin(), trackItr);


    // Work out the predicted charge at each PMT
    for(int iDigit = 0; iDigit < fDigitArray->GetNDigits(); ++iDigit)
      {
        fDigit = (WCSimLikelihoodDigit *)fDigitArray->GetDigit(iDigit);
        fGotTrackParameters = false;
        //std::cout << "chargeObserved = " << chargeObserved << std::endl;
        // This is where we need to work out the digit to vertex parameters
        this->GetTrackParameters(*trackItr);

        // Work out the predicted mean charge, and the likelihood that this digitizes to the measured charge
        predictedCharges[iDigit] += this->ChargeExpectation(*trackItr);

     }
  }

  for(Int_t jDigit =0; jDigit < fDigitArray->GetNDigits(); ++jDigit)
  {
    fDigit = (WCSimLikelihoodDigit *)fDigitArray->GetDigit(jDigit);
    Double_t chargeObserved = fDigit->GetQ();
    //if((jDigit % 100) == 0){ std::cout << "jDigit  = " << jDigit << "/" << fDigitArray->GetNDigits() << "  Observed charge = " << chargeObserved << "   Expected charge = " << predictedCharges[jDigit] << "    ... So -2lnL now = " << Charge2LnL << std::endl;  }
    Charge2LnL += fDigitizer->GetMinus2LnL( predictedCharges[jDigit], chargeObserved);
    X = fDigit->GetX();
    Y = fDigit->GetY();
    Z = fDigit->GetZ();
    Qobs = fDigit->GetQ();
    Qpred = predictedCharges[jDigit];//this->DigitizerExpectation(chargeExpected);
    t->Fill();

  }
  std::cout << "Final -2LnL = " << Charge2LnL;

  //std::cout << "Writing file to " << str.c_str() << std::endl;
  file->cd();
  t->Write();
  file->Close();
  if(file) delete file;
  return Charge2LnL;
}


///////////////////////////////////////////////////////////////////////////
// Charge at a given PMT is Poisson distributed.  This calculates the
// expected mean by summing direct and indirect light contributions.
///////////////////////////////////////////////////////////////////////////
double WCSimChargeLikelihood::ChargeExpectation(WCSimLikelihoodTrack * myTrack)
{

    // std::cout << "*** WCSimChargeLikelihood::ChargeExpectation() *** Calculating the total expected mean charge at the PMT" << std::endl;
    double muDir=0, muIndir=0;
    if(fGotTrackParameters)
    {
      muDir = this->GetMuDirect(myTrack);
      muIndir = this->GetMuIndirect(myTrack);
 //     if(muDir) std::cout << " Direct charge expectation = " << muDir << " and indirect expectation = " << muIndir << "    total = " << muDir+muIndir << std::endl;
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

double WCSimChargeLikelihood::GetMuDirect(WCSimLikelihoodTrack * myTrack )
{
  //std::cout << "*** WCSimChargeLikelihood::GetMuDirect() *** Calculating the direct light contribution to the expected charge" << std::endl;
    double muDirect = 0.0 ;
  if(!fGotTrackParameters) this->GetTrackParameters(myTrack);
  if(fGotTrackParameters)
  {
    double lightFlux = this->GetLightFlux(myTrack);
//    std::cout << "Energy = " << fEnergy << "  and light flux = " << lightFlux <<std::endl;
    //std::cout << " **** CALCULATING COEFFICIENTS ***** " << std::endl;
    fTuner->CalculateCoefficients(myTrack, fDigit);

    std::vector<Double_t> integrals;
    if(!fCalculateIntegrals)
    {
        integrals = fTuner->LookupChIntegrals(myTrack, fDigit);
    }
    else
    {
        integrals = fTuner->CalculateChIntegrals(myTrack, fDigit);
    }

    //std::cout << "i0 = " << i0 << "  i1 = " << i1 << "   i2 = " << i2 << std::endl
    //          << "fCoeffsCh = " << fCoeffsCh[0] << "," << fCoeffsCh[1] << "," << fCoeffsCh[2] << std::endl;

    muDirect = lightFlux * ( integrals[0] * fCoeffsCh[0] + integrals[1] * fCoeffsCh[1] + integrals[2] * fCoeffsCh[2] );
    if(muDirect < 0)
    {
        //std::cout << "muDir is NEGATIVE! " << "  i0 = " << integrals[0] << "   i1 = " << integrals[1] << "   i2 = " << integrals[2] << "    fCoeffsInd = " << fCoeffsInd[0] << "," << fCoeffsInd[1] << "," << fCoeffsInd[2] << std::endl;
        muDirect = 0;
    }


  }
  else std::cout << "Error: did not get track parameters first, returning 0" << std::endl;
    //std::cout << "muDirect = " << muDirect << std::endl;
  return muDirect;

}

///////////////////////////////////////////////////////////////////////////
// Calculate the expected contribution to the Poisson mean from indirect
// (eg. scatter, reflected) Cherenkov light
///////////////////////////////////////////////////////////////////////////
double WCSimChargeLikelihood::GetMuIndirect(WCSimLikelihoodTrack * myTrack)
{

  //std::cout << "*** WCSimChargeLikelihood::GetMuIndirect() *** Calculating the indirect light contribution to the expected charge" << std::endl;
  double muIndirect = 0.0;
  if(!fGotTrackParameters) this->GetTrackParameters(myTrack);
  if(fGotTrackParameters)
  {
    double lightFlux = this->GetLightFlux(myTrack);
	std::vector<Double_t> integrals;
    if(!fCalculateIntegrals)
    {
		integrals = fTuner->LookupIndIntegrals(myTrack);
    }
    else
    {
		integrals = fTuner->CalculateIndIntegrals(myTrack);
    }
    muIndirect = lightFlux * ( integrals[0] * fCoeffsInd[0] + integrals[1] * fCoeffsInd[1] + integrals[2] * fCoeffsInd[2] );
    if(muIndirect < 0)
    {
        //std::cout << "IT'S NEGATIVE! " << "  i0 = " << integrals[0] << "   i1 = " << integrals[1] << "   i2 = " << integrals[2] << "    fCoeffsInd = " << fCoeffsInd[0] << "," << fCoeffsInd[1] << "," << fCoeffsInd[2] << std::endl;
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
double WCSimChargeLikelihood::GetLightFlux(WCSimLikelihoodTrack * myTrack)
{
  //std::cout << "*** WCSimChargeLikelihood::GetLightFlux() *** Getting the light flux as a function of the track KE" << std::endl;

   Double_t energy = myTrack->GetE();
   Double_t factorFromG = 1/(4.0*TMath::Pi());
  // We normalize g to 1, but it's angular we've swept away a phi term, so need 1/(2pi)

   Double_t factorFromSolidAngle = 1; //2*TMath::Pi();
  // There's an overall normalization of 2pi neglected because the solid angle subtended by a cone is 2pi(1 - cosTheta)
  // I've absorbed this into SolidAngle now

   // Double_t nPhotons = 341.726*fEnergy;
//   Double_t nPhotons = 36.9427* - 3356.71;
  Double_t nPhotons = 17.91 + 39.61*energy;
  if(nPhotons < 0) return 0;
  // From a linear fit to nPhotons as a function of energy (it's very linear)
  //
  // Calculating these is a waste of time because they're so simple, but this is where the factors come from

  return factorFromG * factorFromSolidAngle * nPhotons;
}



void WCSimChargeLikelihood::GetTrackParameters(WCSimLikelihoodTrack * myTrack)
{
  // Get the coefficients to multiply the integrals by
  std::vector<Double_t> coeffs = fTuner->CalculateCoefficientsVector(myTrack, fDigit);

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
    Bool_t tempFlag = fCalculateIntegrals;
    fCalculateIntegrals = true;
    this->Calc2LnL();
    // Charge expectation = integral (flux * rho * solid angle * transmission * acceptance * angular emission * ds)
    fCalculateIntegrals = tempFlag;
    return 0;
}