#include <iostream>
#include <cmath>
#include <sstream>
#include <string>
#include <exception>

#include "WCSimChargePredictor.hh"
#include "WCSimDigitizerLikelihood.hh"
#include "WCSimEmissionProfiles.hh"
#include "WCSimRecoDigit.hh"
#include "WCSimLikelihoodDigit.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimLikelihoodTrackBase.hh"
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
#include "TStopwatch.h"
#include "TTree.h"


///////////////////////////////////////////////////////////////////////////
/// Constructor
WCSimChargePredictor::WCSimChargePredictor()
{
  fTuner = NULL;
  fDigitArray = NULL;
  fDigit = NULL;
  fNumCalculations = 0;
  fTuner = NULL;
  fEmissionProfiles=  NULL;
}

//
WCSimChargePredictor::WCSimChargePredictor(WCSimLikelihoodDigitArray * myDigitArray, WCSimEmissionProfiles * myProfiles)
{
    // std::cout << "*** WCSimChargeLikelihood::WCSimChargeLikelihood() *** Created charge likelihood object" << std::endl;
    fTuner = NULL;
    fDigitArray = NULL;
    fDigit = NULL;
    this->Initialize( myDigitArray, myProfiles );
    fNumCalculations = 0;
}



void WCSimChargePredictor::Initialize( WCSimLikelihoodDigitArray * myDigitArray, WCSimEmissionProfiles * myEmissionProfiles)
{
    // std::cout << "*** WCSimChargeLikelihood::Initialize() *** Initializing charge likelihood with tuned values" << std::endl;

    fCoeffsCh.push_back(-999999);
    fCoeffsCh.push_back(-999999);
    fCoeffsCh.push_back(-999999);

    fCoeffsInd.push_back(-999999);
    fCoeffsInd.push_back(-999999);
    fCoeffsInd.push_back(-999999);

    fGotTrackParameters = -1; //false
    fDigitArray = myDigitArray;
    fEmissionProfiles = myEmissionProfiles;
    fTuner = new WCSimLikelihoodTuner(fDigitArray, myEmissionProfiles);

    return;
}


///////////////////////////////////////////////////////////////////////////
/// Copy constructor
WCSimChargePredictor::WCSimChargePredictor(const WCSimChargePredictor &other)
{
  fTuner = NULL;
  fDigitArray = NULL;
  fDigit = NULL;
  this->Initialize( other.fDigitArray, other.fEmissionProfiles );
}


///////////////////////////////////////////////////////////////////////////
/// Copy constructor
WCSimChargePredictor& WCSimChargePredictor::operator=(const WCSimChargePredictor &rhs)
{
  fTuner      = rhs.fTuner;;
  fDigitArray = rhs.fDigitArray;;
  fDigit      = rhs.fDigit;;
  this->Initialize( rhs.fDigitArray, rhs.fEmissionProfiles );
  return *this;
}


///////////////////////////////////////////////////////////////////////////
/// Destructor
WCSimChargePredictor::~WCSimChargePredictor()
{
    if(fTuner != NULL)
    {
      std::cout << "Delete fTuner" << std::endl;
      delete fTuner;
      fTuner = NULL;
    }
    std::cout << "Done!" << std::endl;
}



void WCSimChargePredictor::ClearTracks()
{
  // std::cout << " *** WCSimChargeLikelihood::ClearTracks() *** " << std::endl;
  // This erases the vector of track pointers.  It doesn't delete the tracks that are pointed to.
  // This should be handled by the process that creates the tracks.
  fTracks.clear();
}



void WCSimChargePredictor::SetTracks( std::vector<WCSimLikelihoodTrackBase*> myTracks )
{
  // std::cout << " *** WCSimChargeLikelihood::SetTracks() *** " << std::endl;
  fTracks = myTracks;
  // fTracks.at(0)->Print();
  fGotTrackParameters = -1; //false
}

void WCSimChargePredictor::AddTrack( WCSimLikelihoodTrackBase * myTrack )
{
  // std::cout << " *** WCSimChargeLikelihood::AddTrack() *** " << std::endl;
  fTracks.push_back(myTrack);
  fGotTrackParameters = -1; //false
  return;
}

void WCSimChargePredictor::UpdateDigitArray( WCSimLikelihoodDigitArray * myDigitArray)
{
  // std::cout << " *** WCSimChargeLikelihood::UpdateDigitArray() *** " << std::endl;
  fDigitArray = myDigitArray;
  fTuner->UpdateDigitArray(fDigitArray);
  return;
}



/**
 * Calculate the predicted mean number of photons incident on a given PMT
 */
double WCSimChargePredictor::GetPredictedCharge(WCSimLikelihoodDigit *myDigit)
{
  // std::cout << "*** WCSimChargeLikelihood::Calc2LnL() *** Calculating the charge log likelihood" << std::endl;
  unsigned int trackNum;
  Double_t predictedCharge = 0;

  // Work out how many photons we expect at each PMT from each track
  for(unsigned int iTrack = 0; iTrack < fTracks.size(); iTrack++)
  {
    trackNum = iTrack;

    // Work out the predicted number of photons at each PMT
		fDigit = myDigit;
		fGotTrackParameters = -1; //false

		// Get integral coefficients
		this->GetTrackParameters(iTrack);

		// Work out the predicted mean charge at this PMT
		predictedCharge = this->ChargeExpectation(iTrack);
  }

  return predictedCharge;
}

double WCSimChargePredictor::GetPredictedCharge(Int_t digitNum)
{
	WCSimLikelihoodDigit * myDigit = fDigitArray->GetDigit(digitNum);
	return GetPredictedCharge(myDigit);
}


/**
 * Charge at a given PMT is sampled from some distribution by the PMT digitiser.
 * This calculates the expected mean by summing direct and indirect light contributions.
*/
double WCSimChargePredictor::ChargeExpectation(Int_t trackIndex)
{

    // std::cout << "*** WCSimChargeLikelihood::ChargeExpectation() *** Calculating the total expected mean charge at the PMT" << std::endl;
    double muDir=0, muIndir=0;

    if(fGotTrackParameters == trackIndex)
    {
      muDir   = this->GetMuDirect(trackIndex);
      muIndir = this->GetMuIndirect(trackIndex);
      // if(muDir) { std::cout << " Direct charge expectation = " << muDir
      //                       << " and indirect expectation = "  << muIndir
      //                       << "    total = " << muDir+muIndir << std::endl; }
    }
    else
    {
        std::cerr << "Error: did not get track parameters first" << std::endl;
        exit(EXIT_FAILURE);
    }

    // Sometimes the prediction is very slightly negative because of our Taylor expansion, so we'll take zero in these cases:
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

/**
 * This function does the same thing as ChargeExpectation but for a specified PMT.
 * This is necessary for the time likelihood (and must remain public).
*/
double WCSimChargePredictor::ChargeExpectation(Int_t trackIndex, WCSimLikelihoodDigit * myDigit)
{
    // std::cout << "*** WCSimChargeLikelihood::ChargeExpectation() *** Calculating the total expected mean charge at the specified PMT" << std::endl;

  fDigit = myDigit;
  if(fGotTrackParameters != trackIndex) this->GetTrackParameters(trackIndex);
  double charge_expectation =  ChargeExpectation(trackIndex);
  // std::cout << "TubeID = " << fDigit->GetTubeId()
  //           << ", predicted mu = " << charge_expectation
  //           << ", registered charge = " << fDigit->GetQ() << std::endl;

  return charge_expectation;
}


/**
 * Calculate the expected contribution to the Poisson mean from direct
 * Cherenkov light
*/
double WCSimChargePredictor::GetMuDirect(Int_t trackIndex)
{
  //std::cout << "*** WCSimChargeLikelihood::GetMuDirect() *** Calculating the direct light contribution to the expected charge" << std::endl;
  double muDirect = 0.0 ;

  // Get the coefficients to multiply the integral
  if(fGotTrackParameters != trackIndex) { this->GetTrackParameters(trackIndex); }

  if(fGotTrackParameters == trackIndex)
  {
    double lightFlux = this->GetLightFlux(trackIndex);
    std::vector<Double_t> integrals
      = fTuner->GetChIntegrals(fTracks[trackIndex], fDigit);

    // std::cout << "Energy = " << fEnergy << "  and light flux = " << lightFlux <<std::endl;
    // std::cout << " **** CALCULATING COEFFICIENTS ***** " << std::endl;
    // fTuner->CalculateCoefficients(fTracks[trackIndex], fDigit);

//   std::cout << "Digit = " << fDigit->GetTubeId() << "   charge = " << fDigit->GetQ() << "   integrals[0] = " << integrals[0] << "  integrals[1] = " << integrals[1] << "   integrals[2] = " << integrals[2] << std::endl
//           << "fCoeffsCh = " << fCoeffsCh[0] << "," << fCoeffsCh[1] << "," << fCoeffsCh[2] << std::endl;

    muDirect = lightFlux * ( integrals[0] * fCoeffsCh[0] + integrals[1] * fCoeffsCh[1] + integrals[2] * fCoeffsCh[2] );
    if(muDirect < 0)
    {
        // std::cout << "muDir is NEGATIVE! " << "  i0 = " << integrals[0] << "   i1 = " << integrals[1] << "   i2 = " << integrals[2] << "    fCoeffsInd = " << fCoeffsInd[0] << "," << fCoeffsInd[1] << "," << fCoeffsInd[2] << std::endl;
        muDirect = 0;
    }
  }
  else
  {
      std::cerr << "Error: GetMuDirect did not get track parameters first, aborting" << std::endl;
      exit(EXIT_FAILURE);
  }

  // std::cout << "TubeID = " << fDigit->GetTubeId()
  //           << ", muDirect = " << muDirect << std::endl;
  return muDirect;
}

/**
 * Calculate the expected contribution to the Poisson mean from indirect
 *  (eg. scattered, reflected) Cherenkov light
*/
double WCSimChargePredictor::GetMuIndirect(Int_t trackIndex)
{
  return 0.01;
  //std::cout << "*** WCSimChargeLikelihood::GetMuIndirect() *** Calculating the indirect light contribution to the expected charge" << std::endl;
  double muIndirect = 0.0;
  if(fGotTrackParameters != trackIndex) { this->GetTrackParameters(trackIndex); }
  if(fGotTrackParameters == trackIndex) //FIXME: no longer necessary?
  {
    double lightFlux = this->GetLightFlux(trackIndex);
  	std::vector<Double_t> integrals = fTuner->GetIndIntegrals(fTracks[trackIndex]);
    
    muIndirect = lightFlux * ( integrals[0] * fCoeffsInd[0] + integrals[1] * fCoeffsInd[1] + integrals[2] * fCoeffsInd[2] );
    if(muIndirect < 0.01)
    {
        //std::cout << "IT'S NEGATIVE! " << "  i0 = " << integrals[0] << "   i1 = " << integrals[1] << "   i2 = " << integrals[2] << "    fCoeffsInd = " << fCoeffsInd[0] << "," << fCoeffsInd[1] << "," << fCoeffsInd[2] << std::endl;
        muIndirect = 0.01;
    }
  }
  else
  {
      std::cerr << "Error: GetMuIndirect did not get track parameters first, aborting" << std::endl;
      exit(EXIT_FAILURE);
  }

  // std::cout << "TubeID = " << fDigit->GetTubeId()
  //           << " Mu indirect = " << muIndirect << std::endl;
  return muIndirect;
}

/**
 * Total number of Cherenkov photons emitted by a particle is a function of
 * its energy and particle type.  This function calculates it (just for muons at the moment)
 * @TODO Extend this properly to different particle types
 */
double WCSimChargePredictor::GetLightFlux(Int_t trackIndex)
{
  return fTuner->GetLightFlux(fTracks.at(trackIndex));
}


/**
 * @TODO Check we're not wasting time by calling this too often
 */
void WCSimChargePredictor::GetTrackParameters(Int_t trackIndex)
{
  // Get the coefficients to multiply the integrals by
  std::vector<Double_t> coeffs
    = fTuner->CalculateCoefficientsVector(fTracks[trackIndex], fDigit);

  fCoeffsCh[0] = coeffs[0];
  fCoeffsCh[1] = coeffs[1];
  fCoeffsCh[2] = coeffs[2];

  fCoeffsInd[0] = coeffs[3];
  fCoeffsInd[1] = coeffs[4];
  fCoeffsInd[2] = coeffs[5];
  // std::cout << "*** WCSimShcargeLikelihood::GetTrackParameters() *** " << std::endl
  //           << "fCoeffsCh = " << fCoeffsCh[0] << "," << fCoeffsCh[1] << "," << fCoeffsCh[2] << std::endl
  //           << "fCoeffsInd = " << fCoeffsInd[0] << "," << fCoeffsInd[1] << "," << fCoeffsInd[2] << std::endl;

  fGotTrackParameters = trackIndex;
  return;

}

/**
 *  DEBUGGING:
 * Calculate the likelihood for a given PMT hit Q and a known track, using the
 * values of the individual functions rather than tabulated integrals.
 * We won't use this in the fit, but it's useful for validation
*/
//FIXME: this won't end well - with all the integrating for each digit?...
Double_t WCSimChargePredictor::CalculateExactLikelihood(WCSimLikelihoodDigit *myDigit)
{
    // Remember previous setting
    Bool_t tempFlag = fTuner->GetCalculateIntegrals();

    // Force it to calculate the integrals just this once
    fTuner->SetCalculateIntegrals(true);
    Double_t likelihood = this->GetPredictedCharge(myDigit);

    // Restore previous setting
    fTuner->SetCalculateIntegrals(tempFlag);

    // Charge expectation = integral (flux * rho * solid angle * transmission * acceptance * angular emission * ds)
    return likelihood;
}
