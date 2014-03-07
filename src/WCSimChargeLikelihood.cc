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
/// Constructor
WCSimChargeLikelihood::WCSimChargeLikelihood(WCSimLikelihoodDigitArray * myDigitArray)
{
    // std::cout << "*** WCSimChargeLikelihood::WCSimChargeLikelihood() *** Created charge likelihood object" << std::endl;
    fTuner = NULL;
    fDigitArray = NULL;
    fDigit = NULL;
    fDigitizer = NULL;
    this->Initialize( myDigitArray );
}



void WCSimChargeLikelihood::Initialize( WCSimLikelihoodDigitArray * myDigitArray)
{
    // std::cout << "*** WCSimChargeLikelihood::Initialize() *** Initializing charge likelihood with tuned values" << std::endl;

    fCoeffsCh.push_back(-999999);
    fCoeffsCh.push_back(-999999);
    fCoeffsCh.push_back(-999999);

    fCoeffsInd.push_back(-999999);
    fCoeffsInd.push_back(-999999);
    fCoeffsInd.push_back(-999999);

    fGotTrackParameters = false;
    fDigitArray = myDigitArray;

    fTuner = new WCSimLikelihoodTuner(fDigitArray);
    
    fDigitizer = new WCSimDigitizerLikelihood();

    return;
}



///////////////////////////////////////////////////////////////////////////
/// Destructor
WCSimChargeLikelihood::~WCSimChargeLikelihood()
{
    std::cout << "Delete fTuner" << std::endl;
    if(fTuner != NULL) delete fTuner;
    std::cout << "Delete fDigitizer" << std::endl;
    if(fDigitizer != NULL) delete fDigitizer;
    std::cout << "Done!" << std::endl;
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
  // fTracks.at(0)->Print();
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



/**
 * Calculate the (log) likelihood of the measured PMT hit
 * for a given set of track parameters.  First we have to calculate the
 * expected charge given the parameters, then it's use the digitiser to
 * work out the likelihood of the signal we measured
 */
double WCSimChargeLikelihood::Calc2LnL()
{

  // std::cout << "fDigiPDF = " << fDigiPDF << "   " << fCoeffsCh.at(2) << std::endl;
  // fDigiPDF->Draw("COLZ");
  // std::cout << "*** WCSimChargeLikelihood::Calc2LnL() *** Calculating the charge log likelihood" << std::endl;
  Double_t Charge2LnL = 0.0;
  Int_t trackNum      = 0;

  // A tree for diagnostic purposes
  Double_t X, Y, Z, Qobs, Qpred, Like;
  std::string str("likelihoodDoubleBins.root");
  // str->Form("likelihood_%f_%f.root", fEnergy, fTrack->GetZ());
  TFile * file = new TFile(str.c_str(),"RECREATE");
  TTree * t = new TTree("likelihood","likelihood");
  t->Branch("X",&X,"X/D");
  t->Branch("Y",&Y,"Y/D");
  t->Branch("Z",&Z,"Z/D");
  t->Branch("Qobs",&Qobs,"Qobs/D");
  t->Branch("Qpred",&Qpred,"Qpred/D");
  t->Branch("minus2LnL",&Like,"Like/D");
  t->Branch("trackNum",&trackNum,"TrackNum/I");


  Double_t * predictedCharges = new Double_t[fDigitArray->GetNDigits()];
  for(Int_t i = 0; i < fDigitArray->GetNDigits(); ++ i)
  {
    predictedCharges[i] = 0.0;
  }

  // Work out how many photons we expect at each PMT from each track
  for( std::vector<WCSimLikelihoodTrack *>::iterator trackItr = fTracks.begin(); trackItr != fTracks.end(); ++trackItr)
  {
    trackNum = std::distance(fTracks.begin(), trackItr);

    // Work out the predicted number of photons at each PMT
    for(int iDigit = 0; iDigit < fDigitArray->GetNDigits(); ++iDigit)
      {
        fDigit = (WCSimLikelihoodDigit *)fDigitArray->GetDigit(iDigit);
        fGotTrackParameters = false;

        // Get integral coefficients
        this->GetTrackParameters(*trackItr);

        // Work out the predicted mean charge at this PMT
        predictedCharges[iDigit] += this->ChargeExpectation(*trackItr);
     }
  }

  // Now work out the probability of the measured charge at each PMT given the prediction
  for(Int_t jDigit =0; jDigit < fDigitArray->GetNDigits(); ++jDigit)
  {
    fDigit                  = (WCSimLikelihoodDigit *)fDigitArray->GetDigit(jDigit);
    Double_t chargeObserved = fDigit->GetQ();
    Double_t minus2LnL      = fDigitizer->GetMinus2LnL( predictedCharges[jDigit], chargeObserved);
    Charge2LnL += minus2LnL;

    // Print diagnotics
    // if((jDigit % 1) == 0){ std::cout << "jDigit  = " << jDigit << "/" << fDigitArray->GetNDigits() << "  Observed charge = " << chargeObserved << "   Expected charge = " << predictedCharges[jDigit] << "    ... So -2lnL now = " << Charge2LnL << std::endl;  }

    // Write the diagnostic tree
    X     = fDigit->GetX();
    Y     = fDigit->GetY();
    Z     = fDigit->GetZ();
    Qobs  = fDigit->GetQ();
    Qpred = predictedCharges[jDigit];//this->DigitizerExpectation(chargeExpected);
    Like  = minus2LnL;
    t->Fill();
  }
  std::cout << "Final -2LnL = " << Charge2LnL;

  // Save diagnostic tree to file
  // std::cout << "Writing file to " << str.c_str() << std::endl;
  file->cd();
  t->Write();
  file->Close();
  if(file) delete file;

  // Clean up pointers
  delete [] predictedCharges;

  return Charge2LnL;
}


/**
 * Charge at a given PMT is sampled from some distribution by the PMT digitiser.
 * This calculates the expected mean by summing direct and indirect light contributions.
*/
double WCSimChargeLikelihood::ChargeExpectation(WCSimLikelihoodTrack * myTrack)
{

    // std::cout << "*** WCSimChargeLikelihood::ChargeExpectation() *** Calculating the total expected mean charge at the PMT" << std::endl;
    double muDir=0, muIndir=0;

    if(fGotTrackParameters)
    {
      muDir   = this->GetMuDirect(myTrack);
      muIndir = this->GetMuIndirect(myTrack);
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
 * Calculate the expected contribution to the Poisson mean from direct
 * Cherenkov light
*/
double WCSimChargeLikelihood::GetMuDirect(WCSimLikelihoodTrack * myTrack )
{
  //std::cout << "*** WCSimChargeLikelihood::GetMuDirect() *** Calculating the direct light contribution to the expected charge" << std::endl;
  double muDirect = 0.0 ;

  // Get the coefficients to multiply the integral
  if(!fGotTrackParameters) { this->GetTrackParameters(myTrack); }

  if(fGotTrackParameters)
  {
    double lightFlux = this->GetLightFlux(myTrack);
    std::vector<Double_t> integrals = fTuner->GetChIntegrals(myTrack, fDigit);

    // std::cout << "Energy = " << fEnergy << "  and light flux = " << lightFlux <<std::endl;
    // std::cout << " **** CALCULATING COEFFICIENTS ***** " << std::endl;
    // fTuner->CalculateCoefficients(myTrack, fDigit);

    // std::cout << "i0 = " << i0 << "  i1 = " << i1 << "   i2 = " << i2 << std::endl
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

  // std::cout << "muDirect = " << muDirect << std::endl;
  return muDirect;
}

/**
 * Calculate the expected contribution to the Poisson mean from indirect
 *  (eg. scattered, reflected) Cherenkov light
*/
double WCSimChargeLikelihood::GetMuIndirect(WCSimLikelihoodTrack * myTrack)
{

  //std::cout << "*** WCSimChargeLikelihood::GetMuIndirect() *** Calculating the indirect light contribution to the expected charge" << std::endl;
  double muIndirect = 0.0;
  if(!fGotTrackParameters) this->GetTrackParameters(myTrack);
  if(fGotTrackParameters)
  {
    double lightFlux = this->GetLightFlux(myTrack);
  	std::vector<Double_t> integrals = fTuner->GetIndIntegrals(myTrack);
    
    muIndirect = lightFlux * ( integrals[0] * fCoeffsInd[0] + integrals[1] * fCoeffsInd[1] + integrals[2] * fCoeffsInd[2] );
    if(muIndirect < 0)
    {
        //std::cout << "IT'S NEGATIVE! " << "  i0 = " << integrals[0] << "   i1 = " << integrals[1] << "   i2 = " << integrals[2] << "    fCoeffsInd = " << fCoeffsInd[0] << "," << fCoeffsInd[1] << "," << fCoeffsInd[2] << std::endl;
        muIndirect = 0;
    }
//    std::cout << "Mu indirect = " << muIndirect << std::endl;
  }
  else
  {
      std::cerr << "Error: GetMuIndirect did not get track parameters first, aborting" << std::endl;
      exit(EXIT_FAILURE);
  }
  return muIndirect;
}

/**
 * Total number of Cherenkov photons emitted by a particle is a function of
 * its energy and particle type.  This function calculates it (just for muons at the moment)
 * @TODO Extend this properly to different particle types
 */
double WCSimChargeLikelihood::GetLightFlux(WCSimLikelihoodTrack * myTrack)
{
  //std::cout << "*** WCSimChargeLikelihood::GetLightFlux() *** Getting the light flux as a function of the track KE" << std::endl;

  Double_t energy      = myTrack->GetE();
  Double_t factorFromG = 1/(4.0*TMath::Pi()); // We want solid angle fraction

  Double_t factorFromSolidAngle = 1; //2*TMath::Pi();
  // There was an overall normalization of 2pi that used to be neglected
  // because the solid angle subtended by a cone is 2pi(1 - cosTheta)
  // I've absorbed this into the SolidAngle now

  // Double_t nPhotons = 341.726*fEnergy;    // <-- EARLIER PARAMETRISATION
  // Double_t nPhotons = 36.9427* - 3356.71; // <-- EARLIER PARAMETRISATION

  // From a linear fit to nPhotons as a function of energy (it's very linear)
  Double_t nPhotons = 17.91 + 39.61*energy;

  // Things get sketchy at really low energies
  if(nPhotons < 0) return 0;

  return factorFromG * factorFromSolidAngle * nPhotons;
}


/**
 * @TODO Check we're not wasting time by calling this too often
 */
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
  // std::cout << "*** WCSimShcargeLikelihood::GetTrackParameters() *** " << std::endl
  //           << "fCoeffsCh = " << fCoeffsCh[0] << "," << fCoeffsCh[1] << "," << fCoeffsCh[2] << std::endl
  //           << "fCoeffsInd = " << fCoeffsInd[0] << "," << fCoeffsInd[1] << "," << fCoeffsInd[2] << std::endl;

  fGotTrackParameters = true;
  return;

}

/**
 *  DEBUGGING:
 * Calculate the likelihood for a given PMT hit Q and a known track, using the
 * values of the individual functions rather than tabulated integrals.
 * We won't use this in the fit, but it's useful for validation
*/
Double_t WCSimChargeLikelihood::CalculateExactLikelihood()
{
    // Remember previous setting
    Bool_t tempFlag = fTuner->GetCalculateIntegrals();

    // Force it to calculate the integrals just this once
    fTuner->SetCalculateIntegrals(true);
    Double_t likelihood = this->Calc2LnL();

    // Restore previous setting
    fTuner->SetCalculateIntegrals(tempFlag);

    // Charge expectation = integral (flux * rho * solid angle * transmission * acceptance * angular emission * ds)
    return likelihood;
}
