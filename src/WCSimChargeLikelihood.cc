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
#include "TStopwatch.h"
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
    fNumCalculations = 0;
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

    fGotTrackParameters = -1; //false
    fDigitArray = myDigitArray;

    fTuner = new WCSimLikelihoodTuner(fDigitArray);
    fDigitizer = new WCSimDigitizerLikelihood();

    return;
}


///////////////////////////////////////////////////////////////////////////
/// Copy constructor
WCSimChargeLikelihood::WCSimChargeLikelihood(const WCSimChargeLikelihood &other)
{
  fTuner = NULL;
  fDigitArray = NULL;
  fDigit = NULL;
  fDigitizer = NULL;
  this->Initialize( other.fDigitArray );
}


///////////////////////////////////////////////////////////////////////////
/// Copy constructor
WCSimChargeLikelihood& WCSimChargeLikelihood::operator=(const WCSimChargeLikelihood &rhs)
{
  fTuner      = rhs.fTuner;;
  fDigitArray = rhs.fDigitArray;;
  fDigit      = rhs.fDigit;;
  fDigitizer  = rhs.fDigitizer;;
  this->Initialize( rhs.fDigitArray );
  return *this;
}


///////////////////////////////////////////////////////////////////////////
/// Destructor
WCSimChargeLikelihood::~WCSimChargeLikelihood()
{
    if(fTuner != NULL)
    {
      std::cout << "Delete fTuner" << std::endl;
      delete fTuner;
      fTuner = NULL;
    }
    if(fDigitizer != NULL)
    {
      std::cout << "Delete fDigitizer" << std::endl;
      delete fDigitizer;
      fDigitizer = NULL;
    }
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
  fGotTrackParameters = -1; //false
}

void WCSimChargeLikelihood::AddTrack( WCSimLikelihoodTrack * myTrack )
{
  // std::cout << " *** WCSimChargeLikelihood::AddTrack() *** " << std::endl;
  fTracks.push_back(myTrack);
  fGotTrackParameters = -1; //false
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
double WCSimChargeLikelihood::Calc2LnL(WCSimLikelihoodDigit *myDigit)
{
  // std::cout << "*** WCSimChargeLikelihood::Calc2LnL() *** Calculating the charge log likelihood" << std::endl;
  Double_t Charge2LnL = 0.0;
  Int_t trackNum      = 0;

  // A tree for diagnostic purposes
  //Double_t X, Y, Z, Qobs, Qpred, Like;
  //std::string str("likelihoodDoubleBins.root");
  // str->Form("likelihood_%f_%f.root", fEnergy, fTrack->GetZ());
  //TFile * file = new TFile(str.c_str(),"RECREATE");
  //TTree * t = new TTree("likelihood","likelihood");
  //t->Branch("X",&X,"X/D");
  //t->Branch("Y",&Y,"Y/D");
  //t->Branch("Z",&Z,"Z/D");
  //t->Branch("Qobs",&Qobs,"Qobs/D");
  //t->Branch("Qpred",&Qpred,"Qpred/D");
  //t->Branch("minus2LnL",&Like,"Like/D");
  //t->Branch("trackNum",&trackNum,"TrackNum/I");


  //Double_t * predictedCharges = new Double_t[fDigitArray->GetNDigits()];
  //for(Int_t i = 0; i < fDigitArray->GetNDigits(); ++ i)
  //{
    //predictedCharges[i] = 0.0;
  //}
  Double_t predictedCharge = 0;

  // Work out how many photons we expect at each PMT from each track
  for(unsigned int iTrack = 0; iTrack < fTracks.size(); iTrack++)
  {
    trackNum = iTrack;

    // Work out the predicted number of photons at each PMT
    //for(int iDigit = 0; iDigit < fDigitArray->GetNDigits(); ++iDigit)
      //{
        //fDigit = (WCSimLikelihoodDigit *)fDigitArray->GetDigit(iDigit);
        fDigit = myDigit;
        fGotTrackParameters = -1; //false

        // Get integral coefficients
        this->GetTrackParameters(iTrack);

        // Work out the predicted mean charge at this PMT
        //predictedCharges[iDigit] += this->ChargeExpectation(iTrack);
        predictedCharge = this->ChargeExpectation(iTrack);
     //}
  }

  // Now work out the probability of the measured charge at each PMT given the prediction
  //for(Int_t jDigit =0; jDigit < fDigitArray->GetNDigits(); ++jDigit)
  //{
    //fDigit                  = (WCSimLikelihoodDigit *)fDigitArray->GetDigit(jDigit);
    Double_t chargeObserved = fDigit->GetQ();
    //Double_t minus2LnL      = fDigitizer->GetMinus2LnL( predictedCharges[jDigit], chargeObserved);
    Double_t minus2LnL = fDigitizer->GetMinus2LnL( predictedCharge, chargeObserved);
    Charge2LnL += minus2LnL;

    // Print diagnotics
    // { std::cout << "Digit  = " << fDigit->GetTubeId() << "/" << fDigitArray->GetNDigits() << "  Observed charge = " << chargeObserved << "   Expected charge = " << predictedCharge << "    ... So -2lnL now = " << Charge2LnL << std::endl;  }

    //// Write the diagnostic tree
    //X     = fDigit->GetX();
    //Y     = fDigit->GetY();
    //Z     = fDigit->GetZ();
    //Qobs  = fDigit->GetQ();
    //Qpred = predictedCharges[jDigit];//this->DigitizerExpectation(chargeExpected);
    //Like  = minus2LnL;
    //t->Fill();
  //}
  // std::cout << "Final -2LnL = " << Charge2LnL;

  // Save diagnostic tree to file
  // std::cout << "Writing file to " << str.c_str() << std::endl;
  //file->cd();
  //t->Write();
  //file->Close();
  //if(file) delete file;

  // Clean up pointers
  //delete [] predictedCharges;

  return Charge2LnL;
}


/**
 * Charge at a given PMT is sampled from some distribution by the PMT digitiser.
 * This calculates the expected mean by summing direct and indirect light contributions.
*/
double WCSimChargeLikelihood::ChargeExpectation(Int_t trackIndex)
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
double WCSimChargeLikelihood::ChargeExpectation(Int_t trackIndex, WCSimLikelihoodDigit * myDigit)
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
double WCSimChargeLikelihood::GetMuDirect(Int_t trackIndex)
{
  //std::cout << "*** WCSimChargeLikelihood::GetMuDirect() *** Calculating the direct light contribution to the expected charge" << std::endl;
  double muDirect = 0.0 ;

  // Get the coefficients to multiply the integral
  if(fGotTrackParameters != trackIndex) { this->GetTrackParameters(trackIndex); }

  if(fGotTrackParameters == trackIndex) //FIXME: no longer necessary?
  {
    double lightFlux = this->GetLightFlux(trackIndex);
    std::vector<Double_t> integrals
      = fTuner->GetChIntegrals(fTracks[trackIndex], fDigit);

    // std::cout << "Energy = " << fEnergy << "  and light flux = " << lightFlux <<std::endl;
    // std::cout << " **** CALCULATING COEFFICIENTS ***** " << std::endl;
    // fTuner->CalculateCoefficients(fTracks[trackIndex], fDigit);

    // std::cout << "integrals[0] = " << integrals[0] << "  integrals[1] = " << integrals[1] << "   integrals[2] = " << integrals[2] << std::endl
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
double WCSimChargeLikelihood::GetMuIndirect(Int_t trackIndex)
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
double WCSimChargeLikelihood::GetLightFlux(Int_t trackIndex)
{
  return fTuner->GetLightFlux(fTracks.at(trackIndex));
}


/**
 * @TODO Check we're not wasting time by calling this too often
 */
void WCSimChargeLikelihood::GetTrackParameters(Int_t trackIndex)
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
Double_t WCSimChargeLikelihood::CalculateExactLikelihood(WCSimLikelihoodDigit *myDigit)
{
    // Remember previous setting
    Bool_t tempFlag = fTuner->GetCalculateIntegrals();

    // Force it to calculate the integrals just this once
    fTuner->SetCalculateIntegrals(true);
    Double_t likelihood = this->Calc2LnL(myDigit);

    // Restore previous setting
    fTuner->SetCalculateIntegrals(tempFlag);

    // Charge expectation = integral (flux * rho * solid angle * transmission * acceptance * angular emission * ds)
    return likelihood;
}
