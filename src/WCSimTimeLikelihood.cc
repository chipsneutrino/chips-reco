#include "WCSimTimeLikelihood.hh"
#include "WCSimRecoDigit.hh"
#include "WCSimLikelihoodDigit.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimLikelihoodTrack.hh"

#include "TCollection.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"

#include <iostream>
#include <cmath>

///////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////
WCSimTimeLikelihood::WCSimTimeLikelihood(WCSimLikelihoodDigitArray * myDigitArray, WCSimLikelihoodTrack * myTrack, WCSimChargeLikelihood *myChargeLikelihood )
{   
  //std::cout << "*** WCSimTimeLikelihood::WCSimTimeLikelihood() *** Created time likelihood object" << std::endl;
  this->Initialize( myDigitArray, myTrack, myChargeLikelihood  );
}

///////////////////////////////////////////////////////////////////////////
// Set starting values. TODO: Eventually WCSimLikelihoodTuner will be used
// to work these out, but for now they're hardcoded
///////////////////////////////////////////////////////////////////////////
void WCSimTimeLikelihood::Initialize( WCSimLikelihoodDigitArray * myDigitArray, WCSimLikelihoodTrack * myTrack, WCSimChargeLikelihood *myChargeLikelihood )
{
  //std::cout << "*** WCSimTimeLikelihood::Initialize() *** Initializing time likelihood" << std::endl; 

  fTrack = myTrack;
  fGotTrackParameters = false;
  fDigitArray = myDigitArray;
  fChargeLikelihood = myChargeLikelihood;

  //fEnergyInterval = 0.1;    // steps by which to increment variables 
  //fR0Interval = 10;        // in the lookup table from one bin
  //fCosTheta0Interval = 0.1; // to the next

  //for(int i = 0; i < 100; ++i)
  //  { for(int j = 0; j < 1000; ++j)
  //    { for(int k = 0; k < 100; ++k)
  //      { fChIntegral[i][j][k] = 0.05; }
  //    }
  //    fIndIntegral[1] = 0.06;
  //    fIndIntegral[2] = 0.07;
  //  }

  fEnergy = fTrack->GetE();
  
  //get the external variables from a file
  this->GetExternalVariables( "./config/timeParams.root" );

  //defining our likelihood function with 7 parameters for x from 0 to 1e4
  fLikelihoodFunction = new TF1("fLikelihoodFunction", fFunctionForm, -1e4, 1e4, 7);
  
  //defining our charge dependent parameter function for 9 charge bins
  fChargeParameterFunction = new TF1("fChargeParameterFunction", "pol4", 0., 9.);

//XXX: debugging
//fDebugFile = new TFile("debug_file.root", "recreate");
//fDebugHist9 = new TH1F("fDebugHist9", "corrected time in bin 9", 800, 0., 1930.);

  return;
}   

///////////////////////////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////////////////////////
WCSimTimeLikelihood::~WCSimTimeLikelihood()
{
  delete fLikelihoodFunction;
  delete fChargeParameterFunction;

//XXX: debug code
//fDebugFile->Close();
//delete fDebugFile;

}


///////////////////////////////////////////////////////////////////////////
// Calculate the (log) likelihood of the measured PMT hit time
// for a given set of track parameters.  
/////////////////////////////////////////////////////////////////////////////
Double_t WCSimTimeLikelihood::Calc2LnL()
{
  //std::cout << "*** WCSimTimeLikelihood::Calc2LnL() *** Calculating the time log likelihood" << std::endl; 

  Double_t Time2LnL = 0.;
  for(int iDigit = 0; iDigit < fDigitArray->GetNDigits(); ++iDigit)
  {	
    fDigit = (WCSimLikelihoodDigit *)fDigitArray->GetDigit(iDigit);
    fGotTrackParameters = false;
    Double_t timeObserved = fDigit->GetT();
    // This is where we need to work out the digit to vertex parameters
    //this->GetTrackParameters();

    //here apply the time correction for TOF effects
    Double_t timeCorrected = this->CorrectedTime(timeObserved);

    double localLnL = 0;

    double likelihood = this->TimeLikelihood(timeCorrected);
    if ( TMath::IsNaN(likelihood) ) {
      //std::cout << "For some reason likelihood is not a number, so we take 1e-100 instead\n";
      //likelihood = 1e-100;
      localLnL = 2e3;
    }              
    else if (likelihood <= 0.) {
      //std::cout << "For some reason likelihood is not positive, so we take 1e-100 instead\n";
      //likelihood = 1e-100;
      localLnL = 2e3;
    }              
    else {         
      localLnL = -2.0 * TMath::Log( likelihood );
    }              
    //double localLnL = 1. - likelihood ;
    Time2LnL += localLnL;
    //std::cout << "Digit: " << iDigit << ", adding: " << localLnL
    //        << ", so Time2LnL = " << Time2LnL << std::endl;

  }

//XXX: debug code
//fDebugHist9->Write();

  return Time2LnL;
}


///////////////////////////////////////////////////////////////////////////
// TODO: explain what we're doing here
///////////////////////////////////////////////////////////////////////////
Double_t WCSimTimeLikelihood::CorrectedTime( Double_t primaryTime )
{
  //std::cout << "*** WCSimTimeLikelihood::CorrectedTime() *** Calculating time corrected for TOF effects" << std::endl; 
  //if(!fGotTrackParameters) this->GetTrackParameters();

  //TODO: all the track parameters could be initialised at the very beginning

  //vertex positions
  double vtxX = fTrack->GetX();
  double vtxY = fTrack->GetY();
  double vtxZ = fTrack->GetZ();
  
  //track direction and midpoint position
  double vtxPhi = fTrack->GetPhi();
  double vtxTheta = fTrack->GetTheta();
  double midpoint = fTrackMidpoint;

  //hopefully got my trigonometrics right
  double midX = vtxX + (TMath::Sin(vtxTheta)*TMath::Cos(vtxPhi) * midpoint);
  double midY = vtxY + (TMath::Sin(vtxTheta)*TMath::Sin(vtxPhi) * midpoint);
  double midZ = vtxZ + (TMath::Cos(vtxTheta) * midpoint);

  //pmt position
  double pmtX = fDigit->GetX();
  double pmtY = fDigit->GetY();
  double pmtZ = fDigit->GetZ();

  //distance from midpoint to current pmt
  double distX = pmtX - midX;
  double distY = pmtY - midY;
  double distZ = pmtZ - midZ;
  //0.01 to convert from cm to m
  double distance = 0.01* TMath::Sqrt( distX*distX + distY*distY + distZ*distZ );

  //tof correction

  //speed of light
  double n_water = 1.33;  //TODO: should be read from somewhere (also spectral dependent)
  double c = TMath::C();
  double c_water = c/n_water; //speed in m/s

  double light_tof = (1e9* distance / c_water); //1e9 to convert from s to ns

  double mid_distance = 0.01* midpoint; //0.01 to convert from cm to m
  double particle_tof = (1e9* mid_distance / c); //1e9 to convert from s to ns

  //total tof correction
  double tof_corr = light_tof + particle_tof;

  double timeOffset = 950.; //ns

  //total correction includes true track time
  double trackTime = fTrack->GetT();
  double correctedTime = primaryTime - tof_corr - timeOffset - trackTime;

  //std::cout << "Primary time: " << primaryTime << ", corrected time: " << correctedTime << std::endl;
  return correctedTime;
}


///////////////////////////////////////////////////////////////////////////
// TODO: explain what we're doing here
///////////////////////////////////////////////////////////////////////////
Double_t WCSimTimeLikelihood::TimeLikelihood( WCSimLikelihoodDigit* myDigit, Double_t correctedTime )
{
  fDigit = myDigit;
  return this->TimeLikelihood(correctedTime);
}


Double_t WCSimTimeLikelihood::TimeLikelihood( Double_t correctedTime )
{
  //std::cout << "*** WCSimTimeLikelihood::TimeLikelihood() *** Calculating the likelihood of the observed hit time" << std::endl; 
  //if(!fGotTrackParameters) this->GetTrackParameters();

  double charge = this->GetPredictedCharge();
  //std::cout << "Predicted charge is: " << charge << std::endl;
  //if charge is essentially 0, there won't be any hits, so the likelihood's 0 
  if (charge < 1e-10) {
    //std::cout << "Charge is 0, so likelihood too\n";
    return 0;
  }

  int ibin;
  for (ibin = 0; ibin < fNumChargeCuts; ibin++) {
    if (charge >= fChargeCuts[ibin] && charge < fChargeCuts[ibin+1]) break;
  }
  double chargeBin[1];
  chargeBin[0] = ibin + 0.5; //evaluate these parameters at the bin center
  //std::cout << "\b, charge bin is: " << ibin << std::endl;

//XXX: debug code
//if (ibin == 8) {
//  fDebugHist9->Fill(correctedTime);
//}

  double timeFuncParams[7];  //each bin has a function with seven parameters
  for (int ipar = 0; ipar < 7; ipar++) { //loop over parameters
    //evaluate each one from the fits at the bin center
    timeFuncParams[ipar] = fChargeParameterFunction->EvalPar(chargeBin, fSecondaryParameters[ipar]);
  }

  double x[1]; x[0] = correctedTime;
  double likelihood = fLikelihoodFunction->EvalPar(x, timeFuncParams);

  //std::cout << "Returning likelihood: " << likelihood << std::endl;
  return likelihood;
}


///////////////////////////////////////////////////////////////////////////
//// TODO: explain what we're doing here
/////////////////////////////////////////////////////////////////////////////
void WCSimTimeLikelihood::GetExternalVariables( const char *fName )
{
  TFile *inFile = new TFile(fName);

  TTree *inTree = (TTree*) inFile->Get("timeParams");
  if (!inTree) {
    std::cerr << "Problems reading TFile with external variables!\n";
    return;
  }

  double readEnergy;
  inTree->SetBranchAddress("energy", &readEnergy);
  inTree->SetBranchAddress("midpoint", &fTrackMidpoint);
  inTree->SetBranchAddress("params", &fSecondaryParameters);
  inTree->SetBranchAddress("nqcuts", &fSizeChargeCuts);
  inTree->SetBranchAddress("qcuts", &fChargeCuts);
  
  for (int ientry = 0; ientry < inTree->GetEntries(); ientry++) {
    inTree->GetEntry(ientry);
    fNumChargeCuts = fSizeChargeCuts - 1;
    
    //TODO: the energy will be probably binned somehow
    //  the binning should be defined somewhere!
    double binSize = 250; //MeV
    if ( TMath::Abs(readEnergy - fEnergy) < binSize ) break;
  }

  //now we should have the variables for out energy
  //  written in the correct places
  //XXX: or the last value from the loop...

  inFile->Close();
  delete inFile;
  //delete inTree;
}


///////////////////////////////////////////////////////////////////////////
//// TODO: explain what we're doing here
/////////////////////////////////////////////////////////////////////////////
Double_t WCSimTimeLikelihood::GetPredictedCharge()
{
  //TODO: until the charge likelihood part works, we'll take the recorded charge
  double registeredCharge = fDigit->GetQ();
  return registeredCharge;

  //TODO: this is how it should look like!
  double predictedCharge = fChargeLikelihood->ChargeExpectation(fTrack);
  return predictedCharge;
}


///////////////////////////////////////////////////////////////////////////
// TODO: explain what we're doing here
///////////////////////////////////////////////////////////////////////////
void WCSimTimeLikelihood::GetTrackParameters()
{
  fGotTrackParameters = true;
  return;
}


///////////////////////////////////////////////////////////////////////////
// TODO: explain what we're doing here
// The real body of the corrected time likelihood function
// This is a sum of a gaussian peak and a convolution of
//  the exponential curve with a gaussian resolution
///////////////////////////////////////////////////////////////////////////
Double_t WCSimTimeLikelihood::fFunctionForm(double *x, double *p)
{
	// Parameters:
	// 0 - Gaus1 normalisation
	// 1 - Gaus1 mean
	// 2 - Gaus1 sigma
	// 3 - Gaus2Expo normalisation
	// 4 - Gaus2Expo Gaus mean
	// 5 - Gaus2Expo Gaus sigma
	// 6 - Gaus2Expo Expo decay

	// Gaussian
	double gaus1 = TMath::Gaus(x[0],p[1],p[2],true);

	// Expo part of GausExpo convolution
	double ggePart1 = p[3]*TMath::Exp((pow((p[5]/p[6]),2)/2)-((x[0]-p[4])/p[6]));
	double erfArg = (1/TMath::Sqrt(2))*((p[5]/p[6])-((x[0]-p[4])/p[5]));
	double ggePart2 = TMath::Erfc(erfArg);

	return p[0]*(gaus1 + ggePart1*ggePart2);	
}
