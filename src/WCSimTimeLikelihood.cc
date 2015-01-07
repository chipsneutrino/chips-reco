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
WCSimTimeLikelihood::WCSimTimeLikelihood(WCSimLikelihoodDigitArray * myDigitArray)
{   
  //std::cout << "*** WCSimTimeLikelihood::WCSimTimeLikelihood() *** Created time likelihood object" << std::endl;
  fDigitArray = NULL;
  fDigit = NULL;
  this->Initialize( myDigitArray);
}



void WCSimTimeLikelihood::Initialize( WCSimLikelihoodDigitArray * myDigitArray)
{
  //std::cout << "*** WCSimTimeLikelihood::Initialize() *** Initializing time likelihood" << std::endl; 

  fGotTrackParameters = -1; //false
  fDigitArray = myDigitArray;

  //defining our likelihood function with 7 parameters
  fLikelihoodFunction = new TF1("fLikelihoodFunction", fFunctionForm, -1e4, 1e4, 7);
  
  //defining our charge dependent parameter function
  fChargeParameterFunction = new TF1("fChargeParameterFunction", "pol4", 0., 1e3);

  //defining our energy dependent parameter
  fEnergyParameterFunction = new TF1("fEnergyParameterFunction", "pol3", 0., 1e5);

  //Here was code doing function normalisation, but got moved to GetTrackParameters

//XXX: debugging
//fDebugFile = new TFile("debug_file.root", "recreate");
//fDebugHist9 = new TH1F("fDebugHist9", "corrected time in bin 9", 800, 0., 1930.);

  return;
}   



///////////////////////////////////////////////////////////////////////////
// Destructor
WCSimTimeLikelihood::~WCSimTimeLikelihood()
{
  //if (fNormFuncParams) delete[] fNormFuncParams;
  delete fLikelihoodFunction;
  delete fChargeParameterFunction;
  delete fEnergyParameterFunction;

//XXX: debug code
//fDebugFile->Close();
//delete fDebugFile;

}



void WCSimTimeLikelihood::ClearTracks()
{
  // std::cout << " *** WCSimTimeLikelihood::ClearTracks() *** " << std::endl;
  // This erases the vector of track pointers.  It doesn't delete the tracks that are pointed to.
  // This should be handled by the process that creates the tracks.
  fTracks.clear();
}

void WCSimTimeLikelihood::SetTracks( std::vector<WCSimLikelihoodTrack*> myTracks )
{
  // std::cout << " *** WCSimTimeLikelihood::SetTracks() *** " << std::endl;
  fTracks = myTracks;
  fGotTrackParameters = -1; //false
}

void WCSimTimeLikelihood::AddTrack( WCSimLikelihoodTrack * myTrack )
{
  // std::cout << " *** WCSimTimeLikelihood::AddTrack() *** " << std::endl;
  fTracks.push_back(myTrack);
  fGotTrackParameters = -1; //false
  return;
}

void WCSimTimeLikelihood::UpdateDigitArray( WCSimLikelihoodDigitArray * myDigitArray)
{
  // std::cout << " *** WCSimTimeLikelihood::UpdateDigitArray() *** " << std::endl;
  fDigitArray = myDigitArray;
  return;
}



/**
 * Calculate the (-2 log) likelihood of the measured hit time
 * at a given PMT for a given charge predictions fitting
 * (previously set) track parameters
 */
Double_t WCSimTimeLikelihood::Calc2LnL(WCSimLikelihoodDigit *myDigit, std::vector<Double_t> predictedCharges )
{
  //std::cout << "*** WCSimTimeLikelihood::Calc2LnL() *** Calculating the time log likelihood" << std::endl; 

  //TODO: somewhere here should be a check of the number of tracks
  //      in the track array, followed by appropriate measures
  //      However, we can now only calculate the likelihood for
  //      one track, so we will only pay attention to the first one

  if (predictedCharges.size() < fTracks.size()){
    std::cout << "Error, not enough charge predictions for all tracks!\n";
    return 0;
  }
  if (fTracks.size() < 1) {
    std::cout << "Error, no tracks in the track vector!\n";
    return 0;
  }

  Double_t time2LnL = 0.;
  //For(int iDigit = 0; iDigit < fDigitArray->GetNDigits(); ++iDigit)
  //{	
    //-2 Log likelihood for individual digit
    //double localLnL = 0;

    fDigit = myDigit;
    Double_t timeObserved = fDigit->GetT();

    //check if there is a hit
    if (timeObserved <= 0.) {
      //if time is 0, there's no hit - read next digit
      return 0;
    }              

    //here apply the time correction for TOF effects
    //FIXME: not just first track (index 0)
    Double_t timeCorrected = this->CorrectedTime(0, timeObserved);

    //FIXME: not just first track (index 0)
    double likelihood = this->TimeLikelihood(0, timeCorrected, predictedCharges);

    if ( TMath::IsNaN(likelihood) ) {
      std::cout << "For digit " << fDigit->GetTubeId() << " likelihood is not a number!\n";
      time2LnL = 2e3;
    }              
    else if (likelihood <= 0.) {
      //std::cout << "For some reason likelihood is not positive\n";
      time2LnL = 2e3;
    }              
    else {         
      time2LnL = -2.0 * TMath::Log( likelihood );
      //if (likelihood > 1.0) std::cout << "For digit " << iDigit << " likelihood is " << likelihood << "\n";
    }              

    //time2LnL += localLnL;
    //std::cout << "Digit: " << iDigit << ", adding: " << localLnL
    //        << ", so time2LnL = " << time2LnL << std::endl;

  //}

//XXX: debug code
//fDebugHist9->Write();

  return time2LnL;
}


///////////////////////////////////////////////////////////////////////////
// TODO: explain what we're doing here in new comment style
///////////////////////////////////////////////////////////////////////////
Double_t WCSimTimeLikelihood::CorrectedTime( Int_t trackIndex, Double_t primaryTime )
{
  //std::cout << "*** WCSimTimeLikelihood::CorrectedTime() *** Calculating time corrected for TOF effects" << std::endl; 

  double correctedTime = 0.;

  if(fGotTrackParameters != trackIndex) this->GetTrackParameters(trackIndex);
  if(fGotTrackParameters == trackIndex) //FIXME: not necessary anymore?
  {
    //TODO: all the track parameters could be initialised at the very beginning
    //vertex positions
    double vtxX = fTracks[trackIndex]->GetX();
    double vtxY = fTracks[trackIndex]->GetY();
    double vtxZ = fTracks[trackIndex]->GetZ();
    
    //track direction and midpoint position
    double vtxPhi = fTracks[trackIndex]->GetPhi();
    double vtxTheta = fTracks[trackIndex]->GetTheta();
    double midpoint = fTrackMidpoint;

    //hopefully got my trigonometry right
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
    double trackTime = fTracks[trackIndex]->GetT();
    correctedTime = primaryTime - tof_corr - timeOffset - trackTime;
  }

  //std::cout << "Primary time: " << primaryTime << ", corrected time: " << correctedTime << std::endl;
  return correctedTime;
}


///////////////////////////////////////////////////////////////////////////
// TODO: explain what we're doing here in new comment style
///////////////////////////////////////////////////////////////////////////
Double_t WCSimTimeLikelihood::TimeLikelihood( Int_t trackIndex, WCSimLikelihoodDigit* myDigit, Double_t correctedTime, std::vector<Double_t> predictedCharges )
{
  fDigit = myDigit;
  return this->TimeLikelihood(trackIndex, correctedTime, predictedCharges);
}


Double_t WCSimTimeLikelihood::TimeLikelihood( Int_t trackIndex, Double_t correctedTime, std::vector<Double_t> predictedCharges )
{
  //std::cout << "*** WCSimTimeLikelihood::TimeLikelihood() *** Calculating the likelihood of the observed hit time" << std::endl; 

  //TODO: somehow handle the fact that we're given a whole vector of tracks
  //      but can currently only compute the likelihood for a single one...

  //FIXME: only gets charge for first track in vector
  double charge = predictedCharges[0];
  //find out in which charge bin does this value lie
  int ibin;
  for (ibin = 0; ibin < fNumChargeCuts; ibin++) {
    if (charge >= fChargeCuts[ibin] && charge < fChargeCuts[ibin+1]) break;
  }
  double chargeBin[1]; chargeBin[0] = ibin;
  //std::cout << "\b, charge bin is: " << ibin << std::endl;

//XXX: debug code
//if (ibin == 8) {
//  fDebugHist9->Fill(correctedTime);
//}

  double timeFuncParams[7];  //each bin has a function with seven parameters
  //the first (norm) parameter is already known after normalisation
  timeFuncParams[0] = fNormFuncParams[ibin];
  //loop from 1 to use the normalised first parameter
  for (int ipar = 1; ipar < 7; ipar++) { //loop over parameters
    //evaluate each one from the fits at the bin center
    timeFuncParams[ipar] = fChargeParameterFunction->EvalPar(chargeBin, fSecondaryParameters[ipar]);
  }

  double x[1]; x[0] = correctedTime;
  double likelihood = fLikelihoodFunction->EvalPar(x, timeFuncParams);

  if ( TMath::IsNaN(likelihood) ) {
    std::cout << "For digit " << fDigit->GetTubeId() << " likelihood is NaN.";
    //std::cout << " Norm paramter is " << timeFuncParams[0] << "\n";
    std::cout << " Predicted charge is " << charge << ", charge bin " << ibin << "\n";
  }
  if ( likelihood > 1.0 ) {
    std::cout << "For digit " << fDigit->GetTubeId() << " likelihood is " << likelihood;
    std::cout << " Norm paramter is " << timeFuncParams[0] << "\n";
  }

  //std::cout << "Returning likelihood: " << likelihood << std::endl;
  return likelihood;
}


///////////////////////////////////////////////////////////////////////////
//// TODO: explain what we're doing here in new comment style
/////////////////////////////////////////////////////////////////////////////
void WCSimTimeLikelihood::GetExternalVariables(const char *fName )
{
  TFile *inFile = new TFile(fName);

  TTree *inTree = (TTree*) inFile->Get("timeParams2");
  if (!inTree) {
    std::cerr << "Problems reading TFile with external variables!\n";
    return;
  }

  inTree->SetBranchAddress("midpointcoeffs", fTrackMidpointCoeffs);
  //inTree->SetBranchAddress("params", fSecondaryParameters);
  inTree->SetBranchAddress("tertiaryParams", fTertiaryParameters);
  inTree->SetBranchAddress("nqcuts", &fSizeChargeCuts);
  inTree->SetBranchAddress("qcuts", fChargeCuts);

  inTree->GetEntry(0); //everything should be in the first entry
  fNumChargeCuts = fSizeChargeCuts - 1;

  //let's see
  //inTree->Show();

  inFile->Close();
  delete inFile;
  //delete inTree;
}



///////////////////////////////////////////////////////////////////////////
// TODO: explain what we're doing here in new comment style
///////////////////////////////////////////////////////////////////////////
void WCSimTimeLikelihood::GetTrackParameters(Int_t trackIndex)
{
  //std::cout << "*** WCSimTimeLikelihood::GetTrackParameters() *** " << std::endl;

  //get the external variables from a file
  //right now it doesn't require track information, but in the future
  //  it might eg. need to know the particle type
  std::string path = getenv("WCSIMANAHOME");
  path.append("/config/timeParams.root");
  this->GetExternalVariables(path.c_str());

  fEnergy = fTracks[trackIndex]->GetE();

  //calculate the track midpoint by evaluating a linear function
  TF1 *midpointFunc = new TF1("midpointFunc", "pol1", 0, 1e5);
  midpointFunc->SetParameters(fTrackMidpointCoeffs);
  fTrackMidpoint = midpointFunc->Eval(fEnergy);
  delete midpointFunc;
  
  //Code moved from Initialize function because it needs ext. variables

  //Evaluate the secondary parameters (from the tertiary ones)
  double x[1]; x[0] = fEnergy;
  for (int i = 0; i < 7; i++) {
    for (int j = 0; j < 5; j++) {
      fSecondaryParameters[i][j] =
        fEnergyParameterFunction->EvalPar(x, fTertiaryParameters[i][j]);
    }
  }

  //Here we should normalise the likelihood function for each bin center
  //  and store the norm factor parameters for later
  //(That way we do the integrals for each charge bin instead of each PMT)
  fNormFuncParams = new Double_t[fNumChargeCuts];
  for (int ibin = 0; ibin < fNumChargeCuts; ibin++) {

    double chargeBin[1]; chargeBin[0] = ibin;
    double timeFuncParams[7];  //each bin has a function with seven parameters
    for (int ipar = 0; ipar < 7; ipar++) { //loop over parameters
      //evaluate each one from the fits at the bin center
      timeFuncParams[ipar] = fChargeParameterFunction->EvalPar(chargeBin, fSecondaryParameters[ipar]);
    }

    //set parameters for current charge bin and integrate the function
    double integral = fLikelihoodFunction->Integral(-1e3, 1e4, timeFuncParams);
    //std::cout << "Integral from " << a << " to " << b << " for bin " << ibin << " equals " << integral << "\n";

    //if integral is nan => alert?
    //if integral > 1. => alert?

    //scale the norm factor accordingly and save it in the array
    fNormFuncParams[ibin] = timeFuncParams[0] / integral;

    timeFuncParams[0] = fNormFuncParams[ibin];
    //integral = fLikelihoodFunction->Integral(a, b, timeFuncParams);
    //std::cout << "After normalising: " << integral << "\n";
  }

  fGotTrackParameters = trackIndex;
  return;
}


///////////////////////////////////////////////////////////////////////////
// TODO: explain what we're doing here in new comment style
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
