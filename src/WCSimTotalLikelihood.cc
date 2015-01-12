#include "WCSimChargePredictor.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimLikelihoodTrack.hh"
#include "WCSimTimeLikelihood.hh"
#include "WCSimTotalLikelihood.hh"
#include "WCSimAnalysisConfig.hh"

#include <iostream>
#include <numeric>
#include <vector>

#ifndef REFLEX_DICTIONARY
ClassImp(WCSimTotalLikelihood)
#endif

///////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////
WCSimTotalLikelihood::WCSimTotalLikelihood( WCSimLikelihoodDigitArray * myLikelihoodDigitArray ) : 
  fLikelihoodDigitArray(myLikelihoodDigitArray),
  fTimeLikelihood(myLikelihoodDigitArray),
	fDigitizerLikelihood()
{  
  fChargeLikelihoodVector.push_back(
      WCSimChargePredictor(fLikelihoodDigitArray) );
}

///////////////////////////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////////////////////////
WCSimTotalLikelihood::~WCSimTotalLikelihood()
{
  fChargeLikelihoodVector.clear();
}

///////////////////////////////////////////////////////////////////////////
//  Add some tracks
///////////////////////////////////////////////////////////////////////////
void WCSimTotalLikelihood::SetTracks( std::vector<WCSimLikelihoodTrack> &myTracks)
{
  this->ResetTracks();
  fTracks = myTracks;

  std::vector<WCSimLikelihoodTrack *> myTrackPtrs;
  std::vector<WCSimLikelihoodTrack>::iterator trackIter;
  //FIXME This is a really ugly solution to the new idea
  //  of one Charge Likelihood object for one track...
  int i = 0;
  for(trackIter = fTracks.begin(); trackIter != fTracks.end(); ++trackIter)
  {
    std::cout << "TrackIter -> " << std::endl;
    myTrackPtrs.push_back(&(*trackIter)); 
    if (i > 0) {
      fChargeLikelihoodVector.push_back(
          WCSimChargePredictor(fLikelihoodDigitArray) );
    }
    fChargeLikelihoodVector[i].AddTrack(&(*trackIter));
    i++;
  }
  fTimeLikelihood.SetTracks(myTrackPtrs);

  return;
}

///////////////////////////////////////////////////////////////////////////
//  Clear current tracks and remove from the charge and time likelihoods
///////////////////////////////////////////////////////////////////////////
void WCSimTotalLikelihood::ResetTracks()
{
  for(unsigned int i = 0; i < fChargeLikelihoodVector.size(); i++) {
    fChargeLikelihoodVector[i].ClearTracks();
  }

  fTimeLikelihood.ClearTracks();
  fTracks.clear();
}

//////////////////////////////////////////////////////////////////////////
//  Calculate -2 Ln(likelihood) for the two components and combine them
/////////////////////////////////////////////////////////////////////////
Double_t WCSimTotalLikelihood::Calc2LnL()
{
  // nb. at the moment we're just adding these together
  // may have to account for correlations somewhere
  Double_t minus2LnL = 0; 

  for(int iDigit = 0; iDigit < fLikelihoodDigitArray->GetNDigits(); ++iDigit) {
    minus2LnL += Calc2LnL(iDigit);
  } //for iDigit

  // std::cout << "-2 ln(Likelihood) = " << minus2LnL << std::endl;
  return minus2LnL;
}

/////////////////////////////////////////////////////////////////////////
//  Calculate the total predicted charge for this digit
/////////////////////////////////////////////////////////////////////////
Double_t WCSimTotalLikelihood::CalcPredictedCharge(unsigned int iDigit)
{
  std::vector<Double_t> predictedCharges = CalcPredictedCharges(iDigit);
  Double_t totalCharge = 0.0;
  for(unsigned int iTrack = 0; iTrack < predictedCharges.size(); ++iTrack)
  {
    totalCharge += predictedCharges.at(iTrack);
  }
  return totalCharge;
}

/////////////////////////////////////////////////////////////////////////
//  Calculate the predicted charge for this digit from each track
/////////////////////////////////////////////////////////////////////////
std::vector<Double_t> WCSimTotalLikelihood::CalcPredictedCharges(unsigned int iDigit)
{
  std::vector<Double_t> predictedCharges;
  for (unsigned int iTrack = 0; iTrack < fTracks.size(); iTrack++) {
    predictedCharges.push_back( fChargeLikelihoodVector.at(iTrack).GetPredictedCharge(iDigit) );
  }
  return predictedCharges;
}

Double_t WCSimTotalLikelihood::Calc2LnL(int iDigit)
{
  WCSimLikelihoodDigit *digit = fLikelihoodDigitArray->GetDigit(iDigit);
  std::vector<Double_t> predictedCharges = CalcPredictedCharges(iDigit);
  double totalCharge = std::accumulate(predictedCharges.begin(), predictedCharges.end(), 0.0);
  double minus2LnL = fDigitizerLikelihood.GetMinus2LnL(totalCharge, digit->GetQ());
  // std::cout << "Recorded charge = " << digit->GetQ() << " and predicted charge = " << totalCharge << " so adds " << minus2LnL << " to -2LnL" << std::endl;

  if(WCSimAnalysisConfig::Instance()->GetUseChargeAndTime())
  {
    minus2LnL += fTimeLikelihood.Calc2LnL(digit, predictedCharges);
  }

  return minus2LnL;
}

void WCSimTotalLikelihood::SetLikelihoodDigitArray(
		WCSimLikelihoodDigitArray* likelihoodDigitArray) {
	fLikelihoodDigitArray = likelihoodDigitArray;

  fChargeLikelihoodVector.clear();
  std::vector<WCSimLikelihoodTrack> tmpTracks = fTracks;
  fChargeLikelihoodVector.push_back(WCSimChargePredictor(fLikelihoodDigitArray));
	fTimeLikelihood = WCSimTimeLikelihood(fLikelihoodDigitArray);
  SetTracks(tmpTracks);
  return;
}
