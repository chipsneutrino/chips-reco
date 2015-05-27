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

  ClearVectors();

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
  ClearVectors();
  fTracks.clear();
}

//////////////////////////////////////////////////////////////////////////
//  Calculate -2 Ln(likelihood) for the two components and combine them
/////////////////////////////////////////////////////////////////////////
Double_t WCSimTotalLikelihood::Calc2LnL()
{
  // nb. at the moment we're just adding these together
  // may have to account for correlations somewhere
  ClearVectors();
  Double_t minus2LnL = 0; 

  for(int iDigit = 0; iDigit < fLikelihoodDigitArray->GetNDigits(); ++iDigit) {
    double digit2LnL = Calc2LnL(iDigit);
    if(digit2LnL < 0) { return -99999; }
    else
    {
		  minus2LnL += Calc2LnL(iDigit);
	  }
  } //for iDigit

  fSetVectors = true;
  std::cout << "-2 ln(Likelihood) = " << minus2LnL << std::endl;
  return minus2LnL;
}

/////////////////////////////////////////////////////////////////////////
//  Calculate the total predicted charge for this digit
/////////////////////////////////////////////////////////////////////////
Double_t WCSimTotalLikelihood::CalcPredictedCharge(unsigned int iDigit)
{
  std::vector<Double_t> predictedCharges = CalcPredictedCharges(iDigit);
  double totalCharge = std::accumulate(predictedCharges.begin(), predictedCharges.end(), 0.0);
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
  // std::cout << "iDigit = " << iDigit << std::endl;
  WCSimLikelihoodDigit *digit = fLikelihoodDigitArray->GetDigit(iDigit);
  // if(digit->GetQ() == 0) { return 0; }
  std::vector<Double_t> predictedCharges = CalcPredictedCharges(iDigit);
  // std::cout << "Size of predicted charge vector = " << predictedCharges.size() << std::endl;
  //  for(unsigned int i = 0; i<predictedCharges.size(); ++i)
  //  {
  //    std::cout << "predictedCharges[" << i << "] = " << predictedCharges.at(i) << std::endl;
  //  }
  double totalCharge = CalcPredictedCharge(iDigit);
  if( TMath::IsNaN( totalCharge ) )
  {
    return -99999;
  }


  double minus2LnL = fDigitizerLikelihood.GetMinus2LnL(totalCharge, digit->GetQ());
  // std::cout << "Recorded charge = " << digit->GetQ() << " and predicted charge = " << totalCharge << " so adds " << minus2LnL << " to -2LnL" << std::endl;

  if(WCSimAnalysisConfig::Instance()->GetUseChargeAndTime())
  {
    minus2LnL += fTimeLikelihood.Calc2LnL(digit, predictedCharges);
  }
  

  fMeasuredCharges.at(iDigit) = digit->GetQ();
  fPredictedCharges.at(iDigit) = totalCharge;
  fTotal2LnL.at(iDigit) = minus2LnL;

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
  ClearVectors();
  return;
}

void WCSimTotalLikelihood::ClearVectors() {
	  fSetVectors = false;
	  fMeasuredCharges.clear();
	  fPredictedCharges.clear();
	  fTotal2LnL.clear();
	  fMeasuredCharges.resize(fLikelihoodDigitArray->GetNDigits());
	  fPredictedCharges.resize(fLikelihoodDigitArray->GetNDigits());
	  fTotal2LnL.resize(fLikelihoodDigitArray->GetNDigits());
}

std::vector<double> WCSimTotalLikelihood::GetMeasuredChargeVector() const {
	if( !fSetVectors )
	{
		std::cerr << "Error: You're asking for the measured charge vector, but this hasn't been set" << std::endl;
		std::cerr << "       Have you reset the tracks or likelihood digit array since calculating the likelihood?" << std::endl;
		assert(fSetVectors);
	}
	return fMeasuredCharges;
}

std::vector<double> WCSimTotalLikelihood::GetPredictedChargeVector() const {
	if( !fSetVectors )
	{
		std::cerr << "Error: You're asking for the predicted charge vector, but this hasn't been set" << std::endl;
		std::cerr << "       Have you reset the tracks or likelihood digit array since calculating the likelihood?" << std::endl;
		assert(fSetVectors);
	}
	return fPredictedCharges;
}

std::vector<double> WCSimTotalLikelihood::GetTotal2LnLVector() const {
	if( !fSetVectors )
	{
		std::cerr << "Error: You're asking for the total 2LnL vector, but this hasn't been set" << std::endl;
		std::cerr << "       Have you reset the tracks or likelihood digit array since calculating the likelihood?" << std::endl;
		assert(fSetVectors);
	}
	return fTotal2LnL;
}
