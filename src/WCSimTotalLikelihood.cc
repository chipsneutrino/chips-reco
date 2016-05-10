#include "WCSimChargePredictor.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimLikelihoodTrack.hh"
#include "WCSimTimeLikelihood3.hh"
#include "WCSimTotalLikelihood.hh"
#include "WCSimAnalysisConfig.hh"
#include "TMath.h"

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
	fDigitizerLikelihood()
{  
  fEmissionProfileManager = new WCSimEmissionProfileManager();
//  fChargeLikelihoodVector.push_back(
//      WCSimChargePredictor(fLikelihoodDigitArray, fEmissionProfileManager) );
  fTimeLikelihood = new WCSimTimeLikelihood3(fLikelihoodDigitArray, fEmissionProfileManager);

  ClearVectors();

}

///////////////////////////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////////////////////////
WCSimTotalLikelihood::~WCSimTotalLikelihood()
{
  fChargeLikelihoodVector.clear();
  delete fEmissionProfileManager;
  delete fTimeLikelihood;
}

///////////////////////////////////////////////////////////////////////////
//  Add some tracks
///////////////////////////////////////////////////////////////////////////
void WCSimTotalLikelihood::SetTracks( std::vector<WCSimLikelihoodTrackBase*> &myTracks)
{
  this->ResetTracks();
  fTracks = myTracks;
  fChargeLikelihoodVector.reserve(fTracks.size());

  std::vector<WCSimLikelihoodTrackBase*> myTrackPtrs;
  std::vector<WCSimLikelihoodTrackBase*>::iterator trackIter;
  //FIXME This is a really ugly solution to the new idea
  //  of one Charge Likelihood object for one track...
  int i = 0;
  for(trackIter = fTracks.begin(); trackIter != fTracks.end(); ++trackIter)
  {
    // std::cout << "TrackIter -> " << std::endl;
//    if (i > 0) {
    fChargeLikelihoodVector.push_back(WCSimChargePredictor(fLikelihoodDigitArray, fEmissionProfileManager));
//    }
    fChargeLikelihoodVector[i].AddTrack(*trackIter);
    i++;
  }
  fTimeLikelihood->SetTracks(fTracks);
  fEmissionProfileManager->CacheAtLeast(fTracks.size());

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
  fChargeLikelihoodVector.clear();

  fTimeLikelihood->ClearTracks();
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
    if(TMath::IsNaN(digit2LnL))
    {
      std::cerr << iDigit << " was NaN in Calc2LnL!!" << std::endl;
      return -99999;
    }
    if(digit2LnL < 0) { 
        return -99999; }
    else
    {
	  minus2LnL += digit2LnL;
	}
  } //for iDigit

  double chargeLnL = 0.0;
  double timeLnL = 0.0;
  for(int iDigit = 0; iDigit < fLikelihoodDigitArray->GetNDigits(); ++iDigit)
  {
    if(fCharge2LnL.at(iDigit) > 0) { chargeLnL += fCharge2LnL.at(iDigit); }
    if(fTime2LnL.at(iDigit) > 0) { timeLnL += fTime2LnL.at(iDigit); }
  }
  std::cout << "Time component = " << timeLnL << " and charge component = " << chargeLnL << " so total = " << minus2LnL << std::endl;
  fSetVectors = true;
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

  WCSimLikelihoodDigit *digit = fLikelihoodDigitArray->GetDigit(iDigit);
  double minus2LnL = 0.0;
  double totalCharge = -999.9;
    // if(digit->GetQ() == 0) { return 0; }
  // std::cout << "Size of predicted charge vector = " << predictedCharges.size() << std::endl;
  //  for(unsigned int i = 0; i<predictedCharges.size(); ++i)
  //  {
  //    std::cout << "predictedCharges[" << i << "] = " << predictedCharges.at(i) << std::endl;
  //  }
  double timePart = 0, chargePart = 0;

  if( WCSimAnalysisConfig::Instance()->GetUseCharge())
  {
	  totalCharge = CalcPredictedCharge(iDigit);
	  if( TMath::IsNaN( totalCharge ) )
	  {
		return -99999;
	  }
	  chargePart = fDigitizerLikelihood.GetMinus2LnL(totalCharge, digit->GetQ());
	  minus2LnL += chargePart;
	  fPredictedCharges.at(iDigit) = totalCharge;
	  fCharge2LnL.at(iDigit) = chargePart;
  }
  else
  {
	  fPredictedCharges.at(iDigit) = -999.9;
	  fCharge2LnL.at(iDigit) = -999.9;
  }

  if( WCSimAnalysisConfig::Instance()->GetUseTime())
  {
	  timePart = 0.0;
	  timePart = fTimeLikelihood->Calc2LnL(iDigit);
	  minus2LnL += timePart;
	  fTime2LnL.at(iDigit) = timePart;
  }
  else
  {
	  fTime2LnL.at(iDigit) = -999.9;
  }

  // if( (digit->GetTubeId() == 1500) && digit->GetQ() > 10)
  // {
     //std::cout << "Digit " << digit->GetTubeId() << " has recorded charge = " << digit->GetQ() << " and predicted charge before effiencies = " << totalCharge << " so charge adds " << chargePart << " to -2LnL and time adds " << timePart << " to minus2LnL " << minus2LnL << std::endl;
  // }

  fMeasuredCharges.at(iDigit) = digit->GetQ();
  fTotal2LnL.at(iDigit) = minus2LnL;
  return minus2LnL;
}

void WCSimTotalLikelihood::SetLikelihoodDigitArray(
		WCSimLikelihoodDigitArray* likelihoodDigitArray) {
    fLikelihoodDigitArray = likelihoodDigitArray;

  fChargeLikelihoodVector.clear();
  fChargeLikelihoodVector.resize(0);
  std::vector<WCSimLikelihoodTrackBase*> tmpTracks = fTracks;
  fChargeLikelihoodVector.push_back(WCSimChargePredictor(fLikelihoodDigitArray, fEmissionProfileManager));
  if(fTimeLikelihood != NULL)
  {
    delete fTimeLikelihood;
    fTimeLikelihood = new WCSimTimeLikelihood3(fLikelihoodDigitArray, fEmissionProfileManager);
  }
  SetTracks(tmpTracks);
  ClearVectors();
  return;
}

void WCSimTotalLikelihood::ClearVectors() {
	  fSetVectors = false;
	  fMeasuredCharges.clear();
	  fPredictedCharges.clear();
	  fTotal2LnL.clear();
	  fCharge2LnL.clear();
	  fTime2LnL.clear();
	  fMeasuredCharges.resize(fLikelihoodDigitArray->GetNDigits());
	  fPredictedCharges.resize(fLikelihoodDigitArray->GetNDigits());
	  fTotal2LnL.resize(fLikelihoodDigitArray->GetNDigits());
	  fCharge2LnL.resize(fLikelihoodDigitArray->GetNDigits());
	  fTime2LnL.resize(fLikelihoodDigitArray->GetNDigits());
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

std::vector<double> WCSimTotalLikelihood::GetPredictedTimeVector() const {
	if( !fSetVectors )
	{
		std::cerr << "Error: You're asking for the predicted time vector, but this hasn't been set" << std::endl;
		std::cerr << "       Have you reset the tracks or likelihood digit array since calculating the likelihood?" << std::endl;
		assert(fSetVectors);
	}
  if(! WCSimAnalysisConfig::Instance()->GetUseTime())
  {
    std::cerr << "You're asking for the predicted time vector, but the time component of the likelihood is switched off" << std::endl;
    std::cerr << "Returning a vector of zeroes" << std::endl;
    std::vector<double> zeroes;
    for(int i  = 0; i < fLikelihoodDigitArray->GetNDigits(); ++i)
    {
      zeroes.push_back(0.0);
    }
    return zeroes;
  }
	return fTimeLikelihood->GetAllPredictedTimes();
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

std::vector<double> WCSimTotalLikelihood::GetCharge2LnLVector() const {
	if( !fSetVectors )
	{
		std::cerr << "Error: You're asking for the charge 2LnL vector, but this hasn't been set" << std::endl;
		std::cerr << "       Have you reset the tracks or likelihood digit array since calculating the likelihood?" << std::endl;
		assert(fSetVectors);
	}
	return fCharge2LnL;
}

std::vector<double> WCSimTotalLikelihood::GetTime2LnLVector() const {
	if( !fSetVectors )
	{
		std::cerr << "Error: You're asking for the time 2LnL vector, but this hasn't been set" << std::endl;
		std::cerr << "       Have you reset the tracks or likelihood digit array since calculating the likelihood?" << std::endl;
		assert(fSetVectors);
	}
	return fTime2LnL;
}


WCSimEmissionProfileManager * WCSimTotalLikelihood::GetEmissionProfileManager()
{
  return fEmissionProfileManager;
}
