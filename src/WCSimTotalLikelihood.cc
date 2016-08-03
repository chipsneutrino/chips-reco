#include "WCSimChargePredictor.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimLikelihoodTrack.hh"
#include "WCSimTimeLikelihood3.hh"
#include "WCSimTotalLikelihood.hh"
#include "WCSimAnalysisConfig.hh"
#include "TMath.h"

#include <algorithm>
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
  fTimeScaleFactor = 1.0;

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
    std::cout << "Setting track " << std::endl;
    (*trackIter)->Print();
    fChargeLikelihoodVector.push_back(WCSimChargePredictor(fLikelihoodDigitArray, fEmissionProfileManager));
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
  
  std::vector<std::vector<double> > chargePredictions = CalcPredictedCharges();
  CalcChargeLikelihoods(chargePredictions);
  CalcTimeLikelihoods(chargePredictions);


  double minus2LnL       = 0.0;  // Total -2LnL
  double chargeComponent = 0.0;  // Total charge -2LnL
  double timeComponent   = 0.0;  // Total time -2LnL
  double capSubtracted   = 0.0;  // If (charge+time) exceeds a cap we cut it off
                                 // - this keeps track of the extra 2LnL we removed
                                 // so we make make sure everything adds up

  for(int iDigit = 0; iDigit < fLikelihoodDigitArray->GetNDigits(); ++iDigit) 
  {
    // Likelihood components for this PMT
    double charge2LnL = fCharge2LnL[iDigit];
    double time2LnL   = fTime2LnL[iDigit];
    double total2LnL  = charge2LnL + time2LnL;

    // Add to the running totals, capping the likelihood per digit if necessary
    chargeComponent += charge2LnL;
    timeComponent   += time2LnL;
    if(total2LnL > WCSimTimeLikelihood3::fMaximum2LnL)
    {
        capSubtracted += (total2LnL - WCSimTimeLikelihood3::fMaximum2LnL);
        total2LnL = WCSimTimeLikelihood3::fMaximum2LnL;
    }
    fTotal2LnL[iDigit] = total2LnL;

    if(TMath::IsNaN(total2LnL))
    {
      std::cerr << "iDigit = " << iDigit << ", and charge = " << charge2LnL << " and time = " << time2LnL << std::endl;
      std::cerr << iDigit << " was NaN in Calc2LnL!!" << std::endl;
      return -99999;
    }
    if(total2LnL < 0) 
    { 
        return -99999; 
    }
    else
    {
	  minus2LnL += total2LnL;
	}
  } //for iDigit

  std::cout << "Time component = " << timeComponent << " and charge component = " << chargeComponent << " so total = " << minus2LnL <<  "  (cap of " << WCSimTimeLikelihood3::fMaximum2LnL << " meant " << capSubtracted << " overflow was subtracted)" << std::endl;
  fSetVectors = true;
  return minus2LnL;
}

void WCSimTotalLikelihood::CalcChargeLikelihoods( const std::vector<std::vector<double> >& chargePredictions )
{
  // Charge component of the likelihoods
  // ===================================
  
  // The integral lookup tables mean it's not time-consuming to flip between
  // multiple tracks and work this out one digit at a time 
  for(int iDigit = 0; iDigit < fLikelihoodDigitArray->GetNDigits(); ++iDigit)
  {
    fMeasuredCharges.at(iDigit) = fLikelihoodDigitArray->GetDigit(iDigit)->GetQ();
    if( WCSimAnalysisConfig::Instance()->GetUseCharge())
    {

      // Work out the predicted charge for this PMT
      double totalCharge = SumTotalCharges(chargePredictions[iDigit]);
      assert(!TMath::IsNaN( totalCharge ) );

      // Deal with whether or not the PMT was hit (not used for now)
      double hitUnhit2LnL = 0.0;
      fHit2LnL.at(iDigit) = hitUnhit2LnL;

      // Run the prediction and measured charge through the digitiser
      WCSimLikelihoodDigit* digit  = fLikelihoodDigitArray->GetDigit(iDigit);
      double chargePart            = fDigitizerLikelihood.GetMinus2LnL(totalCharge, digit->GetQ());
      fPredictedCharges.at(iDigit) = totalCharge;
      fCharge2LnL.at(iDigit)       = chargePart;
    }
    else
    {
        fPredictedCharges.at(iDigit) = -999.9;
        fCharge2LnL.at(iDigit)       = 0.0;
    }
  }
}

void WCSimTotalLikelihood::CalcTimeLikelihoods( const std::vector<std::vector<double> >& chargePredictions )
{
    // Time component of the likelihoods
    // =================================
    if( WCSimAnalysisConfig::Instance()->GetUseTime())
    {
        // Get the individual -2LnL time components from each track
        // Vector is vec[digit][track]
        // Address them directly without .at() because it's faster and they should be guaranteed to be filled
        std::vector<std::vector<double> > timeParts = fTimeLikelihood->Calc2LnLComponentsAllDigits();

        // If there is more than one track we need to add likelihoods weighted by predicted charge
        if(timeParts.size() > 1)
        {
            size_t nTracks = fTracks.size();
            
            // Add up the total predicted charges
            for(int iDigit = 0; iDigit < fLikelihoodDigitArray->GetNDigits(); ++iDigit)
            {
                double total = 0;
                for(size_t iTrack = 0; iTrack < nTracks; ++iTrack)
                {
                    // Add up all the nonzero predicted charges (should be all of them)
                    total += chargePredictions[iDigit][iTrack] * (chargePredictions[iDigit][iTrack] > 0);
                }
                if(total > 0)
                {
                    // Convert the log-likelihood back into probabilities, add the two tracks
                    // together weighted by their predicted charges,  and take the log again
                    double timeProb = 0.0;
                    for(size_t iTrack = 0; iTrack < timeParts[iDigit].size(); ++iTrack)
                    {
                       double trackWeight = chargePredictions[iDigit][iTrack] / total;
                       timeProb += exp(-0.5 * timeParts[iDigit][iTrack]) * trackWeight;
                    }
                    double timePart = -2.0 * TMath::Log(timeProb) * fTimeScaleFactor;
                    if(TMath::IsNaN(timePart))
                    {
                        std::cerr << "NaN for " << iDigit << "  timePart = " << timePart << std::endl;
                        for(size_t iTrack = 0; iTrack < timeParts[iDigit].size(); ++iTrack)
                        {
                           double trackWeight = chargePredictions[iDigit][iTrack] / total;
                           std::cout << "  track " << iTrack << " weight = " << chargePredictions[iDigit][iTrack] << " / " << total
                                     << " and -2LnL = " << timeParts[iDigit][iTrack] << " so prob = " << exp(-0.5*timeParts[iDigit][iTrack]) * trackWeight << std::endl;
                        }
                        assert(!TMath::IsNaN(timePart));
                    }
                    fTime2LnL[iDigit] = timePart;
                }
                else
                {
                    fTime2LnL[iDigit] = 0.0;
                }
                if( fTime2LnL[iDigit] < 0 )
                { 
                    if(fTime2LnL[iDigit] < -1e-6)
                    { 
                        std::cerr << "Warning: negative -2LnL of " << fTime2LnL[iDigit] << " for digit " << iDigit << " - setting to zero" << std::endl; 
                    }
                    fTime2LnL[iDigit] = 0; 
                }
            }
        }
        else // Weighting is expensive so only do it if there are >1 tracks
        {
            for(int iDigit = 0; iDigit < fLikelihoodDigitArray->GetNDigits(); ++iDigit)
            {
                fTime2LnL[iDigit] = timeParts[0][iDigit];
            }   
        }
    }
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

std::vector<std::vector<Double_t> > WCSimTotalLikelihood::CalcPredictedCharges()
{
    std::vector<std::vector<Double_t> > predictions(fLikelihoodDigitArray->GetNDigits(), std::vector<double>(fTracks.size()));
    for(int iDigit = 0; iDigit < fLikelihoodDigitArray->GetNDigits(); ++iDigit)
    {
       predictions[iDigit] = CalcPredictedCharges(iDigit); 
    }
    return predictions;
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
    fTimeScaleFactor = 1.;
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
      fHit2LnL.clear();
	  fMeasuredCharges.resize(fLikelihoodDigitArray->GetNDigits());
	  fPredictedCharges.resize(fLikelihoodDigitArray->GetNDigits());
	  fTotal2LnL.resize(fLikelihoodDigitArray->GetNDigits(), 0.0);
	  fCharge2LnL.resize(fLikelihoodDigitArray->GetNDigits(), 0.0);
	  fTime2LnL.resize(fLikelihoodDigitArray->GetNDigits(), 0.0);
	  fHit2LnL.resize(fLikelihoodDigitArray->GetNDigits(), 0.0);
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

std::vector<double> WCSimTotalLikelihood::GetHit2LnLVector() const {
    if( !fSetVectors )
    {
		std::cerr << "Error: You're asking for the hit 2LnL vector, but this hasn't been set" << std::endl;
		std::cerr << "       Have you reset the tracks or likelihood digit array since calculating the likelihood?" << std::endl;
		assert(fSetVectors);
	}
    return fHit2LnL;
}

WCSimEmissionProfileManager * WCSimTotalLikelihood::GetEmissionProfileManager()
{
  return fEmissionProfileManager;
}

void WCSimTotalLikelihood::SetTimeScaleFactor(double scaleFactor)
{
    fTimeScaleFactor = scaleFactor;
    return;
}

double WCSimTotalLikelihood::GetLastCharge2LnL() const
{
    return std::accumulate(fCharge2LnL.begin(), fCharge2LnL.end(), 0.0);
}

double WCSimTotalLikelihood::GetLastTime2LnL() const
{
    return std::accumulate(fTime2LnL.begin(), fTime2LnL.end(), 0.0);
}

double WCSimTotalLikelihood::GetLastHit2LnL() const
{
    return std::accumulate(fHit2LnL.begin(), fHit2LnL.end(), 0.0);
}

double WCSimTotalLikelihood::GetLastTotal2LnL() const
{
    return std::accumulate(fTotal2LnL.begin(), fTotal2LnL.end(), 0.0);
}

double WCSimTotalLikelihood::GetHit2LnL(double predictedCharge)
{
    return 0;
    return -2.0*TMath::Log(1-TMath::Exp(-predictedCharge));
}

double WCSimTotalLikelihood::GetUnhit2LnL(double predictedCharge)
{
    return 2*predictedCharge;
}

// Add up a vector of predicted charges for each PMT, and apply a floor if it's too small
double WCSimTotalLikelihood::SumTotalCharges(const std::vector<double>& totalCharges) const
{
    double totalCharge = std::accumulate(totalCharges.begin(), totalCharges.end(), 0.0);
    if(totalCharge < WCSimChargePredictor::GetMinAllowed()){ totalCharge = WCSimChargePredictor::GetMinAllowed(); }
    return totalCharge;
}
