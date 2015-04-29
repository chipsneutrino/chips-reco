/**
 * \class WCSimTotalLikelihood
 * This class wraps around WCSimChargeLikelihood and
 * WCSimTimeLikelihood to calculate the total combined
 * likelihood to measure the recorded PMT signals given
 * a set of hypothesised tracks.
 */
#ifndef WCSIMTOTALLIKELIHOOD_H
#define WCSIMTOTALLIKELIHOOD_H

#include "WCSimChargePredictor.hh"
#include "WCSimTimeLikelihood.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimLikelihoodTrack.hh"

#include "TObject.h"
#include <vector>

class WCSimTotalLikelihood : public TObject
{
  public:
      /**
       * Constructor
       * @param myLikelihoodDigitArray Set of PMT responses to calculate the likelihood for
       */
      WCSimTotalLikelihood( WCSimLikelihoodDigitArray * myLikelihoodDigitArray );

      virtual ~WCSimTotalLikelihood();

      /**
       * Specify the tracks used to calculate the likelihood
       * @param myTracks Vector of all the tracks to consider
       */
      void SetTracks(std::vector<WCSimLikelihoodTrack> &myTracks);

      /**
       * Clear the vector of tracks being considered
       */
      void ResetTracks();

      /**
       * Calculate the combined charge and time likelihood
       * @return Total -2 log(likelihood) from charge and time components
       */
      Double_t Calc2LnL(int iDigit);

      /**
       * Calculate the combined charge and time likelihood
       * @return Total -2 log(likelihood) from charge and time components
       */
      Double_t Calc2LnL();


      /**
       * Calculate the predicted charge at a given PMT
       * @return The total of the expected mean number of photons, for all tracks
       * @param nDigit The number of the PMT in the LikelihoodDigitArray
       */
      Double_t CalcPredictedCharge(unsigned int iDigit);
      
      /**
       * Calculate the predicted charge at a given PMT
       * @return The expected mean number of photons, one entry per track
       * @param nDigit The number of the PMT in the LikelihoodDigitArray
       */
      std::vector<Double_t> CalcPredictedCharges(unsigned int iDigit);

      void SetLikelihoodDigitArray(WCSimLikelihoodDigitArray * likelihoodDigitArray);

      std::vector<double> GetMeasuredChargeVector() const;
      std::vector<double> GetPredictedChargeVector() const;
      std::vector<double> GetTotal2LnLVector() const;

 
  protected:
  private:
      void ClearVectors();

      WCSimLikelihoodDigitArray * fLikelihoodDigitArray; ///< Event to build likelihood for
      std::vector<WCSimChargePredictor> fChargeLikelihoodVector; ///< Charge component of likelihood calculation
      WCSimTimeLikelihood fTimeLikelihood; ///< Time component of likelihood calculation
      WCSimDigitizerLikelihood fDigitizerLikelihood;
      std::vector<WCSimLikelihoodTrack> fTracks; ///< Tracks to consider when calculating the likelihood

      bool fSetVectors;
      std::vector<double> fMeasuredCharges;
      std::vector<double> fPredictedCharges;
      std::vector<double> fTotal2LnL;


	ClassDef(WCSimTotalLikelihood,1)
		
};

#endif
