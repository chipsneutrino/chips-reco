/**
 * \class WCSimTotalLikelihood
 * This class wraps around WCSimChargeLikelihood and
 * WCSimTimeLikelihood to calculate the total combined
 * likelihood to measure the recorded PMT signals given
 * a set of hypothesised tracks.
 */
#ifndef WCSIMTOTALLIKELIHOOD_H
#define WCSIMTOTALLIKELIHOOD_H

#include "WCSimChargeLikelihood.hh"
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
      void SetTracks(std::vector<WCSimLikelihoodTrack> myTracks);

      /**
       * Clear the vector of tracks being considered
       */
      void ResetTracks();

      /**
       * Calculate the combined charge and time likelihood
       * @return Total -2 log(likelihood) from charge and time components
       */
      Double_t Calc2LnL();
 
  protected:
  private:
      WCSimLikelihoodDigitArray * fLikelihoodDigitArray; ///< Event to build likelihood for
      WCSimChargeLikelihood fChargeLikelihood; ///< Charge component of likelihood calculation
       WCSimTimeLikelihood fTimeLikelihood; ///< Time component of likelihood calculation
      std::vector<WCSimLikelihoodTrack> fTracks; ///< Tracks to consider when calculating the likelihood

	ClassDef(WCSimTotalLikelihood,1)
		
};

#endif
