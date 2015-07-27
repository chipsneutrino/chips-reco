/*
 * WCSimTimeLikelihood2.hh
 *
 *  Created on: 8 Jul 2015
 *      Author: andy
 */

#ifndef INCLUDE_WCSIMTIMELIKELIHOOD2_HH_
#define INCLUDE_WCSIMTIMELIKELIHOOD2_HH_

#include "TString.h"
#include <map>
class WCSimEmissionProfiles;
class WCSimPMTManager;
class WCSimLikelihoodTrackBase;
class WCSimTimePredictor;
class WCSimLikelihoodDigitArray;

class WCSimTimeLikelihood2 {
public:
	WCSimTimeLikelihood2();
	WCSimTimeLikelihood2(WCSimLikelihoodDigitArray * myDigitArray, WCSimEmissionProfiles * myEmissionProfiles);
	virtual ~WCSimTimeLikelihood2();

    /**
     * Specify the tracks used to calculate the likelihood
     * @param myTracks Vector of all the tracks to consider
     */
    void SetTracks(std::vector<WCSimLikelihoodTrackBase*> &myTracks);

    /**
     * Clear the vector of tracks being considered
     */
    void ResetTracks();
    void ClearTracks();

    /**
     * Calculate the combined charge and time likelihood
     * @return Total -2 log(likelihood) from charge and time components
     */
    Double_t Calc2LnL(const unsigned int &iDigit);

    /**
     * Calculate the combined charge and time likelihood
     * @return Total -2 log(likelihood) from charge and time components
     */
    Double_t Calc2LnL();

    std::vector<Double_t> GetAllPredictedTimes();

    double GetPMTTimeResolution(const unsigned int &iDigit);

    double GetGaussianMinusTwoLnL(const double &x, const double &mean, const double &sigma);

private:
    double GetPMTTimeResolution(WCSimLikelihoodDigit * myDigit);
    void GetPredictedTimes();
    bool IsGoodDigit(WCSimLikelihoodDigit * myDigit);
  	WCSimTimePredictor * fTimePredictor;
	  std::vector<WCSimLikelihoodTrackBase*> fTracks; ///< Tracks to consider when calculating the likelihood
    WCSimLikelihoodDigitArray * fLikelihoodDigitArray; ///< Event to build likelihood for
    WCSimPMTManager * fPMTManager;
    WCSimEmissionProfiles * fEmissionProfiles;
    std::vector<double> fAllPreds;
    std::map<TString, double> fPMTTimeConstantMap;
    ClassDef(WCSimTimeLikelihood2,0)
};

#endif /* INCLUDE_WCSIMTIMELIKELIHOOD2_HH_ */
