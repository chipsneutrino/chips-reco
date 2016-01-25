/*
 * WCSimTimeLikelihood3.hh
 *
 *  Created on: 8 Jul 2015
 *      Author: andy
 */

#ifndef INCLUDE_WCSIMTIMELIKELIHOOD3_HH_
#define INCLUDE_WCSIMTIMELIKELIHOOD3_HH_

#include "TString.h"
#include <map>
class WCSimEmissionProfileManager;
class WCSimPMTManager;
class WCSimLikelihoodDigit;
class WCSimLikelihoodTrackBase;
class WCSimTimePredictor;
class WCSimLikelihoodDigitArray;
class TGraph;
class TGraphErrors;
class TGraphAsymmErrors;
class TH2D;
class TH1D; 


double ConvertMeanTimeToMeanFirstTime(double * x, double * par);
double ConvertMeanTimeToMeanFirstTime(const double &t, const double &mean, const double &sigma, const double &nPhotons);
double IntegrateMeanTimeToMeanFirstTime(double * x, double * par);
double IntegrateMeanTimeToMeanFirstTime( const double &t, const double &mean, const double &sigma, const double &nPhotons);

class WCSimTimeLikelihood3 {
public:
	WCSimTimeLikelihood3();
	WCSimTimeLikelihood3(WCSimLikelihoodDigitArray * myDigitArray, WCSimEmissionProfileManager * myEmissionProfileManager);
	virtual ~WCSimTimeLikelihood3();

    /**
     * Specify the tracks used to calculate the likelihood
     * @param myTracks Vector of all the tracks to consider
     */
    void SetTracks(std::vector<WCSimLikelihoodTrackBase*> &myTracks);

    /**
     * Clear the vector of tracks being considered - two aliases for the same thing
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

    /**
     * Get the vector of the mean predicted time for each PMT
     * @return Vector of the mean of the predicted hit times of each PMT
     */
    std::vector<Double_t> GetAllPredictedTimes();

    /**
     * @input Tube ID of the desired PMT
     * @return Time resolution of a given PMT, taking account of its type and hit charge
     */
    double GetPMTTimeResolution(const unsigned int &iDigit);

    /**
     * Given a track and a PMT hit, work out the mean and RMS predicted arrival time of photons
     * from the track at the PMT by stepping through a coarse emission profile, calculating the
     * arrival time, and weighting by the probability of emitting a photons in each step
     * @param myTrack The track whose contribution is being calculated
     * @param myDigit The hit PMT
     * @return A pair contianing the weighted <mean, rms> of the predicted hit times
     */
    std::pair<double, double> GetArrivalTimeMeanSigma(WCSimLikelihoodTrackBase * myTrack, WCSimLikelihoodDigit * myDigit);

    void SavePlot();

private:
    /**
     * @param LikelihoodDigit for the desired PMT
     * @return Time resolution of a given PMT, taking account of its type and hit charge
     */
    double GetPMTTimeResolution(WCSimLikelihoodDigit * myDigit);

    /**
     * Check if the LikelihoodDigit is one we can use for calculating a predicted time
     * Currently just requires a non-zero charge
     * @param LikelihoodDigit for the desired PMT
     * @return True if digit is good
     */
    bool IsGoodDigit(WCSimLikelihoodDigit * myDigit);


	std::vector<WCSimLikelihoodTrackBase*> fTracks; ///< Tracks to consider when calculating the likelihood
    WCSimLikelihoodDigitArray * fLikelihoodDigitArray; ///< Event to build likelihood for
    WCSimPMTManager * fPMTManager; ///< PMT manager so we can look up the PMT type to work out its time resolution
    WCSimEmissionProfileManager * fEmissionProfileManager; ///< To get the coarse emission profiles for a given track to probability weight arrival times
    std::vector<double> fAllPreds; ///< Vector of the mean predicted arrival time for each PMT
    std::map<TString, double> fPMTTimeConstantMap; ///< Map that stores the time constant used to calculate the time resolution for each unique PMT name
	TGraphErrors * fDistVsPred;
	TGraphAsymmErrors * fTMinusTPredAll;
	TGraphAsymmErrors * fTMinusTPredSource;
	TGraph * fDistVsProb;
	TGraph * fDistVsSpreadProb;
  TGraphErrors * fTMinusTPredVsQ;
  TH2D * fHitQVsRMS;
  TH1D * fTMinusTPredHist;
  TH1D * fTMinusTPredSharpHist;
    ClassDef(WCSimTimeLikelihood3,0)

  float fSpeedOfParticle;
};
double ConvertMeanTimeToMeanFirstTime(double * x, double * par);

#endif /* INCLUDE_WCSIMTIMELIKELIHOOD3_HH_ */
