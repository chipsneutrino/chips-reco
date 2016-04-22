/*
 * WCSimTimeLikelihood3.hh
 *
 *  Created on: 8 Jul 2015
 *      Author: andy
 */

#ifndef INCLUDE_WCSIMTIMELIKELIHOOD3_HH_
#define INCLUDE_WCSIMTIMELIKELIHOOD3_HH_

#include "TF1.h"
#include "TString.h"
#include "WCSimTrackParameterEnums.hh"
#include <map>
#include <vector>
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
class TH2F;
class TH1F; 


class LogFinder{
public: 
    LogFinder(double n, double predMean, double predRMS, double t, double reso) : 
        fN(n), fPredMean(predMean), fPredRMS(predRMS), fT(t), fReso(reso){};
   double operator()(const double x) {
     return this->Eval(x);
   } 
 
   double Derivative(const double x) const
   {
       return DiffLogGaussFirstArrivalGradient(x);
   }
   double Eval(const double x) {
     return DiffLogGaussFirstArrival(x);
   }

   bool ApproxEqual(double x, double y, double epsilon){ return fabs(x-y) < epsilon; }
   double LogGaussian(const double x) const;
   double LogGaussianGradient(const double x) const;
   double LogFirstArrival(const double x) const;
   double LogFirstArrivalGradient(const double x) const;
   double DiffLogGaussFirstArrival(const double x) const;
   double DiffLogGaussFirstArrivalGradient(const double x) const;
private:
   double fN;
   double fPredMean;
   double fPredRMS;
   double fT;
   double fReso;
};

/**
 * Container holding the parameters used by the time likelihood
 * to hold the predicted hit time, its width and probability weighting
 */
class TimePrediction{
public:
	TimePrediction(double mean, double rms, double prob) : fMean(mean), fRMS(rms), fProb(prob){};
	inline double GetMean() const{ return fMean; }
	inline double GetRMS()  const{ return fRMS;  }
	inline double GetProb() const{ return fProb; }

private:
	double fMean; ///< Predicted mean photon arrival time in ns
	double fRMS;  ///< RMS of predicted mean arrival times in ns
	double fProb; ///< Total probability weight to emit towards the PMT
};

/**
 * Container holding the parameters used when stepping along the track
 * to tell us which bins of the time likelihood to look up and where we
 * are relative to the PMT
 */
class StepParameters{
public:
	StepParameters(double cosTheta = 0, int sBin = 0, double deltaCosTheta = 0, double magToPMT = 0) :
		fCosTheta(cosTheta), fSBin(sBin), fDeltaCosTheta(deltaCosTheta), fMagToPMT(magToPMT){};
	inline double GetCosTheta() 	 const { return fCosTheta; }
	inline double GetDeltaCosTheta() const { return fDeltaCosTheta; }
	inline double GetMagToPMT() 	 const { return fMagToPMT; }
	inline int    GetSBin() 		 const { return fSBin; }
	bool operator < (const StepParameters& other) const{ return fCosTheta < other.fCosTheta; }
private:
	double fCosTheta;      ///< The value of cosTheta needed to hit the PMT at this step
	int fSBin;        	   ///< Which s bin of the emission profile we are in
	double fDeltaCosTheta; ///< The width of the cosTheta bin we need to emit in to hit the PMT
	double fMagToPMT;      ///< The distance to the PMT in cm
};



class WCSimTimeLikelihood3 {
public:
    static const double fLog2;
    static const double fSqrtPi;
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
     TimePrediction PredictArrivalTime(WCSimLikelihoodTrackBase * myTrack, WCSimLikelihoodDigit * myDigit);

private:

    /**
     * @brief Work out the speed at which the charge particle travels, either by looking it up
     * or using a speed that was provided at run-time
     * @param myTrack The particle track whose speed we want
     */
    void SetParticleSpeed(WCSimLikelihoodTrackBase * myTrack);

    /**
     * Work out the average refractive index and hence the speed of light in the water, either by looking it up
     * or using a speed that was provided at run-time
     *
     * @param myDigit The PMT we're interested in - we might need to to quantum efficiency weighting
     * @return The effective refractive index
     */
    double GetRefractiveIndex(WCSimLikelihoodDigit * myDigit);

    void GetFourNearestHists(WCSimLikelihoodTrackBase * myTrack);
    std::vector<TH2F*> fSCosThetaVec;
    std::vector<TH1F*> fSVec;
    std::vector<double> fEnergiesVec;
    double fLastE;
    TrackType::Type fLastType;

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


    double FindCrossing(const double n, const double predMean, const double predRMS, const double t, const double reso);
    std::vector<double> FindCrossings(const double n, const double predMean, const double predRMS, const double t, const double reso);
    double FindOverlap(const double n, const double predMean, const double predRMS, const double t, const double reso);
    double ApproximateIntegral(const double n, const double predMean, const double predRMS, const double t, const double reso);


    double ConvertMeanTimeToMeanFirstTime(const double &t, const double &mean, const double &sigma, const double &nPhotons);
    double ConvertMeanTimeToMeanFirstTime(double * x, double * par);
    double IntegrateGaussian(const double mean, const double sigma, const double from, const double to);
    double IntegrateGaussianFromMinusInfinity(const double mean, const double sigma, const double to);
    double IntegrateGaussianToInfinity(const double mean, const double sigma, const double from);
    double IntegrateMeanTimeToMeanFirstTime( const double &t, const double &mean, const double &sigma, const double &nPhotons);
    double IntegrateMeanTimeToMeanFirstTime(double * x, double * par);
    double MinimumOfArrivalAndResolution(const double &t, const double &meanArr, const double rmsArr, const double pmtTime, const double &pmtRes, const double &nPhotons);
    double MinimumOfArrivalAndResolutionTF1(double * x, double * par);


	std::vector<WCSimLikelihoodTrackBase*> fTracks; ///< Tracks to consider when calculating the likelihood
    WCSimLikelihoodDigitArray * fLikelihoodDigitArray; ///< Event to build likelihood for
    WCSimPMTManager * fPMTManager; ///< PMT manager so we can look up the PMT type to work out its time resolution
    WCSimEmissionProfileManager * fEmissionProfileManager; ///< To get the coarse emission profiles for a given track to probability weight arrival times
    std::vector<double> fAllPreds; ///< Vector of the mean predicted arrival time for each PMT
    std::map<TString, double> fPMTTimeConstantMap; ///< Map that stores the time constant used to calculate the time resolution for each unique PMT name

	float fSpeedOfParticle;
    static const double fMinimumLikelihood;
    static const double fMaximumLnL;
private:
    ClassDef(WCSimTimeLikelihood3,0)
};
double ConvertMeanTimeToMeanFirstTime(double * x, double * par);

#endif /* INCLUDE_WCSIMTIMELIKELIHOOD3_HH_ */
