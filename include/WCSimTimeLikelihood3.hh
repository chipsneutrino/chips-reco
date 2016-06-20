/*
 * WCSimTimeLikelihood3.hh
 *
 *  Created on: 8 Jul 2015
 *      Author: andy
 */

#ifndef INCLUDE_WCSIMTIMELIKELIHOOD3_HH_
#define INCLUDE_WCSIMTIMELIKELIHOOD3_HH_

#include "TF1.h"
#include "TRandom3.h"
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


/**
 * \class Logfinder
 * \brief Class to perform calculations with the difference between the logs of
 * the first photon arrival and PMT hit time PDFs
 *
 * When calculating the overlap between the two PDFs for the arrival time
 * of the first photon and the PMT hit time smeared by its resolution, we
 * need to find where these PDFs intersect.  At this point their logs are equal.
 * This class calculates those logs and the difference between them (and also the
 * derivative of the logs and difference) in a way that can be used by ROOT's
 * root-finding algorithms
 */
class LogFinder{
public: 
	/**
	 * Constructor
	 * @param n Number of photons hitting PMT
	 * @param predMean Mean of predicted arrival time distribution
	 * @param predRMS RMS of predicted arrival time distribution
	 * @param t Measured PMT hit time
	 * @param reso PMT time resolution
	 */
    LogFinder(double n, double predMean, double predRMS, double t, double reso) : 
        fN(n), fPredMean(predMean), fPredRMS(predRMS), fT(t), fReso(reso){};

   double operator()(const double x) {
     return this->Eval(x);
   } 
 
   /**
    * Calculate derivative of the difference between the logs of the first arrival time
    * PDF and the PMT hit time PDF
    * @param x The time at which to calculate the derivative
    * @return The derivative of the difference between the logs of the first arrival time
    * and PMT hit time PDFs
    */
   double Derivative(const double x) const
   {
       return DiffLogGaussFirstArrivalGradient(x);
   }

   /**
    * Calculate the difference between the logs of the first arrival time
    * PDF and the PMT hit time PDF
    * @param x The time at which to calculate the derivative
    * @return The difference between the logs of the first arrival time and PMT hit time PDFs
    */
   double Eval(const double x) {
     return DiffLogGaussFirstArrival(x);
   }

   /**
    * Determine whether two numbers are equal to within a tolerance epsilon
    * @param x The first number
    * @param y The second number
    * @param epsilon The absolute tolerance
    * @return True if the numbers are equal within the tolerance (i.e., the magnitude of
    * their difference is less than epsilon)
    */
   bool ApproxEqual(double x, double y, double epsilon){ return fabs(x-y) < epsilon; }

   /**
    * Calculate the log of a Gaussian PDF described by the PMT hit time (fT) and its resolution (fReso)
    * @param x The value at which to calculate the log
    * @return The log of a unit-normalised Gaussian PDF evaluated at x
    */
   double LogGaussian(const double x) const;

   /**
    * Calculate the gradient of the log of a Gaussian PDF described by the hit time (fT)
    * and its resolution (fReso)
    * @param x The value at which to calculate the gradient
    * @return Gradient of the log of a unit-normalised Gaussian evaluated at x
    */
   double LogGaussianGradient(const double x) const;

   /**
    * Calculate the log of the PDF describing the arrival time of the first photon, given fN photons
    * with Gaussian-distributed arrival times, centred on fPredMean with width fPredReso
    * @param x The value at which to calculate the log
    * @return Log of the PDF for the arrival time of the first photon, evaluated at x
    */
   double LogFirstArrival(const double x) const;

   /**
    * Calculate the gradient of the log of the PDF describing the arrival time of the first photon,
    * given fN photons with Gaussian-distributed arrival times, centred on fPredMean with width fPredReso
    * @param x The value at which to calculate the gradient
    * @return Gradient of the log of the PDF for the arrival time of the first photon, evaluated at x
    */
   double LogFirstArrivalGradient(const double x) const;

   /**
    * Calculate the difference between the PDFs describing the measured PMT hit time and the arrival
    * time of the first photon
    * @param x The value at which to evaluate the log
    * @return Log(PMT time resolution Gaussian PDF) - Log(First arrival PDF) at x
    */
   double DiffLogGaussFirstArrival(const double x) const;

   /**
    * Calculate the gradient of the  difference between the PDFs describing the measured PMT hit
    * time and the arrival time of the first photon
    * @param x The value at which to evaluate the gradient
    * @return Gradient of Log(PMT time resolution Gaussian PDF) - Log(First arrival PDF) at x
    */
   double DiffLogGaussFirstArrivalGradient(const double x) const;
private:

   double fN; 		 ///< Number of photoelectrons detected by the PMT
   double fPredMean; ///< Predicted mean arrival time of photons at the PMT
   double fPredRMS;  ///< Predicted RMS of arrival time of photons at the PMT
   double fT; 		 ///< Hit time recorded by the PMT
   double fReso; 	 ///< PMT time resolution
};

/**
 * \class TimePrediction
 *
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
 * \class StepParameters
 *
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


/**
 * \class WCSimTimeLikelihood3
 *
 * \brief Class to calculate the expected arrival time of photons at the PMT, compare this with
 * the measured time, and calculate a likelihood based on this comparison
 */
class WCSimTimeLikelihood3 {
public:
    static const double fLog2;   ///< Log (base-e) of 2 (used a lot in calculations here)
    static const double fSqrtPi; ///< Square root of pi (used a lot in calculations here)


	WCSimTimeLikelihood3();
	WCSimTimeLikelihood3(WCSimLikelihoodDigitArray * myDigitArray, WCSimEmissionProfileManager * myEmissionProfileManager);
	virtual ~WCSimTimeLikelihood3();

    /**
     * Specify the tracks used to calculate the likelihood
     * @param myTracks Vector of all the tracks to consider
     */
    void SetTracks(std::vector<WCSimLikelihoodTrackBase*> &myTracks);

    /**
     * Clear the vector of tracks being considered - alias for ClearTracks()
     */
    void ResetTracks();

    /**
     * Clear the vector of tracks being considered - alias for ResetTracks()
     */
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

    /**
     * Load the emission profile histograms corresponding to the four energies closest to the energy of myTrack
     * for which we have profiles
     * @param myTrack That track for whose energy we want the nearby histograms
     */
    void GetFourNearestHists(WCSimLikelihoodTrackBase * myTrack);

    std::vector<TH2F*> fSCosThetaVec; ///< Vector of the 4 nearest (cosTheta,s) time emission profile histograms
    std::vector<TH1F*> fSVec;		  ///< Vector of the 4 nearest time emission profile histograms in s
    std::vector<double> fEnergiesVec; ///< Vector of the 4 nearest energies for which we have emission profiles
    double fLastE;					  ///< Energy of the last track for which we loaded nearby profile histograms
    TrackType::Type fLastType;		  ///< Particle type of the last track for which we loaded nearby profile histograms

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

    /**
     * Find a crossing point between the PMT hit time and photon first arrival time PDFs
     * @param n Number of photoelectrons detected by the PMT
     * @param predMean Predicted mean arrival time of photons at the PMT
     * @param predRMS Predicted RMS of the arrival time of photons at the PMT
     * @param t Measured hit time recorded by the PMT
     * @param reso Time resolution of the PMT
     * @return A time at which the two PDFs intersect
     */
    double FindCrossing(const double n, const double predMean, const double predRMS, const double t, const double reso);

    /**
     * Find multiple (usually all) crossing points between the PMT hit time and photon
     * first arrival time PDFs
     * @param n Number of photoelectrons detected by the PMT
     * @param predMean Predicted mean arrival time of photons at the PMT
     * @param predRMS Predicted RMS of the arrival time of photons at the PMT
     * @param t Measured hit time recorded by the PMT
     * @param reso Time resolution of the PMT
     * @return A vector of times at which the two PDFs intersect
     */
    std::vector<double> FindCrossings(const double n, const double predMean, const double predRMS, const double t, const double reso);

    /**
     * Find the area of overlap between the PMT hit time and photon first arrival time PDFs
     * @param n Number of photoelectrons detected by the PMT
     * @param predMean Predicted mean arrival time of photons at the PMT
     * @param predRMS Predicted RMS of the arrival time of photons at the PMT
     * @param t Measured hit time recorded by the PMT
     * @param reso Time resolution of the PMT
     * @return The area of the overlapping region between the PMT hit time and photon first arrival time PDFs
     */
    double FindOverlap(const double n, const double predMean, const double predRMS, const double t, const double reso);

    /**
     * Find an approximation to the area of overlap between the PMT hit time and photon
     * first arrival time PDFs, using logarithms and the trapezium rule in case doing it
     * normally evaluates to zero because of floating point arithmetic
     * @param n Number of photoelectrons detected by the PMT
     * @param predMean Predicted mean arrival time of photons at the PMT
     * @param predRMS Predicted RMS of the arrival time of photons at the PMT
     * @param t Measured hit time recorded by the PMT
     * @param reso Time resolution of the PMT
     * @return The approximate area of the overlapping region between the PMT hit time and photon first arrival time PDFs
     */
    double ApproximateIntegral(const double n, const double predMean, const double predRMS, const double t, const double reso);

    /**
     * Calculate the value of the PDF that the first of n photons will arrive at
     * a certain time, assuming that each photon's arrival time is a Gaussian
     * with the same mean and RMS
     * @param x Time at which to evaluate the PDF
     * @param mean Predicted mean arrival time of photons at the PMT
     * @param sigma Predicted RMS of the arrival time of photons at the PMT
     * @param nPhotons Number of photoelectrons detected by the PMT
     * @return The PDF that the first of n photons arrives before a certain time, evaluated at x
     */
    double ConvertMeanTimeToMeanFirstTime(const double &x, const double &mean, const double &sigma, const double &nPhotons);

    /**
     * Alternative signature for ConvertMeanTimeToMeanFirstTime so that it can
     * be used with a TF1
     * @param x Array where x[0] is the time at which to evaluate the PDF
     * @param par Array when par[0] is the predicted mean arrival time, par[1] the RMS, and par[2] the number of photoelectrons
     * @return The PDF that the first of n photons arrives before a certain time, evaluated at x[0]
     */
    double ConvertMeanTimeToMeanFirstTime(double * x, double * par);
    double ConvertMeanTimeToMeanFirstTimeTF1(double * x, double * par) { return ConvertMeanTimeToMeanFirstTime(x, par); };


    /**
     * Integrate a unit-normalised Gaussian that describes the PDF for
     * the time measured by a PMT with a known mean and width between two values
     * @param mean The mean of the Gaussian
     * @param sigma The width of the Gaussian
     * @param from The value of x to integrate from
     * @param to The value of x to integrate up to
     * @return The integral of the Gaussian between the two values
     */
    double IntegrateGaussian(const double mean, const double sigma, const double from, const double to);

    /**
     * Integrate a unit-normalised Gaussian with a known mean and width, from
     * minus infinity up to some value of x
     * @param mean The mean of the Gaussian
     * @param sigma The width of the Gaussian
     * @param to The value of x to integrate up to
     * @return The integral of the Gaussian from -infinity to the chosen value
     */
    double IntegrateGaussianFromMinusInfinity(const double mean, const double sigma, const double to);

    /**
     * Integrate a unit-normalised Gaussian with a known mean and width, from
     * some value of x up to infinity
     * @param mean The mean of the Gaussian
     * @param sigma The width of the Gaussian
     * @param to The value of x to integrate from
     * @return The integral of the Gaussian from the chosen value to infinity
     */
    double IntegrateGaussianToInfinity(const double mean, const double sigma, const double from);

    /**
     * Integrate the PDF for the arrival time of the first photon from
     * -infinity up to some value
     * @param t The value up to which the integral should be evaluated
     * @param mean The mean arrival time of each individual photon
     * @param sigma The RMS arrival time of each individual photon
     * @param nPhotons The total number of photons
     * @return The integral of the PDF for the first arrival time from -infinity to t
     */
    double IntegrateMeanTimeToMeanFirstTime( const double &t, const double &mean, const double &sigma, const double &nPhotons);

    /**
     * Alternative signature for IntegrateMeanTimeToMeanFirstTime so that it can
     * be used with a TF1
     * @param x Array where x[0] is the value up to which the integral should be evaluated
     * @param par Array where par[0] is the mean and par[1] the RMS of the individual photons' arrival times,
     * and par[2] is the total number of photons
     * @return The integral of the PDF for the first arrival time from -infinity to x[0]
     */
    double IntegrateMeanTimeToMeanFirstTime(double * x, double * par);

    /**
     * Find and return the minimum of the PDFs for the arrival time of the first
     * photon and the time recorded by the PMT smeared by its resolution
     * @param t The time at which to find the minimum
     * @param meanArr The mean arrival time of each photon
     * @param rmsArr The RMS arrival time of each photon
     * @param pmtTime The hit time recorded by the PMT
     * @param pmtRes The time resolution of the PMT
     * @param nPhotons The toal number of photoelectrons recorded by the PMT
     * @return The minimum of the PDFs for the first photon's arrival time and the PMT's recorded time
     */
    double MinimumOfArrivalAndResolution(const double &t, const double &meanArr, const double rmsArr, const double pmtTime, const double &pmtRes, const double &nPhotons);

    /**
     * Alternative signature for MinimumOfArrivalAndResolution so that
     * it can be used with a TF1
     * @param x Array where x[0] is the time at which to find the minimum
     * @param par Array where par[0] and par[1] are the predicted mean and RMS arrival times of each photon,
     * par[2] and par[3] are the time recorded by the PMT and its resolution, and par[3] is the number of
     * detected photoelectrons
     * @return The minimum of the PDFs for the first photon's arrival time and the PMT's recorded time
     */
    double MinimumOfArrivalAndResolutionTF1(double * x, double * par);

    double GetScatteredTimeLikelihood(WCSimLikelihoodDigit * myDigit);


	std::vector<WCSimLikelihoodTrackBase*> fTracks; 		///< Tracks to consider when calculating the likelihood
    WCSimLikelihoodDigitArray * fLikelihoodDigitArray; 		///< Event to build likelihood for
    WCSimPMTManager * fPMTManager; 							///< PMT manager so we can look up the PMT type to work out its time resolution
    WCSimEmissionProfileManager * fEmissionProfileManager;  ///< To get the coarse emission profiles for a given track to probability weight arrival times
    std::vector<double> fAllPreds; 							///< Vector of the mean predicted arrival time for each PMT
    std::map<TString, double> fPMTTimeConstantMap; 			///< Map that stores the time constant used to calculate the time resolution for each unique PMT name

	float fSpeedOfParticle; 				///< The speed at which the charged particle is assumed to trave, divided by c
    static const double fMinimumLikelihood; ///< The minimum likelihood we're allowed to return before cutting off
    static const double fMaximumLnL; 		///< The maximum log-likelihood we're allowed to return before cutting off
private:
    ClassDef(WCSimTimeLikelihood3,0)
};

#endif /* INCLUDE_WCSIMTIMELIKELIHOOD3_HH_ */
