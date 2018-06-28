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
#include "WCSimTimeLikelihood3.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimLikelihoodTrackBase.hh"
#include "WCSimEmissionProfileManager.hh"

#include "TObject.h"
#include <vector>

class WCSimTotalLikelihood: public TObject {
	public:
		/**
		 * Constructor
		 * @param myLikelihoodDigitArray Set of PMT responses to calculate the likelihood for
		 */
		WCSimTotalLikelihood(WCSimLikelihoodDigitArray * myLikelihoodDigitArray);

		virtual ~WCSimTotalLikelihood();

		/**
		 * Specify the tracks used to calculate the likelihood
		 * @param myTracks Vector of all the tracks to consider
		 */
		void SetTracks(std::vector<WCSimLikelihoodTrackBase*> &myTracks);

		/**
		 * Clear the vector of tracks being considered
		 */
		void ResetTracks();

		/**
		 * Calculate the combined charge and time likelihood
		 * @return Total -2 log(likelihood) from charge and time components
		 */
		Double_t Calc2LnL();

		/**
		 * Calculate the predicted charge for all PMTs
		 * @return Vector of the expected mean number of photons from each track in the form vec[digit][track]
		 */
		std::vector<std::vector<Double_t> > CalcPredictedCharges();

		/**
		 * Calculate the predicted charge at a given PMT
		 * @return The expected mean number of photons, one entry per track
		 * @param nDigit The number of the PMT in the LikelihoodDigitArray
		 */
		std::vector<Double_t> CalcPredictedCharges(unsigned int iDigit);

		/**
		 * @brief Fill the vector fCharge2LnL with charge likelihoods for each PMT
		 */
		void CalcChargeLikelihoods(const std::vector<std::vector<double> >& chargePredictions);

		/**
		 * @brief
		 */
		void CalcTimeLikelihoods(const std::vector<std::vector<double> >& chargePredictions);

		/**
		 * @brief Given some predicted charge, calculate the probability that a PMT was actually hit
		 *
		 * @param predictedCharge The predicted mean number of p.e.
		 *
		 * @return Probability that the number of p.e. measured was non-zero
		 */
		Double_t GetHit2LnL(double predictedCharge);

		/**
		 * @brief Given some predicted charge, calcualte the probability that a PMT was not hit
		 *
		 * @param predictedCharge The predicted mean number of p.e.
		 *
		 * @return Probability that the number of p.e. measured was zero
		 */
		Double_t GetUnhit2LnL(double predictedCharge);

		void SetLikelihoodDigitArray(WCSimLikelihoodDigitArray * likelihoodDigitArray);

		void SetTimeScaleFactor(double scaleFactor = 1.0);
		inline double GetTimeScaleFactor() {
			return fTimeScaleFactor;
		}

		std::vector<double> GetMeasuredChargeVector() const;
		std::vector<double> GetPredictedChargeVector() const;
		std::vector<double> GetTotal2LnLVector() const;
		std::vector<double> GetHit2LnLVector() const; // Component detailing whether or not a PMT was hit
		std::vector<double> GetCharge2LnLVector() const; // Leigh: Get the charge component only
		std::vector<double> GetTime2LnLVector() const; // Leigh: Get the time component only
		std::vector<double> GetPredictedTimeVector() const;
		double GetLastCutoff2LnL() const;
		double GetLastHit2LnL() const;
		double GetLastTime2LnL() const;
		double GetLastCharge2LnL() const;
		double GetLastTotal2LnL() const;

		WCSimEmissionProfileManager * GetEmissionProfileManager();

	protected:
	private:
		void ClearVectors();
		double SumTotalCharges(const std::vector<double>& totalCharges) const;

		WCSimLikelihoodDigitArray * fLikelihoodDigitArray; ///< Event to build likelihood for
		std::vector<WCSimChargePredictor> fChargeLikelihoodVector; ///< Charge component of likelihood calculation
		WCSimTimeLikelihood3 * fTimeLikelihood; ///< Time component of likelihood calculation
		WCSimDigitizerLikelihood fDigitizerLikelihood;
		std::vector<WCSimLikelihoodTrackBase*> fTracks; ///< Tracks to consider when calculating the likelihood
		double fTimeScaleFactor; ///< Optional factor to scale the time likelihood by to alter its weight compared to charge

		bool fSetVectors;
		std::vector<double> fMeasuredCharges;
		std::vector<double> fPredictedCharges;
		std::vector<double> fTotal2LnL;
		std::vector<double> fCharge2LnL; // Leigh: Store the charge part only
		std::vector<double> fTime2LnL; // Leigh: Store the time part only
		std::vector<double> fHit2LnL; // Likelihood that a hit PMT was hit or an unhit PMT was unhit
		std::vector<double> fCutoff2LnL; // Likelihood thrown away by the maximum cap
		WCSimEmissionProfileManager * fEmissionProfileManager;

		ClassDef(WCSimTotalLikelihood,1)

};

#endif
