/**
 * \class WCSimLikelihoodFitter
 * This class is used to vary the hypothesized tracks
 * and perform the minimization of the resulting negative
 * log likelihood to reconstruct the event.
 */

#pragma once

#include "WCSimFitterParameters.hh"
#include "WCSimFitterPlots.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimLikelihoodTrackBase.hh"
#include "WCSimRecoDigit.hh"
#include "WCSimRootEvent.hh"
#include "WCSimTotalLikelihood.hh"
#include "WCSimFitterTrackParMap.hh"
#include "WCSimTrackParameterEnums.hh"
#include "WCSimLikelihoodFitter.hh"

#include <exception>
#include <vector>
#include <map>

class WCSimFitterTree;
class WCSimFitterConfig;

class WCSimCosmicFitter : public WCSimLikelihoodFitter
{
public:
	/**
		 * Constructor
		 * @param Event object from WCSim to reconstruct
		 */
	WCSimCosmicFitter(WCSimFitterConfig *config);
	virtual ~WCSimCosmicFitter();
	void RunFits();

	//        static bool RingSort(const std::pair<WCSimRecoRing*,double> &a, const std::pair<WCSimRecoRing*,double> &b);

protected:
	void SeedEvent();

	/**
		 * @brief Run the minimiser over the current LikelihoodDigitArray, with the current track configuration
		 *
		 * @param minAlgorithm Name of fitting algorithm to use - has to be one of the options for ROOT::Math::Factory::CreateMinimizer("Minuit2", minAlgorithm);
		 */
	void Fit(const char *minAlgorithm = "Simplex");

	/**
		 * @brief Minimise using the current track configuration and return the best-fit -2LnL
		 *
		 * @param minAlgorithm Algorithm to use: passed to ROOT::Math::Factory::CreateMinimizer("Minuit2", minAlgorithm);
		 *
		 * @return Best-fit -2LnL
		 */
	double FitAndGetLikelihood(const char *minAlgorithm = "Simplex"); ///< Does the fit and returns the best -2LnL from this minimization

	/**
		 * @brief Minimise using the current track configuration, but only moving the vertex parallel to the current track direction
		 */
	void FitAlongTrack();

	/**
		 * @brief Fit a single event
		 *
		 * @param iEvent Number of the event to fit
		 */
	void FitEventNumber(Int_t iEvent);

	// A special version to only vary one of the tracks at a time.
	Double_t WrapFuncAlongSingleTrack(const Double_t *x);

	//        Double_t fMinimum; ///< Value of -2 log(likelihood) at the best-fit point
	//        Double_t fMinimumTimeComponent; ///< Value of -2 log(time likelihood) at best-fit point
	//        Double_t fMinimumChargeComponent; ///< Value of -2 log(charge likelihood) at best-fit point

	//        Bool_t fFailed; ///< Set to true if the fitter fails, to flag the bad event

	//        Bool_t fIsFirstCall; ///< Flags whether this is the first time the minimizer had calculated a likelihood (to print the seed)
	//        Int_t fStatus; ///< Minimizer convergence status

	//        WCSimFitterPlots * fFitterPlots; ///< Responsible for making quick plots of results
	//        WCSimFitterTree * fFitterTree; ///< Responsible for saving fit information
	//        WCSimFitterConfig * fFitterConfig; ///< For configuring the fitter options
	unsigned int fSingleTrackToFit; ///< WrapFuncAlongSingleTrack only wants one track to vary at a time. This one.

	ClassDef(WCSimCosmicFitter, 1);
};
