/**
 Simon McCoy * \class WCSimLikelihoodFitter
 * This class is used to vary the hypothesized tracks
 * and perform the minimization of the resulting negative
 * log likelihood to reconstruct the event.
 */

#ifndef WCSIMLIKELIHOODFITTER_H
#define WCSIMLIKELIHOODFITTER_H

#include "WCSimFitterParameters.hh"
#include "WCSimFitterPlots.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimLikelihoodTrackBase.hh"
#include "WCSimOutputTree.hh"
#include "WCSimRecoDigit.hh"
#include "WCSimRecoSummary.hh"
#include "WCSimRootEvent.hh"
#include "WCSimTotalLikelihood.hh"
#include "WCSimFitterTrackParMap.hh"
#include "WCSimTrackParameterEnums.hh"

#include <exception>
#include <vector>
#include <map>

//class WCSimChargeLikelihood;
class WCSimOutputTree;
class WCSimFitterConfig;

struct FitterArgIsNaN : public std::exception
{
	const char *what() const throw()
	{
		return "Caught NaN argument to fitter";
	}
};

class WCSimLikelihoodFitter
{
public:
	/**
		 * Constructor
		 * @param Event object from WCSim to reconstruct
		 */
	WCSimLikelihoodFitter(WCSimFitterConfig *config);
	virtual ~WCSimLikelihoodFitter();
	void SetFitterPlots(WCSimFitterPlots *fitterPlots);
	void SetOutputTree(WCSimOutputTree *fitterTree);
	void RunFits();
	void RunSurfaces();
	void RunLikelihoodPlots();

	static bool RingSort(const std::pair<WCSimRecoRing *, double> &a, const std::pair<WCSimRecoRing *, double> &b);

protected:
	EventHeader BuildEventHeader();
	WCSimRecoSummary BuildRecoSummary();
	WCSimHitComparison BuildHitComparison();
	PidInfo BuildPidInfo();
	TruthInfo BuildTruthInfo();
	std::string GetRecoType();

	void SetFitInfoStage();

	void SeedEvent();
	void FixVertex();
	void FreeVertex();
	void FixEnergy();
	void FreeEnergy();
	void FixDirection();
	void FreeDirection();
	void FixConversionLength();
	void FreeConversionLength();
	void FixTime();
	void FreeTime();

	void FitEnergy();			///< Wrapper to fit the energy - calls Fit() currently
	void FitEnergyGridSearch(); ///< Alternative way to fit the energy - with a grid search.  Not used currently
	void FitVertex();			///< Wrapper to call Fit() in case we ever want to change the vertex method
	void FitTime();				///< Wrappter to fit the time - just calls Fit() at the moment

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
		 * @brief Special minimisation routine for fitting pi0 events
		 *
		 * @param minAlgorithm Algorithm to use: passed to ROOT::Math::Factory::CreateMinimizer("Minuit2", minAlgorithm);
		 */
	void FitPiZero(const char *minAlgorithm = "Simplex");

	/**
		 * @brief Minimise using the current track configuration, but only moving the vertex parallel to the current track direction
		 */
	void FitAlongTrack();

	/**
		 * @brief Alternative minimiser that uses Metropolis Hastings MCMC to sample the parameter space
		 *
		 * @param nTries Number of steps to take around the parameter space, calculating likelihood each time
		 */
	void MetropolisHastings(const int nTries = 500);

	/**
		 * @brief Alternative MCMC-based mimimiser that only moves the track vertex along its current direction
		 *
		 * @param nTries Number of steps to take around the parameter space, calculating likelihood each time
		 */
	void MetropolisHastingsAlongTrack(const int nTries = 500);

	/**
		 * @brief Called by FitEventNumber if we're doing a pi0 fit; uses a different sequence of fixing/freeing track parameters
		 */
	void FitPiZeroEvent();

	/**
		 * @brief Fit a single event
		 *
		 * @param iEvent Number of the event to fit
		 */
	void FitEventNumber(Int_t iEvent);

	/**
		 * @brief Build a LikelihoodDigitArray and reset the fitter to allow the new event to be fitted
		 *
		 * @param iEvent Number of the event to setup
		 */
	void SetEvent(Int_t iEvent);

	/**
		 * @brief Clear everything ready to fit a new event
		 */
	void ResetEvent();

	/**
		 * @brief Called after doing a minimisation to update the best-fit tracks to the latest ones
		 */
	void UpdateBestFits();

	/**
		 * Work out how many free parameters Minuit needs to run over
		 * @return Number of parameters to use in Minuit
		 */
	UInt_t GetNPars();

	/**
		 * Wrapper to WCSimTotalLikelihood::Calc2LnL so that
		 * we can use it as Minuit's FCN
		 * @param x Array of parameters used to construct the track objects for the fit
		 * @return -2 log(likelihood) for this set of tracks
		 */
	Double_t WrapFunc(const Double_t *x);

	/**
		 * @brief Wrapper to WCSimTotalLikelihood::Calc2LnL so we can use it as Minuit's FCN to minimise
		 *
		 * @param x Array of parameters used to construct the track objects for the fit
		 *
		 * @return -2 log(likelihood) for this set of tracks
		 */
	Double_t WrapFuncAlongTrack(const Double_t *x);

	/**
		 * @brief Special wrapper to WCSimTotalLikelihood::Calc2LnL for pi0 events that can handle the mass constraint
		 *
		 * @param x Array of parameters used to construct the track objects for the fit
		 *
		 * @return -2 log(likelihood) for this set of tracks
		 */
	Double_t WrapFuncPiZero(const Double_t *x);

	/**
		 * @brief Penalise fit vertices outside the detector or outside the allowed fit region, depending on how
		 *        far outside they are (just returning a constant big number gives Minuit NaN trouble)
		 *
		 * @param tracksToFit The tracks for which the likelihood is being calculated
		 *
		 * @return Penalty term for out of bounds vertices
		 */
	Double_t GetPenalty(const std::vector<WCSimLikelihoodTrackBase *> &tracksToFit);

	/**
		 * Get the minimum value of -2 log(likelihood) returned by the fit
		 * @return Minimum value of -2 log(likelihood)
		 */
	Double_t GetMinimum();

	/**
		 * Get the convergence status of the minimizer.  Zero
		 * means it converged with no error, everything else is an error code
		 */
	Int_t GetStatus();

	/**
		 * Run the previous WCSimAnalysis fitter that uses a Hough transform
		 * to seed our fit in roughly the right part of parameter space
		 * @param myReco WCSimReco object used to run the old fitter
		 */
	void SeedParams();

	/**
		 * Get a vector containing the track(s) at the best-fit point
		 * @return Best-fit track(s) for this event
		 */
	std::vector<WCSimLikelihoodTrackBase *> GetBestFit();

	/**
		 * As a temporary measure to get results out of the fit, we'll perform
		 * a grid search to determine the energy, instead of turning Minuit loose on it
		 * (until we have a working method that interpolates for a smooth surface)
		 * @param best2LnL Will be updated with the best-fit -2Ln(likelihood)
		 * @param bestEnergies Will be updated with the best-fit energies for each independent track
		 */
	void PerformEnergyGridSearch(Double_t &best2LnL, std::vector<Double_t> &bestEnergies);

	/**
		 * If Minuit encounters a problem (e.g. LnL doesn't change with the track parameters) then
		 * the next iteration's track parameters can be NaN - we want to catch this so it doesn't
		 * break everything
		 * @param x The track parameters being used
		 */
	void CheckTrackParametersForNaN(const Double_t *x) const;

	/**
		 * If Minuit encounters a problem (e.g. LnL doesn't change with the track parameters) then
		 * the next iteration's track parameters can be NaN - we want to catch this so it doesn't
		 * break everything.  This version is for when we're not using the full set of track parameter
		 * (eg we're fitting along the track) so you have to tell it the size of the array instead of
		 * guessing it from the number of track parameters
		 * @param x The track parameters being used
		 * @param sizeofX Number of parameters (size of the x array)
		 */
	void CheckTrackParametersForNaN(const Double_t *x, const unsigned int sizeofX) const;

	/**
		 * @brief Fill the reco and reco-true plots in WCSimFitterPlots
		 */
	void FillPlots();

	/**
		 * @brief Fill the FitTree in WCSimOutputTree
		 */
	void FillTree();

	/**
		 * @brief Sweep out the likelihood as a function of one fit parameter around the best-fit
		 *
		 * @param trackPar Track number and track parameters to sweep out
		 */
	void Make1DSurface(std::pair<unsigned int, FitterParameterType::Type> trackPar);

	/**
		 * @brief Sweep out the likelihood as a function of two fit parameters around the best-fit
		 *
		 * @param trackPar Two sets of track number and parameter to sweep out
		 */
	void Make2DSurface(
		std::pair<std::pair<unsigned int, FitterParameterType::Type>,
				  std::pair<unsigned int, FitterParameterType::Type>>
			trackPar);

	/**
		 * @brief Sweep out a 2D surface of the displacement in the vertex in distance and time, adjusting parallel to the best-fit track only
		 */
	void Make2DSurfaceAlongTrack();

	/**
		 * @brief Find out whether one of the MC truth tracks would have escaped the detector
		 *
		 * @param iTrack Track number to check
		 *
		 * @return True if the track ends outside the detector, according to the emission profile
		 */
	Bool_t GetTrueTrackEscapes(unsigned int iTrack) const;

	/**
		 * @brief Find out whether one of the best-fit tracks would have escaped the detector
		 *
		 * @param iTrack Track number to check
		 *
		 * @return True if the track ends outside the detector, according to the emission profile
		 */
	Bool_t GetFitTrackEscapes(unsigned int iTrack) const;

	/**
		 * @brief Find out whether a given would have escaped the detector
		 *
		 * @param track Track object to check
		 *
		 * @return True if the track ends outside the detector, according to the emission profile
		 */
	Bool_t GetTrackEscapes(WCSimLikelihoodTrackBase *track) const;

	/**
		 * @brief Check if the second decay photon from a pi0 has its energy constrained so the the invariant mass is that of a pi0
		 *
		 * @return True if the constraint is being applied
		 */
	Bool_t GetUsePiZeroMassConstraint() const;

	/**
		 * @brief Work out the energy required for the second pi0 decay photon to give an invariant mass equal to that of a pi0
		 *
		 * @param x The set of track parameters passed to the minimiser
		 *
		 * @return The required energy of the second tracl.  Note this is prevented from going outside the bounds provided to the fitter.
		 */
	Double_t GetPiZeroSecondTrackEnergy(const Double_t *x);

	/**
		 * @brief Check if a point is both inside the detector and inside the bounds set using the FitterInterface
		 *
		 * @param iTrack Number of the track whose bounds we're checking against
		 * @param x X coordinate of the point, in cm
		 * @param y Y coordinate of the point, in cm
		 * @param z Z coordinate of the point, in cm
		 *
		 * @return True if the particle is inside all these bounds
		 */
	Bool_t IsInsideAllowedRegion(const Int_t &iTrack, const Double_t &x, const Double_t &y, const Double_t &z);

	/**
		 * @brief Check if a point is outside the detector or outside the bounds set using the FitterInterface
		 *
		 * @param iTrack Number of the track whose bounds we're checking against
		 * @param x X coordinate of the point, in cm
		 * @param y Y coordinate of the point, in cm
		 * @param z Z coordinate of the point, in cm
		 *
		 * @return True if the particle is outside any of these bounds
		 */
	Bool_t IsOutsideAllowedRegion(const Int_t &iTrack, const Double_t &x, const Double_t &y, const Double_t &z);

	/**
		 * @brief Adjust the track vertex back inside the detector by moving it along its direction
		 *
		 * @param iTrack Number of the track to be moved
		 * @param seedX Vertex X position, in cm (will be updated)
		 * @param seedY Vertex Y position, in cm (will be updated)
		 * @param seedZ Vertex Z position, in cm (will be updated)
		 * @param dirX Track direction x component
		 * @param dirY Track direction y component
		 * @param dirZ Track direction z component
		 */
	void MoveBackInside(const Int_t iTrack, Double_t &seedX, Double_t &seedY, Double_t &seedZ, const Double_t &dirX,
						const Double_t &dirY, const Double_t &dirZ);

	WCSimTotalLikelihood *fTotalLikelihood;							///< Class used to calculate the total (combined charge and time) likelihood that we minimize
	WCSimRootEvent *fRootEvent;										///< Simulated event from WCSim to be reconstructed
	WCSimLikelihoodDigitArray fLikelihoodDigitArrayObj;				///< Array of PMT responses
	WCSimLikelihoodDigitArray *fLikelihoodDigitArray;				///< Array of PMT responses
	std::vector<TrackType::Type> fTypes;							///< Particle type to fit for
	std::map<Int_t, Int_t> fParMap;									///< Map relating the number of tracks (key) to the number of Minuit parameters (value)
	std::vector<WCSimLikelihoodTrackBase *> fBestFit;				///< The best-fit track or tracks
	std::vector<WCSimLikelihoodTrackBase *> *fTrueLikelihoodTracks; ///< The MC-truth tracks for this event

	/**
		 * Used to rescale all the track parameters from the range 0 to 1
		 * (which was useful for the minimizer) to physical values used to
		 * make track objects
		 * This was for testing a different minimizer algorithm and
		 * isn't currently used
		 * @param x Track vertex x
		 * @param y Track vertex y
		 * @param z Track vertex z
		 * @param t Track start time
		 * @param th Track angle to z axix
		 * @param phi Track azimuthal angle
		 * @param E Track energy
		 * @param type Particle type
		 * @return A track object created with the rescaled values
		 */
	WCSimLikelihoodTrackBase *RescaleParams(Double_t x, Double_t y, Double_t z, Double_t t, Double_t th,
											Double_t phi, Double_t E, TrackType::Type type);

	/**
		 * @brief Checks that the event has some hits before trying to fit it
		 *
		 * @return True if there are some hits to try and reconstruct
		 */
	Bool_t CanFitEvent() const;

	Int_t fEvent;					  ///< Number of the event currently being fitted
	Double_t fMinimum;				  ///< Value of -2 log(likelihood) at the best-fit point
	Double_t fMinimumTimeComponent;	  ///< Value of -2 log(time likelihood) at best-fit point
	Double_t fMinimumChargeComponent; ///< Value of -2 log(charge likelihood) at best-fit point
	Double_t fMinimumHitComponent;	  ///< Value of -2 log(hit likelihood) at best-fit point
	Double_t fMinimumCutoffComponent; ///< Amount of -2LnL lost due to cutoff at best-fit

	Bool_t fFailed; ///< Set to true if the fitter fails, to flag the bad event

	Bool_t fIsFirstCall; ///< Flags whether this is the first time the minimizer had calculated a likelihood (to print the seed)
	Int_t fStatus;		 ///< Minimizer convergence status

	WCSimFitterPlots *fFitterPlots;			   ///< Responsible for making quick plots of results
	WCSimOutputTree *fOutputTree;			   ///< Responsible for saving fit information
	WCSimFitterTrackParMap fFitterTrackParMap; ///< Stores and looks-up track parameters in an array
	WCSimFitterConfig *fFitterConfig;		   ///< For configuring the fitter options

	int fCalls; ///< How many times the likelihood has been calculated for this event
	ClassDef(WCSimLikelihoodFitter, 1);
};

#endif // WCSIMLIKELIHOODFITTER_H
