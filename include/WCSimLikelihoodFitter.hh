/**
 * \class WCSimLikelihoodFitter
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
#include "WCSimReco.hh"
#include "WCSimRecoDigit.hh"
#include "WCSimRootEvent.hh"
#include "WCSimTotalLikelihood.hh"
#include "WCSimFitterTrackParMap.hh"
#include "WCSimTrackParameterEnums.hh"

#include <vector>
#include <map>

//class WCSimChargeLikelihood;
class WCSimFitterTree;
class WCSimFitterConfig;


class WCSimLikelihoodFitter
{
    public:
        /**
         * Constructor
         * @param Event object from WCSim to reconstruct
         */
		WCSimLikelihoodFitter(WCSimFitterConfig * config);
        virtual ~WCSimLikelihoodFitter();
        void SetFitterPlots(WCSimFitterPlots * fitterPlots);
        void SetFitterTree(WCSimFitterTree * fitterTree);
        void RunFits();
        void RunSurfaces();

        static bool RingSort(const std::pair<WCSimRecoRing*,double> &a, const std::pair<WCSimRecoRing*,double> &b);


    protected:


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

        void FitEnergy(); ///< Wrapper to fit the energy - calls Fit() currently
        void FitEnergyGridSearch(); ///< Alternative way to fit the energy - with a grid search.  Not used currently
        void FitVertex(); ///< Wrapper to call Fit() in case we ever want to change the vertex method
        void FitTime();
        void Fit(const char * minAlgorithm = "Simplex");
        double FitAndGetLikelihood(const char * minAlgorithm = "Simplex"); ///< Does the fit and returns the best -2LnL from this minimization
        void FitPiZero(const char * minAlgorithm = "Simplex");
        void FitAlongTrack();
        void MetropolisHastings(const int nTries = 500);
        void MetropolisHastingsAlongTrack(const int nTries = 500);

        void FitPiZeroEvent();

        void FitEventNumber(Int_t iEvent);
        void SetEvent(Int_t iEvent);
        void ResetEvent();

        /**
         * Perform the minimization
         */
        void Minimize2LnL();
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
        Double_t WrapFunc(const Double_t * x);
        Double_t WrapFuncAlongTrack(const Double_t * x);
        Double_t WrapFuncPiZero(const Double_t * x);

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
        std::vector<WCSimLikelihoodTrackBase*> GetBestFit();

        /**
         * Create a track object using the seed parameters
         * @return Track object at the seed value
         */
        WCSimLikelihoodTrackBase* GetSeedParams();

        /**
         * As a temporary measure to get results out of the fit, we'll perform
         * a grid search to determine the energy, instead of turning Minuit loose on it
         * (until we have a working method that interpolates for a smooth surface)
         * @param best2LnL Will be updated with the best-fit -2Ln(likelihood)
         * @param bestEnergies Will be updated with the best-fit energies for each independent track
         */
        void PerformEnergyGridSearch(Double_t &best2LnL, std::vector<Double_t> &bestEnergies);


        void FillPlots();
        void FillTree();
        void FillHitComparison();

        void Make1DSurface(std::pair<unsigned int, FitterParameterType::Type> trackPar);
        void Make2DSurface(std::pair<std::pair<unsigned int, FitterParameterType::Type>, std::pair<unsigned int, FitterParameterType::Type> > trackPar);

        Bool_t GetTrueTrackEscapes(unsigned int iTrack) const;
        Bool_t GetFitTrackEscapes( unsigned int iTrack) const;
        Bool_t GetTrackEscapes(WCSimLikelihoodTrackBase * track) const;
        Bool_t GetUsePiZeroMassConstraint() const;

        Double_t GetPiZeroSecondTrackEnergy(const Double_t * x);


        WCSimTotalLikelihood * fTotalLikelihood; ///< Class used to calculate the total (combined charge and time) likelihood that we minimize
        WCSimRootEvent * fRootEvent; ///< Simulated event from WCSim to be reconstructed
        WCSimLikelihoodDigitArray * fLikelihoodDigitArray; ///< Array of PMT responses
        std::vector<TrackType::Type> fTypes; ///< Particle type to fit for
        std::map<Int_t, Int_t> fParMap; ///< Map relating the number of tracks (key) to the number of Minuit parameters (value)
        std::vector<WCSimLikelihoodTrackBase*> fBestFit;  ///< The best-fit track or tracks
        std::vector<WCSimLikelihoodTrackBase*> * fTrueLikelihoodTracks; ///< The MC-truth tracks for this event

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
        WCSimLikelihoodTrackBase * RescaleParams(Double_t x, Double_t y, Double_t z, Double_t t,
                                           Double_t th, Double_t phi, Double_t E,
                                           TrackType::Type type);
        

        Bool_t fUseHoughFitterForSeed;

        Int_t fEvent; ///< Number of the event currently being fitted
        Double_t fMinimum; ///< Value of -2 log(likelihood) at the best-fit point

        Bool_t fIsFirstCall; ///< Flags whether this is the first time the minimizer had calculated a likelihood (to print the seed)
        Int_t fStatus; ///< Minimizer convergence status


        WCSimFitterPlots * fFitterPlots;
        WCSimFitterTree * fFitterTree;
        WCSimFitterTrackParMap fFitterTrackParMap;
        WCSimFitterConfig * fFitterConfig;



        int fCalls;
    ClassDef(WCSimLikelihoodFitter,1);
};

#endif // WCSIMLIKELIHOODFITTER_H
