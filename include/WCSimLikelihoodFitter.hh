/**
 * \class WCSimLikelihoodFitter
 * This class is used to vary the hypothesized tracks
 * and perform the minimization of the resulting negative
 * log likelihood to reconstruct the event.
 */

#ifndef WCSIMLIKELIHOODFITTER_H
#define WCSIMLIKELIHOODFITTER_H

#include "WCSimRootEvent.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimReco.hh"
#include "WCSimRecoDigit.hh"
#include "WCSimLikelihoodTrack.hh"
#include "WCSimTotalLikelihood.hh"
#include <vector>
#include <map>

//class WCSimChargeLikelihood;

class WCSimLikelihoodFitter
{
    public:
        /**
         * Constructor
         * @param Event object from WCSim to reconstruct
         */
        WCSimLikelihoodFitter( WCSimRootEvent*);
        virtual ~WCSimLikelihoodFitter();

        /**
         * Perform the minimization
         * @param nTracks Number of tracks to include in the fit
         */
        void Minimize2LnL(Int_t nTracks);

        /**
         * Work out how many free parameters Minuit needs to run over
         * @param nTracks Number of tracks to include in the fit
         * @return Number of parameters to use in Minuit
         */
        Int_t GetNPars(Int_t nTracks);

        /**
         * Wrapper to WCSimTotalLikelihood::Calc2LnL so that
         * we can use it as Minuit's FCN
         * @param x Array of parameters used to construct the track objects for the fit
         * @return -2 log(likelihood) for this set of trakcs
         */
        Double_t WrapFunc(const Double_t * x);

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
        void SeedParams( WCSimReco * myReco );

        /**
         * Get a vector containing the track(s) at the best-fit point
         * @return Best-fit track(s) for this event
         */
        std::vector<WCSimLikelihoodTrack> GetBestFit();

        /**
         * Create a track object using the seed parameters
         * @return Track object at the seed value
         */
        WCSimLikelihoodTrack GetSeedParams();

    protected:
    private:
        WCSimTotalLikelihood * fTotalLikelihood; ///< Class used to calculate the total (combined charge and time) likelihood that we minimize
        WCSimRootEvent * fRootEvent; ///< Simulated event from WCSim to be reconstructed
        WCSimLikelihoodDigitArray * fLikelihoodDigitArray; ///< Array of PMT responses
        WCSimLikelihoodTrack::TrackType fType; ///< Particle type to fit for
        std::map<Int_t, Int_t> fParMap; ///< Map relating the number of tracks (key) to the number of Minuit parameters (value)
        std::vector<WCSimLikelihoodTrack> fBestFit;  ///< The best-fit track or tracks

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
        WCSimLikelihoodTrack RescaleParams(Double_t x, Double_t y, Double_t z, Double_t t,
                                           Double_t th, Double_t phi, Double_t E,
                                           WCSimLikelihoodTrack::TrackType type);
        

        Double_t fMinimum; ///< Value of -2 log(likelihood) at the best-fit point
        Double_t fSeedVtxX; ///< Seed vertex x coordinate
        Double_t fSeedVtxY; ///< Seed vertex y coordinate
        Double_t fSeedVtxZ; ///< Seed vertex z coordinate
        Double_t fSeedTheta; ///< Seed angle to the z axis
        Double_t fSeedPhi; ///< Seed azimuthal angle
        Bool_t fIsFirstCall; ///< Flags whether this is the first time the minimizer had calculated a likelihood (to print the seed)
        Int_t fStatus; ///< Minimizer convergence status
	
    ClassDef(WCSimLikelihoodFitter,1);
};

#endif // WCSIMLIKELIHOODFITTER_H
