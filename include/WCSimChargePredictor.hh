/**
 * \class WCSimChargePredictor
 *
 * This class is used to calculate the predicted
 * mean charge at a PMT due to a given track hypothesis
 */
#ifndef WCSIMCHARGEPREDICTOR_HH
#define WCSIMCHARGEPREDICTOR_HH

#include <vector>

#include "WCSimDigitizerLikelihood.hh"
#include "WCSimRecoDigit.hh"
#include "WCSimLikelihoodDigit.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimLikelihoodTrackBase.hh"
#include "WCSimLikelihoodTuner.hh"
#include "WCSimRecoEvent.hh"

#include "TFile.h"
#include "TH2D.h"
#include "TMath.h"
#include "TTree.h"
class WCSimEmissionProfilesManager;

class WCSimChargePredictor
{
    public:
        /**
         * Constructor
         * @param myDigitArray The PMT responses for a single event
         */
        WCSimChargePredictor( WCSimLikelihoodDigitArray * myDigitArray, WCSimEmissionProfileManager * myEmissionProfileManager);

        //ROOT requires a default ctor to generate dictionary
        //for a vector of charge likelihood objects - do not use
        WCSimChargePredictor();

        ///Copy constructor
        WCSimChargePredictor(const WCSimChargePredictor &other);
        ///Assignment operator
        WCSimChargePredictor& operator= (const WCSimChargePredictor &rhs);

        virtual ~WCSimChargePredictor();

        /**
         * Add another track to calculate the likelihood for several particles at once
         * @param myTrack Track object to add
         */
        void AddTrack( WCSimLikelihoodTrackBase * myTrack);

        /**
         * Set all the tracks that contribute to the likelihood at once
         * @param myTrack Vector of all the track objects to consider
         */
        void SetTracks( std::vector<WCSimLikelihoodTrackBase*> myTrack );

        /// Remove all the tracks currently loaded
        void ClearTracks();

        /**
         * Replace the current event with hits from a different one
         * @param myDigitArray New array of PMT responses
         */
        void UpdateDigitArray( WCSimLikelihoodDigitArray * myDigitArray);

        /**
         * Calculate the predicted mean number of photons at a PMT given the hypothesised track
         * @param myDigit Current PMT to give result for
         * @return Predicted mean number of photons at this PMT
         */
        Double_t GetPredictedCharge(WCSimLikelihoodDigit *myDigit);

        /**
         * Calculate the predicted mean number of photons at a PMT given the hypothesised track
         * @param digitNum ID number of the desired PMT
         * @return -2 log(likelihood)
         */
        Double_t GetPredictedCharge(Int_t digitNum);

        /**
         * For debugging: force the integrals to be calculated iteratively
         * Makes no permanent changes to how the likelihood is calculated in the future
         * @return -2 log(likelihood) with the integrals calculated analytically, not looked-up
         */
        //FIXME: this won't end well
        Double_t CalculateExactLikelihood(WCSimLikelihoodDigit *myDigit); // Calc2LnL can (optionally) use tabulated integrals, this always calculates them by hand

        /**
         * Calculate the predicted mean number of photons arriving at the PMT specified by
         * WCSimChargeExpectation::fDigit due to an individual particle track
         * @param trackIndex Index of the charged particle track in fTracks vector
         * @return Predicted number of photons at the PMT originating from the track
         */
        Double_t ChargeExpectation(Int_t trackIndex);

        /**
         * Calculate the predicted mean number of photons arriving at the specified PMT
         * due to an individual particle track. For use by the time likelihood
         * @param trackIndex Index of the charged particle track in fTracks vector
         * @param myDigit The PMT in question
         * @return Predicted number of photons at the PMT originating from the track
         */
        Double_t ChargeExpectation(Int_t trackIndex, WCSimLikelihoodDigit *myDigit);

        /**
         * Getter for total number of photons emitted by track over its whole length
         * @param trackIndex Index of the charged particle track in fTracks vector
         * @return Average total number of photons emitted
         */
        Double_t GetLightFlux(Int_t trackIndex);

        /**
         * Get the predicted number of indirect (scattered, reflected etc.) photons
         * striking the PMT in WCSimChargeLikelihood::fDigit originating from a given track
         * @param trackIndex Index of the charged particle track in fTracks vector
         * @return Predicted number of indirect photons hitting the PMT
         */
        Double_t GetMuIndirect(Int_t trackIndex);

        /**
         * Get the predicted number of direct (unscattered Cherenkov) photons
         * striking the PMT in WCSimChargeLikelihood::fDigit originating from a given track
         * @param trackIndex Index of the charged particle track in fTracks vector
         * @return Predicted number of direct photons hitting the PMT
         */
        Double_t GetMuDirect(Int_t trackIndex);

        /**
         * Set the coefficients of the quadratic expansion to several geometric
         * factors that pre-multiply the integrals along the track length
         * @param trackIndex Index of the charged particle track in fTracks vector
         */
        void GetTrackParameters(Int_t trackIndex);

        static double GetMinAllowed() { return 1e-6; }

    protected:
    private:
        /**
         * Called by constructor.  Initializes vectors to hold the track integrals, and
         * the WCSimLikelihoodTuner and WCSimDigitizerLikelihood member variables
         * @param myDigitArray PMT responses for this event
         */
        void Initialize( WCSimLikelihoodDigitArray * myDigitArray, WCSimEmissionProfileManager * myEmissionProfileManager);

        /** Parameters from a parabolic fit to solid angle * transmission * PMT acceptance
         * as a function of distance travelled by particle.  This is for direct Cherenkov light,
         */
        std::vector<Double_t> fCoeffsCh;

        /** Parameters from a parabolic fit to solid angle * transmission * PMT acceptance
         * as a function of distance travelled by particle.  This is for indirect
         * scattered/reflected light,
         */
        std::vector<Double_t> fCoeffsInd;


        /*
         VARIABLE NAMES FOLLOW arXiv.0902.2222v2
        */
        Double_t fEnergy;    ///< For integral lookup tables: particle energy
        Double_t fR0;        ///< For integral lookup tables: vertex to PMT distance
        Double_t fCosTheta0; ///< For integral lookup tables: angle to PMT as viewed from vertex
        Double_t fEta;       ///< For integral lookup tables: angle of incidence of emitted light at the PMT

        Double_t fRadius;    ///< For scattering table: distance from centre of tank to the source
        Double_t fAngle;     ///< For scattering table: angle between source, centre of tank, and PMT
        Double_t fTheta;     ///< For scattering table: angle between source direction and ray to PMT (same as theta0)
        Double_t fPhi;       ///< For scattering table: Angle between plane containing tank centre, PMT and source, and plane containing the track and the tank centre

      // The track and event parameters for which we calculate the likelihood
        std::vector<WCSimLikelihoodTrackBase *>   fTracks;     ///< Vector of simultaneous tracks contributing to the likelihood
        WCSimLikelihoodDigitArray           * fDigitArray; ///< Response for all the detector's PMTs for this event
        WCSimLikelihoodDigit                * fDigit;      ///< The PMT being considered
        WCSimLikelihoodTuner                * fTuner;      ///< Class to calculate effects of absorption, geometry, etc. and perform integrals
        WCSimEmissionProfileManager * fEmissionProfileManager;

      // The fitted functions are defined using various variables that relate the track to the
      // PMT hit in question.  I calculate these in GetTrackParameters() and set a flag when
      // this has been done

      /// Flag to check if GetTrackParameters() has been called to calculate
      /// the integral coefficients. Shows index of the track in fTracks
      /// vector or is negative if track parameters have not been calculated
      Int_t fGotTrackParameters; 
	  Int_t fNumCalculations;



};

#endif // WCSIMCHARGEPREDICTOR_HH
