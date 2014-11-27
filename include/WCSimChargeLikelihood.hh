/**
 * \class WCSimChargeLikelihood
 *
 * This class is used to calculate the contribution to the
 * total log-likelihood due to the charges measured at
 * each PMT (including if a PMT measures no charge)
 */
#ifndef WCSIMCHARGELIKELIHOOD_H
#define WCSIMCHARGELIKELIHOOD_H

#include <vector>

#include "WCSimDigitizerLikelihood.hh"
#include "WCSimRecoDigit.hh"
#include "WCSimLikelihoodDigit.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimLikelihoodTrack.hh"
#include "WCSimLikelihoodTuner.hh"
#include "WCSimRecoEvent.hh"

#include "TFile.h"
#include "TH2D.h"
#include "TMath.h"
#include "TTree.h"

class WCSimChargeLikelihood
{
    public:
        /**
         * Constructor
         * @param myDigitArray The PMT responses for a single event
         */
        WCSimChargeLikelihood( WCSimLikelihoodDigitArray * myDigitArray);

        virtual ~WCSimChargeLikelihood();

        /**
         * Add another track to calculate the likelihood for several particles at once
         * @param myTrack Track object to add
         */
        void AddTrack( WCSimLikelihoodTrack * myTrack);

        /**
         * Set all the tracks that contribute to the likelihood at once
         * @param myTrack Vector of all the track objects to consider
         */
        void SetTracks( std::vector<WCSimLikelihoodTrack*> myTrack );

        /// Remove all the tracks currently loaded
        void ClearTracks();

        /**
         * Replace the current event with hits from a different one
         * @param myDigitArray New array of PMT responses
         */
        void UpdateDigitArray( WCSimLikelihoodDigitArray * myDigitArray);

        /**
         * Calculate -2 log(likelihood) for the current PMT response, given the current set of hypothesised tracks
         * @return -2 log(likelihood)
         */
        Double_t Calc2LnL();

        /**
         * For debugging: force the integrals to be calculated iteratively
         * Makes no permanent changes to how the likelihood is calculated in the future
         * @return -2 log(likelihood) with the integrals calculated analytically, not looked-up
         */
        Double_t CalculateExactLikelihood(); // Calc2LnL can (optionally) use tabulated integrals, this always calculates them by hand

        /**
         * Calculate the predicted mean number of photons arriving at the PMT specified by
         * WCSimChargeExpectation::fDigit due to an individual particle track
         * @param myTrack The charged particle track
         * @return Predicted number of photons at the PMT originating from the track
         */
        Double_t ChargeExpectation(WCSimLikelihoodTrack * myTrack);

        /**
         * Calculate the predicted mean number of photons arriving at the specified PMT
         * due to an individual particle track. For use by the time likelihood
         * @param myTrack The charged particle track
         * @param myDigit The PMT in question
         * @return Predicted number of photons at the PMT originating from the track
         */
        Double_t DigitChargeExpectation(WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit);

        /**
         * Getter for total number of photons emitted by track over its whole length
         * @param myTrack The charged particle track
         * @return Average total number of photons emitted
         */
        Double_t GetLightFlux(WCSimLikelihoodTrack * myTrack);

        /**
         * Get the predicted number of indirect (scattered, reflected etc.) photons
         * striking the PMT in WCSimChargeLikelihood::fDigit originating from a given track
         * @param myTrack The charged particle track
         * @return Predicted number of indirect photons hitting the PMT
         */
        Double_t GetMuIndirect(WCSimLikelihoodTrack * myTrack);

        /**
         * Get the predicted number of direct (unscattered Cherenkov) photons
         * striking the PMT in WCSimChargeLikelihood::fDigit originating from a given track
         * @param myTrack The charged particle track
         * @return Predicted number of direct photons hitting the PMT
         */
        Double_t GetMuDirect(WCSimLikelihoodTrack * myTrack);

        /**
         * Set the coefficients of the quadratic expansion to several geometric
         * factors that pre-multiply the integrals along the track length
         * @param myTrack The charged particle track
         */
        void GetTrackParameters(WCSimLikelihoodTrack * myTrack);

    protected:
    private:
        /**
         * Called by constructor.  Initializes vectors to hold the track integrals, and
         * the WCSimLikelihoodTuner and WCSimDigitizerLikelihood member variables
         * @param myDigitArray PMT responses for this event
         */
        void Initialize( WCSimLikelihoodDigitArray * myDigitArray);

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
        std::vector<WCSimLikelihoodTrack *>   fTracks;     ///< Vector of simultaneous tracks contributing to the likelihood
        WCSimLikelihoodDigitArray           * fDigitArray; ///< Response for all the detector's PMTs for this event
        WCSimLikelihoodDigit                * fDigit;      ///< The PMT being considered
        WCSimLikelihoodTuner                * fTuner;      ///< Class to calculate effects of absorption, geometry, etc. and perform integrals
        WCSimDigitizerLikelihood            * fDigitizer;  ///< Class to calculate the likelihood of getting the measured PMT hit given the predicted value

      // The fitted functions are defined using various variables that relate the track to the
      // PMT hit in question.  I calculate these in GetTrackParameters() and set a flag when
      // this has been done
      Bool_t fGotTrackParameters; ///< Flag to check if GetTrackParameters() has been called to calculate the integral coefficients
      Int_t fNumCalculations;



};

#endif // WCSIMCHARGELIKELIHOOD_H
