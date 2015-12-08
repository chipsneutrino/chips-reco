/**
 * \class WCSimLikelihoodTuner
 * This class is used to calculate a number of lower level
 * functions that contribute to the charge likelihood.
 * Ultimately these combine to give three integrals and
 * three corresponding coefficients to premultiply them
 * for direct and indirect light, which WCSimChargeLikelihood
 * uses to predict the mean number of photons at a PMT.
 */

#ifndef WCSIMLIKELIHOODTUNER_H
#define WCSIMLIKELIHOODTUNER_H

#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TObjArray.h"
#include "TString.h"
#include "TTree.h"

#include "WCSimEmissionProfileManager.hh"
#include "WCSimIntegralLookup.hh"
#include "WCSimLikelihoodTrackBase.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimRootGeom.hh"
#include "WCSimTrackParameterEnums.hh"
#include "WCSimPMTManager.hh"

class WCSimLikelihoodTunerCache
{
public:
  WCSimLikelihoodTunerCache();
  WCSimLikelihoodTunerCache( double const &s, WCSimLikelihoodDigit const * digit, WCSimLikelihoodTrackBase const * track);
  bool Contains(double const &s, WCSimLikelihoodDigit const * digit, WCSimLikelihoodTrackBase const * track);

  TVector3 GetEmissionPos();
  TVector3 GetToPMT();
  double GetAngleToPMT();
  double GetCosAngleToPMT();
  double GetDistanceToPMT();
  WCSimLikelihoodDigit const * GetDigit();
  double GetS();

private:
  WCSimLikelihoodDigit const * fDigit;
  WCSimLikelihoodTrackBase const * fTrack;
  double fS;
  TVector3 fEmissionPos;
  TVector3 fToPMT;
  double fAngleToPMT;
  double fCosAngleToPMT;
  double fDistanceToPMT;

};


class WCSimLikelihoodTuner
{
    public:
        /// Default constructor for an unknown detector type.  Need to UpdateDigitArray before using
        WCSimLikelihoodTuner();

        /**
         * Constructor using the detector settings from a given event
         * @param myDigitArray PMT response for an event
         */
        WCSimLikelihoodTuner(WCSimLikelihoodDigitArray * myDigitArray, WCSimEmissionProfileManager * myEmissionProfileManager);

        /**
         * Called by both constructors to set up internal variables
         * @todo Teach it to dynamically recognise whether we have a muon/electron and load the appropriate file
         */
        void Initialize();
        virtual ~WCSimLikelihoodTuner();

        /**
         * Load in settings from a new event
         * @param myDigitArray PMT response from an event
         */
        void UpdateDigitArray(WCSimLikelihoodDigitArray * myDigitArray);
   
        /**
         * Probability a photon survives the distance from where it was emitted to the PMT
         * @param s Distance from track vertex at which the photon was emitted
         * @param myTrack Track that emitted the photon
         * @param myDigit PMT the photon is travelling to
         * @return Probability the photon survives to the PMT without being absorbed
         */
        Double_t TransmissionFunction();

        /**
         * Probability that the PMT registers a hit as a function of the incoming photon angle
         * @param s Distance from track vertex at which the photon was emitted
         * @param myTrack Track that emitted the photon
         * @param myDigit PMT the photon is travelling to
         * @return Probability that the PMT registers a hit due to the photon
         */
        Double_t Efficiency(Double_t s, WCSimLikelihoodTrackBase * myTrack, WCSimLikelihoodDigit * myDigit);

        /**
         * Apply the wavelength-weighted relative quantum efficiency
         * @param myTrack Track that emitted the photon
         */
        Double_t QuantumEfficiency(const double &s, WCSimLikelihoodTrackBase * myTrack, WCSimLikelihoodDigit * myDigit);

        /**
         * Factor for the solid angle of the PMT as viewed from where the photon was emitted
         * @param s Distance from track vertex at which the photon was emitted
         * @param myTrack Track that emitted the photon
         * @param myDigit PMT the photon is travelling to
         * @return Total (not fractional) solid angle of the photon
         */
        Double_t SolidAngleFraction(Double_t s, WCSimLikelihoodTrackBase * myTrack, WCSimLikelihoodDigit * myDigit);

        /**
         * Relates the amount of direct and indirect light as a function of several geometric factors
         * Currently just a constant
         * @todo Implement this fully
         * @param s Distance from track vertex
         * @return
         */
        Double_t ScatteringTable(Double_t s);

        /**
         * Get efficiency * transmission * solid angle for direct and indirect photons
         * @param s Distance from track vertex at which the photon was emitted
         * @param myTrack Track that emitted the photon
         * @param myDigit PMT the photon is travelling to
         * @return Vector of efficiency * transmission * solid angle: [0] is direct and [1] is indirect
         */
        std::vector<Double_t> CalculateJ( Double_t s, WCSimLikelihoodTrackBase * myTrack, WCSimLikelihoodDigit * myDigit);
        
        /**
         * How far along the track should we integrate to, before the particle escapes a cylindrical detector?
         * @param myTrack Track object to calculate it for
         * @return Distance in the track direction where the integral should stop
         */
        Double_t CalculateCylinderCutoff(WCSimLikelihoodTrackBase * myTrack);

        /**
         * How far along the track should we integrate to, before the particle escapes a mailbox detector?
         * @param myTrack Track object to calculate it for
         * @return Distance in the track direction where the integral should stop
         */
        Double_t CalculateMailBoxCutoff(WCSimLikelihoodTrackBase * myTrack);

        /**
         * Wrapper to calculate the cutoff for the appropriate geometry
         * @param myTrack Track object to calculate it for
         */
        void CalculateCutoff(WCSimLikelihoodTrackBase * myTrack);

        // DEBUGGING:
        // By default the tuner will use a config file to determine whether to calculate or lookup the integrals
        // This funtion can be used to override it
        /**
         * DEBUGGING
         * By default the tuner uses a config file to determine whether to calculate or look up the integrals
         * This can be used to override that
         * @param calc True to force integrals to be calculated, false to force lookup table
         */
        void SetCalculateIntegrals(Bool_t calc);

        Bool_t GetCalculateIntegrals() const;

        /**
         * Wrapper to get the direct integrals.  Will either calculate or look them up depending on fCalculateIntegrals
         * @param myTrack Track to get the integrals for
         * @param myDigit PMT being considered
         * @param sPower Power of s under the integral (0,1,2)
         * @return Integral for a chosen power of s
         */
        Double_t GetChIntegrals(WCSimLikelihoodTrackBase * myTrack, WCSimLikelihoodDigit * myDigit, int sPower); 

        /**
         * Wrapper to get the direct integrals.  Will either calculate or look them up depending on fCalculateIntegrals
         * @param myTrack Track to get the integrals for
         * @param myDigit PMT being considered
         * @return Integrals for all three powers of s
         */
        std::vector<Double_t> GetChIntegrals(WCSimLikelihoodTrackBase * myTrack, WCSimLikelihoodDigit * myDigit); 


        /**
         * Wrapper to get the indirect integrals.  Will either calculate or look them up depending on fCalculateIntegrals
         * @param myTrack Track to get the integrals for
         * @param myDigit PMT being considered
         * @param sPower Power of s under the integral (0,1,2)
         * @return Integral for a chosen power of s
         */
        Double_t GetIndIntegrals(WCSimLikelihoodTrackBase * myTrack, Int_t sPower); 

        /**
         * Wrapper to get the direct integrals.  Will either calculate or look them up depending on fCalculateIntegrals
         * @param myTrack Track to get the integrals for
         * @param myDigit PMT being considered
         * @return Integrals for all three powers of s
         */
        std::vector<Double_t> GetIndIntegrals(WCSimLikelihoodTrackBase * myTrack); 
        
        // Get the integrals using a lookup table
        /// Get the direct integrals by looking them up, for a given power of s
        Double_t LookupChIntegrals(WCSimLikelihoodTrackBase * myTrack, WCSimLikelihoodDigit * myDigit, int sPower); 

        /// Get the direct integrals by looking them up, for all three powers of s
        std::vector<Double_t> LookupChIntegrals(WCSimLikelihoodTrackBase * myTrack, WCSimLikelihoodDigit * myDigit); 

        /// Get the indirect integrals by looking them up, for a given power of s
        Double_t LookupIndIntegrals(WCSimLikelihoodTrackBase * myTrack, Int_t sPower);

        /// Get the indirect integrals by looking them up, for all three powers of s
        std::vector<Double_t> LookupIndIntegrals(WCSimLikelihoodTrackBase * myTrack); 

        /// Get the direct integrals by calculating them numerically, for a given power of s
        Double_t CalculateChIntegrals(WCSimLikelihoodTrackBase * myTrack, WCSimLikelihoodDigit * myDigit, int i);

        /// Get the direct integrals by calculating them numerically, for all powers of s
        std::vector<Double_t> CalculateChIntegrals(WCSimLikelihoodTrackBase * myTrack, WCSimLikelihoodDigit * myDigit);

        /// Get the indirect integrals by calculating them numerically, for a given power of s
        Double_t CalculateIndIntegrals(WCSimLikelihoodTrackBase * myTrack, int i);

        /// Get the indirect integrals by calculating them numerically, for all powers of s
        std::vector<Double_t> CalculateIndIntegrals(WCSimLikelihoodTrackBase * myTrack);



        /**
         * Find the components of the quadratic expansion we use to multiply the integrals
         * @param s Distance from track vertex at which the photon was emitted
         * @param myTrack Track that emitted the photon
         * @param myDigit PMT the photon is travelling to
         */
        void CalculateCoefficients(WCSimLikelihoodTrackBase * myTrack, WCSimLikelihoodDigit * myDigit);

        /**
         * As with WCSimLikelihoodTuner::CalculateCoefficients but returns the coefficients too
         * rather than just setting the member variables
         * @return Vector of the quadratic expansion coefficients: [0-2] direct, [3-5] indirect
         */
        std::vector<Double_t> CalculateCoefficientsVector(WCSimLikelihoodTrackBase * myTrack, WCSimLikelihoodDigit * myDigit);
   
        Double_t GetLightFlux(WCSimLikelihoodTrackBase * myTrack);        

        double GetTrackLengthForPercentile(WCSimLikelihoodTrackBase * myTrack, const double &percentile);

        double GetCutoff(WCSimLikelihoodTrackBase * myTrack);

    protected:
    private:

        /**
         * Check is a position is inside or outside the detector
         * @param pos Position to check
         * @return True if pos is outside the detector, false if inside
         */
        Bool_t IsOutsideDetector(const TVector3 &pos);

        Double_t fDirCoeffs[3]; ///< Quadratic expansion coefficients to multiply the track integrals, for direct light
        Double_t fIndCoeffs[3]; ///< Quadratic expansion coefficients to multiply the track integrals, for indirect light
		    Double_t fAverageQE;    ///< Quantum efficiency averaged over all emitted wavelenghts (not used any more)

        Double_t fExtent[3];    ///< Dimensions of the detector: (x,y,z) for mailbox, (r,r,z) for cylinder
        Bool_t fConstrainExtent; ///< Should we stop integrating when we reach the edge of the detector?
        double fMaxDistance;  ///< The furthest a photon can travel to hit a PMT (i.e. distance from the bottom left to top right of the detector)
        WCSimLikelihoodDigitArray::GeomType_t fGeomType; ///< Detector geometry type


        // Work out where to cut off integrals, and cache the last one
        WCSimLikelihoodTrackBase * fLastCutoff; ///< The last track we calculated a cutoff for (for cacheing)
        Double_t fCutoffIntegral; ///< Distance along the particle's trajectory where it leaves the detector

	  	  Bool_t   fCalculateIntegrals;  ///< True if we should calculate integrals numerically, false to look them up in a table
      
	  	  TrackType::Type fIntegralParticleType; ///< The particle type whose table we've already loaded
	  	  WCSimEmissionProfileManager * fEmissionProfileManager; ///< The emission profile handler
	  	  WCSimPMTManager * fPMTManager;

        void MakeCache(const double &s , WCSimLikelihoodTrackBase * myTrack, WCSimLikelihoodDigit * myDigit);
        WCSimLikelihoodTunerCache fCache;




};

#endif // WCSIMLIKELIHOODTUNER_H
