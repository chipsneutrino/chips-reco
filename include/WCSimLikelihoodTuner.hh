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

#include "WCSimEmissionProfiles.hh"
#include "WCSimLikelihoodTrack.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimRootGeom.hh"



class WCSimLikelihoodTuner
{
    public:
        /// Default constructor for an unknown detector type.  Need to UpdateDigitArray before using
        WCSimLikelihoodTuner();

        /**
         * Constructor using the detector settings from a given event
         * @param myDigitArray PMT response for an event
         */
        WCSimLikelihoodTuner(WCSimLikelihoodDigitArray * myDigitArray);

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
         * Get the emission profiles (fraction of photons emitted in a given direction at a given
         * distance along the track) for a specific particle type
         * @param myTrack Track object to load profiles for
         */
        void LoadEmissionProfiles(WCSimLikelihoodTrack * myTrack);


        /**
         * Probability a photon survives the distance from where it was emitted to the PMT
         * @param s Distance from track vertex at which the photon was emitted
         * @param myTrack Track that emitted the photon
         * @param myDigit PMT the photon is travelling to
         * @return Probability the photon survives to the PMT without being absorbed
         */
        Double_t TransmissionFunction(Double_t s, WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit);

        /**
         * Probability that the PMT registers a hit as a function of the incoming photon angle
         * @param s Distance from track vertex at which the photon was emitted
         * @param myTrack Track that emitted the photon
         * @param myDigit PMT the photon is travelling to
         * @return Probability that the PMT registers a hit due to the photon
         */
        Double_t Efficiency(Double_t s, WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit);

        // Double_t QuantumEfficiency(WCSimLikelihoodTrack * myTrack);

        /**
         * Factor for the solid angle of the PMT as viewed from where the photon was emitted
         * @param s Distance from track vertex at which the photon was emitted
         * @param myTrack Track that emitted the photon
         * @param myDigit PMT the photon is travelling to
         * @return Total (not fractional) solid angle of the photon
         */
        Double_t SolidAngle(Double_t s, WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit);

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
        std::vector<Double_t> CalculateJ( Double_t s, WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit);
        
        /**
         * How far along the track should we integrate to, before the particle escapes a cylindrical detector?
         * @param vtx Track vertex
         * @param dir Track direction unit vector
         * @return Distance in the track direction where the integral should stop
         */
        Double_t CalculateCylinderCutoff(const TVector3 &vtx, const TVector3 &dir);

        /**
         * How far along the track should we integrate to, before the particle escapes a mailbox detector?
         * @param vtx Track vertex
         * @param dir Track direction unit vector
         * @return Distance in the track direction where the integral should stop
         */
        Double_t CalculateMailBoxCutoff(const TVector3 &vtx, const TVector3 &dir);

        /**
         * Wrapper to calculate the cutoff for the appropriate geometry
         * @param myTrack Track object to calculate it for
         */
        void CalculateCutoff(WCSimLikelihoodTrack * myTrack);

        /**
         * Load the lookup table containing the track integrals
         * @param myTrack Track to load the table for
         */
	      void LoadTabulatedIntegrals(WCSimLikelihoodTrack * myTrack);

	    // Bin lookups:
	    Int_t GetEBin(Double_t energy);
        Int_t GetSBin(Double_t sMax);
        Int_t GetESBin(Double_t energy, Double_t sMax);
        Int_t GetIntegralBin(Double_t energy, Double_t R0, Double_t cosTheta0, Double_t s);
       
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
        Double_t GetChIntegrals(WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit, int sPower); 

        /**
         * Wrapper to get the direct integrals.  Will either calculate or look them up depending on fCalculateIntegrals
         * @param myTrack Track to get the integrals for
         * @param myDigit PMT being considered
         * @return Integrals for all three powers of s
         */
        std::vector<Double_t> GetChIntegrals(WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit); 


        /**
         * Wrapper to get the indirect integrals.  Will either calculate or look them up depending on fCalculateIntegrals
         * @param myTrack Track to get the integrals for
         * @param myDigit PMT being considered
         * @param sPower Power of s under the integral (0,1,2)
         * @return Integral for a chosen power of s
         */
        Double_t GetIndIntegrals(WCSimLikelihoodTrack * myTrack, Int_t sPower); 

        /**
         * Wrapper to get the direct integrals.  Will either calculate or look them up depending on fCalculateIntegrals
         * @param myTrack Track to get the integrals for
         * @param myDigit PMT being considered
         * @return Integrals for all three powers of s
         */
        std::vector<Double_t> GetIndIntegrals(WCSimLikelihoodTrack * myTrack); 
        
        // Get the integrals using a lookup table
        /// Get the direct integrals by looking them up, for a given power of s
        Double_t LookupChIntegrals(WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit, int sPower); 

        /// Get the direct integrals by looking them up, for all three powers of s
        std::vector<Double_t> LookupChIntegrals(WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit); 

        /// Get the indirect integrals by looking them up, for a given power of s
        Double_t LookupIndIntegrals(WCSimLikelihoodTrack * myTrack, Int_t sPower);

        /// Get the indirect integrals by looking them up, for all three powers of s
        std::vector<Double_t> LookupIndIntegrals(WCSimLikelihoodTrack * myTrack); 

        /// Get the direct integrals by calculating them numerically, for a given power of s
        Double_t CalculateChIntegrals(WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit, int i);

        /// Get the direct integrals by calculating them numerically, for all powers of s
        std::vector<Double_t> CalculateChIntegrals(WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit);

        /// Get the indirect integrals by calculating them numerically, for a given power of s
        Double_t CalculateIndIntegrals(WCSimLikelihoodTrack * myTrack, int i);

        /// Get the indirect integrals by calculating them numerically, for all powers of s
        std::vector<Double_t> CalculateIndIntegrals(WCSimLikelihoodTrack * myTrack);



        /**
         * Find the components of the quadratic expansion we use to multiply the integrals
         * @param s Distance from track vertex at which the photon was emitted
         * @param myTrack Track that emitted the photon
         * @param myDigit PMT the photon is travelling to
         */
        void CalculateCoefficients(WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit);

        /**
         * As with WCSimLikelihoodTuner::CalculateCoefficients but returns the coefficients too
         * rather than just setting the member variables
         * @return Vector of the quadratic expansion coefficients: [0-2] direct, [3-5] indirect
         */
        std::vector<Double_t> CalculateCoefficientsVector(WCSimLikelihoodTrack * myTrack, WCSimLikelihoodDigit * myDigit);
   
        /**
         * Produce a table of direct integrals that can later be loaded and used to look them up
         * @param myType Particle type
         * @param filename Where to save the tables
         */
        void TabulateIndirectIntegrals( WCSimLikelihoodTrack::TrackType myType, TString filename);

        /**
         * Produce a table of indirect integrals that can later be loaded and used to look them up
         * @param myType Particle type
         * @param filename Where to save the tables
         */
        void TabulateDirectIntegrals(WCSimLikelihoodTrack::TrackType myType, TString filename);

        /**
         * Produce tables of direct and indirect integrals that can later be loaded and used to look them up
         * @param myType Particle type
         * @param filename Where to save the tables
         */
        void TabulateIntegrals(WCSimLikelihoodTrack::TrackType myType, TString filename);

        Double_t GetLightFlux(WCSimLikelihoodTrack * myTrack);        
    
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
        WCSimLikelihoodDigitArray::GeomType_t fGeomType; ///< Detector geometry type


        // Work out where to cut off integrals, and cache the last one
        WCSimLikelihoodTrack * fLastCutoff; ///< The last track we calculated a cutoff for (for cacheing)
        Double_t fCutoffIntegral; ///< Distance along the particle's trajectory where it leaves the detector

        // Read in from the emission profile files

	  
	  	// These are used to speed up reading in the pre-computed integrals
	  	// by ensuring that if the appropriate file is already open, and the appropriate
	  	// tree has already been loaded, we don't waste time getting them again
        // The way integrals are loaded has changed to just use a binary file not a tree
        // - some of these should be safe to delete now
	  	Bool_t   fCalculateIntegrals;  ///< True if we should calculate integrals numerically, false to look them up in a table
	  	TFile  * fRhoIntegralFile;     ///< File containing the integrals of the 1D emission profile
	  	TFile  * fRhoGIntegralFile;    ///< File containing the integrals of the 2D emission profile
	  	TTree  * fRhoIntegralTree;     ///< Not used any more
	  	TTree  * fRhoGIntegralTree;    ///< Not used any more
        Int_t    fIntegralEnergyBin;   ///< Not used any more
        Double_t fIntegralSMax;        ///< Not used any more
      
      
      // The binning scheme for the table of indirect integrals      
      Int_t    fNBinsRho;  ///< Total number of bins in all indirect integrals at all energies
      Int_t    fNSBinsRho; ///< How many bins of s are there for each energy in the indirect table
      Int_t    fNEBinsRho; ///< How different energies are there tables for
      Double_t fSMinRho;   ///< Lowest distance along track (should be 0)
      Double_t fSMaxRho;   ///< Highest distance along track
      Double_t fEMinRho;   ///< Lowest tabulated energy
      Double_t fEMaxRho;   ///< Highest tabulated energy


      // The binning scheme for the table of direct integrals
      // Co-ordinates are defined in WCSimLikelihoodTuner.cc and follow the MiniBooNE paper
      Int_t    fNBinsRhoG;      ///< Total number of bins in all direct integrals at all energies
      Int_t    fNR0Bins;        ///< Number of bins in R0
      Int_t    fNCosTheta0Bins; ///< Number of bins in cosTheta0
	  Int_t    fNSBins;         ///< Number of bins in distance along track
	  Int_t    fNEBins;         ///< Number of track energy bins
	  Double_t fR0Min;          ///< Minimum value of R0 (should be 0)
	  Double_t fR0Max;          ///< Maximum value of R0
	  Double_t fEMin;           ///< Minimum tabulated energy
	  Double_t fEMax;           ///< Maximum tabulated energy
      Double_t fSMin;           ///< Lowest distance along track (should be 0)
	  Double_t fSMax;           ///< Highest distance along track
      Double_t fCosTheta0Min;   ///< Minimum value of cosTheta0 (should be -1)
      Double_t fCosTheta0Max;   ///< Maximum value of cosTheta0 (should be +1)
      
      // The table of integrals
	  Double_t * fRhoIntegrals;  ///< Table of indirect integrals flattened to a 1D array
	  Double_t * fRhoGIntegrals; ///< Table of direct integrals flattened to a 1D array

	  WCSimLikelihoodTrack::TrackType fIntegralParticleType; ///< The particle type whose table we've already loaded
    WCSimEmissionProfiles fEmissionProfiles; ///< The emission profile handler 	  
};

#endif // WCSIMLIKELIHOODTUNER_H
