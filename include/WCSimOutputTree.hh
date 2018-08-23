#ifndef WCSIMOUTPUTTREE_HH
#define WCSIMOUTPUTTREE_HH

/*
 * WCSimOutputTree.hh
 *
 *  Created on: 23 Jan 2015
 *      Author: ajperch
 */

#include "WCSimLikelihoodTrackBase.hh"
#include "WCSimLikelihoodTrack.hh"
#include "WCSimHitComparison.hh"
#include "WCSimLikelihoodTrackFactory.hh"
#include "WCSimRecoRing.hh"
#include "WCSimFitterConfig.hh"

#include "TObject.h"
#include "TBranch.h"

#include <string>
#include <vector>

class TFile;
class TString;
class TTree;
class WCSimRootGeom;
class WCSimLikelihoodRecoEvent;
class WCSimLikelihoodDigitArray;
class WCSimTrueEvent;
class WCSimRootEvent;
class WCSimRecoSummary;

// TODO: Update the output tree to make it simpler + include the PID tree by default
// TODO: Add a print() function to each of the outputTree classes

/*
 * \class EventHeader
 * Holds metadata type data for the given event
 */
class EventHeader: public TObject {
	public:
		EventHeader();
		EventHeader(std::string inputFile, int inputEvent, bool failed, std::string recoType);
		EventHeader(const EventHeader& other);
		EventHeader& operator=(const EventHeader& rhs);
		~EventHeader();

		void Print();

		std::string GetUID() {
			return fUID;
		}
		std::string GetInputFile() {
			return fInputFile;
		}
		int GetEventNum() {
			return fInputEvent;
		}
		bool IsFailed() {
			return fFailed;
		}
		std::string GetRecoType() {
			return fRecoType;
		}

	private:
		void BuildUID();

		std::string fUID;			///< Unique Identification Number for a file/event combination
		std::string fInputFile;		///< Path to input WCSim file
		int fInputEvent;			///< Event number within fInputFile
		bool fFailed;				///< True if the event failed
		std::string fRecoType;		///< Type of fit conducted in std::string form

		ClassDef(EventHeader,1)
};

/*
 * \class ParameterInfo
 * Holds the WCSimParameters used in the reconstruction
 */
class ParameterInfo: public TObject {
	public:
		ParameterInfo();
		~ParameterInfo();

		void Print();

		// Slicer parameters
		Double_t GetSlicerClusterDistance() {
			return fSlicerClusterDistance;
		}

		UInt_t GetSlicerMinSize() {
			return fSlicerMinSize;
		}

		Double_t GetSlicerChargeCut() {
			return fSlicerChargeCut;
		}

		Double_t GetSlicerTimeCut() {
			return fSlicerTimeCut;
		}

		Bool_t IterateSlicing() {
			return fIterateSlicing;
		}

		// Veto slicing parameters
		Double_t GetVetoClusterDistance() {
			return fVetoClusterDistance;
		}

		UInt_t GetVetoMinSize() {
			return fVetoMinSize;
		}

		Double_t GetVetoMinChargeCut() {
			return fVetoMinChargeCut;
		}

		Double_t GetVetoMaxChargeCut() {
			return fVetoMaxChargeCut;
		}

		Double_t GetVetoTimeCut() {
			return fVetoTimeCut;
		}

		// Integral parameters
		Bool_t CalculateIntegrals() {
			return fCalculateIntegrals;
		}

		Bool_t TruncateIntegrals() {
			return fTruncateIntegrals;
		}

		Bool_t ConstrainExtent() {
			return fConstrainExtent;
		}

		// Likelihood tuner parameters
		Bool_t UseTransmission() {
			return fUseTransmission;
		}

		Bool_t UseAngularEfficiency() {
			return fUseAngularEfficiency;
		}

		Bool_t UseGlassCathodeReflection() {
			return fUseGlassCathodeReflection;
		}

		Bool_t UseScatteringTable() {
			return fUseScatteringTable;
		}

		Bool_t UseNewAngularEfficiency() {
			return fUseNewAngularEfficiency;
		}

		Bool_t UseTrackFit() {
			return fUseTrackFit;
		}

		// Fitter parameters
		Bool_t UseTime() {
			return fUseTime;
		}

		Bool_t UseCharge() {
			return fUseCharge;
		}

		Bool_t EqualiseChargeAndTime() {
			return fEqualiseChargeAndTime;
		}

		Bool_t SaveWCSimRootEvent() {
			return fSaveWCSimRootEvent;
		}

		std::string GetDigiType() {
			return fDigiType;
		}

		Bool_t SaveSeedInfo() {
			return fSaveSeedInfo;
		}

		Bool_t SaveStageInfo() {
			return fSaveStageInfo;
		}

		Bool_t SaveHitComparison() {
			return fSaveHitComparison;
		}

		Bool_t SaveParameterInfo() {
			return fSaveParameters;
		}

		// Speed of light parameters
		Bool_t UseCustomParticleSpeed() {
			return fUseCustomParticleSpeed;
		}

		Bool_t UseCustomSpeedOfLight() {
			return fUseCustomSpeedOfLight;
		}

		Bool_t UseFittedSpeedOfLight() {
			return fUseFittedSpeedOfLight;
		}

		Double_t GetCustomParticleSpeed() {
			return fCustomParticleSpeed;
		}

		Double_t GetCustomSpeedOfLight() {
			return fCustomSpeedOfLight;
		}

		Double_t GetFittedSpeedOfLight() {
			return fFittedSpeedOfLight;
		}

		// Not currently used parameters
		Bool_t GetUseSimpleTimeResolution() {
			return fUseSimpleTimeResolution;
		}

		Bool_t GetUseSimpleTimeSlew() {
			return fUseSimpleTimeSlew;
		}

		Bool_t GetUseSimpleRefractiveIndex() {
			return fUseSimpleRefractiveIndex;
		}

	private:
		// Slicer parameters
		Double_t fSlicerClusterDistance; 	///< Max distance in cm between hits in the slices
		UInt_t fSlicerMinSize; 				///< Minimum number of hits for a slice to be formed
		Double_t fSlicerChargeCut; 			///< Only consider digits above the charge threshold
		Double_t fSlicerTimeCut; 			///< Maximum gap allowed between hits when ordered by time in order to be associated with the previous hit.
		Bool_t fIterateSlicing; 			///< Iterate the charge cut to make as many slices as we have tracks

		// Veto slicing parameters
		Double_t fVetoClusterDistance; 		///< Max distance in cm between hits in the slices
		UInt_t fVetoMinSize; 				///< Minimum number of hits for a slice to be formed
		Double_t fVetoMinChargeCut; 		///< Only consider digits above the charge threshold, set initially
		Double_t fVetoMaxChargeCut; 		///< at the min value, and iterates up to the max value.
		Double_t fVetoTimeCut; 				///< Maximum gap allowed between hits when ordered by time in order to be associated with the previous hit.

		// Integral parameters
		Bool_t fCalculateIntegrals;     	///< True if charge likelihood should calculate integrals, false to look them up
		Bool_t fTruncateIntegrals; 			///< True if charge likelihood should use full lookup tables for integrals, false to use one that doesn't cut off
		Bool_t fConstrainExtent;            ///< True if integrals should cut off when the particle leaves the detector

		// Likelihood tuner parameters
		Bool_t fUseTransmission;            ///< True if we should account for absorption of photons in the water
		Bool_t fUseAngularEfficiency;       ///< True if we should account for the PMT efficiency as a function of angle
		Bool_t fUseGlassCathodeReflection;  ///< True if we should account for photons being reflected off the PMT glass
		Bool_t fUseScatteringTable;         ///< True if we should use the scattering table, false for flat 1% chance
		Bool_t fUseNewAngularEfficiency;         	///< True if we should use the new solid angle and angular efficiency
		Bool_t fUseTrackFit;         				///< True if we should use the muon track photon emission fit instead of the emission profiles

		// Fitter parameters
		Bool_t fUseTime;                    ///< True if we should include timing information in the likelihood
		Bool_t fUseCharge;                  ///< True if we should include charge information in the likelihood
		Bool_t fEqualiseChargeAndTime; 		///< After we've seeded the vertex, energy and time, we can scale the time likelihood to weight it the same as the charge
		Bool_t fSaveWCSimRootEvent; 		///< Whether to save a full copy of the WCSimRootEvent fitted, or just a link to the original file
		std::string fDigiType;            	///< Name of digitiser type to use in WCSimDigitizerLikelihood::DigiType_t
		Bool_t fSaveSeedInfo;			   	///< True if we should save the SeedInfo in the output file
		Bool_t fSaveStageInfo;			   	///< True if we should save the StageInfo in the output file
		Bool_t fSaveHitComparison;		   	///< True if we should save the HitComparison in the output file
		Bool_t fSaveParameters;				///< True if we should save the ParameterInfo in the output file

		// Speed of light parameters
		Bool_t fUseCustomParticleSpeed; 	///< Normally we assume particles travel at c - this allows us to switch and set it manually
		Bool_t fUseCustomSpeedOfLight; 		///< Normally we assume light travels at c/(average n) - this allows us to switch and set it manually
		Bool_t fUseFittedSpeedOfLight; 		///< Normally we assume light travels at c/(average n) - this uses a fitted speed instead
		Double_t fCustomParticleSpeed;      ///< The speed of the propagating particle as a fraction of c
		Double_t fCustomSpeedOfLight;       ///< The speed of the propagating particle as a fraction of c
		Double_t fFittedSpeedOfLight;       ///< The effective speed of light from our fit as a fraction of c

		// Not currently used parameters
		Bool_t fUseSimpleTimeResolution;	///<
		Bool_t fUseSimpleTimeSlew;			///<
		Bool_t fUseSimpleRefractiveIndex;	///<

		ClassDef(ParameterInfo,1)
};

/*
 * \class TruthInfo
 * Very similar to WCSimTruthSummary for WCSim, used for comparison with WCSimRecoSummary
 */
class TruthInfo: public TObject {
	public:
		TruthInfo();
		TruthInfo(int type, int beamPDG, float beamEnergy, int nPrimaries, std::vector<int> primaryPDGs,
				std::vector<double> primaryEnergies, std::vector<TVector3> primaryDirs);
		TruthInfo(const TruthInfo& other);
		TruthInfo& operator=(const TruthInfo& rhs);
		~TruthInfo();

		void Print();

		void SetVtxTime(float t);
		void SetVtx(float x, float y, float z);
		void SetBeamDir(float x, float y, float z);

		int GetType() {
			return fType;
		}
		float GetBeamE() {
			return fBeamEnergy;
		}
		int GetBeamPDG() {
			return fBeamPDG;
		}

		int GetPrimaryPDG(int p) {
			if (p < fNPrimaries) {
				return fPrimaryPDGs[p]; // Index starts at 0
			} else {
				std::cerr << "GetPrimaryPDG(index) out of range [0..." << fNPrimaries - 1 << "]" << std::endl;
				return -999;
			}
		}
		double GetPrimaryEnergy(int p) {
			if (p < fNPrimaries) {
				return fPrimaryEnergies[p]; // Index starts at 0
			} else {
				std::cerr << "GetPrimaryEnergies(index) out of range [0..." << fNPrimaries - 1 << "]" << std::endl;
				return -999;
			}
		} // Index starts at 0
		TVector3 GetPrimaryDir(int p) {
			if (p < fNPrimaries) {
				return fPrimaryDirs[p]; // Index starts at 0
			} else {
				std::cerr << "GetPrimaryDir(index) out of range [0..." << fNPrimaries - 1 << "]" << std::endl;
				return TVector3(-999, -999, -999);
			}
		} // Index starts at 0

		int GetNPrimaries() {
			return fNPrimaries;
		}
		std::vector<int> GetPrimaryPDGs() {
			return fPrimaryPDGs;
		}
		std::vector<double> GetPrimaryEnergies() {
			return fPrimaryEnergies;
		}
		std::vector<TVector3> GetPrimaryDirs() {
			return fPrimaryDirs;
		}

		bool IsCC() {
			return fIsCC;
		}
		bool IsNC() {
			return fIsNC;
		}
		bool IsQE() {
			return fIsQE;
		}
		bool IsRes() {
			return fIsRes;
		}
		bool IsDIS() {
			return fIsDIS;
		}
		bool IsCoherent() {
			return fIsCoherent;
		}
		bool IsNueElectronElastic() {
			return fIsNueElectronElastic;
		}
		bool IsInverseMuonDecay() {
			return fIsInverseMuonDecay;
		}
		bool IsOther() {
			return fIsOther;
		}

		float GetVtxTime() const {
			return fVtxTime;
		}
		float GetVtxX() const {
			return fVtxX;
		}
		float GetVtxY() const {
			return fVtxY;
		}
		float GetVtxZ() const {
			return fVtxZ;
		}

		float GetBeamDirX() const {
			return fBeamDirX;
		}
		float GetBeamDirY() const {
			return fBeamDirY;
		}
		float GetBeamDirZ() const {
			return fBeamDirZ;
		}

	private:
		int fType;  							///< Interaction type code from WCSimTruthSummary
		int fBeamPDG;  							///< PDG code of the beam particle (usually a neutrino)
		float fBeamEnergy; 						///< Energy of the incident beam particle

		int fNPrimaries;
		std::vector<int> fPrimaryPDGs; 			///< Vector of primary PDG's
		std::vector<double> fPrimaryEnergies; 	///< Vector of primary energies
		std::vector<TVector3> fPrimaryDirs; 	///< Vector of primary directions <TVector3>

		bool fIsCC;								///< True if event is Charged-Current
		bool fIsNC;								///< True if event is Neutral-Current
		bool fIsQE;								///< True if event is Quasielastic
		bool fIsRes;							///< True if event is Resonant
		bool fIsDIS;							///< True if event is Deep Inelastic Scattering
		bool fIsCoherent;						///< True if event is Coherent
		bool fIsNueElectronElastic;				///< True if event is electron elastic
		bool fIsInverseMuonDecay;				///< True is event is inverse muon decay
		bool fIsOther; 							///< A CC event that doesn't fall into any of the above categories
												// sometimes there isn't a nuance code for the type of event Genie has made

		float fVtxTime;							///< Neutrino interaction vertex time
		float fVtxX;							///< Neutrino interaction vertex x-position
		float fVtxY;							///< Neutrino interaction vertex y-position
		float fVtxZ;							///< Neutrino interaction vertex z-position
		float fBeamDirX;						///< Neutrino interaction vertex x-direction
		float fBeamDirY;						///< Neutrino interaction vertex y-direction
		float fBeamDirZ;						///< Neutrino interaction vertex z-direction

		ClassDef(TruthInfo,1)
};

/*
 * \class PidInfo
 * This class holds variables destined to be used in various PIDs
 * Combines the old HitInfo and RecoInfo classes
 */
class PidInfo: public TObject {
	public:
		PidInfo();
		PidInfo(const PidInfo& other);
		PidInfo& operator=(const PidInfo& rhs);
		~PidInfo();

		void Print();

		// Old HitInfo Functions
		void SetHitInfo(bool veto, int NVetoHits, int NHits, int NHitsUpstream, int NHitsInBottom, int NHitsInTop,
				int NHitsAboveMid, float totalQ, float totalQUpstream, float totalQInBottom, float totalQInTop,
				float totalQAboveMid);

		bool GetVeto() const {
			return fVeto;
		}

		int GetNVetoHits() const {
			return fNVetoHits;
		}

		int GetNHits() const {
			return fNHits;
		}
		int GetNHitsUpstream() const {
			return fNHitsUpstream;
		}
		int GetNHitsDownstream() const {
			return fNHitsDownstream;
		}
		int GetNHitsInBottom() const {
			return fNHitsInBottom;
		}
		int GetNHitsInTop() const {
			return fNHitsInTop;
		}
		int GetNHitsAboveMid() const {
			return fNHitsAboveMid;
		}
		int GetNHitsBelowMid() const {
			return fNHitsBelowMid;
		}

		float GetFracHitsUpstream() const {
			return fFracHitsUpstream;
		}
		float GetFracHitsDownstream() const {
			return fFracHitsDownstream;
		}
		float GetFracHitsInBottom() const {
			return fFracHitsInBottom;
		}
		float GetFracHitsInTop() const {
			return fFracHitsInTop;
		}
		float GetFracHitsAboveMid() const {
			return fFracHitsAboveMid;
		}
		float GetFracHitsBelowMid() const {
			return fFracHitsBelowMid;
		}

		float GetTotalQ() const {
			return fTotalQ;
		}
		float GetTotalQUpstream() const {
			return fTotalQUpstream;
		}
		float GetTotalQDownstream() const {
			return fTotalQDownstream;
		}
		float GetTotalQInBottom() const {
			return fTotalQInBottom;
		}
		float GetTotalQInTop() const {
			return fTotalQInTop;
		}
		float GetTotalQAboveMid() const {
			return fTotalQAboveMid;
		}
		float GetTotalQBelowMid() const {
			return fTotalQBelowMid;
		}

		float GetFracQUpstream() const {
			return fFracQUpstream;
		}
		float GetFracQDownstream() const {
			return fFracQDownstream;
		}
		float GetFracQInBottom() const {
			return fFracQInBottom;
		}
		float GetFracQInTop() const {
			return fFracQDownstream;
		}
		float GetFracQAboveMid() const {
			return fFracQAboveMid;
		}
		float GetFracQBelowMid() const {
			return fFracQBelowMid;
		}

		// Old RecoInfo Functions
		void SetLikelihoods(float total, float hit, float charge, float time);

		void SetTotalQ(float in, float out, float hole);
		void SetNHits(int in, int out, int hole);
		void SetPredQ(float in, float out, float hole);

		void SetEnergy(float E);
		void SetVtxTime(float t);
		void SetVtx(float x, float y, float z);
		void SetEnd(float x, float y, float z);
		void SetDir(float x, float y, float z);

		void SetEscapes(bool escapes);

		float GetTotal2LnL() const {
			return fTotal2LnL;
		}
		float GetHit2LnL() const {
			return fHit2LnL;
		}
		float GetCharge2LnL() const {
			return fCharge2LnL;
		}
		float GetTime2LnL() const {
			return fTime2LnL;
		}
		float GetCutoff2LnL() const {
			return fCutoff2LnL;
		}

		float GetTotalQInRing() const {
			return fTotalQInRing;
		}
		float GetTotalQOutsideRing() const {
			return fTotalQOutsideRing;
		}
		float GetTotalQInRingHole() const {
			return fTotalQInRingHole;
		}

		int GetNHitsInRing() {
			return fNHitsInRing;
		}
		int GetNHitsOutsideRing() {
			return fNHitsOutsideRing;
		}
		int GetNHitsInRingHole() {
			return fNHitsInRingHole;
		}

		float GetPredQInRing() const {
			return fPredQInRing;
		}
		float GetPredQOutsideRing() const {
			return fPredQOutsideRing;
		}
		float GetPredQInRingHole() const {
			return fPredQInRingHole;
		}

		float GetFracTotalQInRing() const {
			return fFracTotalQInRing;
		}
		float GetFracTotalQOutsideRing() const {
			return fFracTotalQOutsideRing;
		}
		float GetFracTotalQInRingHole() const {
			return fFracTotalQInRingHole;
		}

		float GetFracNHitsInRing() {
			return fFracNHitsInRing;
		}
		float GetFracNHitsOutsideRing() {
			return fFracNHitsOutsideRing;
		}
		float GetFracNHitsInRingHole() {
			return fFracNHitsInRingHole;
		}

		float GetFracPredQInRing() const {
			return fFracPredQInRing;
		}
		float GetFracPredQOutsideRing() const {
			return fFracPredQOutsideRing;
		}
		float GetFracPredQInRingHole() const {
			return fFracPredQInRingHole;
		}

		float GetPredictedOverTotalCharge() const {
			return fPredictedOverTotalCharge;
		}

		float GetEnergy() const {
			return fEnergy;
		}
		float GetVtxTime() const {
			return fVtxTime;
		}
		float GetVtxX() const {
			return fVtxX;
		}
		float GetVtxY() const {
			return fVtxY;
		}
		float GetVtxZ() const {
			return fVtxZ;
		}
		float GetVtxRho() const {
			return fVtxRho;
		}

		float GetEndX() const {
			return fEndX;
		}
		float GetEndY() const {
			return fEndY;
		}
		float GetEndZ() const {
			return fEndZ;
		}
		float GetEndRho() const {
			return fEndRho;
		}

		float GetDirX() const {
			return fDirX;
		}
		float GetDirY() const {
			return fDirY;
		}
		float GetDirZ() const {
			return fDirZ;
		}

		bool Escapes() const {
			return fEscapes;
		}

	private:
		// I think this class should not hold variables that can also be found in RecoSummary
		// That class can deal with multiple tracks while this one cannot

		// Old HitInfo Variables
		bool fVeto;					///< True if the reconstructed track exits the detector

		int fNVetoHits;             ///< Number of PMTs hit in the veto

		int fNHits;					///< Number of PMTs hit in the event
		int fNHitsUpstream;			///< Number of PMTs hit in the upstream(x<0) region of the detector
		int fNHitsDownstream;		///< Number of PMTs hit in the downstream(x>=0) region of the detector
		int fNHitsInBottom;			///< Number of PMTs hit in the bottom endcap
		int fNHitsInTop;			///< Number of PMTs hit in the top endcap
		int fNHitsAboveMid;         ///< Number of PMTs hit with Z>0
		int fNHitsBelowMid;         ///< Number of PMTs hit with Z<0

		float fFracHitsUpstream;	///< fFracHitsUpstream = fNHitsUpstream / fNHits
		float fFracHitsDownstream;	///< fFracHitsDownstream = fNHitsDownstream / fNHits
		float fFracHitsInBottom;	///< fFracHitsInBottom = fNHitsInBottom / fNHits
		float fFracHitsInTop;		///< fFracHitsInTop = fNHitsInTop / fNHits
		float fFracHitsAboveMid;    ///< fFracHitsAboveMid = fNHitsAboveMid / fNHits
		float fFracHitsBelowMid;    ///< fFracHitsBelowMid = fNHitfNHitsBelowMidsAboveMid / fNHits

		float fTotalQ;				///< Total collected charge on all PMTs in the event
		float fTotalQUpstream;		///< Total collected charge on PMTs in the upstream(x<0) region of the detector
		float fTotalQDownstream;	///< Total collected charge on PMTs in the downstream(x>=0) region of the detector
		float fTotalQInBottom;		///< Total collected charge on PMTs in the bottom endcap
		float fTotalQInTop;			///< Total collected charge on PMTs in the top endcap
		float fTotalQAboveMid;      ///< Total collected charge on PMTs with Z>0
		float fTotalQBelowMid;      ///< Total collected charge on PMTs with Z<0

		float fFracQUpstream;		///< fFracQUpstream = fTotalQUpstream / fTotalQ
		float fFracQDownstream;		///< fFracQDownstream = fTotalQDownstream / fTotalQ
		float fFracQInBottom;		///< fFracQInBottom = fTotalQInBottom / fTotalQ
		float fFracQInTop;			///< fFracQInTop = fTotalQInTop / fTotalQ
		float fFracQAboveMid;       ///< fFracQAboveMid = fTotalQAboveMid / fTotalQ
		float fFracQBelowMid;       ///< fFracQBelowMid = fTotalQBelowMid / fTotalQ

		// Old RecoInfo variables
		float fTotal2LnL;				///< Final total -2LnL
		float fHit2LnL;					///< Final hit component of total -2LnL
		float fCharge2LnL;				///< Final charge component of total -2LnL
		float fTime2LnL;				///< Final time component of total -2LnL
		float fCutoff2LnL;				///< Final cutoff componenet of total -2LnL

		float fTotalQInRing;			///< Total charge within the reconstructed ring
		float fTotalQOutsideRing;		///< Total charge outside the reconstructed ring
		float fTotalQInRingHole;		///< Total charge within the hole of the reconstructed ring

		int fNHitsInRing;				///< Number of hits within the reconstructed ring
		int fNHitsOutsideRing;			///< Number of hits outside the reconstructed ring
		int fNHitsInRingHole;			///< Number of hits within the hole of the reconstructed ring

		float fPredQInRing;				///< Predicted charge inside the reconstructed ring
		float fPredQOutsideRing;		///< Predicted charge outside the reconstructed ring
		float fPredQInRingHole;			///< Predicted charge within the hole of the reconstructed ring

		float fFracTotalQInRing;		///< Fraction of actual total charge inside the reconstructed ring
		float fFracTotalQOutsideRing;	///< Fraction of actual total charge outside the reconstructed ring
		float fFracTotalQInRingHole;	///< Fraction of actual total charge within the hole of the reconstructed ring

		float fFracNHitsInRing;			///< Fraction of total hits inside the reconstructed ring
		float fFracNHitsOutsideRing;	///< Fraction of total hits outside the reconstructed ring
		float fFracNHitsInRingHole;		///< Fraction of total hits within the hole of the reconstructed ring

		float fFracPredQInRing;			///< Fraction of predicted charge inside the reconstructed ring
		float fFracPredQOutsideRing;	///< Fraction of predicted charge outside the reconstructed ring
		float fFracPredQInRingHole;		///< Fraction of predicted charge within the hole of the reconstructed ring

		float fPredictedOverTotalCharge;		///< Total predicted charge over total actual charge

		float fVtxTime;					///< Reconstructed neutrino interaction vertex time
		float fEnergy;					///< Reconstructed primary particle energy
		float fVtxX;					///< Reconstructed neutrino interaction x-position
		float fVtxY;					///< Reconstructed neutrino interaction y-position
		float fVtxZ;					///< Reconstructed neutrino interaction z-position
		float fVtxRho;					///< Reconstructed neutrino interaction radial position

		float fEndX;					///< Reconstructed track end x-position
		float fEndY;					///< Reconstructed track end x-position
		float fEndZ;					///< Reconstructed track end x-position
		float fEndRho;					///< Reconstructed track end radial position

		float fDirX;					///< Reconstructed track x-direction
		float fDirY;					///< Reconstructed track y-direction
		float fDirZ;					///< Reconstructed track z-direction

		bool fEscapes;					///< True if the track escapes the detector

		ClassDef(PidInfo,1)
};

/*
 * \class SeedInfo
 * This class holds info relating to the Hough transform ring finder and resulting seed given to the main fit
 */
class SeedInfo: public TObject {
	public:
		SeedInfo();
		SeedInfo(std::vector<WCSimLikelihoodTrackBase*> tracks, int slices, std::vector<WCSimRecoRing*> ringVec,
				std::vector<double> ringTime);
		SeedInfo(const SeedInfo& other);
		SeedInfo& operator=(const SeedInfo& rhs);
		~SeedInfo();

		void Print();

		int GetNumTracks() {
			return int(fTracks.size());
		}

		std::vector<WCSimLikelihoodTrackBase*> GetTracks() {
			return fTracks;
		}
		WCSimLikelihoodTrackBase* GetTrack(int p) {
			if (p >= 0 && (size_t) p < fTracks.size()) {
				return fTracks[p]; // Index starts at 0
			} else {
				std::cerr << "GetTrack(index) out of range [0..." << fTracks.size() - 1 << "]" << std::endl;
				return NULL;
				//return WCSimLikelihoodTrackFactory::MakeTrack(TrackType::Unknown, -999, -999, -999, -999, -999, -999, -999, -999);
			}
		}

		int GetNumSlices() {
			return fNSlices;
		}
		int GetNumRings() {
			if (fRingHeight.size() == fRingTime.size()) {
				return fRingHeight.size();
			} else {
				std::cout << "fRingVec and fRingTime are not the same size!!!" << std::endl;
				return -999;
			}
		}

		std::vector<double> GetRingHeights() {
			return fRingHeight;
		}
		double GetRingHeight(int p) {
			if (p >= 0 && (size_t) p < fRingHeight.size()) {
				return fRingHeight[p];
			} else {
				std::cerr << "GetRing(index) out of range [0..." << fRingHeight.size() - 1 << "]" << std::endl;
				return -999;
			}
		}

		std::vector<double> GetRingTimes() {
			return fRingTime;
		}
		double GetRingTime(int p) {
			if (p >= 0 && (size_t) p < fRingTime.size()) {
				return fRingTime[p];
			} else {
				std::cerr << "GetRingTime(index) out of range [0..." << fRingTime.size() - 1 << "]" << std::endl;
				return -999;
			}
		}

		std::vector<double> GetRingAngles() {
			return fRingAngle;
		}
		double GetRingAngle(int p) {
			if (p >= 0 && (size_t) p < fRingAngle.size()) {
				return fRingAngle[p];
			} else {
				std::cerr << "GetRingAngle(index) out of range [0..." << fRingAngle.size() - 1 << "]" << std::endl;
				return -999;
			}
		}

		std::vector<TVector3> GetRingVtxs() {
			return fRingVtx;
		}
		TVector3 GetRingVtx(int p) {
			if (p >= 0 && (size_t) p < fRingVtx.size()) {
				return fRingVtx[p];
			} else {
				std::cerr << "GetRingVtx(index) out of range [0..." << fRingVtx.size() - 1 << "]" << std::endl;
				return TVector3(-999, -999, -999);
			}
		}

		std::vector<TVector3> GetRingDirs() {
			return fRingDir;
		}
		TVector3 GetRingDir(int p) {
			if (p >= 0 && (size_t) p < fRingDir.size()) {
				return fRingDir[p];
			} else {
				std::cerr << "GetRingDir(index) out of range [0..." << fRingDir.size() - 1 << "]" << std::endl;
				return TVector3(-999, -999, -999);
			}
		}

	private:
		std::vector<WCSimLikelihoodTrackBase*> fTracks; ///< Stores a vector of the tracks actually passed to the fitter.

		int fNSlices; 									///< Number of slices produced by the seeder.
		std::vector<double> fRingHeight; 				///< Stores the ring heights produces by the seeder.
		std::vector<double> fRingTime; ///< Stores the vertex times of the rings. These are the same for all rings belonging to the same slice.
		std::vector<double> fRingAngle; 				///< Stores the ring angles
		std::vector<TVector3> fRingVtx; 				///< Stores TVector3 objects containing the ring vtx positions
		std::vector<TVector3> fRingDir; 				///< Stores Tvector3 objects containing the ring dirs.

		ClassDef(SeedInfo,1)
};

/*
 * \class StageInfo
 * This class holds info relating to all the PMTs and tracks at all stages of the fit
 * WARNING: Causes vastly increased file size, only use if required
 */
class StageInfo: public TObject {
	public:
		StageInfo();
		StageInfo(const StageInfo& other);
		StageInfo& operator=(const StageInfo& rhs);
		~StageInfo();

		void SetStageInfo(int stageNCalls, std::vector<WCSimHitPrediction> stagePreds,
				std::vector<WCSimLikelihoodTrackBase*> stageTracks);

		void Print();

		int GetNumStages() {
			return int(fStageNcalls.size());
		}
		int GetNumTracks() {
			return int(fStageTracks[0].size());
		}

		std::vector<int> GetNCallsVec() {
			return fStageNcalls;
		}
		int GetNCalls(int stage) {
			if (stage >= 0 && (size_t) stage < fStageNcalls.size()) {
				return fStageNcalls[stage]; // Index starts at 0
			} else {
				std::cerr << "fStageNcalls(index) out of range [0..." << fStageNcalls.size() - 1 << "]" << std::endl;
				return -999;
			}
		}

		std::vector<WCSimHitPrediction> GetStagePreds(int stage) {
			if (stage >= 0 && (size_t) stage < fStagePreds.size()) {
				return fStagePreds[stage]; // Index starts at 0
			} else {
				std::cerr << "GetStagePreds(index) out of range [0..." << fStagePreds.size() - 1 << "]" << std::endl;
				std::vector<WCSimHitPrediction> error;
				WCSimHitPrediction bestFit(-999, -999, -999, -999, -999, -999);
				error.push_back(bestFit);
				return error;
			}
		}

		std::vector<WCSimLikelihoodTrackBase*> GetStageTracks(int stage) {
			if (stage >= 0 && (size_t) stage < fStageTracks.size()) {
				return fStageTracks[stage]; // Index starts at 0
			} else {
				std::cerr << "GetStageTracks(index) out of range [0..." << fStagePreds.size() - 1 << "]" << std::endl;
				std::vector<WCSimLikelihoodTrackBase*> error;
				error.push_back(
						WCSimLikelihoodTrackFactory::MakeTrack(TrackType::Unknown, -999, -999, -999, -999, -999, -999,
								-999, -999));
				return error;
			}
		}

		// TODO: Add functions to get total predicted charge etc... Add up the likelihoods in fStagePreds to give total across the stage!

	private:
		std::vector<int> fStageNcalls; 						///< Number of total function calls at the end of the stage
		std::vector<std::vector<WCSimHitPrediction> > fStagePreds; ///< A vector of WCSimHitPrediction's for every PMT at the end of every stage
		std::vector<std::vector<WCSimLikelihoodTrackBase*> > fStageTracks; ///< A vector of the fitted track at the end of every stage

		ClassDef(StageInfo,1)
};

/*
 * \class WCSimOutputTree
 * Deals with the filling and saving of the output tree from the reconstruction
 */
class WCSimOutputTree: public TObject {
	public:
		WCSimOutputTree();
		virtual ~WCSimOutputTree();

		TString FormTime();
		std::string FitType(WCSimFitterConfig * config);
		TString InitOutputTree(TString fileName, bool modifyInput, WCSimFitterConfig * fitterConfig);

		void GetExistingTree();
		void MakeNewTree();
		void MakeNewFit(std::string type);

		void SaveTree();
		void SaveFile();

		void FillSeedInfo(SeedInfo& seedInfo);

		void FillStageInfo(int stageNCalls, std::vector<WCSimHitPrediction> stagePreds,
				std::vector<WCSimLikelihoodTrackBase*> stageTracks);

		void FillHitComparison(WCSimHitComparison& hitComp);

		void Fill(EventHeader& eventHead, TruthInfo& truthInfo, WCSimRecoSummary& summ, PidInfo& pidInfo);

		void SetSaveFileName(TString saveName);
		TString GetSaveFileName() const;

		void SetModifyFile(bool modify);

		void SetInputFile(TString inputFile);
		TString GetInputFileName() const;

	private:
		void DeletePointersIfExist();

		TString fInputFile;	///< The input file name, included in eventHeader, but for ease of use store seperately

		TFile * fSaveFile;						///<
		TString fSaveFileName;					///<
		TTree * fResultsTree;					///< The main results tree

		TTree * fGeoTree;						///< The geometry tree, holds geometry and PMT info
		WCSimRootGeom * fGeometry;				///< The WCSim geometry object that is then saved in the geometry tree

		TTree * fParameterTree;					///< The parameter tree, holds ParameterInfo
		ParameterInfo * fParameterInfo;			///< Holds the WCSimParameters for the reconstruction used

		bool fModifyFile;
		EventHeader * fEventHeader;				///< Holds info about the input event
		TruthInfo * fTruthInfo;					///< Akin to WCSimTruthSummary, holds the truth information
		WCSimRecoSummary * fRecoSummary;		///< Result of the fit, for comparison with the TruthInfo
		PidInfo * fPidInfo;			///< Holds variables to be used in the PID, combines the old HitInfo and RecoInfo
		SeedInfo * fSeedInfo;					///< Holds info on the output from the Hough transform ring finder
		StageInfo * fStageInfo;	///< Holds info on how the fit parameters and minus2LnL evolve over the stages of the fit
		WCSimHitComparison * fHitComparison;	///< Holds the final HitComparison for all the PMTs

		TBranch * fEventHeader_branch;
		TBranch * fTruthInfo_branch;
		TBranch * fRecoSummary_branch;
		TBranch * fPidInfo_branch;
		TBranch * fSeedInfo_branch;
		TBranch * fStageInfo_branch;
		TBranch * fHitComparison_branch;

		ClassDef(WCSimOutputTree, 3)
};

#endif /* WCSIMOUTPUTTREE_HH */
