#ifndef WCSIMPARAMETERS_HH
#define WCSIMPARAMETERS_HH

#include "TObject.h"

class WCSimParameters: public TObject {
	public:
		static WCSimParameters* Instance();

		static Double_t SpeedOfLight();
		static Double_t CherenkovAngle();
		static Double_t TimeResolution(Double_t Q);
		static Double_t WCSimTimeResolution(Double_t Q, Double_t timeConst);
		static Double_t TimeSlew(Double_t Q);
		static Double_t RefractiveIndex(Double_t r);

		static Double_t ThetaC();     // Cherenkov Angle
		static Double_t CosThetaC();  // Cosine of Cherenkov Angle

		static void PrintParameters();
		void RunPrintParameters();

		void SetUseTimeOnly(Bool_t doIt = true);
		void SetUseChargeOnly(Bool_t doIt = true);
		void SetUseChargeAndTime(Bool_t doIt = true);

		Double_t GetTimeResolution(Double_t Q);
		Double_t GetTimeSlew(Double_t Q);
		Double_t GetRefractiveIndex(Double_t r);

		Double_t GetSimpleTimeResolution(Double_t Q);
		Double_t GetSimpleTimeSlew() {
			return 0.0;
		}

		Double_t GetSimpleRefractiveIndex() {
			return 1.33;
		}

		// Get and Set methods

		// Slicer parameters
		void SetSlicerClusterDistance(Double_t val) {
			fSlicerClusterDistance = val;
		}
		Double_t GetSlicerClusterDistance() {
			return fSlicerClusterDistance;
		}

		void SetSlicerMinSize(UInt_t val) {
			fSlicerMinSize = val;
		}
		UInt_t GetSlicerMinSize() {
			return fSlicerMinSize;
		}

		void SetSlicerChargeCut(Double_t val) {
			fSlicerChargeCut = val;
		}
		Double_t GetSlicerChargeCut() {
			return fSlicerChargeCut;
		}

		void SetSlicerTimeCut(Double_t val) {
			fSlicerTimeCut = val;
		}
		Double_t GetSlicerTimeCut() {
			return fSlicerTimeCut;
		}

		void SetIterateSlicing(Double_t val) {
			fIterateSlicing = val;
		}
		Double_t GetIterateSlicing() {
			return fIterateSlicing;
		}

		// Veto slicing parameters
		void SetVetoClusterDistance(Double_t val) {
			fVetoClusterDistance = val;
		}
		Double_t GetVetoClusterDistance() {
			return fVetoClusterDistance;
		}

		void SetVetoMinSize(UInt_t val) {
			fVetoMinSize = val;
		}
		UInt_t GetVetoMinSize() {
			return fVetoMinSize;
		}

		void SetVetoMinChargeCut(Double_t val) {
			fVetoMinChargeCut = val;
		}
		Double_t GetVetoMinChargeCut() {
			return fVetoMinChargeCut;
		}

		void SetVetoMaxChargeCut(Double_t val) {
			fVetoMaxChargeCut = val;
		}
		Double_t GetVetoMaxChargeCut() {
			return fVetoMaxChargeCut;
		}

		void SetVetoTimeCut(Double_t val) {
			fVetoTimeCut = val;
		}
		Double_t GetVetoTimeCut() {
			return fVetoTimeCut;
		}

		// Integral parameters
		void SetCalculateIntegrals(Bool_t val) {
			fCalculateIntegrals = val;
		}
		Bool_t CalculateIntegrals() {
			return fCalculateIntegrals;
		}

		void SetTruncateIntegrals(Bool_t val) {
			fTruncateIntegrals = val;
		}
		Bool_t TruncateIntegrals() {
			return fTruncateIntegrals;
		}

		void SetConstrainExtent(Bool_t val) {
			fConstrainExtent = val;
		}
		Bool_t ConstrainExtent() {
			return fConstrainExtent;
		}

		// Likelihood tuner parameters
		void SetUseTransmission(Bool_t val) {
			fUseTransmission = val;
		}
		Bool_t UseTransmission() {
			return fUseTransmission;
		}

		void SetUseAngularEfficiency(Bool_t val) {
			fUseAngularEfficiency = val;
		}
		Bool_t UseAngularEfficiency() {
			return fUseAngularEfficiency;
		}

		void SetUseGlassCathodeReflection(Bool_t val) {
			fUseGlassCathodeReflection = val;
		}
		Bool_t UseGlassCathodeReflection() {
			return fUseGlassCathodeReflection;
		}

		void SetUseScatteringTable(Bool_t val) {
			fUseScatteringTable = val;
		}
		Bool_t UseScatteringTable() {
			return fUseScatteringTable;
		}

		// Fitter parameters
		void SetUseTime(Bool_t val) {
			fUseTime = val;
		}
		Bool_t UseTime() {
			return fUseTime;
		}

		void SetUseCharge(Bool_t val) {
			fUseCharge = val;
		}
		Bool_t UseCharge() {
			return fUseCharge;
		}

		void SetEqualiseChargeAndTime(Bool_t val) {
			fEqualiseChargeAndTime = val;
		}
		Bool_t EqualiseChargeAndTime() {
			return fEqualiseChargeAndTime;
		}

		void SetSaveWCSimRootEvent(Bool_t val) {
			fSaveWCSimRootEvent = val;
		}
		Bool_t SaveWCSimRootEvent() {
			return fSaveWCSimRootEvent;
		}

		void SetDigiType(std::string val) {
			fDigiType = val;
		}
		std::string GetDigiType() {
			return fDigiType;
		}

		void SetSaveSeedInfo(Bool_t val) {
			fSaveSeedInfo = val;
		}
		Bool_t SaveSeedInfo() {
			return fSaveSeedInfo;
		}

		void SetSaveStageInfo(Bool_t val) {
			fSaveStageInfo = val;
		}
		Bool_t SaveStageInfo() {
			return fSaveStageInfo;
		}

		void SetSaveHitComparison(Bool_t val) {
			fSaveHitComparison = val;
		}
		Bool_t SaveHitComparison() {
			return fSaveHitComparison;
		}

		// Speed of light parameters
		void SetUseCustomParticleSpeed(Bool_t val) {
			fUseCustomParticleSpeed = val;
		}
		Bool_t UseCustomParticleSpeed() {
			return fUseCustomParticleSpeed;
		}

		void SetUseCustomSpeedOfLight(Bool_t val) {
			fUseCustomSpeedOfLight = val;
		}
		Bool_t UseCustomSpeedOfLight() {
			return fUseCustomSpeedOfLight;
		}

		void SetUseFittedSpeedOfLight(Bool_t val) {
			fUseFittedSpeedOfLight = val;
		}
		Bool_t UseFittedSpeedOfLight() {
			return fUseFittedSpeedOfLight;
		}

		void SetCustomParticleSpeed(Double_t val) {
			fCustomParticleSpeed = val;
		}
		Double_t GetCustomParticleSpeed() {
			return fCustomParticleSpeed;
		}

		void SetCustomSpeedOfLight(Double_t val) {
			fCustomSpeedOfLight = val;
		}
		Double_t GetCustomSpeedOfLight() {
			return fCustomSpeedOfLight;
		}

		void SetFittedSpeedOfLight(Double_t val) {
			fFittedSpeedOfLight = val;
		}
		Double_t GetFittedSpeedOfLight() {
			return fFittedSpeedOfLight;
		}

		// Not currently used parameters
		void SetUseSimpleTimeResolution(Bool_t val) {
			fUseSimpleTimeResolution = val;
		}
		Bool_t GetUseSimpleTimeResolution() {
			return fUseSimpleTimeResolution;
		}

		void SetUseSimpleTimeSlew(Bool_t val) {
			fUseSimpleTimeSlew = val;
		}
		Bool_t GetUseSimpleTimeSlew() {
			return fUseSimpleTimeSlew;
		}

		void SetUseSimpleRefractiveIndex(Bool_t val) {
			fUseSimpleRefractiveIndex = val;
		}
		Bool_t GetUseSimpleRefractiveIndex() {
			return fUseSimpleRefractiveIndex;
		}

	private:
		WCSimParameters();
		~WCSimParameters();

		// Slicer parameters
		Double_t fSlicerClusterDistance; 	///< Max distance in cm between hits in the slices
		UInt_t fSlicerMinSize; 				///< Minimum number of hits for a slice to be formed
		Double_t fSlicerChargeCut; 			///< Only consider digits above the charge threshold
		Double_t fSlicerTimeCut; 			///< Maximum gap allowed between hits when ordered by time in order to be associated with the previous hit.
		Double_t fIterateSlicing; 			///< Iterate the charge cut to make as many slices as we have tracks

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

		// Fitter parameters
		Bool_t fUseTime;                    ///< True if we should include timing information in the likelihood
		Bool_t fUseCharge;                  ///< True if we should include charge information in the likelihood
		Bool_t fEqualiseChargeAndTime; 		///< After we've seeded the vertex, energy and time, we can scale the time likelihood to weight it the same as the charge
		Bool_t fSaveWCSimRootEvent; 		///< Whether to save a full copy of the WCSimRootEvent fitted, or just a link to the original file
		std::string fDigiType;            	///< Name of digitizer type to use in WCSimDigitizerLikelihood::DigiType_t
		Bool_t fSaveSeedInfo;			   	///< True if we should save the SeedInfo in the output file
		Bool_t fSaveStageInfo;			   	///< True if we should save the StageInfo in the output file
		Bool_t fSaveHitComparison;		   	///< True if we should save the HitComparison in the output file

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

		ClassDef(WCSimParameters,0)
};

#endif

