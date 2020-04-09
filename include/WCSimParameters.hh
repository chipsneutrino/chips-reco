#pragma once

#include "TObject.h"
#include <string>
#include <map>

class WCSimParameters : public TObject
{
public:
	static WCSimParameters *Instance();

	/**
		 * Set the path of the configuration file to be loaded
		 * @param config Path of a text file containing all the configuration parameters
		 */
	void SetConfigFile(const char *config);

	static Double_t SpeedOfLight();
	static Double_t CherenkovAngle();
	static Double_t TimeResolution(Double_t Q);
	static Double_t WCSimTimeResolution(Double_t Q, Double_t timeConst);
	static Double_t TimeSlew(Double_t Q);
	static Double_t RefractiveIndex(Double_t r);

	static Double_t ThetaC();	 // Cherenkov Angle
	static Double_t CosThetaC(); // Cosine of Cherenkov Angle

	static void PrintParameters();
	void RunPrintParameters();

	void SetUseTimeOnly(Bool_t doIt = true);
	void SetUseChargeOnly(Bool_t doIt = true);
	void SetUseChargeAndTime(Bool_t doIt = true);

	Double_t GetTimeResolution(Double_t Q);
	Double_t GetTimeSlew(Double_t Q);
	Double_t GetRefractiveIndex(Double_t r);

	Double_t GetSimpleTimeResolution(Double_t Q);
	Double_t GetSimpleTimeSlew()
	{
		return 0.0;
	}

	Double_t GetSimpleRefractiveIndex()
	{
		return 1.33;
	}

	// Get and Set methods

	// Slicer parameters
	void SetSlicerClusterDistance(Double_t val)
	{
		fSlicerClusterDistance = val;
	}
	Double_t GetSlicerClusterDistance()
	{
		return fSlicerClusterDistance;
	}

	void SetSlicerMinSize(UInt_t val)
	{
		fSlicerMinSize = val;
	}
	UInt_t GetSlicerMinSize()
	{
		return fSlicerMinSize;
	}

	void SetSlicerChargeCut(Double_t val)
	{
		fSlicerChargeCut = val;
	}
	Double_t GetSlicerChargeCut()
	{
		return fSlicerChargeCut;
	}

	void SetSlicerTimeCut(Double_t val)
	{
		fSlicerTimeCut = val;
	}
	Double_t GetSlicerTimeCut()
	{
		return fSlicerTimeCut;
	}

	void SetIterateSlicing(Bool_t val)
	{
		fIterateSlicing = val;
	}
	Bool_t IterateSlicing()
	{
		return fIterateSlicing;
	}

	// Veto slicing parameters
	void SetVetoClusterDistance(Double_t val)
	{
		fVetoClusterDistance = val;
	}
	Double_t GetVetoClusterDistance()
	{
		return fVetoClusterDistance;
	}

	void SetVetoMinSize(UInt_t val)
	{
		fVetoMinSize = val;
	}
	UInt_t GetVetoMinSize()
	{
		return fVetoMinSize;
	}

	void SetVetoMinChargeCut(Double_t val)
	{
		fVetoMinChargeCut = val;
	}
	Double_t GetVetoMinChargeCut()
	{
		return fVetoMinChargeCut;
	}

	void SetVetoMaxChargeCut(Double_t val)
	{
		fVetoMaxChargeCut = val;
	}
	Double_t GetVetoMaxChargeCut()
	{
		return fVetoMaxChargeCut;
	}

	void SetVetoTimeCut(Double_t val)
	{
		fVetoTimeCut = val;
	}
	Double_t GetVetoTimeCut()
	{
		return fVetoTimeCut;
	}

	// Integral parameters
	void SetCalculateIntegrals(Bool_t val)
	{
		fCalculateIntegrals = val;
	}
	Bool_t CalculateIntegrals()
	{
		return fCalculateIntegrals;
	}

	void SetTruncateIntegrals(Bool_t val)
	{
		fTruncateIntegrals = val;
	}
	Bool_t TruncateIntegrals()
	{
		return fTruncateIntegrals;
	}

	void SetConstrainExtent(Bool_t val)
	{
		fConstrainExtent = val;
	}
	Bool_t ConstrainExtent()
	{
		return fConstrainExtent;
	}

	// Likelihood tuner parameters
	void SetUseTransmission(Bool_t val)
	{
		fUseTransmission = val;
	}
	Bool_t UseTransmission()
	{
		return fUseTransmission;
	}

	void SetUseAngularEfficiency(Bool_t val)
	{
		fUseAngularEfficiency = val;
	}
	Bool_t UseAngularEfficiency()
	{
		return fUseAngularEfficiency;
	}

	void SetUseGlassCathodeReflection(Bool_t val)
	{
		fUseGlassCathodeReflection = val;
	}
	Bool_t UseGlassCathodeReflection()
	{
		return fUseGlassCathodeReflection;
	}

	void SetUseScatteringTable(Bool_t val)
	{
		fUseScatteringTable = val;
	}
	Bool_t UseScatteringTable()
	{
		return fUseScatteringTable;
	}

	void SetUseNewAngularEfficiency(Bool_t val)
	{
		fUseNewAngularEfficiency = val;
	}
	Bool_t UseNewAngularEfficiency()
	{
		return fUseNewAngularEfficiency;
	}

	void SetUseTrackFit(Bool_t val)
	{
		fUseTrackFit = val;
	}
	Bool_t UseTrackFit()
	{
		return fUseTrackFit;
	}

	void SetQScaling(Double_t val)
	{
		fQScaling = val;
	}
	Double_t GetQScaling()
	{
		return fQScaling;
	}

	// Fitter parameters
	void SetUseTime(Bool_t val)
	{
		fUseTime = val;
	}
	Bool_t UseTime()
	{
		return fUseTime;
	}

	void SetUseCharge(Bool_t val)
	{
		fUseCharge = val;
	}
	Bool_t UseCharge()
	{
		return fUseCharge;
	}

	void SetEqualiseChargeAndTime(Bool_t val)
	{
		fEqualiseChargeAndTime = val;
	}
	Bool_t EqualiseChargeAndTime()
	{
		return fEqualiseChargeAndTime;
	}

	void SetSaveWCSimRootEvent(Bool_t val)
	{
		fSaveWCSimRootEvent = val;
	}
	Bool_t SaveWCSimRootEvent()
	{
		return fSaveWCSimRootEvent;
	}

	void SetDigiType(std::string val)
	{
		fDigiType = val;
	}
	std::string GetDigiType()
	{
		return fDigiType;
	}

	void SetSaveSeedInfo(Bool_t val)
	{
		fSaveSeedInfo = val;
	}
	Bool_t SaveSeedInfo()
	{
		return fSaveSeedInfo;
	}

	void SetSaveStageInfo(Bool_t val)
	{
		fSaveStageInfo = val;
	}
	Bool_t SaveStageInfo()
	{
		return fSaveStageInfo;
	}

	void SetSaveHitComparison(Bool_t val)
	{
		fSaveHitComparison = val;
	}
	Bool_t SaveHitComparison()
	{
		return fSaveHitComparison;
	}

	void SetSaveParameterInfo(Bool_t val)
	{
		fSaveParameters = val;
	}
	Bool_t SaveParameterInfo()
	{
		return fSaveParameters;
	}

	// Speed of light parameters
	void SetUseCustomParticleSpeed(Bool_t val)
	{
		fUseCustomParticleSpeed = val;
	}
	Bool_t UseCustomParticleSpeed()
	{
		return fUseCustomParticleSpeed;
	}

	void SetUseCustomSpeedOfLight(Bool_t val)
	{
		fUseCustomSpeedOfLight = val;
	}
	Bool_t UseCustomSpeedOfLight()
	{
		return fUseCustomSpeedOfLight;
	}

	void SetUseFittedSpeedOfLight(Bool_t val)
	{
		fUseFittedSpeedOfLight = val;
	}
	Bool_t UseFittedSpeedOfLight()
	{
		return fUseFittedSpeedOfLight;
	}

	void SetCustomParticleSpeed(Double_t val)
	{
		fCustomParticleSpeed = val;
	}
	Double_t GetCustomParticleSpeed()
	{
		return fCustomParticleSpeed;
	}

	void SetCustomSpeedOfLight(Double_t val)
	{
		fCustomSpeedOfLight = val;
	}
	Double_t GetCustomSpeedOfLight()
	{
		return fCustomSpeedOfLight;
	}

	void SetFittedSpeedOfLight(Double_t val)
	{
		fFittedSpeedOfLight = val;
	}
	Double_t GetFittedSpeedOfLight()
	{
		return fFittedSpeedOfLight;
	}

	// Fitter plots parameters
	void SetPlotVtxXMax(Double_t val)
	{
		fPlotVtxXMax = val;
	}
	Double_t GetPlotVtxXMax()
	{
		return fPlotVtxXMax;
	}

	void SetPlotVtxYMax(Double_t val)
	{
		fPlotVtxYMax = val;
	}
	Double_t GetPlotVtxYMax()
	{
		return fPlotVtxYMax;
	}

	void SetPlotVtxZMax(Double_t val)
	{
		fPlotVtxZMax = val;
	}
	Double_t GetPlotVtxZMax()
	{
		return fPlotVtxZMax;
	}

	void SetPlotVtxTMax(Double_t val)
	{
		fPlotVtxTMax = val;
	}
	Double_t GetPlotVtxTMax()
	{
		return fPlotVtxTMax;
	}

	void SetPlotThetaMax(Double_t val)
	{
		fPlotThetaMax = val;
	}
	Double_t GetPlotThetaMax()
	{
		return fPlotThetaMax;
	}

	void SetPlotPhiMax(Double_t val)
	{
		fPlotPhiMax = val;
	}
	Double_t GetPlotPhiMax()
	{
		return fPlotPhiMax;
	}

	void SetPlotEnergyMax(Double_t val)
	{
		fPlotEnergyMax = val;
	}
	Double_t GetPlotEnergyMax()
	{
		return fPlotEnergyMax;
	}

	Double_t GetPlotMax(std::string param);

	// Ignore paramters

	void SetIgnoreLowCut(Double_t val)
	{
		fIgnoreLowCut = val;
	}
	Double_t GetIgnoreLowCut()
	{
		return fIgnoreLowCut;
	}

	void SetIgnoreHighCut(Double_t val)
	{
		fIgnoreHighCut = val;
	}
	Double_t GetIgnoreHighCut()
	{
		return fIgnoreHighCut;
	}

	// Not currently used parameters
	void SetUseSimpleTimeResolution(Bool_t val)
	{
		fUseSimpleTimeResolution = val;
	}
	Bool_t GetUseSimpleTimeResolution()
	{
		return fUseSimpleTimeResolution;
	}

	void SetUseSimpleTimeSlew(Bool_t val)
	{
		fUseSimpleTimeSlew = val;
	}
	Bool_t GetUseSimpleTimeSlew()
	{
		return fUseSimpleTimeSlew;
	}

	void SetUseSimpleRefractiveIndex(Bool_t val)
	{
		fUseSimpleRefractiveIndex = val;
	}
	Bool_t GetUseSimpleRefractiveIndex()
	{
		return fUseSimpleRefractiveIndex;
	}

private:
	WCSimParameters();
	~WCSimParameters();

	/// Read the configuration text file specified in WCSimAnalysisConfig::fConfName
	void LoadConfig();

	/**
		 * Strip out comments preceded by // or # from the text being parsed
		 * @param str String to remove comments from
		 */
	void IgnoreComments(std::string &str);

	/**
		 * Check whether a line could be trying to set a value
		 * @param str Line to check
		 * @return True if str if of the form "TEXT1 = TEXT2" (i.e. an = sign and text either side)
		 */
	const Bool_t IsGoodLine(std::string str);

	/**
		 * Check for strings containing only whitespace
		 * @param str String to check
		 * @return True if string is blank or just whitespace
		 */
	const Bool_t IsBlankLine(std::string str);

	/**
		 * Extract key and value from configuration file line
		 * @param lhs String to extract the key to
		 * @param rhs String to extract the value to
		 * @param str String in the format "key = value"
		 */
	void ExtractPair(std::string &lhs, std::string &rhs, std::string str);

	/**
		 * Check if a line is valid, then extract the key and value from it and add
		 * them to WCSimAnalysisConfig::fMap
		 * @param str String to parse
		 * @param lineNum Line number in the file
		 */
	void ParseLine(std::string str, Int_t lineNum);

	/**
		 * Add a key and value to WCSimAnalysisConfig::fMap
		 * @param lhs Key string
		 * @param rhs Value string
		 */
	void AddToMap(std::string lhs, std::string rhs);

	/// Iterate through the map and use it to set all the variables it specifies
	void SetFromMap();

	std::string fConfName;					 ///< Path of a text file containing all the configuration parameters
	std::map<std::string, std::string> fMap; ///< Map to store key names and their values to set

	// Slicer parameters
	Double_t fSlicerClusterDistance; ///< Max distance in cm between hits in the slices
	UInt_t fSlicerMinSize;			 ///< Minimum number of hits for a slice to be formed
	Double_t fSlicerChargeCut;		 ///< Only consider digits above the charge threshold
	Double_t fSlicerTimeCut;		 ///< Maximum gap allowed between hits when ordered by time in order to be associated with the previous hit.
	Bool_t fIterateSlicing;			 ///< Iterate the charge cut to make as many slices as we have tracks

	// Veto slicing parameters
	Double_t fVetoClusterDistance; ///< Max distance in cm between hits in the slices
	UInt_t fVetoMinSize;		   ///< Minimum number of hits for a slice to be formed
	Double_t fVetoMinChargeCut;	   ///< Only consider digits above the charge threshold, set initially
	Double_t fVetoMaxChargeCut;	   ///< at the min value, and iterates up to the max value.
	Double_t fVetoTimeCut;		   ///< Maximum gap allowed between hits when ordered by time in order to be associated with the previous hit.

	// Integral parameters
	Bool_t fCalculateIntegrals; ///< True if charge likelihood should calculate integrals, false to look them up
	Bool_t fTruncateIntegrals;	///< True if charge likelihood should use full lookup tables for integrals, false to use one that doesn't cut off
	Bool_t fConstrainExtent;	///< True if integrals should cut off when the particle leaves the detector

	// Likelihood tuner parameters
	Bool_t fUseTransmission;		   ///< True if we should account for absorption of photons in the water
	Bool_t fUseAngularEfficiency;	   ///< True if we should account for the PMT efficiency as a function of angle
	Bool_t fUseGlassCathodeReflection; ///< True if we should account for photons being reflected off the PMT glass
	Bool_t fUseScatteringTable;		   ///< True if we should use the scattering table, false for flat 1% chance
	Bool_t fUseNewAngularEfficiency;   ///< True if we should use the new solid angle and angular efficiency
	Bool_t fUseTrackFit;			   ///< True if we should use the muon track photon emission fit instead of the emission profiles
	Double_t fQScaling;				   ///< Scaling for the predicted charge

	// Fitter parameters
	Bool_t fUseTime;			   ///< True if we should include timing information in the likelihood
	Bool_t fUseCharge;			   ///< True if we should include charge information in the likelihood
	Bool_t fEqualiseChargeAndTime; ///< After we've seeded the vertex, energy and time, we can scale the time likelihood to weight it the same as the charge
	Bool_t fSaveWCSimRootEvent;	   ///< Whether to save a full copy of the WCSimRootEvent fitted, or just a link to the original file
	std::string fDigiType;		   ///< Name of digitizer type to use in WCSimDigitizerLikelihood::DigiType_t
	Bool_t fSaveSeedInfo;		   ///< True if we should save the SeedInfo in the output file
	Bool_t fSaveStageInfo;		   ///< True if we should save the StageInfo in the output file
	Bool_t fSaveHitComparison;	   ///< True if we should save the HitComparison in the output file
	Bool_t fSaveParameters;		   ///< True if we should save the ParameterInfo in the output file

	// Speed of light parameters
	Bool_t fUseCustomParticleSpeed; ///< Normally we assume particles travel at c - this allows us to switch and set it manually
	Bool_t fUseCustomSpeedOfLight;	///< Normally we assume light travels at c/(average n) - this allows us to switch and set it manually
	Bool_t fUseFittedSpeedOfLight;	///< Normally we assume light travels at c/(average n) - this uses a fitted speed instead
	Double_t fCustomParticleSpeed;	///< The speed of the propagating particle as a fraction of c
	Double_t fCustomSpeedOfLight;	///< The speed of the propagating particle as a fraction of c
	Double_t fFittedSpeedOfLight;	///< The effective speed of light from our fit as a fraction of c

	// Fitter plots parameters
	Double_t fPlotVtxXMax;	 ///< Used to set the extent of the vtxX reco-true range in the fitter plots
	Double_t fPlotVtxYMax;	 ///< Used to set the extent of the vtxY reco-true range in the fitter plots
	Double_t fPlotVtxZMax;	 ///< Used to set the extent of the vtxZ reco-true range in the fitter plots
	Double_t fPlotVtxTMax;	 ///< Used to set the extent of the vtxT reco-true range in the fitter plots
	Double_t fPlotThetaMax;	 ///< Used to set the extent of the theta reco-true range in the fitter plots
	Double_t fPlotPhiMax;	 ///< Used to set the extent of the phi reco-true range in the fitter plots
	Double_t fPlotEnergyMax; ///< Used to set the extent of the energy reco-true range in the fitter plots

	// Ignore parameters
	Double_t fIgnoreLowCut;
	Double_t fIgnoreHighCut;

	// Not currently used parameters
	Bool_t fUseSimpleTimeResolution;  ///<
	Bool_t fUseSimpleTimeSlew;		  ///<
	Bool_t fUseSimpleRefractiveIndex; ///<

	ClassDef(WCSimParameters, 0)
};
