/* 
 * \class WCSimAnalysisConfig
 * This class is used to load several configuration parameters from
 * a text file without requiring the code to be recompile
 * It is implemented as a sigleton, accessed by the Instance() method
 *
 * File:   WCSimAnalysisConfig.hh
 * Author: andy
 *
 * Created on 10 January 2014, 12:24
 */

#ifndef WCSIMANALYSISCONFIG_HH
#define	WCSIMANALYSISCONFIG_HH

#include "TObject.h"
#include <string>
#include <map>

class WCSimAnalysisConfig : public TObject
{
public:
    WCSimAnalysisConfig();
    virtual ~WCSimAnalysisConfig();

    /// Class should be a singleton, so it has an Instance operator to access itself
    static WCSimAnalysisConfig * Instance();

    /**
     * Set the path of the configuration file to be loaded
     * @param config Path of a text file containing all the configuration parameters
     */
    void SetConfigFile(const char * config);

    /// Print everything this class will set
    void Print();
    
    /**
     * Should the charge likelihood cut off the integrals where the particle escapes?
     * We'll use different, smaller lookup tables if not
     * @return True if it should truncate them
     */
    Bool_t GetTruncateIntegrals() const;

    /**
     * Should the charge likelihood calculate or use tables for integrals?
     * @return True if it should calculate them
     */
    Bool_t GetCalculateIntegrals() const;

    /**
     * Should the integrals be cut off once the particle reaches the edge of the detector?
     * @return True if the integrals should be cut off
     */
    Bool_t GetConstrainExtent() const;

    /**
     * Which digitizer method should we use to calculate the likelihood that the observed charge
     * is measured given the predicted charge?
     * @return String containing the name of the digitizer type
     */
    std::string GetDigiType() const;
   
    /**
     * Should we account for absorption of photons in the detector water
     * @return True if we should account for transmission
     */ 
    Bool_t GetUseTransmission() const;

    /**
     * Should we account for the efficiency of the PMT as a function of incident photon angle?
     * @return True if we should account for the angular efficiency
     */
    Bool_t GetUseAngularEfficiency() const;
    
    /**
     *  Should we account for the probability that a photon is reflected at the PMT glass
     *  @return True if photons can be reflected by the PMT
     */
    Bool_t GetUseGlassCathodeReflection() const;

    /**
     * Should we use the charge to help calculate the likelihood?
     * Says nothing about using time
     * @return True if we should use the charge
     */
    Bool_t GetUseCharge() const;

    /**
     * Should we use the time to help calculate the likelihood?
     * Says nothing about using charge
     * @return True if we should use the time
     */
    Bool_t GetUseTime() const;

    /**
     * Should we only use the charge to calculate the likelihood?
     * @return True if we should only use the charge
     */
    Bool_t GetUseChargeOnly() const;

    /**
     * Should we only use the time component to calculate the likelihood?
     * @return True if we should only use the time
     */
    Bool_t GetUseTimeOnly() const;

    /**
     * Should we use the full combined (charge + time) information to calculate
     * the likelihood?
     * @return True if we should use charge and time together
     */
    Bool_t GetUseChargeAndTime() const;
    
    /**
     * Should we use the old Hough transform fitter to seed the start values?
     * @return True if we should seed using the old fitter
     */
    Bool_t GetUseHoughFitterForSeed() const;

    /**
     * Should we account for scattering by reading a scattered light percentage from a file
     * that tabulates it as a function of some geometric parameters, or just use a default
     * flat percentage
     * @return True if we should load it from the table
     */
    Bool_t GetUseScatteringTable() const;

    /**
     * Should we assume the particle propagates at some custom speed, or use the default (c)
     * for working out the time predictions?
     * @return True if we should use the custom value
     */
    Bool_t GetUseCustomParticleSpeed() const;

    /**
     * Should we assume that photons propagate at some custom speed,
     * for working out the time predictions?
     * @return True if we should use the custom value
     */
    Bool_t GetUseCustomSpeedOfLight() const;

    /**
     * @brief Get whether the fitter should save full copy of the WCSimRootEvent fitted, or a link to the original file
     *
     * @return True if the fitter saves a full copy, false if it just links to the file
     */
    Bool_t GetSaveWCSimRootEvent() const;
    /**
     * We can work out the speed of light by averaging the refractive index against the
     * photon wavelength spectrum but this doesn't account for scattering etc. that
     * effectively slows it down (by increasing path length to the PMT
     * So I fitted c/n using upward-going muons and electrons - this option uses that value instead
     * -- Custom speed of light takes precedence over fitted speed if both are set to true -- 
     * @return True if we should use the fitted value
     */
    Bool_t GetUseFittedSpeedOfLight() const;

    /**
     * Normally we assume particles propagate at c when working out the time likelihood
     * But we can use a custom speed defined in the config file: this method returns that
     * speed (as a fraction of c)
     * -- Custom speed of light takes precedence over fitted speed if both are set to true -- 
     * @return Speed particles are assumed to travel at by the time likelihood as a fraction of c
     */
    Double_t GetCustomParticleSpeed() const;

    /**
     * Normally we assume light travels at c/(average n) which we calculate beforehand
     * But we can use a custom speed defined in the config file: this method returns that
     * speed (as a fraction of c)
     * @return Speed particles are assumed to travel at by the time likelihood as a fraction of c
     */
    Double_t GetCustomSpeedOfLight() const;

    /** 
     * Return the effective propagation speed of light needed for the time fitter, produced
     * by fitting it using upward-going muons and electrons.  Refractive index slightly higher
     * than if we just worked it out from the inputs to WCSim because of scattering etc
     * @return Speed photons are assumed to travel at, according to the fit
     */
    Double_t GetFittedSpeedOfLight() const;

    Bool_t GetEqualiseChargeAndTime() const;

    void SetUseTimeOnly(Bool_t doIt = true);
    void SetUseChargeOnly(Bool_t doIt = true);
    void SetUseChargeAndTime(Bool_t doIt = true);
    void SetUseCustomParticleSpeed(Bool_t doIt = true);
    void SetUseCustomSpeedOfLight(Bool_t doIt = true);
    void SetEqualiseChargeAndTime(Bool_t doIt = true);

    /**
     * Set the custom speed at which the time likelihood assumes particles travel
     * @input speed Speed at which particle travels, as a fraction of c
     */
    void SetCustomParticleSpeed(const Double_t &speed);

    /**
     * Set the custom speed at which the time likelihood assumes optical photons travel
     * @input speed Speed at which light travels, as a fraction of c
     */
    void SetCustomSpeedOfLight(const Double_t &speed);

    
    /**
     * @brief Can either save a full copy of the WCSimRootEvent or a TNamed whose name
     * tells you where to find the original ROOT file
     *
     * @param doIt True saves the WCSimRootEvent, false saves the TNamed that links to the original
     */
    void SetSaveWCSimRootEvent(Bool_t doIt = true);




private:
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
    
    std::string   fConfName;                   ///< Path of a text file containing all the configuration parameters
    Bool_t        fCalculateIntegrals;         ///< True if charge likelihood should calculate integrals, false to look them up
    Bool_t        fTruncateIntegrals;          ///< True if charge likelihood should use full lookup tabels for integrals, false to use one that doesn't cut off
    Bool_t        fConstrainExtent;            ///< True if integrals should cut off when the particle leaves the detector
    Bool_t        fUseTransmission;            ///< True if we should account for absorption of photons in the water
    Bool_t        fUseAngularEfficiency;       ///< True if we should account for the PMT efficiency as a function of angle
    Bool_t        fUseGlassCathodeReflection;  ///< True if we should account for photons being reflected off the PMT glass
    Bool_t        fUseTime;                    ///< True if we should include timing information in the likelihood
    Bool_t        fUseCharge;                  ///< True if we should include charge information in the likelihood
    Bool_t 	      fUseHoughFitterForSeed;	     ///< True if we should use the old Hough transform fitter to seed the start values
    Bool_t        fUseScatteringTable;         ///< True if we should use the scattering table, false for flat 1% chance
    Bool_t        fUseCustomParticleSpeed;     ///< Normally we assume particles travel at c - this allows us to switch and set it manually
    Bool_t        fUseCustomSpeedOfLight;      ///< Normally we assume light travels at c/(average n) - this allows us to switch and set it manually
    Bool_t        fUseFittedSpeedOfLight;      ///< Normally we assume light travels at c/(average n) - this uses a fitted speed instead
    Bool_t        fEqualiseChargeAndTime;      ///< After we've seeded the vertex, energy and time, we can scale the time likelihood to weight it the same as the charge
    Bool_t        fSaveWCSimRootEvent;   ///< Whether to save a full copy of the WCSimRootEvent fitted, or just a link to the original file
    Double_t      fCustomParticleSpeed;        ///< The speed of the propagating particle as a fraction of c       
    Double_t      fCustomSpeedOfLight;         ///< The speed of the propagating particle as a fraction of c       
    Double_t      fFittedSpeedOfLight;         ///< The effective speed of light from our fit as a fraction of c
    std::string   fDigiType;                   ///< Name of digitizer type to use in WCSimDigitizerLikelihood::DigiType_t

    std::map<std::string, std::string> fMap;   ///< Map to store key names and their values to set
    
    ClassDef(WCSimAnalysisConfig,0)
};

#endif	/* WCSIMANALYSISCONFIG_HH */

