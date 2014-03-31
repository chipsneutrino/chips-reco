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
    

    const char *  fConfName;                   ///< Path of a text file containing all the configuration parameters
    Bool_t        fCalculateIntegrals;         ///< True if charge likelihood should calculate integrals, false to look them up
    Bool_t        fConstrainExtent;            ///< True if integrals should cut off when the particle leaves the detector
    Bool_t        fUseTransmission;            ///< True if we should account for absorption of photons in the water
    Bool_t        fUseAngularEfficiency;       ///< True if we should account for the PMT efficiency as a function of angle
    Bool_t        fUseGlassCathodeReflection;  ///< True if we should account for photons being reflected off the PMT glass
    std::string   fDigiType;            ///< Name of digitizer type to use in WCSimDigitizerLikelihood::DigiType_t

    std::map<std::string, std::string> fMap; ///< Map to store key names and their values to set
    
    ClassDef(WCSimAnalysisConfig,0)
};

#endif	/* WCSIMANALYSISCONFIG_HH */

