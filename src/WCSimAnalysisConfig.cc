#include "WCSimAnalysisConfig.hh"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

ClassImp(WCSimAnalysisConfig)

/// Static pointer to self (so we can make object a singleton)
static WCSimAnalysisConfig * fgConfig = NULL;

/**
 * Constructor
 * @TODO Make this private so that the class is a singleton
 */
WCSimAnalysisConfig::WCSimAnalysisConfig()
{
    // Initialise the member variables to some sensible
    // values.  Note that these will be overriden in a 
    // few lines time.
    fCalculateIntegrals        = false;
    fTruncateIntegrals         = false;
    fConstrainExtent           = true;
    fDigiType                  = "";
    fUseTransmission           = true;
    fUseAngularEfficiency      = true;
    fUseGlassCathodeReflection = true;
    fUseCharge                 = true;
    fUseTime                   = true; 
    fUseScatteringTable        = true;
    fSaveWCSimRootEvent  = false;

    fUseCustomParticleSpeed    = false;
    fCustomParticleSpeed       = 1.0;
    fUseCustomSpeedOfLight     = false;
    fCustomSpeedOfLight        = 1.0;
    fUseFittedSpeedOfLight     = true;
    fFittedSpeedOfLight        = 0.71864;
    fEqualiseChargeAndTime     = false;

    fSaveSeedInfo			   = false;
    fSaveStageInfo			   = false;
    fSaveHitComparison		   = false;

    // Load the file with the default configuration
    // and update the member variables above as necessary
    fConfName = getenv("WCSIMANAHOME");
    fConfName.append("/config/default.cfg");
    this->SetConfigFile(fConfName.c_str());
    return;
}

void WCSimAnalysisConfig::SetConfigFile(const char* config)
{
    std::cout << "Setting config file" << std::endl;
    fConfName = config;
    this->LoadConfig();
    return;
}

void WCSimAnalysisConfig::Print()
{
    std::cout << "*** WCSimAnalysisConfig::Print() *** " << std::endl;
    std::cout << "fCalculateIntegrals = " << fCalculateIntegrals << std::endl
              << "fTruncateIntegrals  = " << fTruncateIntegrals << std::endl
              << "fConstrainExtent = " << fConstrainExtent << std::endl
              << "fDigiType = " << fDigiType << std::endl
              << "fUseTransmission = " << fUseTransmission << std::endl
              << "fUseAngularEfficiency = " << fUseAngularEfficiency << std::endl
              << "fUseGlassCathodeReflection = " << fUseGlassCathodeReflection << std::endl
              << "fUseTime = " << fUseTime << std::endl
              << "fUseCharge = " << fUseCharge << std::endl
              << "fUseScatteringTable = " << fUseScatteringTable << std::endl
              << "fEqualiseChargeAndTime = " << fEqualiseChargeAndTime << std::endl
    		  << "fSaveSeedInfo = " << fSaveSeedInfo << std::endl
    		  << "fSaveStageInfo = " << fSaveStageInfo << std::endl
    		  << "fSaveHitComparison = " << fSaveHitComparison << std::endl;
    return; 
}

Bool_t WCSimAnalysisConfig::GetCalculateIntegrals() const
{
  return fCalculateIntegrals;
}

Bool_t WCSimAnalysisConfig::GetTruncateIntegrals() const
{
  return fTruncateIntegrals;
}

Bool_t WCSimAnalysisConfig::GetConstrainExtent() const
{
  return fConstrainExtent;
}

std::string WCSimAnalysisConfig::GetDigiType() const
{
  return fDigiType;
}

Bool_t WCSimAnalysisConfig::GetUseTransmission() const
{
    return fUseTransmission;
}

Bool_t WCSimAnalysisConfig::GetUseAngularEfficiency() const
{
    return fUseAngularEfficiency;
}

Bool_t WCSimAnalysisConfig::GetUseGlassCathodeReflection() const
{
    return fUseGlassCathodeReflection;
}

Bool_t WCSimAnalysisConfig::GetUseChargeOnly() const
{
    return (fUseCharge && !fUseTime);
}
Bool_t WCSimAnalysisConfig::GetUseTimeOnly() const
{
    return (fUseTime && !fUseCharge);
}
Bool_t WCSimAnalysisConfig::GetUseChargeAndTime() const
{
    return (fUseCharge && fUseTime);
}

Bool_t WCSimAnalysisConfig::GetUseCharge() const
{
    return fUseCharge;
}
Bool_t WCSimAnalysisConfig::GetUseTime() const
{
    return fUseTime;
}

Bool_t WCSimAnalysisConfig::GetUseScatteringTable() const
{
  return fUseScatteringTable;
}

Bool_t WCSimAnalysisConfig::GetUseCustomParticleSpeed() const
{
	return fUseCustomParticleSpeed;
}

Double_t WCSimAnalysisConfig::GetCustomParticleSpeed() const
{
	return fCustomParticleSpeed;
}

Bool_t WCSimAnalysisConfig::GetUseCustomSpeedOfLight() const
{
	return fUseCustomSpeedOfLight;
}

Double_t WCSimAnalysisConfig::GetCustomSpeedOfLight() const
{
	return fCustomSpeedOfLight;
}

Bool_t WCSimAnalysisConfig::GetUseFittedSpeedOfLight() const
{
	return fUseFittedSpeedOfLight;
}

Double_t WCSimAnalysisConfig::GetFittedSpeedOfLight() const
{
	return fFittedSpeedOfLight;
}

Bool_t WCSimAnalysisConfig::GetEqualiseChargeAndTime() const
{
    return fEqualiseChargeAndTime;
}

Bool_t WCSimAnalysisConfig::GetSaveWCSimRootEvent() const
{
    return fSaveWCSimRootEvent;
}

Bool_t WCSimAnalysisConfig::GetSaveSeedInfo() const
{
    return fSaveSeedInfo;
}

Bool_t WCSimAnalysisConfig::GetSaveStageInfo() const
{
    return fSaveStageInfo;
}

Bool_t WCSimAnalysisConfig::GetSaveHitComparison() const
{
    return fSaveHitComparison;
}

void WCSimAnalysisConfig::LoadConfig()
{
    // Open the file and test opening:
    std::cout << "Loading config file" << std::endl;
    std::ifstream inFile;
    inFile.open(fConfName.c_str());
    
    if( !inFile.is_open() )
    {
        std::cerr << "Error: " << __FILE__ << "  " << __LINE__ << " - could not open " << fConfName << std::endl;
        exit(EXIT_FAILURE);
    }   
    
    // Read the file
    std::string line;
    Int_t lineNum = 0;
    while(getline(inFile, line))
    {
      // std::cout << "Line is: " << line << std::endl;
      // std::cout << "lineNum = " << lineNum << std::endl;
        this->IgnoreComments(line);
        this->ParseLine(line, ++lineNum);
    }
    this->SetFromMap();
    
    return;
}

/// Erases comments beginning with // or #
void WCSimAnalysisConfig::IgnoreComments( std::string &str )
{
    // std::cout << "IgnoreComments" << std::endl;
    if( str.find("//") != str.npos) str.erase(str.find("//"));
    if( str.find("#") != str.npos) str.erase(str.find("#"));
    // std::cout << "Without comments " << str << std::endl;
    return;
}

const Bool_t WCSimAnalysisConfig::IsBlankLine(std::string str)
{
    // std::cout << "IsBlankLink" << "   " << str << std::endl;
    if( str.find_first_not_of(' ') == str.npos) return true;
    else return false;
}

const Bool_t WCSimAnalysisConfig::IsGoodLine( std::string str)
{
    // std::cout << "IsGoodLine: " << str << std::endl;

    Bool_t haveEquals = false;
    Bool_t haveLHS    = false;
    Bool_t haveRHS    = false;

    // Look for an equals sign
    if( str.find("=") == str.npos)
    {
        std::cout << "No \"=\" sign found in string: " << std::endl
                  << str << std::endl;
    }
    else haveEquals = true;
    
    // Look for text on the LHS of = sign:
    std::string tempStr = str;
    tempStr.erase(0, tempStr.find_first_not_of("\t "));
    if( tempStr[0] != '0') { haveLHS = true; }
    
    // Look for text on RHS of = sign:
    tempStr = str;
    for( UInt_t rhs = tempStr.find("=") + 1; rhs < tempStr.length(); ++rhs)
    {
        if(tempStr[rhs] != '\t' && tempStr[rhs] != ' ') haveRHS = true;
    }
    // std::cout << "haveEquals = " << haveEquals << std::endl
    //           << "haveLHS = " << haveLHS << std::endl
    //           << "haveRHS = " << haveRHS << std::endl
    //           << "All = " << (haveLHS && haveRHS && haveEquals) << std::endl;

    return (haveEquals && haveLHS && haveRHS);
}

void WCSimAnalysisConfig::ExtractPair(std::string &lhs, std::string &rhs, std::string str)
{
    // std::cout << "ExtractPair" << std::endl;
    UInt_t splitPos = str.find("=");
    
    // Get left hand side of = sign and strip whitespace
    lhs = str.substr(0, splitPos);
    lhs.erase(std::remove(lhs.begin(), lhs.end(), ' '), lhs.end());
    lhs.erase(std::remove(lhs.begin(), lhs.end(), '\t'), lhs.end());
    
    // And the other side
    rhs = str.substr(splitPos+1);
    rhs.erase(std::remove(rhs.begin(), rhs.end(), ' '), rhs.end());
    rhs.erase(std::remove(rhs.begin(), rhs.end(), '\t'), rhs.end()); 
    // std::cout << "str = " << str << std::endl
    //           << "lhs = " << lhs << std::endl
    //           << "rhs = " << rhs << std::endl;
    return;
}

void WCSimAnalysisConfig::ParseLine(std::string str, Int_t lineNum)
{
    // std::cout << "ParseLine: " << str << std::endl;
    if(!(this->IsGoodLine(str)) ) 
    {
        std::cerr << "Error: line " << lineNum << "has the wrong format" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    std::string lhs, rhs;
    this->ExtractPair(lhs, rhs, str);
    this->AddToMap(lhs, rhs);
    return;
}

void WCSimAnalysisConfig::AddToMap(std::string lhs, std::string rhs)
{
    // std::cout << "AddToMap" << std::endl;
    if( fMap.find(lhs) == fMap.end() ){ fMap[lhs] = rhs; }
    else
    { 
      std::cerr << "Error: map already contains the key: " << lhs << std::endl;
      exit(EXIT_FAILURE);
    }
    std::cout << lhs << "  " << fMap[lhs] << std::endl;
    return;
}

void WCSimAnalysisConfig::SetFromMap()
{
    // std::cout << "SetFromMap" << std::endl;
    std::map<std::string, std::string>::const_iterator itr = fMap.begin();

    // Loop through the map, checking if any of the keys correspond to things we can set
    for( ; itr != fMap.end(); ++itr)
    {
      std::cout << (*itr).first << "   " << (*itr).second << std::endl;
        if((*itr).first.compare("CalculateIntegrals") == 0)
        {
            if((*itr).second.compare("true") == 0 || (*itr).second.compare("1") == 0)
            {
                fCalculateIntegrals = true;
            }
            else if((*itr).second.compare("false") == 0 || (*itr).second.compare("0") == 0)
            {
                fCalculateIntegrals = false;
            }
            else
            {
              std::cerr << "Error: " << (*itr).first << " = " << (*itr).second 
                           << " should equal true/false or 1/0" << std::endl;
              exit(EXIT_FAILURE);
            }
        }
        if((*itr).first.compare("TruncateIntegrals") == 0)
        {
            if((*itr).second.compare("true") == 0 || (*itr).second.compare("1") == 0)
            {
                fTruncateIntegrals = true;
            }
            else if((*itr).second.compare("false") == 0 || (*itr).second.compare("0") == 0)
            {
                fTruncateIntegrals = false;
            }
            else
            {
              std::cerr << "Error: " << (*itr).first << " = " << (*itr).second 
                           << " should equal true/false or 1/0" << std::endl;
              exit(EXIT_FAILURE);
            }
        }
        else if((*itr).first.compare("ConstrainExtent") == 0)
        {
            if((*itr).second.compare("true") == 0 || (*itr).second.compare("1") == 0)
            {
                fConstrainExtent = true;
            }
            else if((*itr).second.compare("false") == 0 || (*itr).second.compare("0") == 0)
            {
                fConstrainExtent = false;
            }
            else
            {
              std::cerr << "Error: " << (*itr).first << " = " << (*itr).second 
                           << " should equal true/false or 1/0" << std::endl;
              exit(EXIT_FAILURE);
            }
        }
        else if((*itr).first.compare("UseTransmission") == 0)
        {
            if((*itr).second.compare("true") == 0 || (*itr).second.compare("1") == 0)
            {
                fUseTransmission = true;
            }
            else if((*itr).second.compare("false") == 0 || (*itr).second.compare("0") == 0)
            {
                fUseTransmission = false;
            }
            else
            {
              std::cerr << "Error: " << (*itr).first << " = " << (*itr).second 
                           << " should equal true/false or 1/0" << std::endl;
              exit(EXIT_FAILURE);
            }
        }
        else if((*itr).first.compare("UseAngularEfficiency") == 0)
        {
            if((*itr).second.compare("true") == 0 || (*itr).second.compare("1") == 0)
            {
                fUseAngularEfficiency = true;
            }
            else if((*itr).second.compare("false") == 0 || (*itr).second.compare("0") == 0)
            {
                fUseAngularEfficiency = false;
            }
            else
            {
              std::cerr << "Error: " << (*itr).first << " = " << (*itr).second 
                           << " should equal true/false or 1/0" << std::endl;
              exit(EXIT_FAILURE);
            }
        }
        else if((*itr).first.compare("UseGlassCathodeReflection") == 0)
        {
            if((*itr).second.compare("true") == 0 || (*itr).second.compare("1") == 0)
            {
                fUseGlassCathodeReflection = true;
            }
            else if((*itr).second.compare("false") == 0 || (*itr).second.compare("0") == 0)
            {
                fUseGlassCathodeReflection = false;
            }
            else
            {
              std::cerr << "Error: " << (*itr).first << " = " << (*itr).second 
                           << " should equal true/false or 1/0" << std::endl;
              exit(EXIT_FAILURE);
            }
        }
        else if((*itr).first.compare("UseCharge") == 0)
        {
            if((*itr).second.compare("true") == 0 || (*itr).second.compare("1") == 0)
            {
                fUseCharge = true;
            }
            else if((*itr).second.compare("false") == 0 || (*itr).second.compare("0") == 0)
            {
                fUseCharge = false;
            }
            else
            {
              std::cerr << "Error: " << (*itr).first << " = " << (*itr).second 
                           << " should equal true/false or 1/0" << std::endl;
              exit(EXIT_FAILURE);
            }
        }
        else if((*itr).first.compare("UseTime") == 0)
        {
            if((*itr).second.compare("true") == 0 || (*itr).second.compare("1") == 0)
            {
                fUseTime = true;
            }
            else if((*itr).second.compare("false") == 0 || (*itr).second.compare("0") == 0)
            {
                fUseTime = false;
            }
            else
            {
              std::cerr << "Error: " << (*itr).first << " = " << (*itr).second 
                           << " should equal true/false or 1/0" << std::endl;
              exit(EXIT_FAILURE);
            }
        }
        else if((*itr).first.compare("DigiType") == 0)
        {
            std::cout << "It's " << ((*itr).second) << std::endl;
            fDigiType = ((*itr).second);
        }
        else if((*itr).first.compare("UseScatteringTable") == 0)
        {
             if((*itr).second.compare("true") == 0 || (*itr).second.compare("1") == 0)
             {
                 fUseScatteringTable = true;
             }
             else if((*itr).second.compare("false") == 0 || (*itr).second.compare("0") == 0)
             {
                 fUseScatteringTable = false;
             }
             else
             {
               std::cerr << "Error: " << (*itr).first << " = " << (*itr).second 
                            << " should equal true/false or 1/0" << std::endl;
               exit(EXIT_FAILURE);
             }
        }
        else if((*itr).first.compare("UseCustomParticleSpeed") == 0)
        {
             if((*itr).second.compare("true") == 0 || (*itr).second.compare("1") == 0)
             {
                 fUseCustomParticleSpeed = true;
             }
             else if((*itr).second.compare("false") == 0 || (*itr).second.compare("0") == 0)
             {
                 fUseCustomParticleSpeed = false;
             }
             else
             {
               std::cerr << "Error: " << (*itr).first << " = " << (*itr).second 
                            << " should equal true/false or 1/0" << std::endl;
               exit(EXIT_FAILURE);
             }
        }
        else if((*itr).first.compare("CustomParticleSpeed") == 0)
        {
             std::stringstream ss(itr->second);
             float speed = 0.0;
             if( !(ss >> speed) )
             {
               std::cerr << "Error: " << (*itr).second << " = " << (*itr).second 
                            << " should be a double" << std::endl;
               exit(EXIT_FAILURE);
             }
             fCustomParticleSpeed = speed;
        }
        else if((*itr).first.compare("UseFittedSpeedOfLight") == 0)
        {
             if((*itr).second.compare("true") == 0 || (*itr).second.compare("1") == 0)
             {
                 fUseFittedSpeedOfLight = true;
             }
             else if((*itr).second.compare("false") == 0 || (*itr).second.compare("0") == 0)
             {
                 fUseFittedSpeedOfLight = false;
             }
             else
             {
               std::cerr << "Error: " << (*itr).first << " = " << (*itr).second 
                            << " should equal true/false or 1/0" << std::endl;
               exit(EXIT_FAILURE);
             }
        }
        else if((*itr).first.compare("FittedSpeedOfLight") == 0)
        {
             std::stringstream ss(itr->second);
             float speed = 0.0;
             if( !(ss >> speed) )
             {
               std::cerr << "Error: " << (*itr).second << " = " << (*itr).second 
                            << " should be a double" << std::endl;
               exit(EXIT_FAILURE);
             }
             fFittedSpeedOfLight = speed;
        }
        else if((*itr).first.compare("EqualiseChargeAndTime") == 0)
        {
             if((*itr).second.compare("true") == 0 || (*itr).second.compare("1") == 0)
             {
                 fEqualiseChargeAndTime = true;
             }
             else if((*itr).second.compare("false") == 0 || (*itr).second.compare("0") == 0)
             {
                 fEqualiseChargeAndTime = false;
             }
             else
             {
               std::cerr << "Error: " << (*itr).first << " = " << (*itr).second 
                            << " should equal true/false or 1/0" << std::endl;
               exit(EXIT_FAILURE);
             }
        }
        else if((*itr).first.compare("SaveWCSimRootEvent") == 0)
        {
             if((*itr).second.compare("true") == 0 || (*itr).second.compare("1") == 0)
             {
                 fSaveWCSimRootEvent = true;
             }
             else if((*itr).second.compare("false") == 0 || (*itr).second.compare("0") == 0)
             {
                 fSaveWCSimRootEvent = false;
             }
             else
             {
               std::cerr << "Error: " << (*itr).first << " = " << (*itr).second 
                            << " should equal true/false or 1/0" << std::endl;
               exit(EXIT_FAILURE);
             }
        }
        else if((*itr).first.compare("SaveSeedInfo") == 0)
        {
             if((*itr).second.compare("true") == 0 || (*itr).second.compare("1") == 0)
             {
                 fSaveSeedInfo = true;
             }
             else if((*itr).second.compare("false") == 0 || (*itr).second.compare("0") == 0)
             {
            	 fSaveSeedInfo = false;
             }
             else
             {
               std::cerr << "Error: " << (*itr).first << " = " << (*itr).second
                            << " should equal true/false or 1/0" << std::endl;
               exit(EXIT_FAILURE);
             }
        }
        else if((*itr).first.compare("SaveStageInfo") == 0)
        {
             if((*itr).second.compare("true") == 0 || (*itr).second.compare("1") == 0)
             {
                 fSaveStageInfo = true;
             }
             else if((*itr).second.compare("false") == 0 || (*itr).second.compare("0") == 0)
             {
            	 fSaveStageInfo = false;
             }
             else
             {
               std::cerr << "Error: " << (*itr).first << " = " << (*itr).second
                            << " should equal true/false or 1/0" << std::endl;
               exit(EXIT_FAILURE);
             }
        }
        else if((*itr).first.compare("SaveHitComparison") == 0)
        {
             if((*itr).second.compare("true") == 0 || (*itr).second.compare("1") == 0)
             {
                 fSaveHitComparison = true;
             }
             else if((*itr).second.compare("false") == 0 || (*itr).second.compare("0") == 0)
             {
            	 fSaveHitComparison = false;
             }
             else
             {
               std::cerr << "Error: " << (*itr).first << " = " << (*itr).second
                            << " should equal true/false or 1/0" << std::endl;
               exit(EXIT_FAILURE);
             }
        }
    }
    fMap.clear();
}


// Destructor
WCSimAnalysisConfig::~WCSimAnalysisConfig() 
{
}

WCSimAnalysisConfig* WCSimAnalysisConfig::Instance()
{
  // std::cout << " Instance ! " << std::endl;
  if( fgConfig == NULL ){ fgConfig = new WCSimAnalysisConfig(); }
  assert(fgConfig);
  if( fgConfig ){ }

  return fgConfig;
}

void WCSimAnalysisConfig::SetUseTimeOnly(Bool_t doIt)
{
  if(doIt)
  {
    std::cout << " *** WCSimAnalysisConfig::SetUseTimeOnly *** " << std::endl;
    fUseTime = true;
    fUseCharge = false;
  }
  return;
}

void WCSimAnalysisConfig::SetUseChargeOnly(Bool_t doIt)
{
  if(doIt)
  {
    std::cout << " *** WCSimAnalysisConfig::SetUseChargeOnly *** " << std::endl;
    fUseCharge = true;
    fUseTime = false;
  }
  return;
}

void WCSimAnalysisConfig::SetUseChargeAndTime(Bool_t doIt)
{
  if(doIt)
  {
    std::cout << " *** WCSimAnalysisConfig::SetUseChargeAndTime *** " << std::endl;
    fUseCharge = true;
    fUseTime = true;
  }
  return;
}

void WCSimAnalysisConfig::SetUseCustomParticleSpeed(Bool_t doIt)
{
  if(doIt)
  {
    std::cout << " *** WCSimAnalysisConfig::SetUseCustomParticleSpeed *** " << std::endl;
    fUseCustomParticleSpeed = doIt;
  }
  return;
}

void WCSimAnalysisConfig::SetUseCustomSpeedOfLight(Bool_t doIt)
{
  if(doIt)
  {
    std::cout << " *** WCSimAnalysisConfig::SetUseCustomSpeedOfLight *** " << std::endl;
    fUseCustomSpeedOfLight = doIt;
  }
  return;
}

void WCSimAnalysisConfig::SetCustomParticleSpeed(const Double_t& speed)
{
  std::cout << " *** WCSimAnalysisConfig::SetCustomParticleSpeed to " << speed << "c ***" << std::endl;
  fCustomParticleSpeed = speed;
  return;
}

void WCSimAnalysisConfig::SetCustomSpeedOfLight(const Double_t& speed)
{
  std::cout << " *** WCSimAnalysisConfig::SetCustomSpeedOfLight to " << speed << "c ***" << std::endl;
  fCustomSpeedOfLight = speed;
  return;
}

void WCSimAnalysisConfig::SetEqualiseChargeAndTime(Bool_t doIt)
{
    std::cout << " *** WCSimAnalysisConfig::SetEqualiseChargeAndTime to " << doIt << std::endl;
    fEqualiseChargeAndTime = doIt;
}

void WCSimAnalysisConfig::SetSaveWCSimRootEvent(Bool_t doIt)
{
    std::cout << " *** WCSimAnalysisConfig::SetSaveWCSimRootEvent to " << doIt << std::endl;
    fSaveWCSimRootEvent = doIt;
}

void WCSimAnalysisConfig::SetSaveSeedInfo(Bool_t doIt)
{
    std::cout << " *** WCSimAnalysisConfig::SetSaveSeedInfo to " << doIt << std::endl;
    fSaveSeedInfo = doIt;
}

void WCSimAnalysisConfig::SetSaveStageInfo(Bool_t doIt)
{
    std::cout << " *** WCSimAnalysisConfig::SetSaveStageInfo to " << doIt << std::endl;
    fSaveStageInfo = doIt;
}

void WCSimAnalysisConfig::SetSaveHitComparison(Bool_t doIt)
{
    std::cout << " *** WCSimAnalysisConfig::SetSaveHitComparison to " << doIt << std::endl;
    fSaveHitComparison = doIt;
}


