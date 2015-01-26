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
    fConstrainExtent           = true;
    fDigiType                  = "";
    fUseTransmission           = true;
    fUseAngularEfficiency      = true;
    fUseGlassCathodeReflection = true;
    fUseCharge                 = true;
    fUseTime                   = true; 
    fUseHoughFitterForSeed			   = true;
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
              << "fConstrainExtent = " << fConstrainExtent << std::endl
              << "fDigiType = " << fDigiType << std::endl
              << "fUseTransmission = " << fUseTransmission << std::endl
              << "fUseAngularEfficiency = " << fUseAngularEfficiency << std::endl
              << "fUseGlassCathodeReflection = " << fUseGlassCathodeReflection << std::endl
              << "fUseTime = " << fUseTime << std::endl
              << "fUseCharge = " << fUseCharge << std::endl
    		  << "fSeedWithHough = " << fUseHoughFitterForSeed << std::endl;
    return; 
}

Bool_t WCSimAnalysisConfig::GetCalculateIntegrals() const
{
  return fCalculateIntegrals;
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

Bool_t WCSimAnalysisConfig::GetUseHoughFitterForSeed() const
{
	return fUseHoughFitterForSeed;
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
        else if((*itr).first.compare("UseHoughFitterForSeed") == 0)
        {
            if((*itr).second.compare("true") == 0 || (*itr).second.compare("1") == 0)
            {
                fUseHoughFitterForSeed = true;
            }
            else if((*itr).second.compare("false") == 0 || (*itr).second.compare("0") == 0)
            {
                fUseHoughFitterForSeed = false;
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
