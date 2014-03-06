#include "WCSimAnalysisConfig.hh"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

ClassImp(WCSimAnalysisConfig)

static WCSimAnalysisConfig * fgConfig = NULL;

// Constructors
WCSimAnalysisConfig::WCSimAnalysisConfig()
{
    std::cout << "The constructor" << std::endl;    
    fCalculateIntegrals = false;
    fConstrainExtent = false;
    fDigiType = "";
    this->SetConfigFile("default.cfg");
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
              << "fDigiType = " << fDigiType << std::endl;
    
}

// Getters
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


void WCSimAnalysisConfig::LoadConfig()
{
    // Open the file and test opening:
    std::cout << "Loading config file" << std::endl;
    std::ifstream inFile;
    inFile.open(fConfName);
    
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
      std::cout << "Line is: " << line << std::endl;
      std::cout << "lineNum = " << lineNum << std::endl;
        this->IgnoreComments(line);
        this->ParseLine(line, ++lineNum);
    }
    std::cout << "Here!" << std::endl;
    this->SetFromMap();
    
    return;
}

// Erases comments beginning with // or #
void WCSimAnalysisConfig::IgnoreComments( std::string &str )
{
    std::cout << "IgnoreComments" << std::endl;
    if( str.find("//") != str.npos) str.erase(str.find("//"));
    if( str.find("#") != str.npos) str.erase(str.find("#"));
    std::cout << "Without comments " << str << std::endl;
    return;
}

const Bool_t WCSimAnalysisConfig::IsBlankLine(std::string str)
{
  std::cout << "IsBlankLink" << "   " << str << std::endl;
    if( str.find_first_not_of(' ') == str.npos) return true;
    else return false;
}

const Bool_t WCSimAnalysisConfig::IsGoodLine( std::string str)
{
  std::cout << "IsGoodLine: " << str << std::endl;

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
    std::cout << "haveEquals = " << haveEquals << std::endl
              << "haveLHS = " << haveLHS << std::endl
              << "haveRHS = " << haveRHS << std::endl
              << "All = " << (haveLHS && haveRHS && haveEquals) << std::endl;

    return (haveEquals && haveLHS && haveRHS);
}

void WCSimAnalysisConfig::ExtractPair(std::string &lhs, std::string &rhs, std::string str)
{
  std::cout << "ExtractPair" << std::endl;
    UInt_t splitPos = str.find("=");
    
    // Get left hand side of = sign and strip whitespace
    lhs = str.substr(0, splitPos);
    lhs.erase(std::remove(lhs.begin(), lhs.end(), ' '), lhs.end());
    lhs.erase(std::remove(lhs.begin(), lhs.end(), '\t'), lhs.end());
    
    // And the other side
    rhs = str.substr(splitPos+1);
    rhs.erase(std::remove(rhs.begin(), rhs.end(), ' '), rhs.end());
    rhs.erase(std::remove(rhs.begin(), rhs.end(), '\t'), rhs.end()); 
    std::cout << "str = " << str << std::endl
              << "lhs = " << lhs << std::endl
              << "rhs = " << rhs << std::endl;
    return;
}

void WCSimAnalysisConfig::ParseLine(std::string str, Int_t lineNum)
{
  std::cout << "ParseLine: " << str << std::endl;
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
  std::cout << "AddToMap" << std::endl;
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
  std::cout << "SetFromMap" << std::endl;
    std::map<std::string, std::string>::const_iterator itr = fMap.begin();
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
        else if((*itr).first.compare("DigiType") == 0)
        {
            std::cout << "It's " << ((*itr).second) << std::endl;
            fDigiType = ((*itr).second);
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

  std::cout << " Instance ! " << std::endl;
  if( fgConfig == NULL ){ fgConfig = new WCSimAnalysisConfig(); }
  assert(fgConfig);
  if( fgConfig ){ }

  return fgConfig;
}
