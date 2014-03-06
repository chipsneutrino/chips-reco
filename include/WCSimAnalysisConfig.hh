/* 
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
    static WCSimAnalysisConfig * Instance();
    void SetConfigFile(const char * config);
    void Print();
    
    Bool_t GetCalculateIntegrals() const;
    Bool_t GetConstrainExtent() const;
    std::string GetDigiType() const;
    
    

private:

    void LoadConfig();
    void IgnoreComments(std::string &str);
    
    const Bool_t IsGoodLine(std::string str);
    const Bool_t IsBlankLine(std::string str);
    void ExtractPair(std::string &lhs, std::string &rhs, std::string str);
    void ParseLine(std::string, Int_t lineNum);
    void AddToMap(std::string lhs, std::string rhs);
    void SetFromMap();
    
    const char *  fConfName;
    Bool_t        fCalculateIntegrals;
    Bool_t        fConstrainExtent;
    std::string   fDigiType;

    std::map<std::string, std::string> fMap;
    
    ClassDef(WCSimAnalysisConfig,0)
};

#endif	/* WCSIMANALYSISCONFIG_HH */

