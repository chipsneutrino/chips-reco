/*
 * WCSimEmissionProfileGenerator.hh
 *
 *  Created on: 28 Mar 2014
 *      Author: andy
 */

#ifndef WCSIMEMISSIONPROFILEGENERATOR_HH_
#define WCSIMEMISSIONPROFILEGENERATOR_HH_

#include "TChain.h"
#include "TString.h"
#include <map>


class WCSimEmissionProfileGenerator
{
public:
    WCSimEmissionProfileGenerator( const char * treeName, Int_t particleType );

    virtual ~WCSimEmissionProfileGenerator();

    void AddFile(Double_t energy, Int_t particleType, const char * filename);
    void SetSaveName( const char * savename);
    std::map<Double_t, TChain *> GetChain();
    TString GetChainSaveName();
    void BuildSaveName();
    void Concatenate();
    void ResetEnergy();
    void MakeFluxHistogram();
    TChain * Next();

private:

    const char * fTreeName;
    const char * fSaveName;
    std::map<Double_t, TString> fNameMap;
    Int_t fParticle;
    TString fName;
    std::map<Double_t, TChain *> fChain;
    std::map<Double_t, TChain *>::const_iterator fChainItr;
    Double_t fEnergy;
    Int_t fType;
    Bool_t fIsFirstCall;

    UInt_t fUNKNOWN;

    ClassDef(WCSimEmissionProfileGenerator,1)
};


#endif /* WCSIMEMISSIONPROFILEGENERATOR_HH_ */

