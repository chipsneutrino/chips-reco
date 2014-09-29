/*
 * WCSimConcatenateEmissionProfiles.hh
 *
 *  Created on: 10 Mar 2014
 *      Author: andy
 */

#ifndef WCSIMCONCATENATEEMISSIONPROFILES_HH_
#define WCSIMCONCATENATEEMISSIONPROFILES_HH_

#include <map>
#include "Rtypes.h"

class TFile;
class TH1D;
class TH2D;
class TObjArray;

class WCSimConcatenateEmissionProfiles
{
public:
    WCSimConcatenateEmissionProfiles(const char * saveName = "emissionProfiles.root");
    virtual ~WCSimConcatenateEmissionProfiles();

    void AddFile(const char * fileName, Double_t energy, Int_t pdgType);
    void BuildArrays();
    void SaveArrays();
    void Run();

private:
    Int_t       fPdgType;
    Bool_t      fPdgTypeSet;

    TObjArray * histArray;
    TObjArray * angHistArray;
    TH1D      * hWhichHisto;
    std::map<Double_t, TH1D*> rhoMap;
    std::map<Double_t, TH2D*> gMap;

    TFile * fSaveFile;
};

#endif /* WCSIMCONCATENATEEMISSIONPROFILES_HH_ */
