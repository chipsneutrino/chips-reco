/*
 * WCSimConcatenateEmissionProfiles.cc
 *
 *  Created on: 10 Mar 2014
 *      Author: andy
 */

#include "WCSimConcatenateEmissionProfiles.hh"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TObjArray.h"
#include "TObject.h"

#include <iostream>
#include <map>
#include <cassert>
#include <cstdlib>

WCSimConcatenateEmissionProfiles::WCSimConcatenateEmissionProfiles( const char * saveName )
{
    fSaveFile    = new TFile(saveName, "RECREATE");


    histArray    = new TObjArray();
    histArray->SetOwner(kTRUE);

    angHistArray = new TObjArray();
    angHistArray->SetOwner(kTRUE);

    hWhichHisto  = new TH1D();
    hWhichHisto->SetNameTitle("hWhichHisto","Histogram containing energy binning scheme");

    fPdgTypeSet  = false;

}

WCSimConcatenateEmissionProfiles::~WCSimConcatenateEmissionProfiles()
{
    delete histArray;
    delete angHistArray;
    delete hWhichHisto;
    delete fSaveFile;

}

void WCSimConcatenateEmissionProfiles::AddFile(const char* fileName,
        Double_t energy, Int_t pdgType)
{
    TH1D * hRho = NULL;
    TH2D * hG   = NULL;

    try
    {
        std::map<Double_t,TH1D*>::const_iterator rhoIt = rhoMap.find(energy);
        if( rhoIt != rhoMap.end() ){ throw 1; }

        std::map<Double_t,TH2D*>::const_iterator gIt = gMap.find(energy);
        if( gIt != gMap.end() ){ throw 2; }

        if(fPdgTypeSet && pdgType != fPdgType)
        {
            throw 3;
            fPdgType    = pdgType;
            fPdgTypeSet = true;
        }

        TFile f(fileName,"READ");
        if(f.IsZombie() ) { throw fileName; }

        f.GetObject("f_hS", hRho);
        if(hRho == NULL){ throw 0; }
        hRho->SetDirectory(0);

        hG          = (TH2D*)f.Get("f_hG");
        if(hG == NULL){ throw 0; }
        hG->SetDirectory(0);

    }
    catch(const char * n)
    {
        std::cout << "Could not open file: " << fileName << std::endl;
        abort();
    }
    catch( const Int_t i)
    {
        if(i == 0) { std::cout << "Null pointer returned when getting histograms"       << std::endl; }
        else if(i == 1) { std::cout << "Already have an entry at this energy in rhoMap" << std::endl; }
        else if(i == 2) { std::cout << "Already have an entry at this energy in gMap"   << std::endl; }
        else if(i == 3) { std::cout << "Pdg code " << pdgType << " differs from other files ("
                                    << fPdgType << ")" << std::endl; }
        else { std::cout << "Unknown exception caught.  This is the default exception message"  << std::endl; }
        abort();
    }

    rhoMap[energy] = hRho;
    gMap[energy]   = hG;
}

void WCSimConcatenateEmissionProfiles::BuildArrays()
{
    assert(gMap.size() == rhoMap.size());
    std::cout << "Size of rhoMap = " << rhoMap.size() << std::endl;

    Double_t * energyBins = new Double_t[rhoMap.size() + 1];

    histArray->Clear();
    histArray->Expand( rhoMap.size() );
    angHistArray->Clear();
    angHistArray->Expand( gMap.size() );

    std::map<Double_t, TH1D*>::const_iterator rhoItr;
    std::map<Double_t, TH1D*>::const_iterator rhoBegin = rhoMap.begin();
    std::map<Double_t, TH1D*>::const_iterator rhoEnd   = rhoMap.end();

    for( rhoItr = rhoBegin ; rhoItr != rhoEnd ; ++rhoItr )
    {

        UInt_t distance = std::distance( rhoBegin, rhoItr );
        std::cout << "distance = " << distance << std::endl;
        histArray->AddAt( static_cast<TObject*>((*rhoItr).second), distance );
        (*rhoItr).second->Draw();
        std::cout << (*rhoItr).first << "  " << (*rhoItr).second << std::endl;
        energyBins[distance] = (*rhoItr).first;
    }


    std::map<Double_t, TH2D*>::const_iterator gItr;
    std::map<Double_t, TH2D*>::const_iterator gBegin = gMap.begin();
    std::map<Double_t, TH2D*>::const_iterator gEnd   = gMap.end();

    for( gItr = gBegin ; gItr != gEnd ; ++gItr )
    {
        UInt_t distance = std::distance( gBegin, gItr );
        angHistArray->AddAt( (TObject*)(*gItr).second, distance );
        assert( energyBins[distance] == (*gItr).first );
    }

    energyBins[0] = 0.0;
    energyBins[rhoMap.size()] = 30000.0;
    hWhichHisto->GetXaxis()->Set( rhoMap.size(), energyBins );

    delete [] energyBins;
}

void WCSimConcatenateEmissionProfiles::SaveArrays()
{
    std::cout << "Size of histarray = " << histArray->GetSize() << std::endl;

    fSaveFile->cd();
    hWhichHisto->Write();
    angHistArray->Write("angHistArray",TObject::kSingleKey);
    histArray->Write("histArray",TObject::kSingleKey);
    fSaveFile->Close();

}

void WCSimConcatenateEmissionProfiles::Run()
{
    BuildArrays();
    SaveArrays();
    return;
}
