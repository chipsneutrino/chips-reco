/*
 * WCSimEmissionProfileGenerator.cc
 *
 *  Created on: 28 Mar 2014
 *      Author: andy
 */

#include "WCSimEmissionProfileGenerator.hh"
#include "WCSimConcatenateEmissionProfiles.hh"
#include "TAxis.h"
#include "TBranch.h"
#include "TChain.h"
#include "TDatabasePDG.h"
#include "TFile.h"
#include "TGraph.h"
#include "TParticlePDG.h"
#include "TString.h"
#include <cassert>
#include <iostream>

#ifndef REFLEX_DICTIONARY
ClassImp(WCSimEmissionProfileGenerator)
#endif

WCSimEmissionProfileGenerator::WCSimEmissionProfileGenerator( const char * treeName, Int_t particleType)
{

    fEnergy = -1;
    fUNKNOWN = (UInt_t)(-1);    
    fTreeName = treeName;
    fName = TString("");
    fSaveName = "emissionProfiles.root";
    fType = particleType;
    fIsFirstCall = true;
    // fChain = NULL;
    // TODO Auto-generated constructor stub
}

void WCSimEmissionProfileGenerator::BuildSaveName()
{
  TDatabasePDG * myPDG = new TDatabasePDG();
  fName.Form("ep_%s_%.0f.root", ((TParticlePDG*)myPDG->GetParticle(fType))->GetName(), fEnergy);
  fNameMap[fEnergy] = fName;
  delete myPDG;
  return;
}

void WCSimEmissionProfileGenerator::ResetEnergy()
{
  // if( fChain != NULL ) 
  // { 
  //   delete fChain; 
  //   fChain = NULL;
  // }
  fName = TString("");
  fEnergy = -1;
  return;
}



WCSimEmissionProfileGenerator::~WCSimEmissionProfileGenerator()
{
    // TODO Auto-generated destructor stub
}

void WCSimEmissionProfileGenerator::AddFile(Double_t energy, Int_t type, const char * fileName)
{
    fEnergy = energy;
    fType = type;
    
    if( fChain.find(energy) == fChain.end() )
    {
      fChain[energy] = new TChain("ntuple");
      BuildSaveName();
    }

    fChain[energy]->Add(fileName);
    return;
}

std::map<Double_t, TChain *> WCSimEmissionProfileGenerator::GetChain()
{
  return fChain;
}

TString WCSimEmissionProfileGenerator::GetChainSaveName()
{
  return fName;
}

void WCSimEmissionProfileGenerator::Concatenate()
{
    WCSimConcatenateEmissionProfiles conc( fSaveName );
    std::map<Double_t, TString>::iterator nameItr = fNameMap.begin();

    for ( ; nameItr != fNameMap.end() ; ++nameItr )
    {
        const char * name = (*nameItr).second;
        Double_t energy = (*nameItr).first;
        std::cout << "Name = " << name << "   energy = " << energy << std::endl;
        conc.AddFile( name, energy, fParticle );
    }
    conc.Run();
    return;
}

void WCSimEmissionProfileGenerator::SetSaveName(const char * savename)
{
    fSaveName = savename;
    return;
}

TChain * WCSimEmissionProfileGenerator::Next()
{
  if(fIsFirstCall)
  {
    fChainItr = fChain.begin();
    fIsFirstCall = false;
    std::cout << "First call" << std::endl;
  }

  TChain * returnMe = 0;
  if( fChainItr != fChain.end() )
  {
    returnMe = (*fChainItr).second;
    fName = fNameMap[(*fChainItr).first];
    std::cout << "Energy = " << (*fChainItr).first << "   Name = " << (*fChainItr).second << std::endl;
  }
  fChainItr++;
  std::cout << "About to return " << returnMe << std::endl;
  return returnMe;
}

void WCSimEmissionProfileGenerator::MakeFluxHistogram()
{
  TFile * fluxFile = new TFile("flux.root","RECREATE");
  TGraph * gr = new TGraph();

  std::map<Double_t, TChain *>::iterator chainItr = fChain.begin();
  for( ; chainItr != fChain.end(); ++chainItr)
  {
    std::cout << "Now processing tree " << std::distance(fChain.begin(), chainItr) << std::endl; 
    TChain * chain = (*chainItr).second;
    Int_t lastEvent = -1;
    Int_t eventID = -999;

    Int_t nPhotons = chain->GetEntries();
    Int_t nEvents  = 100;

    TBranch * br = NULL;
    chain->SetBranchAddress("eventID",&eventID, &br);
    //for(Int_t iEntry = 0; iEntry < chain->GetEntries(); ++iEntry)
    //{
    //  chain->GetEntry(iEntry, 0);
    //  if( eventID != lastEvent ){ nEvents++; }
    //}
    std::cout << "Energy = " << (*chainItr).first << "    nPhotons = " << nPhotons << "     nEvents = " << nEvents << std::endl;
    gr->SetPoint( gr->GetN(), (*chainItr).first, nPhotons/static_cast<Double_t>(nEvents) );
  }

  gr->SetName("gFlux");
  gr->SetTitle("Flux");
  gr->GetXaxis()->SetTitle("Energy / MeV");
  gr->GetYaxis()->SetTitle("Mean number of photons");
  gr->SetLineColor(kBlue);
  gr->Write();
  fluxFile->Write();
  fluxFile->Close();
  delete fluxFile;
}
