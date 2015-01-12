#include <map>
/*
 * RunEmissionProfileGenerator.C
 *
 *  Created on: 28 Mar 2014
 *      Author: andy
 */
void RunEmissionProfileGenerator()
{
    const char * treename = "ntuple";
    Int_t particlePDG = 14;

    /*
    const char * WCSimDir = gSystem->Getenv("WCSIMHOME");
    TString toLoad(WCSimDir);
    toLoad = toLoad + "/libWCSimRoot.so";
    gSystem->Load(toLoad.Data());
    */

    gApplication->ProcessLine(".L WCSimConcatenateEmissionProfiles.cc++");
    gApplication->ProcessLine(".L WCSimEmissionProfileGenerator.cc++");
    WCSimEmissionProfileGenerator * myGenerator = new WCSimEmissionProfileGenerator(treename, particlePDG);
    myGenerator->SetSaveName("emissionProfiles.root");

    // Repeat this for all the energies you need
    ////////////////////////////////////////////
    //                  Energy, PDG ID, path
    myGenerator->AddFile(1250,    14, "/unix/fnu/ajperch/ep_jan5/ep_muon_1250_fast_0_photons.root");
    myGenerator->AddFile(1500,    14, "/unix/fnu/ajperch/ep_jan5/ep_muon_1500_fast_0_photons.root");
    myGenerator->AddFile(1750,    14, "/unix/fnu/ajperch/ep_jan5/ep_muon_1750_fast_0_photons.root");
 
     /* ...
       ...
       ...
    */


     TChain * chain = 0;
     do
     {
       chain = myGenerator->Next();
       if( chain ) { chain->Process("WCSimEmissionProfileSelector.cc+", myGenerator->GetChainSaveName()); }
     }
     while( chain );
     ////////////////////////////////////////////
     
     myGenerator->Concatenate();

  myGenerator->MakeFluxHistogram();
}


