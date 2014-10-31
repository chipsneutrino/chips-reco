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

    gApplication->ProcessLine(".L WCSimConcatenateEmissionProfiles.cc+");
    gApplication->ProcessLine(".L WCSimEmissionProfileGenerator.cc+");
    WCSimEmissionProfileGenerator * myGenerator = new WCSimEmissionProfileGenerator(treename, particlePDG);
    myGenerator->SetSaveName("emissionProfiles.root");

    // Repeat this for all the energies you need
    ////////////////////////////////////////////
    //                  Energy, PDG ID, path
    myGenerator->AddFile(1000,    14, "./out/ep_muon_1000_0_photons.root");
    myGenerator->AddFile(1250,    14, "./out/ep_muon_1250_0_photons.root");
    myGenerator->AddFile(1500,    14, "./out/ep_muon_1500_0_photons.root");
    myGenerator->AddFile(1750,    14, "./out/ep_muon_1750_0_photons.root");
    myGenerator->AddFile(2000,    14, "./out/ep_muon_2000_0_photons.root");
    myGenerator->AddFile(2250,    14, "./out/ep_muon_2250_0_photons.root");
    myGenerator->AddFile(2500,    14, "./out/ep_muon_2500_0_photons.root");
    myGenerator->AddFile(2750,    14, "./out/ep_muon_2750_0_photons.root");
    myGenerator->AddFile(3000,    14, "./out/ep_muon_3000_0_photons.root");
    myGenerator->AddFile(3250,    14, "./out/ep_muon_3250_0_photons.root");
 
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


