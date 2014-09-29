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

    gSystem->Load("../libWCSimRoot.so");
    WCSimEmissionProfileGenerator * myGenerator = new WCSimEmissionProfileGenerator(treename, particlePDG);
    myGenerator->SetSaveName("emissionProfiles.root");

    // Repeat this for all the energies you need
    ////////////////////////////////////////////
    myGenerator->AddFile(1000, 14, "../ep_muon_1000_0_photons.root");
    myGenerator->AddFile(1100, 14, "../ep_muon_1100_0_photons.root");
    myGenerator->AddFile(1200, 14, "../ep_muon_1200_0_photons.root");
    myGenerator->AddFile(1300, 14, "../ep_muon_1300_0_photons.root");
    myGenerator->AddFile(1400, 14, "../ep_muon_1400_0_photons.root");
    myGenerator->AddFile(1500, 14, "../ep_muon_1500_0_photons.root");
    myGenerator->AddFile(1600, 14, "../ep_muon_1600_0_photons.root");
    myGenerator->AddFile(1700, 14, "../ep_muon_1700_0_photons.root");
    myGenerator->AddFile(1800, 14, "../ep_muon_1800_0_photons.root");
    myGenerator->AddFile(1900, 14, "../ep_muon_1900_0_photons.root");
    myGenerator->AddFile(2000, 14, "../ep_muon_2000_0_photons.root");
    myGenerator->AddFile(2100, 14, "../ep_muon_2100_0_photons.root");
    myGenerator->AddFile(2200, 14, "../ep_muon_2200_0_photons.root");
    myGenerator->AddFile(2300, 14, "../ep_muon_2300_0_photons.root");
    myGenerator->AddFile(2400, 14, "../ep_muon_2400_0_photons.root");
    myGenerator->AddFile(2500, 14, "../ep_muon_2500_0_photons.root");
    myGenerator->AddFile(2600, 14, "../ep_muon_2600_0_photons.root");
    myGenerator->AddFile(2700, 14, "../ep_muon_2700_0_photons.root");
    myGenerator->AddFile(2800, 14, "../ep_muon_2800_0_photons.root");
    myGenerator->AddFile(2900, 14, "../ep_muon_2900_0_photons.root");
    myGenerator->AddFile(3000, 14, "../ep_muon_3000_0_photons.root");
    myGenerator->AddFile(3100, 14, "../ep_muon_3100_0_photons.root");
    myGenerator->AddFile(3200, 14, "../ep_muon_3200_0_photons.root");
    myGenerator->AddFile(3300, 14, "../ep_muon_3300_0_photons.root");
    myGenerator->AddFile(3400, 14, "../ep_muon_3400_0_photons.root");
    myGenerator->AddFile(3500, 14, "../ep_muon_3500_0_photons.root");
    myGenerator->AddFile(3600, 14, "../ep_muon_3600_0_photons.root");
    myGenerator->AddFile(3700, 14, "../ep_muon_3700_0_photons.root");
    myGenerator->AddFile(3800, 14, "../ep_muon_3800_0_photons.root");
    myGenerator->AddFile(3900, 14, "../ep_muon_3900_0_photons.root");
 
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


