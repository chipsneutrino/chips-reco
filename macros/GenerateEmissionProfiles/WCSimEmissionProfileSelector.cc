#define WCSimEmissionProfileSelector_cxx
// The class definition in WCSimEmissionProfileSelector.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("WCSimEmissionProfileSelector.C")
// Root > T->Process("WCSimEmissionProfileSelector.C","some options")
// Root > T->Process("WCSimEmissionProfileSelector.C+")
//

#include "WCSimEmissionProfileSelector.hh"
#include <cassert>
#include <TH2D.h>
#include <TH1D.h>
#include <TMath.h>
#include <TStyle.h>
#include <iostream>


void WCSimEmissionProfileSelector::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).
   fSaveName = TString("emissionProfile.root");
   TString option = GetOption();
   if(option.EndsWith(".root")) { fSaveName = option; }

   fNBinsS  = 250;
   fSMin    = 0;
   fSMax    = 2500;
   fNBinsTh = 1000;
   fThMin   = -1.;
   fThMax   = 1.;
   f_hS.SetNameTitle("f_hS","Emission profile in s");
   f_hSTh.SetNameTitle("f_hSTh","Emission profile in (cosTheta,s)");
   f_hSTh.SetBins(fNBinsTh, fThMin, fThMax, fNBinsS, fSMin, fSMax);
   f_hG.SetNameTitle("f_hG","Emission profile in (cosTheta,s) - normalized");
   f_hG.SetBins(fNBinsTh, fThMin, fThMax, fNBinsS, fSMin, fSMax);
   f_hS.SetBins(fNBinsS, fSMin, fSMax);

   MM_TO_CM = 0.1;

   eventID = 0;
   std::cout << "Done!" << std::endl;
}

void WCSimEmissionProfileSelector::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t WCSimEmissionProfileSelector::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either WCSimEmissionProfileSelector::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.


    // std::cout << "Getting entry" << std::endl;
    fChain->GetTree()->GetEntry(entry);
    fOverallEventNum = eventID + fChain->GetTreeNumber() * fMaxTreeEvents;

    vtxX *= MM_TO_CM;
    vtxY *= MM_TO_CM;
    vtxZ *= MM_TO_CM;

    if(entry % 100000 == 0 || entry == fChain->GetEntries() - 1)
    {
        std::cout << "Got entry " << entry << "/" << fChain->GetTree()->GetEntries() << std::endl;
    }

    if(entry == 0){ ProcessFirst(); }

    if(optical && !scattered && parentID == 1)
    {
        // std::cout << "FillComVtx" << std::endl;
        FillComVtx( vtxX, vtxY, vtxZ, energy);
        // std::cout << "FillRawDir" << std::endl;
        FillRawDir( vtxdirX, vtxdirY, vtxdirZ);
        // std::cout << "FillRawPos" << std::endl;
        // Double_t R = TMath::Sqrt(vtxX*vtxX + vtxY*vtxY + vtxZ*vtxZ);

        FillRawPos( vtxZ );
        fNEvents[fOverallEventNum]++;
    }



   return kTRUE;
}

void WCSimEmissionProfileSelector::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void WCSimEmissionProfileSelector::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
    std::cout << "Calculating CoM" << std::endl;
    CalculateCoM();
    std::cout << "Filling profiles" << std::endl;
    FillProfiles();
    std::cout << "Normalising profiles" << std::endl;
    NormaliseProfiles();
    std::cout << "Saving profiles" << std::endl;
    SaveProfiles();

    return;
}

// Calculate the QE by interpolating the same arrays
// as WCSim uses
////////////////////////////////////////////////////
Double_t WCSimEmissionProfileSelector::GetQE( const Double_t& lambda )
{
  Double_t wavelengthArr[20] = { 280., 300., 320., 340., 360., 380., 400., 420., 440., 460., 480., 500., 520., 540., 560., 580., 600., 620., 640., 660.};
  Double_t QE[20] = { 0.00, .0502, .2017, .2933, .3306, .3396, .3320, .3168, .2915, .2655, .2268,  .1971, .1641, .1102, .0727, .0499, .0323, .0178, .0061, 0.00};
  Double_t thisQE;

  if(lambda > 660 || lambda < 280) thisQE = 0;
  else
  {
    Int_t i;
    for(i = 0; i < 19; ++i)
    {   
      if(wavelengthArr[i] < lambda && wavelengthArr[i+1] > lambda)
      break;
    }   
    thisQE = QE[i] + (QE[i+1] - QE[i]) * (lambda - wavelengthArr[i]) / (wavelengthArr[i+1] - wavelengthArr[i]);
  }
  return thisQE;   
}

void WCSimEmissionProfileSelector::FillComVtx(const Double_t& x, const Double_t& y, const Double_t& z, const Double_t& E)
{
    // std::cout << "EventID = " << eventID << "    and fComVtx.size() = " << fComVtxX.size() << std::endl;
    fComVtxX.at(fOverallEventNum) += ( x * E );
    fComVtxY.at(fOverallEventNum) += ( y * E );
    fComVtxZ.at(fOverallEventNum) += ( z * E );
    return;
}

void WCSimEmissionProfileSelector::FillRawDir(const Double_t& x,
        const Double_t& y, const Double_t& z)
{
    (fRawDir.at(fOverallEventNum)).push_back(TVector3(x,y,z));
    return;
}

void WCSimEmissionProfileSelector::FillRawPos(const Double_t& z)
{
    (fRawPos.at(fOverallEventNum)).push_back(z);
    return;
}

void WCSimEmissionProfileSelector::CalculateCoM()
{
    std::vector<Double_t>::iterator vtxItr;

    std::vector<Double_t>::iterator vtxBegin[3] = { fComVtxX.begin(), fComVtxY.begin(), fComVtxZ.begin() };
    std::vector<Double_t>::iterator vtxEnd[3]   = { fComVtxX.end(),   fComVtxY.end(),   fComVtxZ.end()   };

    for( UInt_t coord = 0 ; coord < 3 ; ++coord )
    {
        for( vtxItr = vtxBegin[coord] ; vtxItr != vtxEnd[coord] ; ++ vtxItr )
        {
            UInt_t distance = std::distance( vtxBegin[coord], vtxItr);
            std::cout << "Distance = " << distance << "   fNEvents.at(distance) = " << fNEvents.at(distance) << std::endl;
            assert( distance < fNEvents.size() );
            assert( fNEvents.at(distance) > 0 );

            (*vtxItr) /= fNEvents.at(distance);
            if(coord == 2)
            {
                fComDir.at(distance) = TVector3( fComVtxX.at(distance), fComVtxY.at(distance), fComVtxZ.at(distance) );
            }
        }
    }
}

void WCSimEmissionProfileSelector::FillProfiles()
{
    assert( fRawPos.size() == fRawDir.size() );
    std::vector<std::vector<Double_t> >::iterator eventItr      = fRawPos.begin();
    std::vector<std::vector<TVector3> >::iterator dirEventItr = fRawDir.begin();

    Int_t event = 0;
    for( ; eventItr != fRawPos.end() ; ++eventItr)
    {

        assert( (*eventItr).size() == (*dirEventItr).size() );

        std::vector<Double_t>::const_iterator posItr   = (*eventItr).begin();
        std::vector<TVector3>::const_iterator dirItr = (*dirEventItr).begin();

        for( ; posItr != (*eventItr).end() ; ++posItr )
        {
            f_hS.Fill( *posItr );
            // std::cout << "theta = " << *thetaItr << "    offset = " << fThetaOffset.at(event) << std::endl;
            f_hSTh.Fill( TMath::Cos( (*dirItr).Angle( fComDir.at(event) )), *posItr );
            ++dirItr;
        }
        ++dirEventItr;
        ++event;
    }
    f_hS.Print();
    std::cout << f_hS.Integral();

    return;
}

void WCSimEmissionProfileSelector::NormaliseProfiles()
{
    f_hS.Scale(1.0 / f_hS.Integral("width") );

    Int_t nBinsTh = f_hSTh.GetNbinsX();
    for( Int_t sBin = 1; sBin <= f_hSTh.GetNbinsY(); ++sBin)
    {
      Double_t sumTheta = 0.0;

      for( Int_t thBin = 1; thBin <= nBinsTh; ++thBin)
      {
        sumTheta += f_hSTh.GetBinContent(thBin, sBin);
      }

      for( Int_t thBin = 1; thBin <= nBinsTh; ++thBin)
      {
        if(sumTheta>0) f_hG.SetBinContent(thBin, sBin, f_hSTh.GetBinContent(thBin, sBin) * nBinsTh / sumTheta);
        else f_hG.SetBinContent(thBin, sBin, 0);
      }
    }
    std::cout << f_hG.Integral("width") << std::endl;


}

void WCSimEmissionProfileSelector::SaveProfiles()
{
    TFile f(fSaveName.Data() ,"RECREATE");
    f_hS.Write();
    f_hSTh.Write();
    f_hG.Write();
    f.Close();
    return;
}

Double_t WCSimEmissionProfileSelector::Theta(const Double_t &x, const Double_t& y, const Double_t& z)
{
    Double_t R = TMath::Sqrt(x*x + y*y);
    Double_t theta = TMath::ATan2(R,z);
    return theta;
}

void WCSimEmissionProfileSelector::ProcessFirst()
{
    Int_t nTrees = 1 ;
    if(fChain->InheritsFrom("TChain")){ nTrees = ((TChain*)fChain)->GetNtrees(); }
    std::cout << "Total number of trees = " << nTrees << std::endl;

    fMaxTreeEvents = static_cast<Int_t>( fChain->GetMaximum("eventID") ) + 1;
    fEvents = nTrees * fMaxTreeEvents;
    // std::cout << "Max number of events = " << fEvents << std::endl;
    fNEvents.resize(fEvents,0);
    fRawDir.resize(fEvents);
    fRawPos.resize(fEvents);
    fComVtxX.resize(fEvents,0);
    fComVtxY.resize(fEvents,0);
    fComVtxZ.resize(fEvents,0);
    fComDir.resize(fEvents,TVector3(0,0,1));

    f_hNEvents.SetBins(fEvents, 0.0, fEvents);
    return;
}


