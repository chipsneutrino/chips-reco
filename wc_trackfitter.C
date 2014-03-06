#include <vector>

void wc_trackfitter(Double_t energy)
{
    // Load WCSim libraries
    gSystem->Load("libGeom");
    gSystem->Load("libEve");
    gSystem->Load("libMinuit");
    gSystem->Load("../WCSim/libWCSimRoot.so");
    gSystem->Load("./lib/libWCSimAnalysis.so");
 
    // File to analyse
    TString filename("/unix/fnu/ajperch/WCSim/pions.root");

    WCSimInterface * myInterface = WCSimInterface::Instance();
    myInterface->LoadData(filename.Data());
    int nEvts = myInterface->GetNumEvents();
    std::cout << nEvts << " = num events" << std::endl;

    for(int iEvent = 0; iEvent < 1; ++iEvent)
    {

      WCSimRootEvent * myEvent = myInterface->GetWCSimEvent(iEvent);
      WCSimLikelihoodDigitArray * myLikelihoodDigitArray = new WCSimLikelihoodDigitArray(myEvent); // digitized
      
      WCSimLikelihoodTrack * myTrack = new WCSimLikelihoodTrack(0.,0.,0.,0.,0.1353125395,0, energy, WCSimLikelihoodTrack::ElectronLike);
      WCSimLikelihoodTrack * myTrack2 = new WCSimLikelihoodTrack(0.,0.,0.,0.,-0.1353125395,0, energy, WCSimLikelihoodTrack::ElectronLike);
      std::vector<WCSimLikelihoodTrack*> myTrackArray;
      myTrackArray.push_back(myTrack);
      myTrackArray.push_back(myTrack2);
      
      WCSimChargeLikelihood * myChargeLikelihood = new WCSimChargeLikelihood(myLikelihoodDigitArray);
      myChargeLikelihood->SetTracks(myTrackArray);
//      Double_t minus2LnL = myChargeLikelihood->Calc2LnL();
      
      
      
      delete myTrack;
      delete myChargeLikelihood;
    }


    
}
