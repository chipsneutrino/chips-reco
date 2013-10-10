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
    TString filename("/unix/fnu/ajperch/sampleEvents/WCSimOutput/singleEvent.root");

    WCSimInterface * myInterface = WCSimInterface::Instance();
    myInterface->LoadData(filename.Data());
    int nEvts = myInterface->GetNumEvents();
    std::cout << nEvts << " = num events" << std::endl;

    for(int iEvent = 0; iEvent < 1; ++iEvent)
    {

      WCSimRootEvent * myEvent = myInterface->GetWCSimEvent(iEvent);
      WCSimLikelihoodDigitArray * myLikelihoodDigitArray = new WCSimLikelihoodDigitArray(myEvent); // digitized
      
      WCSimLikelihoodTrack * myTrack = new WCSimLikelihoodTrack(0.,0.,0.,0.,-0.5,0, energy, WCSimLikelihoodTrack::MuonLike);
      WCSimLikelihoodTrack * myTrack2 = new WCSimLikelihoodTrack(0.,0.,0.,0.,0.5,0, energy, WCSimLikelihoodTrack::MuonLike);
      std::vector<WCSimLikelihoodTrack*> myTrackArray;
      myTrackArray.push_back(myTrack2);
      myTrackArray.push_back(myTrack);
      
      WCSimChargeLikelihood * myChargeLikelihood = new WCSimChargeLikelihood(myLikelihoodDigitArray, false);
      myChargeLikelihood->SetTracks(myTrackArray);
      WCSimLikelihoodTuner * myTuner = new WCSimLikelihoodTuner( myLikelihoodDigitArray->GetExtent(0),
                                                                 myLikelihoodDigitArray->GetExtent(1),
                                                                 myLikelihoodDigitArray->GetExtent(2),
                                                                 kFALSE);
      Double_t minus2LnL = myChargeLikelihood->Calc2LnL();
      
      
      
      delete myTrack;
      delete myChargeLikelihood;
    }


    
}
