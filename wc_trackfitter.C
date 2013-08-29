void wc_trackfitter(Double_t energy)
{
    // Load WCSim libraries
    gSystem->Load("libGeom");
    gSystem->Load("libEve");
    gSystem->Load("libMinuit");
    gSystem->Load("../WCSim/libWCSimRoot.so");
    gSystem->Load("./lib/libWCSimAnalysis.so");

    // File to analyse
    TString filename("/unix/fnu/ajperch/sampleEvents/WCSimOutput/muons_1500MeV.root");

    WCSimInterface * myInterface = WCSimInterface::Instance();
    myInterface->LoadData(filename.Data());
    int nEvts = myInterface->GetNumEvents();
    std::cout << nEvts << " = num events" << std::endl;

    for(int iEvent = 0; iEvent < 1; ++iEvent)
    {

      WCSimRootEvent * myEvent = myInterface->GetWCSimEvent(iEvent);
      WCSimLikelihoodDigitArray * myLikelihoodDigitArray = new WCSimLikelihoodDigitArray(myEvent); // digitized
      //WCSimLikelihoodDigitArray * myLikelihoodDigitArray = new WCSimLikelihoodDigitArray(myEvent, kTRUE); // undigitized
      WCSimLikelihoodTrack * myTrack = new WCSimLikelihoodTrack(0.,0.,0.,0.,0,0, energy, WCSimLikelihoodTrack::MuonLike);
      WCSimChargeLikelihood * myChargeLikelihood = new WCSimChargeLikelihood(myLikelihoodDigitArray, myTrack);
      WCSimLikelihoodTuner * myTuner = new WCSimLikelihoodTuner( myLikelihoodDigitArray->GetExtent(0),
                                                                 myLikelihoodDigitArray->GetExtent(1),
                                                                 myLikelihoodDigitArray->GetExtent(2),
                                                                 kFALSE);
      Double_t minus2LnL = myChargeLikelihood->Calc2LnL();
      delete myTrack;
      delete myChargeLikelihood;
    }
    
}
