void wc_trackfitter(Double_t energy)
{
    // Load WCSim libraries
    gSystem->Load("libGeom");
    gSystem->Load("libEve");
    gSystem->Load("libMinuit");
    gSystem->Load("../WCSim/libWCSimRoot.so");
    gSystem->Load("./lib/libWCSimAnalysis.so");

    // File to analyse
    TString filename("/home/ajperch/MINOS/waterCherenkov/mcSamples/muons_various.root");

    WCSimInterface * myInterface = WCSimInterface::Instance();
    myInterface->LoadData(filename.Data());
    int nEvts = myInterface->GetNumEvents();
    std::cout << nEvts << " = num events" << std::endl;
    for(int iEvent = 6; iEvent < 7; ++iEvent)
    {

      WCSimRootEvent * myEvent = myInterface->GetWCSimEvent(iEvent);
      WCSimLikelihoodFitter * myFitter = new WCSimLikelihoodFitter(myEvent);
      WCSimLikelihoodTrack * myTrack = new WCSimLikelihoodTrack(0.,0.,500.,0.,0.0, 0.0, energy, WCSimLikelihoodTrack::MuonLike);
//      if( iEvent == 1 )WCSimLikelihoodTrack * myTrack = new WCSimLikelihoodTrack(0.,0.,0.,0.,0.0, 0.0, 1000, WCSimLikelihoodTrack::MuonLike);
//      if( iEvent == 2 )WCSimLikelihoodTrack * myTrack = new WCSimLikelihoodTrack(0.,0.,0.,0.,0.0, 0.0, 1500, WCSimLikelihoodTrack::MuonLike);
//      if( iEvent == 3 )WCSimLikelihoodTrack * myTrack = new WCSimLikelihoodTrack(0.,0.,0.,0.,0.0, 0.0, 2000, WCSimLikelihoodTrack::MuonLike);
//      if( iEvent == 4 )WCSimLikelihoodTrack * myTrack = new WCSimLikelihoodTrack(0.,0.,0.,0.,0.0, 0.0, 2500, WCSimLikelihoodTrack::MuonLike);
//      if( iEvent == 5 )WCSimLikelihoodTrack * myTrack = new WCSimLikelihoodTrack(0.,0.,0.,0.,0.0, 0.0, 3000, WCSimLikelihoodTrack::MuonLike);
//      if( iEvent == 6 )WCSimLikelihoodTrack * myTrack = new WCSimLikelihoodTrack(0.,0.,500.,0.,0.0, 0.0, 500, WCSimLikelihoodTrack::MuonLike);
//      if( iEvent == 7 )WCSimLikelihoodTrack * myTrack = new WCSimLikelihoodTrack(0.,0.,800.,0.,0.0, 0.0, 500, WCSimLikelihoodTrack::MuonLike);
//      if( iEvent == 8 )WCSimLikelihoodTrack * myTrack = new WCSimLikelihoodTrack(0.,0.,-500.,0.,0.0, 0.0, 500, WCSimLikelihoodTrack::MuonLike);
      myFitter->Calc2LnL(myTrack);
      delete myTrack;
      delete myEvent;
      delete myFitter;
    }
}
