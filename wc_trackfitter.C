void wc_trackfitter()
{
    // Load WCSim libraries
    gSystem->Load("libGeom");
    gSystem->Load("libEve");
    gSystem->Load("libMinuit");
    gSystem->Load("../WCSim/libWCSimRoot.so");
    gSystem->Load("./lib/libWCSimAnalysis.so");

    // File to analyse
    TString filename("/home/ajperch/MINOS/waterCherenkov/mcSamples/numu.root");

    WCSimInterface * myInterface = WCSimInterface::Instance();
    myInterface->LoadData(filename.Data());
    int nEvts = myInterface->GetNumEvents();
    std::cout << nEvts << " = num events" << std::endl;
    for(int iEvent = 10; iEvent < 11; ++iEvent)
    {
      WCSimRootEvent * myEvent = myInterface->GetWCSimEvent(iEvent);
      WCSimLikelihoodFitter * myFitter = new WCSimLikelihoodFitter(myEvent);
      myFitter->Minimize2LnL();
    }
}
