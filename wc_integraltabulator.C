void wc_integraltabulator()
{
    // Load WCSim libraries
    gSystem->Load("libGeom");
    gSystem->Load("libEve");
    gSystem->Load("libMinuit");
    gSystem->Load("../WCSim/libWCSimRoot.so");
    gSystem->Load("./lib/libWCSimAnalysis.so");
    WCSimLikelihoodTuner * myTuner = new WCSimLikelihoodTuner(false);
    WCSimLikelihoodTrack::TrackType myType = WCSimLikelihoodTrack::MuonLike;
    TString filename("tabulatedIntegralsRhoG.bin");
    myTuner->TabulateIntegrals(myType, filename);
    //WCSimLikelihoodTrack * myTrack = new WCSimLikelihoodTrack(0,0,0,0,0,0,1500,WCSimLikelihoodTrack::MuonLike);
    //myTuner->LoadTabulatedIntegrals(myTrack);
}
