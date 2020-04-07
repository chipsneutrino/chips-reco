{

    // Load WCSim libraries
    gSystem->Load("libGeom");
    gSystem->Load("libEve");
    TString libWCSimRoot = TString::Format("%s%s",gSystem->Getenv("WCSIMHOME"), "/libWCSimRoot.so");
    TString libWCSimAnalysis = TString::Format("%s%s",gSystem->Getenv("WCSIMANAHOME"), "/lib/libWCSimAnalysis.so");
    gSystem->Load(libWCSimRoot.Data());
    gSystem->Load(libWCSimAnalysis.Data());
    
    // File to analyse
    TString filename("/home/ajperch/CHIPS/WCSimAnalysis/input_files/muon_1500_xz_6jan_recoTest.root");

    WCSimInterface * myInterface = WCSimInterface::Instance();
    myInterface->LoadData(filename.Data());
    int nEvts = myInterface->GetNumEvents();
    
    //  str->Form("likelihood_%f_%f.root", fEnergy, fTrack->GetZ());
    TFile * myFile = new TFile("../out/testminim.root","RECREATE");
    myFile->cd();
    TTree * t = new TTree("likelihood","likelihood");
    Double_t X, Y, Z, Qundigi, Qdigi;
    
    Int_t iEvent = 0;
    t->Branch("event",&iEvent,"event/I");
    t->Branch("X",&X,"X/D");
    t->Branch("Y",&Y,"Y/D");
    t->Branch("Z",&Z,"Z/D");
    t->Branch("Qundigi",&Qundigi,"Qundigi/D");
    t->Branch("Qdigi",&Qdigi,"Qdigi/D");

    for(iEvent = 0; iEvent < nEvts; ++iEvent)
    {
      std::cout << "Processing event " << iEvent << "/" << nEvts << std::endl;
      WCSimRootEvent * myEvent = myInterface->GetWCSimEvent(iEvent);
      WCSimLikelihoodDigitArray * myLikelihoodDigitArray = new WCSimLikelihoodDigitArray(myEvent, kTRUE);
      WCSimLikelihoodDigitArray * myLikelihoodDigitDigiArray = new WCSimLikelihoodDigitArray(myEvent);
 
      for(int iDigit = 0; iDigit < myLikelihoodDigitDigiArray->GetNDigits(); ++iDigit)
      {
        WCSimLikelihoodDigit * fDigit = (WCSimLikelihoodDigit *)myLikelihoodDigitArray->GetDigit(iDigit);
        X = fDigit->GetX();
        Y = fDigit->GetY();
        Z = fDigit->GetZ();
        
        Qundigi = fDigit->GetQ()/static_cast<Float_t>(nEvts);
        fDigit = myLikelihoodDigitDigiArray->GetDigit(iDigit);
        Qdigi = fDigit->GetQ()/static_cast<Float_t>(nEvts);
    //    std::cout << Qundigi << "   " << Qdigi << std::endl;
        t->Fill();
      }
    }
      t->Write();
      myFile->Write();
      myFile->Close();
}
  
