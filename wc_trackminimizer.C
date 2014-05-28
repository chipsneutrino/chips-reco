#include <vector>

void wc_trackminimizer()
{
    TFile * fTree = new TFile("fitTree_newQE.root","RECREATE");
    TTree * fitTree = new TTree("fitTree", "fitTree");
    Int_t event;
    Double_t fitX, fitY, fitZ, fitDistance;
    Double_t seedX, seedY, seedZ, seedDistance;
    Double_t minus2LnL;
    Int_t fitStatus;
    
    fitTree->Branch("event",&event);
    fitTree->Branch("fitX",&fitX);
    fitTree->Branch("fitY",&fitY);
    fitTree->Branch("fitZ",&fitZ);
    fitTree->Branch("seedX",&seedX);
    fitTree->Branch("seedY",&seedY);
    fitTree->Branch("seedZ",&seedZ);
    fitTree->Branch("minus2LnL",&minus2LnL);
    fitTree->Branch("fitDistance",&fitDistance);
    fitTree->Branch("seedDistance",&seedDistance);
    fitTree->Branch("fitStatus",&fitStatus);


    // Load WCSim libraries
    gSystem->Load("libGeom");
    gSystem->Load("libEve");
    gSystem->Load("libMinuit");
    gSystem->Load("/home/mpf/wcsim/WCSim/libWCSimRoot.so");
    gSystem->Load("/home/mpf/wcsim/WCSimAnalysis/lib/libWCSimAnalysis.so");
 
    // File to analyse
    //TString filename("/unix/fnu/ajperch/WCSim/ep_muons_1500_QE_test.root");
    TString filename("/unix/fnu/ajperch/WCSim/ep_muon_1500_*.root");

    // Resolution histograms
    TH1D * hReso = new TH1D("hReso","Vertex resolution plot", 100, 0, 200);
    
    // Initialize Reconstruction
    // =========================
    WCSimInterface * myInterface = WCSimInterface::Instance();
    myInterface->LoadData(filename.Data());
    WCSimReco* myReco = WCSimRecoFactory::Instance()->MakeReco();

    // Reconstruct Events
    // ==================
    int nEvts = myInterface->GetNumEvents();
    // if(nEvts > 500) { nEvts = 500; }
    std::cout << nEvts << " = num events" << std::endl;

    // The best fit track results
    std::vector<std::vector<WCSimLikelihoodTrack>> bestFits;

    for(int iEvent = 0; iEvent < nEvts; ++iEvent)
    {
      // get first entry
      WCSimInterface::LoadEvent(iEvent);
      std::cerr << "Processing event " << iEvent << "/" << nEvts - 1 << std::endl; 
      WCSimRecoEvent* recoEvent = WCSimInterface::RecoEvent();
      WCSimTrueEvent* trueEvent = WCSimInterface::TrueEvent();

      WCSimRootEvent * myEvent = myInterface->GetWCSimEvent(iEvent);
      WCSimLikelihoodFitter * myFitter = new WCSimLikelihoodFitter(myEvent);
      myFitter->SeedParams( myReco );
      WCSimLikelihoodTrack seedTrack = myFitter->GetSeedParams();
      myFitter->Minimize2LnL(1);
      fitStatus = myFitter->GetStatus();
      std::vector<WCSimLikelihoodTrack> fitTracks = myFitter->GetBestFit();
      bestFits.push_back(fitTracks);
      
      event = iEvent;
      fitX = (fitTracks.at(0)).GetX();
      fitY = (fitTracks.at(0)).GetY();
      fitZ = (fitTracks.at(0)).GetZ();
      minus2LnL = myFitter->GetMinimum();
      fitDistance = TMath::Sqrt( (0.  - fitTracks.at(0).GetX())*(0.  - fitTracks.at(0).GetX())
                        +(0. - fitTracks.at(0).GetY())*(0. - fitTracks.at(0).GetY())
                        +(0.  - fitTracks.at(0).GetZ())*(0.  - fitTracks.at(0).GetZ()));
      seedX = seedTrack.GetX();
      seedY = seedTrack.GetY();
      seedZ = seedTrack.GetZ();
      seedDistance = TMath::Sqrt(  (0. - seedX)*(0. - seedX) + (0 - seedY)*(0 - seedY) 
                                 + (0 - seedTrack.GetZ()) * (0 - seedTrack.GetZ()));


      fitTree->Fill();
      delete myFitter;
    }

    for( UInt_t iEvt = 0 ; iEvt < bestFits.size(); ++iEvt)
    {
      for(UInt_t iTrack = 0 ; iTrack < (bestFits.at(iEvt)).size(); ++iTrack)
      {
        WCSimLikelihoodTrack fTrack = (bestFits.at(iEvt)).at(iTrack);
        Float_t dist = TMath::Sqrt( (0.  - fTrack.GetX())*(0.  - fTrack.GetX())
                        +(0. - fTrack.GetY())*(0. - fTrack.GetY())
                        +(0.  - fTrack.GetZ())*(0.  - fTrack.GetZ()));
        hReso->Fill(dist);
      }
    }
    std::cout << "True values were (0,0,0,0) (0,0) 1500MeV" << std::endl;
    fTree->cd();
    fitTree->Write();


    hReso->Write();
    fTree->Close();

//    TFile seedMarkerFile("seedMarker.root","RECREATE");
//    TMarker seedMarker(seedTrack.GetX(), seedTrack.GetY(), 0);
//    seedMarker.SetMarkerSize(1.8);
//    seedMarker.SetMarkerColor(kGreen-2);
//    seedMarker.SetMarkerStyle(29);
//    seedMarker.Write();
//    seedMarkerFile.Close();
//    TFile fitMarkerFile("fitMarker.root","RECREATE");
//    TMarker fitMarker((fitTracks.at(0)).GetX(), (fitTracks.at(0)).GetY(), 1);
//    fitMarker.SetMarkerSize(1.8);
//    fitMarker.SetMarkerColor(kBlue-2);
//    fitMarker.SetMarkerStyle(29);
//    fitMarker.Write();
//    fitMarkerFile.Close();

    
}
