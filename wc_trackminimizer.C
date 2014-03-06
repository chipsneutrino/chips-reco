#include <vector>

void wc_trackminimizer()
{
    TFile fTree("fitTree500.root","RECREATE");
    TTree fitTree("fitTree", "fitTree");
    Int_t event;
    Double_t fitX, fitY, fitZ, fitDistance;
    Double_t seedX, seedY, seedZ, seedDistance;
    Double_t minus2LnL;
    
    fitTree.Branch("event",&event);
    fitTree.Branch("fitX",&fitX);
    fitTree.Branch("fitY",&fitY);
    fitTree.Branch("fitZ",&fitZ);
    fitTree.Branch("seedX",&seedX);
    fitTree.Branch("seedY",&seedY);
    fitTree.Branch("seedZ",&seedZ);
    fitTree.Branch("minus2LnL",&minus2LnL);
    fitTree.Branch("fitDistance",&fitDistance);
    fitTree.Branch("seedDistance",&seedDistance);


    // Load WCSim libraries
    gSystem->Load("libGeom");
    gSystem->Load("libEve");
    gSystem->Load("libMinuit");
    gSystem->Load("../WCSim/libWCSimRoot.so");
    gSystem->Load("./lib/libWCSimAnalysis.so");
 
    // File to analyse
    TString filename("/unix/fnu/ajperch/WCSim/localfile500.root");

    // Resolution histograms
    TH1D hReso("hReso","Vertex resolution plot", 100, 0, 200);
    
    // Initialize Reconstruction
    // =========================
    WCSimInterface * myInterface = WCSimInterface::Instance();
    myInterface->LoadData(filename.Data());
    WCSimReco* myReco = WCSimRecoFactory::Instance()->MakeReco();

    // Reconstruct Events
    // ==================
    int nEvts = myInterface->GetNumEvents();
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
      std::vector<WCSimLikelihoodTrack> fitTracks = myFitter->GetBestFit();
      bestFits.push_back(fitTracks);
      
      event = iEvent;
      fitX = (fitTracks.at(0)).GetX();
      fitY = (fitTracks.at(0)).GetY();
      fitZ = (fitTracks.at(0)).GetZ();
      minus2LnL = myFitter->GetMinimum();
      fitDistance = TMath::Sqrt( (100.  - fitTracks.at(0).GetX())*(100.  - fitTracks.at(0).GetX())
                        +(-200. - fitTracks.at(0).GetY())*(-200. - fitTracks.at(0).GetY())
                        +(320.  - fitTracks.at(0).GetZ())*(320.  - fitTracks.at(0).GetZ()));
      seedX = seedTrack.GetX();
      seedY = seedTrack.GetY();
      seedZ = seedTrack.GetZ();
      seedDistance = TMath::Sqrt(  (100. - seedX)*(100. - seedX) + (-200 - seedY)*(-200 - seedY) 
                                 + (320 - seedTrack.GetZ()) * (320 - seedTrack.GetZ()));


      fitTree.Fill();
      delete myFitter;
    }

    for( UInt_t iEvt = 0 ; iEvt < bestFits.size(); ++iEvt)
    {
      for(UInt_t iTrack = 0 ; iTrack < (bestFits.at(iEvt)).size(); ++iTrack)
      {
        WCSimLikelihoodTrack fTrack = (bestFits.at(iEvt)).at(iTrack);
        Float_t dist = TMath::Sqrt( (100.  - fTrack.GetX())*(100.  - fTrack.GetX())
                        +(-200. - fTrack.GetY())*(-200. - fTrack.GetY())
                        +(320.  - fTrack.GetZ())*(320.  - fTrack.GetZ()));
        hReso.Fill(dist);
      }



    }
    std::cout << "True values were (100,-200,320,0) (pi/4,0) 1500MeV" << std::endl;
    fTree.cd();
    fitTree.Write();
    fTree.Close();


    TFile fOut("resolution.root","RECREATE");
    hReso.Write();
    fOut.Close();

    TFile seedMarkerFile("seedMarker.root","RECREATE");
    TMarker seedMarker(seedTrack.GetX(), seedTrack.GetY(), 0);
    seedMarker.SetMarkerSize(1.8);
    seedMarker.SetMarkerColor(kGreen-2);
    seedMarker.SetMarkerStyle(29);
    seedMarker.Write();
    seedMarkerFile.Close();
    TFile fitMarkerFile("fitMarker.root","RECREATE");
    TMarker fitMarker((fitTracks.at(0)).GetX(), (fitTracks.at(0)).GetY(), 1);
    fitMarker.SetMarkerSize(1.8);
    fitMarker.SetMarkerColor(kBlue-2);
    fitMarker.SetMarkerStyle(29);
    fitMarker.Write();
    fitMarkerFile.Close();

    
}
