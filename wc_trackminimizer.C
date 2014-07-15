#include <vector>

void wc_trackminimizer(const char * file, Int_t skip, Int_t process)
{
    const char * prefix = gSystem->BaseName(file);

    TString savename(prefix);
    savename.ReplaceAll(".root","");
    savename.Append(TString::Format("_%d_charge.root", skip ));
    std::cout << savename.Data() << std::endl;
    TFile * fTree = new TFile(savename.Data(),"RECREATE");
    TTree * fitTree = new TTree("fitTree", "fitTree");
    Int_t event;
    Double_t fitX, fitY, fitZ, fitDistance;
    Double_t seedX, seedY, seedZ, seedDistance;
    Double_t trueX, trueY, trueZ;
    Double_t minus2LnL;
    Int_t fitStatus;
    
    fitTree->Branch("event",&event);
    fitTree->Branch("fitX",&fitX);
    fitTree->Branch("fitY",&fitY);
    fitTree->Branch("fitZ",&fitZ);
    fitTree->Branch("trueX",&trueX);
    fitTree->Branch("trueY",&trueY);
    fitTree->Branch("trueZ",&trueZ);
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
    gSystem->Load("../WCSim_github/libWCSimRoot.so");
    gSystem->Load("./lib/libWCSimAnalysis.so");
 
    // File to analyse
    //TString filename("/unix/fnu/ajperch/WCSim/macfiles/muon_1500_jul14_5.root
    //");
    TString filename(file);

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

    std::vector<Double_t> trueXvec, trueYvec, trueZvec;

    for(int iEvent = skip; iEvent < nEvts; ++iEvent)
    {
      std::cout << "Skip = " << skip << "   Process = " << process << "    iEvt = " << iEvent << std::endl;
      if (iEvent >= skip + process ) { break; } 
      // get first entry
      WCSimInterface::LoadEvent(iEvent);
      std::cerr << "Processing event " << iEvent << "/" << nEvts - 1 << std::endl; 
      WCSimRecoEvent* recoEvent = WCSimInterface::RecoEvent();
      WCSimTrueEvent* trueEvent = WCSimInterface::TrueEvent();
      trueX = trueEvent->GetVtxX();
      trueY = trueEvent->GetVtxY();
      trueZ = trueEvent->GetVtxZ();
      trueXvec.push_back(trueX);
      trueYvec.push_back(trueY);
      trueZvec.push_back(trueZ);



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
      fitDistance = TMath::Sqrt( (trueX  - fitTracks.at(0).GetX())*(trueX  - fitTracks.at(0).GetX())
                        +(trueY - fitTracks.at(0).GetY())*(trueY - fitTracks.at(0).GetY())
                        +(trueZ  - fitTracks.at(0).GetZ())*(trueZ  - fitTracks.at(0).GetZ()));
      seedX = seedTrack.GetX();
      seedY = seedTrack.GetY();
      seedZ = seedTrack.GetZ();
      seedDistance = TMath::Sqrt(  (trueX - seedX)*(trueX - seedX) + (trueY - seedY)*(trueY - seedY) 
                                 + (trueZ - seedTrack.GetZ()) * (trueZ - seedTrack.GetZ()));


      fitTree->Fill();
      delete myFitter;
    }

    for( UInt_t iEvt = 0 ; iEvt < bestFits.size(); ++iEvt)
    {
      for(UInt_t iTrack = 0 ; iTrack < (bestFits.at(iEvt)).size(); ++iTrack)
      {
        WCSimLikelihoodTrack fTrack = (bestFits.at(iEvt)).at(iTrack);
        Float_t dist = TMath::Sqrt( (trueXvec[iEvt]  - fTrack.GetX())*(trueXvec[iEvt]  - fTrack.GetX())
                        +(trueYvec[iEvt] - fTrack.GetY())*(trueYvec[iEvt] - fTrack.GetY())
                        +(trueZvec[iEvt]  - fTrack.GetZ())*(trueZvec[iEvt]  - fTrack.GetZ()));
        hReso->Fill(dist);
      }
    }
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
