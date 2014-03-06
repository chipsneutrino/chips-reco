void wc_trackfitter_sweep(Double_t energy)
{
    // Load WCSim libraries
    gSystem->Load("libGeom");
    gSystem->Load("libEve");
    gSystem->Load("libMinuit");
    gSystem->Load("../WCSim/libWCSimRoot.so");
    gSystem->Load("./lib/libWCSimAnalysis.so");
    gStyle->SetOptStat(0);

    // File to analyse
    TString filename("/unix/fnu/ajperch/WCSim/localfile.root");

    WCSimInterface * myInterface = WCSimInterface::Instance();
    myInterface->LoadData(filename.Data());
    int nEvts = myInterface->GetNumEvents();
    std::cout << nEvts << " = num events" << std::endl;
    Int_t nXBins = 200;
    Int_t nYBins = 200;
    Int_t xMin = 0;
    Int_t xMax = 200;
    Int_t yMin = -280;
    Int_t yMax = -80;
    Double_t minimum = 1000000000000.0;
    TH2D * chargeXZ = new TH2D("chargeXZ","Charge likelihood for an otherwise correct track, in 10cm (x,z) bins", nXBins,xMin,xMax,nYBins,yMin, yMax);
    TMarker marker;
    for(int iEvent = 0; iEvent < 1; ++iEvent)
    {

      WCSimRootEvent * myEvent = myInterface->GetWCSimEvent(iEvent);
      WCSimLikelihoodDigitArray * myLikelihoodDigitArray = new WCSimLikelihoodDigitArray(myEvent); // digitized
      WCSimChargeLikelihood * myChargeLikelihood = new WCSimChargeLikelihood(myLikelihoodDigitArray, true);
      //WCSimLikelihoodDigitArray * myLikelihoodDigitArray = new WCSimLikelihoodDigitArray(myEvent, kTRUE); // undigitized
       //WCSimLikelihoodTrack * myTrack = new WCSimLikelihoodTrack(0.,0.,0.,0.,0,0, energy, WCSimLikelihoodTrack::MuonLike);
      for(Int_t iXBin = 0; iXBin < nXBins; ++iXBin)
      {
        Double_t x = xMin + iXBin * (xMax-xMin)/static_cast<Float_t>nXBins;
//        Double_t x = 0;
//        Double_t z = 0;
        for(Int_t iYBin = 0; iYBin < nYBins; ++iYBin)
        {
          Double_t y = yMin + iYBin * (yMax - yMin)/static_cast<Float_t>nYBins;
          std::cout << std::endl << "x = " << x << "    y = " << y << std::endl;
          WCSimLikelihoodTrack * myTrack = new WCSimLikelihoodTrack(x,y,320.,0.,0.771177,-0.0245231, energy, WCSimLikelihoodTrack::MuonLike);
          std::vector<WCSimLikelihoodTrack*> myTrackVec;
          myTrackVec.push_back(myTrack);
          myChargeLikelihood->SetTracks(myTrackVec);//WCSimLikelihoodTrack::TrackType myType = WCSimLikelihoodTrack::MuonLike; 
          TStopwatch myStopwatch;
          myStopwatch.Start();
          Double_t minus2LnL = myChargeLikelihood->Calc2LnL();
          myStopwatch.Stop();
          chargeXZ->Fill(x,y,minus2LnL);
          delete myTrack;
          myTrackVec.pop_back();
          if(minus2LnL < minimum)
          {
            minimum = minus2LnL; 
            minX = x;
            minY = y;
          }
        }
      }
      delete myChargeLikelihood;
      delete myLikelihoodDigitArray;
      std::cout << "Calculate likelihoods in " << myStopwatch.Print("u");
      
    }
    std::cout << "Minimum = " << minimum << " at (x,y) = (" << minX << "," << minY << ") cm" << std::endl;
    marker.SetX(minX);
    marker.SetY(minY);
    marker.SetMarkerStyle(29);
    marker.SetMarkerSize(2);
    marker.SetMarkerColor(kAzure);
    
    TFile *f1 = new TFile("xY_calc.root","RECREATE");
    chargeXZ->Write();
    f1->Close();


    TCanvas * c1 = new TCanvas("c1","c1",1280,800);
    chargeXZ->Draw("COLZ");
    marker.Draw();
    c1->SaveAs("xY_calc.png");


  delete c1;

}
