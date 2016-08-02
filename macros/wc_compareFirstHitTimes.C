void wc_compareFirstHitTimes(const char * infile)
{
  // Path to WCSim ROOT file
  // =======================
  TString filename(infile);
  if(filename.CompareTo("") == 0 )
  {
    filename = TString("localfile.root");
  }
  gApplication->ProcessLine(".except");

  // Load libraries
  // ==============
  gSystem->Load("libGeom");
  gSystem->Load("libEve");
  gSystem->Load("libMinuit");

  TString libWCSimRoot = TString::Format("%s%s",gSystem->Getenv("WCSIMHOME"), "/libWCSimRoot.so");
  TString libWCSimAnalysis = TString::Format("%s%s",gSystem->Getenv("WCSIMANAHOME"), "/lib/libWCSimAnalysis.so");
  gSystem->Load(libWCSimRoot.Data());
  gSystem->Load(libWCSimAnalysis.Data());
  gApplication->ProcessLine("#include <vector>");

  // Load Data
  // =========
  WCSimInterface::LoadData(filename.Data());
  WCSimInterface::Instance()->LoadEvent(0);
  WCSimLikelihoodDigitArray * myLDA = WCSimInterface::Instance()->GetWCSimLikelihoodDigitArray(0);
  int positive = 0;
  for(int iDigit = 0; iDigit < myLDA->GetNDigits(); ++iDigit)
  {
    if(myLDA->GetDigit(iDigit)->GetQ() > 0){ positive++; }
  }
  std::cout << "Hits: " << positive << std::endl;
  assert(positive);
  

  WCSimTotalLikelihood * myTotal = new WCSimTotalLikelihood(myLDA);
  std::vector<TH1D> allHists;

  TH1D hAllPMT("hAllPMT","Predicted - true times, all events, all digits;t_{pred} - t_{hit};PMTs", 400,-20,20);
  TH1D hAllPMTA("hAllPMT","Predicted - true times, all events, digits in ring;t_{pred} - t_{hit};PMTs", 400,-20,20);
  hAllPMT.SetLineColor(kRed);
  hAllPMTA.SetLineColor(kBlue);
  hAllPMT.SetFillColor(kRed);
  hAllPMTA.SetFillColor(kBlue);
  hAllPMT.SetFillStyle(3002);
  hAllPMTA.SetFillStyle(3002);
  
  TStopwatch sw2; // To time how long this takes per event
  for(int iEvent = 0; iEvent < WCSimInterface::GetNumEvents(); ++iEvent)
  {
    sw2.Start();
    
    // Declare and format our histograms
    // =================================
    TString   histName1 = TString::Format("hPredictMinusTrue_%d", iEvent);
    TString   histName2 = TString::Format("hPredictMinusTrueAllowed_%d", iEvent);
    TH1D hPMT(histName1.Data(),"Predicted - true times (all digits);t_{pred} - t_{hit} (ns);PMTs", 200, -10, 10);
    TH1D hPMTA(histName2.Data(),"Predicted - true times (digits in ring);t_{pred} - t_{hit} (ns);PMTs", 200, -10, 10);
    hPMT.SetLineColor(kRed);
    hPMTA.SetLineColor(kBlue);
    hPMT.SetFillColor(kRed);
    hPMTA.SetFillColor(kBlue);
    hPMT.SetFillStyle(3002);
    hPMTA.SetFillStyle(3002);

    // Read in the event
    // =================
    WCSimInterface::Instance()->LoadEvent(iEvent);
    WCSimLikelihoodDigitArray * myLDA = WCSimInterface::Instance()->GetWCSimLikelihoodDigitArray(0);
    myTotal->SetLikelihoodDigitArray(myLDA);

    // Set the likelihood to calculate using the truth information
    // ===========================================================
	std::vector<WCSimLikelihoodTrackBase*> * trueTracks = WCSimInterface::Instance()->GetTrueLikelihoodTracks();
	myTotal->SetTracks(*trueTracks);
    trueTracks->at(0)->Print();
    
    // Do we want to fiddle with the speeds?
    // =====================================
    WCSimAnalysisConfig::Instance()->SetUseCustomParticleSpeed(false);
	WCSimAnalysisConfig::Instance()->SetCustomParticleSpeed(1.0);
    WCSimAnalysisConfig::Instance()->SetUseCustomSpeedOfLight(false);
	WCSimAnalysisConfig::Instance()->SetCustomSpeedOfLight(0.7119);

    // Calculate the minimum Z at which a PMT could be hit by a photon emitted at 41 degrees
    // =====================================================================================
    double radius  = myLDA->GetExtent(1);
    double vtxZ    = trueTracks->at(0)->GetZ();
    double thetaCh = 41. * TMath::DegToRad();
    double minHeight = vtxZ + radius / TMath::Tan(thetaCh);

    // Actually do the calculation, the loop through individual digits comparing the times
    // ===================================================================================
    double total2LnL = myTotal->Calc2LnL();
    std::vector<double> pTime = myTotal->GetPredictedTimeVector();

    for(size_t iT = 0; iT < pTime.size(); ++iT)
    {
      double time = pTime.at(iT);
      WCSimLikelihoodDigit * myDigit = myLDA->GetDigit(iT);
      if(time <= 0 || myDigit->GetZ() > 950) { continue; }
      if( time > 0 && myDigit->GetT() > 0){ 
        std::cout << "Digit " << iT << " time = " << myDigit->GetT() << " pred = " << time << std::endl;
      }

      hPMT.Fill(time - myDigit->GetT());
      hAllPMT.Fill(time - myDigit->GetT());

      if( myDigit->GetZ() > minHeight )
      {
          hPMTA.Fill(time - myDigit->GetT());
          hAllPMTA.Fill(time - myDigit->GetT());
      }
    }

    // We'll save them all later
    // =========================
    hPMT.SetDirectory(0);
    hPMTA.SetDirectory(0);
    allHists.push_back(hPMT);
    allHists.push_back(hPMTA);
    std::cout << "PROCESSING EVENT " << iEvent << " took " << sw2.RealTime() << " seconds" << std::endl;
    sw2.Reset();
  }
  hAllPMT.SetDirectory(0);
  hAllPMTA.SetDirectory(0);
  allHists.push_back(hAllPMT);
  allHists.push_back(hAllPMTA);


  // Save and delete the histograms
  // ==============================
  TFile * f = new TFile("fCompareTimes.root","RECREATE");
  for(size_t iH = 0; iH < allHists.size(); ++iH)
  {
    allHists.at(iH).Write();
  }
  allHists.clear();
  delete myTotal;
  f->Close();
  delete f;


  std::cout << "Done!" << std::endl;
}
