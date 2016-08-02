/*
 * This macro performs a fit in order to estimate the speed of a charged
 * particle by setting up a likelihood calculator using the MC true parameters
 * and searching for the particle speed which minimises -2LnL
 * 
 */

// Wrap the total likelihood inside a class so we can use the functor version of TF1
class LikelihoodWrapper {
public:
	LikelihoodWrapper(WCSimTotalLikelihood * totalLikelihood)
	{
		fTotalLikelihood = totalLikelihood;
	}
	double Eval(double * x, double * par)
	{
		WCSimAnalysisConfig::Instance()->SetCustomParticleSpeed(*x);
		return fTotalLikelihood->Calc2LnL();
	}
private:
	WCSimTotalLikelihood * fTotalLikelihood;
};

void wc_estimateParticleSpeed(const char * infile = ""){
  // Path to WCSim ROOT file
  // =======================
  TString filename(infile);
  if(filename.CompareTo("") == 0 )
  {
    filename = TString("localfile.root");
  }
  gApplication->ProcessLine(".except");

  // Load Data
  // =========
  WCSimInterface::LoadData(filename.Data());
  int numEvents = WCSimInterface::Instance()->GetNumEvents();
  numEvents = 5;

  std::vector<double> minima(numEvents);

  // Event Loop:
  // ===========
  for( int iEvent = 0; iEvent < numEvents; ++iEvent)
  {
      // Load this event
	  WCSimInterface::Instance()->LoadEvent(iEvent);
	  WCSimLikelihoodDigitArray * myLDA = WCSimInterface::Instance()->GetWCSimLikelihoodDigitArray(iEvent);

      // Check for some PMT hits
	  int positive = 0;
	  for(int iDigit = 0; iDigit < myLDA->GetNDigits(); ++iDigit)
	  {
		if(myLDA->GetDigit(iDigit)->GetQ() > 0){ positive++; }
	  }
	  std::cout << "Hits: " << positive << std::endl;
	  assert(positive > 0);

	  // Set up a track using the MC truth information:
	  WCSimTotalLikelihood * myTotal = new WCSimTotalLikelihood(myLDA);
	  std::vector<WCSimLikelihoodTrackBase*> * trueTracks = WCSimInterface::Instance()->GetTrueLikelihoodTracks();
	  myTotal->SetTracks(*trueTracks);
	  LikelihoodWrapper * myWrapper = new LikelihoodWrapper(myTotal);

      TF1 * fitFunc = new TF1("fitFunc", myWrapper, &LikelihoodWrapper::Eval, 0, 1, 0, "LikelihoodWrapper", "Eval");
	  double min = fitFunc->GetMinimum(0,1);
	  minima.at(iEvent) = min;

	  delete myTotal;
	  delete fitFunc;
  }

  // Make histogram of all the speeds that give the best likelihood
  // ==============================================================
  TH1D * hMinima = new TH1D("hMinima", "Best-fit particle speeds;Speed of particle / c;Events", 50, 0.0, 1.0);
  std::vector<double>::const_iterator minItr = minima.begin();
  double total = 0.0;
  while(minItr != minima.end())
  {
	  hMinima->Fill(*minItr);
	  ++minItr;
	  total += *minItr;
  }

  std::cout << "Mean best-fit speed = " << total/numEvents;

  hMinima->GetXaxis()->CenterTitle();
  hMinima->GetYaxis()->CenterTitle();
  hMinima->SetLineColor(kAzure);
  hMinima->SetLineWidth(2);
  TCanvas * can = new TCanvas("can","",800, 600);
  hMinima->Draw();
  can->SaveAs("bestFitSpeeds.root");
  can->SaveAs("bestFitSpeeds.C");
  can->SaveAs("bestFitSpeeds.pdf");
  can->SaveAs("bestFitSpeeds.png");

}
