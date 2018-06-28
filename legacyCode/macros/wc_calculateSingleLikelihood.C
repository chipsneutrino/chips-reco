#include <vector>
void wc_calculateSingleLikelihood(const char * infile = "", int start=0){
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


  // Load Data
  // =========
  WCSimInterface::LoadData(filename.Data());
  WCSimInterface::Instance()->LoadEvent(start);
  WCSimLikelihoodDigitArray * myLDA = WCSimInterface::Instance()->GetWCSimLikelihoodDigitArray(start);
  int positive = 0;
  for(int iDigit = 0; iDigit < myLDA->GetNDigits(); ++iDigit)
  {
    if(myLDA->GetDigit(iDigit)->GetQ() > 0){ positive++; }
  }
  std::cout << "Hits: " << positive << std::endl;
  assert(positive);
  
  // Setup a likelihood calculator
  WCSimTotalLikelihood * myTotal = new WCSimTotalLikelihood(myLDA);
  myTotal->SetLikelihoodDigitArray(myLDA);

  // Get the true tracks from the event - you can substitute your own ones here
  std::vector<WCSimLikelihoodTrackBase*> * trueTracks = WCSimInterface::Instance()->GetTrueLikelihoodTracks();
  myTotal->SetTracks(*trueTracks);
    
  // Set any custom options you want
  /*
   * WCSimAnalysisConfig::Instance()->SetUseCustomParticleSpeed(true);
   * WCSimAnalysisConfig::Instance()->SetCustomParticleSpeed(0.8374);
   * WCSimAnalysisConfig::Instance()->SetUseCustomSpeedOfLight(false);
   * WCSimAnalysisConfig::Instance()->SetCustomSpeedOfLight(0.70);
  */
  
  double total2LnL = myTotal->Calc2LnL();
  std::cout << "Total 2LnL = " << total2LnL << std::endl;
  
  std::cout << "Done!" << std::endl;

}
